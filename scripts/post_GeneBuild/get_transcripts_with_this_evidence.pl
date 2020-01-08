#!/usr/bin/env perl


# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


use warnings ;
use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Bio::SeqIO;
use Getopt::Long qw(:config no_ignore_case);

my $file;

my $dbhost = 'ecs2d';
my $dbuser = 'ensro';
my $dbname = 'homo_sapiens_core_9_30';
my $dbpass = undef;

my $path;
my $genetype;


my $protein_id;
my $dna_id;

GetOptions(
	    'protein_id:s'    => \$protein_id,
	    'dna_id:s'        => \$dna_id,
	    'host|dbhost|h:s'        => \$dbhost,
	    'dbname|db|D:s'        => \$dbname,
	    'genetype:s'      => \$genetype,
            'path|cs_version:s'          => \$path,
);

unless ( $protein_id || $dna_id ){
  print STDERR "script to print out all the transcripts with a given id as evidence\n";
 
  print STDERR "Usage: $0 -dbname -dbhost ( -protein_id OR -dna_id ) [ -genetype (optional)]\n"; 
  exit(0);
}

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $dbhost,
					    '-user'   => $dbuser,
					    '-dbname' => $dbname,
					    '-pass'   => $dbpass,
					   );



print STDERR "connected to $dbname : $dbhost\n";

unless ($path){
  $path = $db->assembly_type;
}

print STDERR  "path = $path\n";

if ( $protein_id ){
  print "Transcripts based on $protein_id\n";
  my @transcripts;
 
 TRAN:
  foreach my $tran (  Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->find_transcripts_by_protein_evidence($protein_id,$db,$genetype) ){
      &show_transcript($db,$tran->dbID);
  }
}

if ( $dna_id ){
    print "Transcripts based on $dna_id\n";
    my @transcripts;
  TRAN:
    foreach my $tran ( Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->find_transcripts_by_dna_evidence($dna_id,$db,$genetype) ){
	&show_transcript($db,$tran->dbID);
    }
}

############################################################

sub show_transcript{
  my ($db,$t_id) = @_;
  my $tran = $db->get_TranscriptAdaptor->fetch_by_dbID($t_id);
  my $tran_id   = $tran->stable_id || $tran->dbID;
  my $transl_id;
  if ( $tran->translation ){
    $transl_id = $tran->translation->stable_id || $tran->translation->dbID;
  }
  my $gene_id;
  if ( $tran->stable_id ){
    $gene_id = &get_gene_stable_id($db,$tran->dbID);
  }
  else{
    $gene_id = &get_gene_id($db,$tran->dbID);
  }
  
  print "Gene:$gene_id\tTranscript:$tran_id\tPeptide:$transl_id\n";

  my @exons = @{$tran->get_all_Exons};
  print "contig_id\tcontig_name\texon_id\tstart\tend\tphase\tend_phase\tstrand\tlength\n";
  foreach my $exon (@exons){
    
    # if exon is sticky print each component
    if ( $exon->isa('Bio::EnsEMBL::StickyExon') ){
      
      foreach my $exon_c ( @{$exon->get_all_component_Exons} ){
	my $length = $exon_c->end - $exon_c->start +1;
	print $exon_c->contig->dbID."\t".$exon_c->contig->name."\t".
	  $exon_c->dbID."\t".$exon_c->start."\t".$exon_c->end."\t".$exon_c->phase."\t".
	    $exon_c->end_phase."\t".$exon_c->strand."\t".$length."\n";
	&print_evidence($exon_c);
	print "\n";
      }
    }
    else{
      my $length = $exon->end - $exon->start +1;
      print $exon->contig->dbID."\t".$exon->contig->name."\t".
	$exon->dbID."\t".$exon->start."\t".$exon->end."\t".$exon->phase."\t".
	  $exon->end_phase."\t".$exon->strand."\t".$length."\n";
      &print_evidence($exon);
      print "\n";
    }
  }
}
############################################################

sub print_evidence{
  my $exon = shift;
  my @evidence = @{$exon->get_all_supporting_features};
  if ( @evidence ){
    foreach my $evi ( @evidence ){
      my $length = $evi->end - $evi->start + 1;
      my $hlength = $evi->hend - $evi->hstart + 1;
      print "Evidence: ".$evi->dbID."\t".$evi->contig->dbID."\t".$evi->contig->name."\t".
	$evi->start."-".$evi->end."\t".$evi->phase."\t".
	  $evi->end_phase."\t".$evi->strand."\t".$length."\t".
	    $evi->hstart."-".$evi->hend."\t".$hlength."\t".$evi->hseqname."\n";
    }
  }
  else{
    print "No evidence\n";
  }
}


############################################################

sub get_gene_stable_id{
  my $db = shift;
  my $t_id = shift;
  
  my $q = qq( SELECT gs.stable_id
	      FROM   gene_stable_id gs, transcript t, gene g
	      WHERE  t.transcript_id = $t_id AND
	             t.gene_id       = g.gene_id AND
        	     gs.gene_id      = g.gene_id
	    );
  
  my $sth = $db->prepare($q) || $db->throw("can't prepare: $q");
  my $res = $sth->execute || $db->throw("can't execute: $q");
  
  return  $sth->fetchrow_array;
}

############################################################

sub get_gene_id{
  my $db = shift;
  my $t_id = shift;
  
  my $q = qq( SELECT g.gene_id
	      FROM   transcript t, gene g
	      WHERE  t.transcript_id = $t_id     AND
	             t.gene_id       = g.gene_id
	    );
  
  my $sth = $db->prepare($q) || $db->throw("can't prepare: $q");
  my $res = $sth->execute || $db->throw("can't execute: $q");
  
  return  $sth->fetchrow_array;
}
