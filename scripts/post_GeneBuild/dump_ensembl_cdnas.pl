#!/usr/bin/env perl


# Copyright [1999-2013] Genome Research Ltd. and the EMBL-European Bioinformatics Institute
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


=head1 NAME

=head1 DESCRIPTION

dumps in Fastaa format the cdnas of all the genes in a database specified

=head1 OPTIONS

=cut

use warnings ;
use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::SeqIO;
use Getopt::Long qw(:config no_ignore_case);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils;

my $file = 'ensembl_cdnas';
≈
my $dbhost;
my $dbuser = 'admin' ;
my $dbname;
my $dbpass;
my $dbport;

my $dnadbhost;
my $dnadbuser = 'admin' ;
my $dnadbname;
my $dnadbpass;
my $dnadbport;

my $genetype;


GetOptions(
	    'dbname|db|D:s'    => \$dbname,
	    'host|dbhost|h:s'    => \$dbhost,
	    'pass|dbpass|p:s'    => \$dbpass,
	    'port|dbport|P:s'    => \$dbport,
	    'dnadbname:s' => \$dnadbname,
	    'dnadbhost:s' => \$dnadbhost,
	    'dnadbpass:s'    => \$dnadbpass,
	    'dnadbport:s'    => \$dnadbport,
	    'file:s'  => \$file,
	    'genetype:s'   => \$genetype,
	   );

unless ( $dbname && $dbhost && $dnadbname && $dnadbhost ){
  print STDERR "script to dump all the cdnas from the transcripts in a database\n";
 
  print STDERR "Usage: $0 -dbname -dbhost -dbport -dbpass -dnadbname -dnadbhost -dnadbport -dnadbpass\n";
  print STDERR "Optional: -genetype -file (defaulted to ensembl_cdnas)\n";
  exit(0);
}

my $dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					       '-host'   => $dnadbhost,
					       '-user'   => $dnadbuser,
					       '-dbname' => $dnadbname,
					       '-pass'   => $dnadbpass,
					       '-port'   => $dnadbport,
					      );


my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $dbhost,
					    '-user'   => $dbuser,
					    '-dbname' => $dbname,
					    '-pass'   => $dbpass,
					    '-port'   => $dbport,
					   );


print STDERR "connected to $dbname : $dbhost\n";

open (OUT,">$file") or die("unable to open file $file");

my $seqio = Bio::SeqIO->new('-format' => 'Fasta' , -fh => \*OUT ) ;

my  @genes = @{$db->get_GeneAdaptor->fetch_all} ;


GENE:
foreach my $gene(@genes) {
  
    if ($genetype){
	next GENE unless ( $gene->biotype eq $genetype );
    }


  my $gene_id = $gene->dbID();
  my $chr = $gene->chr_name;

 TRANS:
  foreach my $trans ( @{$gene->get_all_Transcripts} ) {
    my $gene_id = $gene->stable_id || $gene->dbID;
    my $tran_id = $trans->stable_id || $trans->dbID;
    my @evidence = &get_evidence($trans);
    
    my $strand = $trans->start_Exon->strand;
    my ($start,$end);
    my @exons;
    if ( $strand == 1 ){
      @exons = sort {$a->start <=> $b->end} @{$trans->get_all_Exons};
      $start = $exons[0]->start;
      $end   = $exons[$#exons]->end;
    }
    else{
      @exons = sort {$b->start <=> $a->end} @{$trans->get_all_Exons};
      $start = $exons[0]->end;
      $end   = $exons[$#exons]->start;
    }
    
    eval {      
      my $tran_seq = $trans->seq;
      
      my $tseq = $trans->translate();
      if ( $tseq->seq =~ /\*/ ) {
	print STDERR "translation of ".$trans->dbID." has stop codons. Skipping!\n";
	Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Evidence($trans);
	next TRANS;
      }
      $tran_seq->display_id("Gene:$gene_id Transcript:$tran_id");
      $tran_seq->desc("HMM:@evidence Chr:$chr Strand:$strand Start:$start End:$end");
      my $result = $seqio->write_seq($tran_seq);
    };
    if( $@ ) {
      print STDERR "unable to process transcript $tran_id, due to \n$@\n";
    }
  }
}

close (OUT);

sub get_evidence{
  my ($trans) = @_;
  my %evi;
  foreach my $exon (@{$trans->get_all_Exons}){
    foreach my $evidence ( @{$exon->get_all_supporting_features} ){
      $evi{$evidence->hseqname}++;
    }
  }
  return keys %evi;
}
