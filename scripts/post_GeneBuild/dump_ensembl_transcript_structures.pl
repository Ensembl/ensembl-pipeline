#!/usr/bin/env perl


# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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


# script to dump info for the Oxford guys, for the Rat genome project

use warnings ;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::SeqFetcher::OBDAIndexSeqFetcher;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Getopt::Long qw(:config no_ignore_case);
use strict;

my $dbname;
my $dbhost;
my $dbuser = 'ensro';
my $dnadbname;
my $dnadbhost;

my $genetype;

my $file = 'ensembl_transcript_structures';

GetOptions( 'host|dbhost|h:s'       => \$dbhost,
	     'dbname|db|D:s'       => \$dbname,
	     'dnadbname:s'     => \$dnadbname,
	     'dnadbhost:s'     => \$dnadbhost,
	     'genetype:s'     => \$genetype,
	     'file:s'         => \$file,
	   );

unless ( $dbhost && $dbname ){
  print STDERR "Usage: $0 -dbhost -dbname -dnadbhost -dnadbname > & log_file\n";
  print STDERR "Optional: -file -genetype\n";
  exit(0);
}

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $dbhost,
					    -user   => $dbuser,
					    -dbname => $dbname,
					   );

my $dnadb;
if ($dnadbhost && $dnadbname ){
  $dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $dnadbhost,
					      -user   => $dbuser,
					      -dbname => $dnadbname,
					     );
  $db->dnadb($dnadb);
}

  open (OUT,     ">$file" )     or die ("cannot open $file");

print STDERR "connected to $dbname : $dbhost\n";
my $sa = $db->get_SliceAdaptor;

my  @ids = @{$db->get_GeneAdaptor->list_geneIds};

print OUT "gene_id\ttranscript_id\tchr_name\tstrand\texon_coordinates\n";

GENE:
foreach my $gene_id( @ids) {
  #print STDERR "retrieveing gene: ".$gene_id."\n";
  
  # this gives the gene in chromosome coordinates:
  my $gene    = $db->get_GeneAdaptor->fetch_by_dbID( $gene_id, 1 );
  my $geneid = $gene->stable_id || $gene->dbID;

  if ($genetype){
    next GENE unless( $genetype eq $gene->type );
  }

 TRANS:
  foreach my $trans ( @{$gene->get_all_Transcripts} ) {
    
    my $trans_id   = $trans->stable_id || $trans->dbID;
    my @exons     = @{$trans->get_all_Exons};
    my $slice     = $exons[0]->contig;
    my $strand    = $exons[0]->strand;
    my $chr_name  = $slice->chr_name;
    

    if ( $dnadb ){
      eval{
	my $tseq = $trans->translate();
	if ( $tseq->seq =~ /\*/ ) {
	  print STDERR "translation of ".$trans->dbID." has stop codons. Skipping!\n";
	  Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Evidence($trans);
	  next TRANS;
	}
      };
    }
    print OUT $geneid."\t".$trans_id."\t".$chr_name."\t".$strand."\t";

    foreach my $exon ( @exons ){
      print OUT $exon->start."-".$exon->end."\t";
    }
    print OUT "\n";
  }
}

close (OUT);

############################################################
