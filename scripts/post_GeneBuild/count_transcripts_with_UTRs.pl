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


use warnings ;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
# use Bio::EnsEMBL::Pipeline::SeqFetcher::OBDAIndexSeqFetcher;
use Getopt::Long qw(:config no_ignore_case);
use strict;

my $dbname;
my $dbhost;
my $dbuser = 'ensro';
my $dbport; 
my $genetype;

GetOptions( 'host|dbhost|h:s'       => \$dbhost,
	     'dbname|db|D:s'       => \$dbname,
	     'dbport:s'        => \$dbport, 
             'genetype:s'     => \$genetype,
	   );

unless ( $dbhost && $dbname ){
  print STDERR "Usage: $0 -dbhost -dbname > & log_file\n";
  print STDERR "Optional: -genetype\n";
  exit(0);
}

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $dbhost,
					    -user   => $dbuser,
					    -port   => $dbport,
                                            -dbname => $dbname,
					   );



print STDERR "connected to $dbname : $dbhost\n";
my $sa = $db->get_SliceAdaptor;

my  @ids = @{$db->get_GeneAdaptor->list_dbIDs()};

my $only_utr_three = 0;
my $only_utr_five  = 0;
my $both = 0;
my $no_utr = 0;
my $total;

GENE:
foreach my $gene_id( @ids) {
    
  # this gives the gene in chromosome coordinates:
  my $gene    = $db->get_GeneAdaptor->fetch_by_dbID( $gene_id, 1 );
  my $gene_id = $gene->stable_id || $gene->dbID;

  if ($genetype){
    next GENE unless( $genetype eq $gene->type );
  }

 TRANS:
  foreach my $trans ( @{$gene->get_all_Transcripts} ){
    
    $total++;

    my $five_seq  = $trans->five_prime_utr_Feature; #->seq;
    my $three_seq = $trans->three_prime_utr_Feature; # ->seq;

    if ( $five_seq && $three_seq ){
      $both++;
      next TRANS;
    }
    
    if( $five_seq && !$three_seq ){
      $only_utr_five++;
      next TRANS;
    }

    if( !$five_seq && $three_seq ){
      $only_utr_three++;
      next TRANS;
    }

    $no_utr++;
  }
}

print STDERR "total transcripts: ".$total."\n";
print STDERR "no UTRs          : ".$no_utr."\n";
print STDERR "with 5'UTR only  : ".$only_utr_five."\n";
print STDERR "with 3'UTR only  : ".$only_utr_three."\n";
print STDERR "with both UTRs   : ".$both."\n";

############################################################
