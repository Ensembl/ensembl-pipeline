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
use Bio::EnsEMBL::Pipeline::SeqFetcher::OBDAIndexSeqFetcher;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Getopt::Long qw(:config no_ignore_case);
use strict;

my $dbname;
my $dbhost;
my $dbuser = 'ensro';

my $genetype;

GetOptions( 'host|dbhost|h:s'       => \$dbhost,
	     'dbname|db|D:s'       => \$dbname,
	     'genetype:s'     => \$genetype,
	   );

unless ( $dbhost && $dbname ){
  print STDERR "Usage: $0 -dbhost -dbname > & log_file\n";
  print STDERR "Optional: -genetype\n";
  exit(0);
}

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $dbhost,
					    -user   => $dbuser,
					    -dbname => $dbname,
					   );



print STDERR "connected to $dbname : $dbhost\n";
my $sa = $db->get_SliceAdaptor;

my  @ids = @{$db->get_GeneAdaptor->list_geneIds};

my $single_exon = 0;
my $only_frameshifts = 0;

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
    
    my $trans_id   = $trans->stable_id || $trans->dbID;
    my $exon_count = scalar( @{$trans->get_all_Exons} );
    if ($exon_count == 1 ){
      $single_exon++;
    }
    if ( !Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->is_spliced($trans) 
	 &&
	 $exon_count>1 ){
      $only_frameshifts++;
    }
  }
}

print STDERR "single exon transcripts          : ".$single_exon."\n";
print STDERR "transcritps with only frameshifts: ".$only_frameshifts."\n";
print STDERR "total                            : ".( $single_exon + $only_frameshifts )."\n";

############################################################
