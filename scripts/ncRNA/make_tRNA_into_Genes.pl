#!/usr/bin/env perl


# Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(stack_trace);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils ;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonUtils ;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use ncRNA_update_config ;

my $dbhost;
my $dbuser;
my $pass;
my $dbname;
my $dbport;
my $help;
my @species_list;
my $write;
my @dbID;


$| = 1; 

&GetOptions(
	    'pass=s'    => \$pass,
	    'species=s' => \@species_list,
	    'write!'   => \$write,
           );

if(!$pass || $help){
  die ("perl make_tRNA_into_Genes.pl
-pass *(password) 
-species (species list)
-write ( actually write them into the database ) 
Gets tRNA simple features out of the core database specified in ncRNA_update_config and 
writes Genes with tRNA biotype into the write db\n");
  $help = 1;
}

my @speciess;
if (scalar(@species_list)) {
  @species_list = split(/,/,join(',',@species_list));
  foreach my $species (@species_list){
    if ($CONFIG->{$species}){
      push @speciess, $species;
    } else {
    print "Skipping species $species\n";
    }
  }
} else {
  @speciess = (keys %$CONFIG);
}

# might want to check some config first
SPECIES:  foreach my $species (@speciess){
  print "checking species $species\n";
  my @tRNA_Genes;
  # check that all the dbs exist and have finished their run!
  my $sdb = new Bio::EnsEMBL::DBSQL::DBAdaptor
    (
     -host   => $CONFIG->{$species}->{"WRITEHOST"},
     -user   => $WRITEUSER,
     -port   => $CONFIG->{$species}->{"WRITEPORT"},
     -pass   => $pass,
     -dbname => $CONFIG->{$species}->{"WRITENAME"},
    );
  my $core_db = new Bio::EnsEMBL::DBSQL::DBAdaptor
    (
     '-host'   => $CONFIG->{$species}->{"DBHOST"},
     '-user'   => 'ensro',
     '-dbname' => $CONFIG->{$species}->{"DBNAME"},
     '-port'   => $CONFIG->{$species}->{"DBPORT"},
    );
  
  my $ga  = $sdb->get_GeneAdaptor;
  my $sfa = $core_db->get_SimpleFeatureAdaptor;
  foreach my $tRNA ( @{$sfa->fetch_all_by_logic_name('trnascan')} ) {
    my $tRNA_Gene = make_gene($tRNA);
    push @tRNA_Genes , $tRNA_Gene if ( $tRNA_Gene) ;
  }
  print "Storing " . scalar ( @tRNA_Genes ) , " tRNA genes in " . $CONFIG->{$species}->{"WRITENAME"} . 
    "@" . $CONFIG->{$species}->{"WRITEHOST"} ."\n";
  if ( $write ) {
    foreach my $tRNA_Gene ( @tRNA_Genes ) {
      #     print $tRNA_Gene->seq_region_name . " " .
      #$tRNA_Gene->start . " " .
      #  $tRNA_Gene->end . " " .
      #    $tRNA_Gene->strand . " " .
      #      $tRNA_Gene->biotype . " " .
      #	$tRNA_Gene->description . "\n" ;
      eval {
	$ga->store($tRNA_Gene);
      }; if ( $@ ) {
	throw("Problem writing gene \n$@\n");
      }
    }
  } else {
    print " Write protect on, not writing\n";
  }
}

exit;

sub make_gene {
  my ($tRNA) = @_;
  my $exon = create_Exon(
			 $tRNA->start,
			 $tRNA->end,
			 -1,
			 -1,
			 $tRNA->strand,
			 $tRNA->analysis,
			 undef,
			 undef,
			 $tRNA->slice,
			);
  my $tran =  new Bio::EnsEMBL::Transcript(-EXONS => [( $exon )]);
  $tran->analysis($tRNA->analysis);
  $tran->version(1);
  if ( $tRNA->display_label eq  'Pseudo' ) {
    $tran->biotype('tRNA_pseudogene');
  } else {
    $tran->biotype('tRNA');
  }
  my ( $tRNA_gene ) = @{convert_to_genes(($tran),$tRNA->analysis)};
  if ( $tRNA->display_label eq  'Pseudo' ) {
    $tRNA_gene->biotype('tRNA_pseudogene');
  } else {
    $tRNA_gene->biotype('tRNA');
    $tRNA_gene->description("tRNA, anticodon " . $tRNA->display_label) ;
  }
  return $tRNA_gene;
}
