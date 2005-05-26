#! /usr/local/ensembl/bin/perl -w

use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Runnable::ESTDescriminator; 

my @gene_stable_ids = @ARGV;

unless (@gene_stable_ids) {
  print STDERR "Invalid input ids.\n" . 
    "Usage : check_shared_ests.pl <list of gene stable ids>\n";
  die;
}

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(-dbname => 'homo_sapiens_core_26_35',
                                             -user   => 'ensro',
                                             -host   => 'ecs2',
					     -port   => 3364);

my $ga = $db->get_GeneAdaptor;

my @genes;

foreach my $gene_stable_id (@gene_stable_ids){
  my $gene = $ga->fetch_by_stable_id($gene_stable_id);
  die "Gene not found" unless $gene;
  push @genes, $gene;
  print STDERR "Found gene : $gene_stable_id\n";
}

my $estdb = Bio::EnsEMBL::DBSQL::DBAdaptor->new(-dbname => 'vivek_homo_sapiens_23_35_est',
						-user   => 'ensro',
						-host   => 'ecs4',
						-port   => 3353);

my $est_descrim = 
  Bio::EnsEMBL::Pipeline::Runnable::ESTDescriminator->new(
    -genes               => \@genes,
    -est_db              => $estdb,
    -est_coverage_cutoff => 0.8);

$est_descrim->print_shared_ESTs;
