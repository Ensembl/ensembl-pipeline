#! /usr/local/ensembl/bin/perl -w

use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Runnable::ESTDiscriminator; 
use Bio::EnsEMBL::Pipeline::SeqFetcher::OBDAIndexSeqFetcher;

my @gene_stable_ids = @ARGV;

unless (@gene_stable_ids) {
  print STDERR "Invalid input ids.\n" . 
    "Usage : est_discriminator.pl <list of gene stable ids>\n";
  die;
}

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(-dbname => 'homo_sapiens_core_26_35',
                                             -user   => 'ensro',
                                             -host   => 'ecs2',
					     -port   => 3365);

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

my $seq_fetcher = 
  Bio::EnsEMBL::Pipeline::SeqFetcher::OBDAIndexSeqFetcher->new(
    -db     => ['/ecs4/work4/dta/est_frameshifts/homo_sapiens_23_35_ests/index/'],
    -format => 'fasta');

my $est_descrim = 
  Bio::EnsEMBL::Pipeline::Runnable::ESTDescriminator->new(
    -genes               => \@genes,
    -est_db              => $estdb,
    -est_coverage_cutoff => 0.8,
    -sequence_fetcher    => $seq_fetcher);


my $est_vs_gene = $est_descrim->run;


# Print output to screen.
#print STDERR "ESTs matched to their most appropriate gene:\n";
#foreach my $est_id (keys %$est_vs_gene){
#  print STDERR "EST : $est_id\tGene(s) : ";
#  foreach my $gene_id (@{$est_vs_gene->{$est_id}}){
#    print STDERR $gene_id . ' ';
#  }
#  print STDERR "\n";
#}

# Display singly matching EST data

print STDOUT "Single Gene vs EST matches : \n";
foreach my $gene_stable_id (@gene_stable_ids){
  my $single_match_est_ids = $est_descrim->single_match_ESTs($gene_stable_id);
  print STDOUT $gene_stable_id . " : " . scalar @$single_match_est_ids . "\n"
#  print STDOUT $gene_stable_id . " : @$single_match_est_ids\n"
}

# Display multiply matching EST data

print STDOUT "Gene vs multiple EST matches : \n";
foreach my $gene_stable_id (@gene_stable_ids){
  my $multiple_match_est_ids = $est_descrim->multiple_match_ESTs($gene_stable_id);
  print STDOUT $gene_stable_id . " : @$multiple_match_est_ids\n"
}

# Display incorrectly matching EST data

print STDOUT "Incorrect EST matches : \n";
foreach my $gene_stable_id (@gene_stable_ids){
  my $incorrect_match_est_ids = $est_descrim->incorrect_match_ESTs($gene_stable_id);
  print STDOUT $gene_stable_id . " : " . scalar @$incorrect_match_est_ids . "\n"
}

