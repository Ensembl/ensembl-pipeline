#! /usr/local/bin/perl -w

use strict;
use lib 't';
use Test;
use Bio::EnsEMBL::Pipeline::GeneDuplication;
use Bio::EnsEMBL::Pipeline::Runnable::BlastDB;
use Bio::SeqIO;

my $verbose = 1 if @ARGV;

BEGIN { $| = 1; plan test => 7;}

ok(1);

# Must give GeneDuplication a blast database that can be seq fetched from - wu_new or ncbi flavs

ok(my $blastdb 
  = Bio::EnsEMBL::Pipeline::Runnable::BlastDB->new(
      -dbfile               => 't/data/relaxins.fa',
      -molecule_type        => 'DNA',
      -workdir              => '/tmp',
      -copy                 => 1,
      -index_type           => 'wu_new',
      -make_fetchable_index => 1,
      ));

ok($blastdb->run);

my $gene_dupl 
  = Bio::EnsEMBL::Pipeline::GeneDuplication->new(
	 '-blastdb'                => $blastdb,
	 '-hit_coverage'           => 50,
	 '-hit_identity'           => 80,
	 '-regex_query_species'    => 'Hsa',
	 '-regex_outgroup_species' => ['Lca','Ssc'],
	 '-genetic_code'           => 1,
	 '-work_dir'               => '/tmp',
	);

ok($gene_dupl->isa("Bio::EnsEMBL::Pipeline::GeneDuplication"));

my $quick_seqio = Bio::SeqIO->new('-file'   => 't/data/relaxins.fa',
				  '-format' => 'fasta');

my $input_seq = $quick_seqio->next_seq;

ok($input_seq->isa("Bio::Seq"));

$input_seq->display_id('test_sequence');

ok(my $result = $gene_dupl->run($input_seq));

ok(@{$result->matches});

exit 1 unless $verbose;

print "Query sequence       : " . $result->query_id . "\n";
print "Distance method used : " . $result->distance_method . "\n\n";

print join "\t", "match_id", "dN", "dS", "\n";
foreach my $match (@{$result->matches}){
  print join "\t", $match->{id}, $match->{dN}, $match->{dS}, "\n";
}
