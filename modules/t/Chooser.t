#! /usr/local/bin/perl -w

use strict;
use lib 't';
use Test;
use Bio::EnsEMBL::Pipeline::GeneDuplication::Chooser;
use Bio::EnsEMBL::Pipeline::Runnable::BlastDB;
use Bio::SeqIO;

my $verbose = 1 if @ARGV;

BEGIN { $| = 1; plan test => 7;}

ok(1);

my $seqio = Bio::SeqIO->new('-file'   => 't/data/relaxins.fa',
			    '-format' => 'fasta');

my $input_seq = $seqio->next_seq;
$input_seq->display_id('Hsa_testRLN');

ok($input_seq->isa("Bio::Seq"));

my $blastdb 
  = Bio::EnsEMBL::Pipeline::Runnable::BlastDB->new(
      -dbfile     => 't/data/relaxins.fa',
      -index_type => 'wu_new',
      -type       => 'DNA');

ok($blastdb->isa("Bio::EnsEMBL::Pipeline::Runnable::BlastDB"));

ok($blastdb->run);

my $hit_sort_tool 
  = Bio::EnsEMBL::Pipeline::GeneDuplication::Chooser->new(
	-query_seq              => $input_seq,
	-blastdb                => $blastdb,
	-identity_cutoff        => 80,
	-coverage_cutoff        => 50,
	-regex_query_species    => 'Hsa',
	-regex_outgroup_species => ['Lca','Ssc'],
	-work_dir               => '/tmp',
	-genetic_code           => 1);

ok($hit_sort_tool->isa("Bio::EnsEMBL::Pipeline::GeneDuplication::Chooser"));

ok(my $ids = $hit_sort_tool->find_recent_duplications);

ok(@$ids);

exit 1 unless $verbose;

print "Recent paralogues of Gene " . $input_seq->display_id . "\n" if scalar @$ids;

foreach my $id (@$ids) {
  print $id . "\n" unless $id eq $input_seq->display_id;
}

