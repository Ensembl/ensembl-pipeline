#! /usr/local/ensembl/bin/perl -w

use strict;
use lib 't';
use Test;
use Bio::EnsEMBL::Pipeline::GeneDuplication::PAML;
use Bio::SeqIO;

BEGIN { $| = 1; plan test => 5;}

ok(1);

my $verbose = 1 if @ARGV;

my $seqio = Bio::SeqIO->new(-file   =>'t/data/paml_test.fa',
			    -format => 'fasta');

my @seqs;
while (my $seq = $seqio->next_seq){
  push (@seqs, $seq);
}

ok(@seqs);

my $paml 
  = Bio::EnsEMBL::Pipeline::GeneDuplication::PAML->new(
      '-executable'      => 'codeml',
      '-aligned_seqs'    => \@seqs,
      '-runmode'         => '-2',
      '-seqtype'         => '1',
      '-model'           => '0',
      '-nssites'         => '0',
      '-icode'           => 1);

ok($paml->isa("Bio::EnsEMBL::Pipeline::GeneDuplication::PAML"));

ok(my $parser = $paml->run_codeml);

exit 1 unless $verbose;

while (my $result = $parser->next_result()) {
  my @otus = $result->get_seqs();

  my $matrix;

  $matrix = $result->get_MLmatrix();

  printf "The omega ratio for sequences %s vs %s was: %g\n",
  $otus[1]->id, $otus[2]->id, $matrix->[1]->[2]->{omega};
}
