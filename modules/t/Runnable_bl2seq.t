
use strict;
use Bio::EnsEMBL::Pipeline::Tools::Bl2seq;
use Bio::EnsEMBL::Pipeline::Runnable::Bl2seq;
use Bio::SeqIO;
use lib 't';
use Test;

BEGIN { $| = 1; plan test => 6;}

ok(1);

ok(open(SEQ, "<t/data/relaxins.fa"));

my $seqio = Bio::SeqIO->new('-fh'     => \*SEQ,
			    '-format' => 'fasta');

my $seq1 = $seqio->next_seq;
my $seq2 = $seqio->next_seq;

ok($seq1->isa("Bio::Seq"));
ok($seq2->isa("Bio::Seq"));

print $seq1->id . "\n";
print $seq2->id . "\n";

my $runnable_bl2seq 
   = Bio::EnsEMBL::Pipeline::Runnable::Bl2seq->new(
           -seq1      => $seq1,	
	   -seq2      => $seq2,
	   -alntype   => 'blastn',
	   -min_score => 10,
	   -workdir   => '/tmp');

ok($runnable_bl2seq->isa("Bio::EnsEMBL::Pipeline::Runnable::Bl2seq"));

ok($runnable_bl2seq->run);

close SEQ;
