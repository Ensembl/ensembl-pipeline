use lib 't';
use Test;
use strict;

BEGIN { $| = 1; plan test => 6;}

use Bio::EnsEMBL::Pipeline::Runnable::Exonerate;
use Bio::SeqIO;

ok(1);

ok(my $estseqio = new Bio::SeqIO(-file   => 't/data/testest.fa', 
                                 -format => 'fasta'));

my @estseqs;

while (my $seq = $estseqio->next_seq) {
  push(@estseqs,$seq);
}

ok(my $exonerate = Bio::EnsEMBL::Pipeline::Runnable::Exonerate->new(
   -query_seqs  => \@estseqs,
   -query_type  => 'DNA',
   -database    => 't/data/AC099340.fa.masked',
   -target_type => 'DNA',
   -exonerate   => '/usr/local/ensembl/bin/exonerate-0.8.2',
   -verbose     => 1,
   -options     => '--exhaustive FALSE --model est2genome --softmasktarget --score 500 --fsmmemory 800  --saturatethreshold 100 --dnahspthreshold 60 --dnawordlen 14 --showalignment FALSE --showvulgar FALSE --percent 90'
   ));

ok($exonerate->run);

ok(my @transcripts = $exonerate->output());

ok(display(@transcripts));

sub display {
  my @transcripts = @_;

  foreach my $transcript (@transcripts) {
    my $exons = $transcript->get_all_Exons;
    my $dafs = $exons->[0]->get_all_supporting_features;

    print "Match : (hseqname) " . 
    $dafs->[0]->hseqname . "\t" . 
    "(length) " . $transcript->length . "\t" . 
    "(exons) " . scalar @$exons . "\t" . 
    "\n";
  }

  return 1
}

