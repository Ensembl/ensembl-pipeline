use lib 't';
use Test;
use strict;

BEGIN { $| = 1; plan test => 9;
	require "Bio/EnsEMBL/Pipeline/pipeConf.pl";
      }

use Bio::EnsEMBL::Pipeline::Runnable::Exonerate;
use Bio::PrimarySeq;
use Bio::Seq;
use Bio::SeqIO;

ok(1);

open(IN,"<t/data/AC099340.fa.masked");

ok(my $genseqio = new Bio::SeqIO(-fh     => \*IN, -format => 'fasta'));

ok(my $genseq = $genseqio->next_seq);

close(IN);

my @estseqs;

open(IN,"<t/data/testest.fa") || die "Can't open testfile";

ok(my $estseqio = new Bio::SeqIO(-fh => \*IN, -format => 'fasta'));

while (my $seq = $estseqio->next_seq) {
  push(@estseqs,$seq);
}

close(IN);

my $exargs = " -w 14 -t 65 -H 100 -D 15 -m 500 ";

ok(my $exe = $::pipeConf{'bindir'} . '/exonerate');

ok(my $exonerate = Bio::EnsEMBL::Pipeline::Runnable::Exonerate->new (
   -EST       => \@estseqs,
   -GENOMIC   => $genseq,
   -EXONERATE => $exe,
   -ARGS      => $exargs,
   -PRINT     => 1));

my $ungapped = 1;

ok($exonerate->run($ungapped));

ok(my @results = $exonerate->output());

ok(display(@results));

sub display {
  my @results = @_;

  foreach my $pair (@results) {
    print $pair->seqname . "\t" . $pair->start . "\t" . $pair->end . "\t" . $pair->score . "\t" . 
      $pair->strand . "\t" . $pair->hseqname . "\t" . $pair->hstart . "\t" . $pair->hend . "\n";
  }
  return 1;
}

