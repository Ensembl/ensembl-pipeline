use lib 't';
use Test;
use strict;

BEGIN { $| = 1; plan test => 7;}

use Bio::EnsEMBL::Pipeline::Runnable::Est2Genome;
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

foreach my $est (@estseqs) {
  ok(my $e2g = Bio::EnsEMBL::Pipeline::Runnable::Est2Genome->new(-genomic => $genseq,
								 -est     => $est ));

  ok($e2g->run);

  ok(my @results = $e2g->output);

  foreach my $f (@results) {
    print "Feature " . $f->gffstring . "\n";
    
    if ($f->sub_SeqFeature) {
      foreach my $sub ($f->sub_SeqFeature) {
	print "Sub " . $sub->gffstring . "\n";
	if ($sub->sub_SeqFeature) {
	  foreach my $subsub ($sub->sub_SeqFeature) {
	    print "  SubSub " . $sub->gffstring . "\n";
	  }
	}
      }
    }
  }
}
