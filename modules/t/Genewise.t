use lib 't';
use Test;
use strict;

BEGIN { $| = 1; plan test => 12;}

use Bio::EnsEMBL::Pipeline::Runnable::Genewise;
use Bio::SeqIO;

ok(1);
print STDERR "genewise test is being created\n";
my $genfile = 't/data/AC099340.fa.masked';
my $pepfile = 't/data/testpep.fa';

open(IN,"<$genfile") || die "Can't open $genfile";

ok(my $genseqio = new Bio::SeqIO(-fh => \*IN));
ok(my $gseq     = $genseqio->next_seq);
close(IN);

open(IN,"<$pepfile") || die "Can't open $pepfile";

ok(my $pepseqio = new Bio::SeqIO(-fh => \*IN));

while (my $pseq = $pepseqio->next_seq) {
  ok($pseq);
  print STDERR "have peptide sequence\n";
  my $reverse = 0;
  
  if ($pseq->id eq 'Q91VS1') {
    $reverse = 1;
    
  }
  
  ok(my $genewise = Bio::EnsEMBL::Pipeline::Runnable::Genewise->new(-genomic  => $gseq,
								    -protein  => $pseq,
								    -reverse  => $reverse,
								   ));
  
  
  
  ok($genewise->run());
  ok(my @results = $genewise->output());
  
  my @results;

  foreach my $res (@results) {
    print "Feature " . $res->gffstring . "\n";
    if ($res->sub_SeqFeature) {
      foreach my $sub ($res->sub_SeqFeature) {
	print "Sub " . $sub->gffstring . "\n";
      }
    }
  }
}
