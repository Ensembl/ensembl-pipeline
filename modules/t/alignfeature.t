use lib 't';
use strict;
use Test;

BEGIN { $| = 1; plan tests => 8;}

use Bio::EnsEMBL::Pipeline::Runnable::AlignFeature;
use Bio::EnsEMBL::Pipeline::SeqFetcher::FileIndex;
use Bio::PrimarySeq;
use Bio::SeqIO;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::SeqFeature;

ok(1);

my $genfile = 't/data/AC099340.fa.masked';
my $pepfile = 't/data/testpep.fa';
my $estfile = 't/data/testest.fa';

open(IN,"<$genfile") || die "Can't open $genfile";

ok(my $genseqio = new Bio::SeqIO(-fh => \*IN));
ok(my $gseq     = $genseqio->next_seq);

close(IN);

open(IN,"<$estfile") || die "Can't open $estfile";

ok(my $estseqio = new Bio::SeqIO(-fh => \*IN));

my @estseqs;

while (my $seq = $estseqio->next_seq) {
  push(@estseqs,$seq);
}

ok(scalar(@estseqs) == 2);



my $seqfetcher = Bio::EnsEMBL::Pipeline::SeqFetcher::FileIndex->new(-seqfile => $estfile);


my $alignfeature = Bio::EnsEMBL::Pipeline::Runnable::AlignFeature->new(
       '-genomic'    => $gseq,
       '-features'   => \@estfp,
       '-seqfetcher' => $seqfetcher);
 
ok($alignfeature);

ok($alignfeature->run);

ok(my @results = $alignfeature->output);

display(@results);


sub display {
  my @results = @_;
  my @methods = qw( seqname start end strand );

  foreach my $obj (@results) {

      printf STDERR "\n";
      foreach my $method_name (@methods) {
        my $value = $obj->$method_name();
        printf STDERR ("%10s = $value\n", $method_name);
      }
    }
}
sub make_feature {
  my ($start,$end,$strand,$hid,$hstart,$hend,$hstrand) = @_;
  
  my $f1 = new Bio::EnsEMBL::SeqFeature(-start => $start,
					-end   => $end,
					-strand => $strand);
  my $f2 = new Bio::EnsEMBL::SeqFeature(-seqname => $hid,
					-start   => $hstart,
					-end     => $hend,
					-hstrand => $hstrand);
  
  my $fp = new Bio::EnsEMBL::FeaturePair(-feature1 => $f1,
					 -feature2 => $f2);
  
  return $fp;
}
