use lib 't';
use Test;
use strict;

BEGIN { $| = 1; plan test => 7;}

use Bio::EnsEMBL::Pipeline::Runnable::MiniEst2Genome;
use Bio::EnsEMBL::Pipeline::SeqFetcher::FileIndex;
use Bio::Index::Fasta;
use Bio::PrimarySeq;
use Bio::SeqIO;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::SeqFeature;

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


my $f1 = make_feature(57293,57493,1,'AK024759',1,201,1,100);
my $f2 = make_feature(74826,74934,1,'AK024759',202,310,1,100);
my $f3 = make_feature(76226,76310,1,'AK024759',311,395,1,100);
my $f4 = make_feature(81608,81708,1,'AK024759',396,496,1,100);
my $f5 = make_feature(100807,100885,1,'AK024759',497,595,1,100);
my $f6 = make_feature(100994,101078,1,'AK024759',576,660,1,100);

my $f7 = make_feature(12309,12723,-1,'BC029503',17,431,-1,100);
my $f8 = make_feature(13411,13658,-1,'BC029503',432,679,-1,100);
my $f9 = make_feature(13747,13850,-1,'BC029503',680,783,-1,100);
my $f10= make_feature(22570,22704,-1,'BC029503',784,918,-1,100);
my $f11= make_feature(24872,24980,-1,'BC029503',919,1027,-1,100);
my $f12= make_feature(30670,30725,-1,'BC029503',1028,1083,-1,100);


ok(my $seqfetcher = Bio::EnsEMBL::Pipeline::SeqFetcher::FileIndex->new(
  -seqfile => 't/data/testest.fa'));


my @feat_pairs = ($f1,$f2,$f3,$f4,$f5,$f6,$f7,$f8,$f9,$f10,$f11,$f12);

ok(my $miniest2genome = Bio::EnsEMBL::Pipeline::Runnable::MiniEst2Genome->new(
  '-genomic'    => $genseq,
  '-features'   => \@feat_pairs,
  '-seqfetcher' => $seqfetcher));

 
ok($miniest2genome->run);

ok(my @results = $miniest2genome->output);

ok(display(@results));

sub display {
  my @results = @_;
  my @methods = qw( seqname start end strand );

  foreach my $obj (@results) {
    print "Feature " . $obj->gffstring . "\n";
    if ($obj->sub_SeqFeature) {
      foreach my $sub ($obj->sub_SeqFeature) {
	print " Sub " . $sub->gffstring . "\n";
	if ($sub->sub_SeqFeature) {
	  foreach my $subsub($sub->sub_SeqFeature) {
	    print "  Subsub " . $subsub->gffstring . "\n";
	  }
	}
      }
    }
  }
  return 1;
}


sub make_feature {
  my ($start,$end,$strand,$hid,$hstart,$hend,$hstrand) = @_;
  
  my $f1 = new Bio::EnsEMBL::SeqFeature(-start => $start,
					-end   => $end,
					-strand => $strand);
  $f1->seqname('AC099340');
  my $f2 = new Bio::EnsEMBL::SeqFeature(-seqname => $hid,
					-start   => $hstart,
					-end     => $hend,
					-hstrand => $hstrand);
  
  my $fp = new Bio::EnsEMBL::FeaturePair(-feature1 => $f1,
					 -feature2 => $f2);
  
  return $fp;
}
