use lib 't';
use strict;
use Test;

BEGIN { $| = 1; plan test => 8;}

use Bio::EnsEMBL::Pipeline::Runnable::MiniGenewise;
use Bio::SeqIO;


my $genfile = 't/data/AC099340.fa.masked';
my $pepfile = 't/data/testpep.fa';

open(IN,"<$genfile") || die "Can't open $genfile";

ok(my $genseqio = new Bio::SeqIO(-fh => \*IN));
ok(my $gseq     = $genseqio->next_seq);

close(IN);

open(IN,"<$pepfile") || die "Can't open $pepfile";

ok(my $pepseqio = new Bio::SeqIO(-fh => \*IN));

# Let's make some features :

#Sub genomic                     111967  112056          +       .       Q9H7A7  1       30      +       .
#Sub genomic                     113024  113125          +       .       Q9H7A7  31      64      +       .
#Sub genomic                     120362  120739          +       .       Q9H7A7  65      190     +       .
#Sub genomic                     122399  122560          +       .       Q9H7A7  192     245     +       .

my $f1 = make_feature(111967,112056,1,'Q9H7A7',1,30,1);
my $f2 = make_feature(113024,113125,1,'Q9H7A7',31,64,1);
my $f3 = make_feature(120362,120739,1,'Q9H7A7',65,190,1);
my $f4 = make_feature(122399,122560,1,'Q9H7A7',192,245,1);

#  Exon 57017 56959 phase 0
#  Exon 49479 49388 phase 2
#  Exon 41734 41691 phase 1
#  Exon 35879 35639 phase 0
#  Exon 34045 33932 phase 1
#  Exon 31884 31740 phase 1
#  Exon 30725 30670 phase 2
#  Exon 24980 24872 phase 1
#  Exon 22704 22570 phase 2
#  Exon 13850 13747 phase 2
#  Exon 13658 13411 phase 1
#  Exon 12723 12520 phase 0

my $f5  = make_feature(56959,57017,-1,'Q91VS1',1,20,1);
my $f6  = make_feature(49388,49479,-1,'Q91VS1',21,52,1);
my $f7  = make_feature(41691,41734,-1,'Q91VS1',52,67,1);
my $f8  = make_feature(35639,35879,-1,'Q91VS1',55,152,1);
my $f9  = make_feature(33932,34045,-1,'Q91VS1',134,184,1);
my $f10 = make_feature(31740,31884,-1,'Q91VS1',185,233,1);
my $f11 = make_feature(30670,30725,-1,'Q91VS1',223,249,1);
my $f12 = make_feature(24872,24980,-1,'Q91VS1',249,288,1);
my $f13 = make_feature(22570,22704,-1,'Q91VS1',288,333,1);
my $f14 = make_feature(13747,13850,-1,'Q91VS1',329,369,1);
my $f15 = make_feature(13411,13658,-1,'Q91VS1',368,450,1);
my $f16 = make_feature(12520,12723,-1,'Q91VS1',450,518,1);

my @features2 = ($f5,$f6,$f7,$f8,$f9,$f10,$f11,$f12,$f13,$f14,$f15,$f16);

my $minigenewise1 = Bio::EnsEMBL::Pipeline::Runnable::MiniGenewise->new('-genomic'    => $gseq,
									'-protein'    => $pepseqio->next_seq,
									'-features'   => [$f1,$f2,$f3,$f4],
								       );

my $minigenewise2 = Bio::EnsEMBL::Pipeline::Runnable::MiniGenewise->new('-genomic'    => $gseq,
									'-protein'    => $pepseqio->next_seq,
									'-features'   => \@features2,
								       );


ok($minigenewise1);
ok($minigenewise2);

$minigenewise1->run;

ok(1);

$minigenewise2->run;

ok(my @results1 = $minigenewise1->output());
ok(my @results2 = $minigenewise2->output());

foreach my $f (@results1,@results2) {
  print "Feature " . $f->gffstring . "\n";
  
  if ($f->sub_SeqFeature) {
    foreach my $sub ($f->sub_SeqFeature) {
      print "Sub " . $sub->gffstring . "\n";
      if ($sub->sub_SeqFeature) {
	foreach my $subsub ($sub->sub_SeqFeature) {
	  print "  Sub sub " . $subsub->gffstring . "\n";
	}
      }
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
