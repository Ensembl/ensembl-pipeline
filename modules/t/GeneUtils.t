use lib 't';
use strict;
use Test;

BEGIN { $| = 1; plan test => 10;}

use Bio::EnsEMBL::Pipeline::SeqFetcher::FileIndex;
use Bio::EnsEMBL::Pipeline::GeneUtils;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::RawContig;
use Bio::SeqIO;

ok(1);

my $pepfile = '/nfs/acari/michele/cvs/ensembl-trunk/ensembl-pipeline/modules/t/data/testpep.fa';

ok( -e $pepfile);

ok(my $dnafile = 't/data/AC099340.fa.masked');

ok(open(IN,"<$dnafile"));

ok(my $seqio = new Bio::SeqIO(-fh => \*IN, -format => 'fasta'));

ok(my $dnaseq = $seqio->next_seq);

ok(my $contig = new Bio::EnsEMBL::RawContig(-id  => $dnaseq->id,
					    -seq => $dnaseq->seq));

ok(my $analysis = new Bio::EnsEMBL::Analysis(-logic_name => 'gene'));

my $f1  = make_feature(111967,112056,1,'Q9H7A7',1,30,1,163);
my $f2  = make_feature(113024,113125,1,'Q9H7A7',31,64,1,179);
my $f3  = make_feature(120362,120739,1,'Q9H7A7',65,190,1,645);
my $f4  = make_feature(122399,122560,1,'Q9H7A7',192,245,1,290);
my $f5  = make_feature(56959,57017,-1,'Q91VS1',1,20,1,147);
my $f6  = make_feature(49388,49479,-1,'Q91VS1',21,52,1,404);
my $f7  = make_feature(41691,41734,-1,'Q91VS1',52,67,1,79);
my $f8  = make_feature(35639,35879,-1,'Q91VS1',55,152,1,158);
my $f9  = make_feature(33932,34045,-1,'Q91VS1',134,184,1,217);
my $f10 = make_feature(31740,31884,-1,'Q91VS1',185,233,1,99);
my $f11 = make_feature(30670,30725,-1,'Q91VS1',223,249,1,114);
my $f12 = make_feature(24872,24980,-1,'Q91VS1',249,288,1,140);
my $f13 = make_feature(22570,22704,-1,'Q91VS1',288,333,1,193);
my $f14 = make_feature(13747,13850,-1,'Q91VS1',329,369,1,194);
my $f15 = make_feature(13411,13658,-1,'Q91VS1',368,450,1,386);
my $f16 = make_feature(12520,12723,-1,'Q91VS1',450,518,1,291);

ok(my @features = ($f1,$f2,$f3,$f4,$f5,$f6,$f7,$f8,$f9,$f10,$f11,$f12,$f13,$f14,$f15,$f16));

my $exon1 = make_feature(56959,57017,-1);
my $sub1a = make_feature(56973,57017,-1,'Q91VS1',1,15,1,163);
my $sub1b = make_feature(56961,56972,-1,'Q91VS1',17,20,1,163);

$exon1->add_sub_SeqFeature($sub1a);
$exon1->add_sub_SeqFeature($sub1b);

my $exon2 = make_feature(49388,49479,-1);
my $sub2  = make_feature(49389,49478,-1,'Q91VS1',22,51,1,404);

$exon2->add_sub_SeqFeature($sub2);

my $exon3 = make_feature(41691,41734,-1);
my $sub3  = make_feature(41691,41732,-1,'Q91VS1',53,66,1,202);

$exon3->add_sub_SeqFeature($sub3);

my $exon4 = make_feature(35639,35879,-1);
my $sub4  = make_feature(35640,35879,-1,'Q91VS1',67,146,1,315);

$exon4->add_sub_SeqFeature($sub4);

my $exon5 = make_feature(33932,34045,-1);
my $sub5  = make_feature(33933,34043,-1,'Q91VS1',148,184,1,287);

my $tranf = new Bio::EnsEMBL::SeqFeature(-start   => 33933,
					 -end     => 57017,
					 -strand  => -1);

$tranf->add_sub_SeqFeature($exon5,'EXPAND');
$tranf->add_sub_SeqFeature($exon4,'EXPAND');
$tranf->add_sub_SeqFeature($exon3,'EXPAND');
$tranf->add_sub_SeqFeature($exon2,'EXPAND');

ok($tranf->add_sub_SeqFeature($exon1,'EXPAND'));

ok(my ($transcript) = Bio::EnsEMBL::Pipeline::GeneUtils::SeqFeature_to_Transcript($tranf,$contig,$analysis,'',1));

foreach my $exon (@{$transcript->get_all_Exons}) {
	print "Exon " . $exon->start . " " . $exon->end . " " . $exon->strand . "\n";
}

print $transcript->translation->start . "\n";
print $transcript->translation->end   . "\n";
print $transcript->translate->seq     . "\n";

my $valid = Bio::EnsEMBL::Pipeline::GeneUtils::check_strand($transcript);

print "Valid strand $valid\n";

my @tran  = Bio::EnsEMBL::Pipeline::GeneUtils::check_introns($transcript,7000);

ok(scalar(@tran) == 1);

foreach my $tran (@tran) {
  print "Tran " . $tran->translate->seq . "\n";
}

my $compl = Bio::EnsEMBL::Pipeline::GeneUtils::check_low_complexity($transcript,80);

print "Complexity " . $compl . "\n";

my $seqfetch = Bio::EnsEMBL::Pipeline::SeqFetcher::FileIndex->new(
	   -seqfile => $pepfile);

my @ids = @{$seqfetch->list_all_ids};

print "Ids @ids\n";

my $cov = Bio::EnsEMBL::Pipeline::GeneUtils::check_coverage(
     $transcript,50,[$seqfetch]);

print "Coverage $cov\n";


sub make_feature {
  my ($start,$end,$strand,$hid,$hstart,$hend,$hstrand) = @_;
  
  my $seqname = 'AC099340';
  
  my $f1 = new Bio::EnsEMBL::SeqFeature(-seqname  => $seqname,
					-start => $start,
					-end   => $end,
					-strand => $strand);
  
  if (defined($hid)) {
    my $f2 = new Bio::EnsEMBL::SeqFeature(-seqname => $hid,
					  -start   => $hstart,
					  -end     => $hend,
					  -hstrand => $hstrand);
    
    my $fp = new Bio::EnsEMBL::FeaturePair(-feature1 => $f1,
					   -feature2 => $f2);
    
    return $fp;
  } else {
    return $f1;
  }
}
