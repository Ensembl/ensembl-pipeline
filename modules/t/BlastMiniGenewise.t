use lib 't';
use Test;
use strict;

BEGIN { $| = 1; plan test => 8;}

use Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise;
use Bio::EnsEMBL::Pipeline::SeqFetcher::FileIndex;
use Bio::Index::Fasta;
use Bio::PrimarySeq;
use Bio::Seq;
use Bio::SeqIO;

ok(1);

my $genfile = 't/data/AC099340.fa.masked';
my $pepfile = 't/data/testpep.fa';

open(IN,"<$genfile") || die "Can't open $genfile";

ok(my $genseqio   = new Bio::SeqIO(-fh => \*IN));

ok(my $gseq       = $genseqio->next_seq);

ok(my $seqfetcher = Bio::EnsEMBL::Pipeline::SeqFetcher::FileIndex->new(-seqfile => $pepfile));

ok(my $ids        = $seqfetcher->list_all_ids);

ok(my $blastminigenewise = Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise->new(
			       '-genomic'    => $gseq,
			       '-ids'	      => $ids,
			       '-seqfetcher' => $seqfetcher
										   ));


ok($blastminigenewise->run());

ok(my @results = $blastminigenewise->output);

foreach my $res (@results) {
  print "Feature " . $res->gffstring . "\n";
  if ($res->sub_SeqFeature) {
    foreach my $sub ($res->sub_SeqFeature) {
      print "Sub " . $sub->gffstring . "\n";
      if ($sub->sub_SeqFeature) {
	foreach my $subsub ($sub->sub_SeqFeature) {
	  print "   subsub " . $subsub->gffstring . "\n";
	}
      }
    }
  }
}
