use lib 't';
use Test;
use strict;

BEGIN {$| = 1; plan test => 6}

use Bio::EnsEMBL::Pipeline::SeqFetcher::FileIndex;

ok(1);

ok(my $index = new Bio::EnsEMBL::Pipeline::SeqFetcher::FileIndex(
   -seqfile => 't/data/testest.fa'));

ok(my @ids = @{$index->list_all_ids});

ok(scalar(@ids) == 2);

foreach my $id (@ids)  {
  ok(my $seq = $index->get_Seq_by_acc($id));
}

