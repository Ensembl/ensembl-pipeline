use lib 't';

use Test;
use strict;

BEGIN { $| = 1; plan test => 13}; 

use Bio::EnsEMBL::Pipeline::Runnable::BlastMiniEst2Genome;
use Bio::EnsEMBL::Pipeline::SeqFetcher::FileIndex;
use Bio::Index::Fasta;
use Bio::PrimarySeq;
use Bio::SeqIO;

ok(1);

my $genfile = 't/data/AC099340.fa.masked';
my $estfile = 't/data/testest.fa';

open(IN,"<$genfile") || die "Can't open $genfile";

ok(my $genseqio = new Bio::SeqIO(-fh => \*IN));
ok(my $gseq     = $genseqio->next_seq);

my $seqfetcher = Bio::EnsEMBL::Pipeline::SeqFetcher::FileIndex->new (
	-seqfile => $estfile);

ok($seqfetcher);

foreach my $id (@{$seqfetcher->list_all_ids}) {

	ok(my $queryseq = $seqfetcher->get_Seq_by_acc($id));

	my $blastminiest2genome = Bio::EnsEMBL::Pipeline::Runnable::BlastMiniEst2Genome->new(
											 '-genomic'    => $gseq,
                       '-seqfetcher' => $seqfetcher,
                       '-queryseq'   => [$queryseq],
                       '-threshold'  => 1e-2
								     );
	
	ok($blastminiest2genome);
	ok($blastminiest2genome->run);

	ok(my @results = $blastminiest2genome->output());


	foreach my $f (@results) {
		print "Feature " . $f->gffstring . "\n";

		if ($f->sub_SeqFeature) {
			foreach my $sub ($f->sub_SeqFeature) {
				print "Sub " . $sub->gffstring . "\n";
			}
		}
	}
}

ok(1);
