use lib 't';
use strict;
use Test;

BEGIN { $| = 1; plan test => 18;}

use Bio::EnsEMBL::Pipeline::Runnable::MiniGenomewise;
use Bio::PrimarySeq;
use Bio::Seq;
use Bio::SeqIO;

ok(1);

ok(my $seqin = Bio::SeqIO->new( -file => "t/data/human.genomic"));
ok(my $genomic_seq = $seqin->next_seq());

ok(my $transcript = Bio::EnsEMBL::Transcript->new());

ok(my $exon1 = Bio::EnsEMBL::Exon->new);
ok(my	$exon2 = Bio::EnsEMBL::Exon->new);

ok($exon1->start(1787) == 1787);
ok($exon1->end  (1935) == 1935);
ok($exon1->strand (1)  == 1);

ok($transcript->add_Exon($exon1));

ok($exon2->start(2084) == 2084);
ok($exon2->end  (2180)  == 2180);
ok($exon2->strand (1)   == 1);

ok($transcript->add_Exon($exon2));

my $minigenomewise = Bio::EnsEMBL::Pipeline::Runnable::MiniGenomewise->new(
																																					 -genomic     => $genomic_seq,
																																					 -transcripts => [$transcript],
																																					);

ok($minigenomewise);

$minigenomewise->run;

ok(1);

ok(my @results = $minigenomewise->output());

foreach my $trans (@results ) {
	ok($trans->isa('Bio::EnsEMBL::Transcript'));
}

