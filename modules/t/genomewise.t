use lib 't';
use strict;
use Test;

## We start with some black magic to print on failure.
BEGIN { $| = 1; plan test => 38;} 

use Bio::EnsEMBL::Pipeline::Runnable::Genomewise;
use Bio::PrimarySeq;
use Bio::Seq;
use Bio::SeqIO;

ok(1);

ok(my $seqin = Bio::SeqIO->new( -file => "t/data/human.genomic"));
ok(my $gen = $seqin->next_seq());

ok(my $run = Bio::EnsEMBL::Pipeline::Runnable::Genomewise->new);

ok($run->seq($gen));

#exon 1794 1935
#exon 2084 2180

ok(my $t    = Bio::EnsEMBL::Transcript->new());

ok(my $exon  = Bio::EnsEMBL::Exon->new);
ok(my $exon2 = Bio::EnsEMBL::Exon->new);

ok($exon->start(1787) == 1787);
ok($exon->end  (1935) == 1935);
ok($exon->strand (1)  == 1);

ok($exon2->start(2084) == 2084);
ok($exon2->end   (2180) == 2180);
ok($exon2->strand (1)   == 1);

ok($t->add_Exon($exon));
ok($t->add_Exon($exon2));

ok($run->add_Transcript($t));

$run->run;

ok(1);

foreach $t ( $run->output ) {
   ok($t->isa('Bio::EnsEMBL::Transcript'));

   foreach my $e ( @{$t->get_all_Exons} ) {
		 ok($e);
     print "Translation ",$t->translation->start," ",$t->translation->start_Exon,"\n";
     print "Translation ",$t->translation->end," ",$t->translation->end_Exon,"\n";
   }
}


ok(my $seqin2 = Bio::SeqIO->new( -file => "t/data/genomewise2.seq"));

ok(my $gen2 = $seqin2->next_seq());

ok(my $run2 = Bio::EnsEMBL::Pipeline::Runnable::Genomewise->new);

ok($run2->seq($gen2));


ok(my $t2 = Bio::EnsEMBL::Transcript->new());

ok(my $exon3 = Bio::EnsEMBL::Exon->new);
ok(my $exon4 = Bio::EnsEMBL::Exon->new);

ok($exon3->start(1001) == 1001);
ok($exon3->end  (1422) == 1422);
ok($exon3->strand (1)  == 1);

ok($exon4->start(1482) == 1482);
ok($exon4->end  (1869) == 1869);
ok($exon4->strand (1)  == 1);

ok($t2->add_Exon($exon3));
ok($t2->add_Exon($exon4));

ok($run2->add_Transcript($t2));

$run2->run;

ok(1);

