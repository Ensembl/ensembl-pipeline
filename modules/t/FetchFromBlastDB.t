use lib 't';
use strict;
use Test;

BEGIN { $| = 1; plan test => 17;}

use Bio::EnsEMBL::Pipeline::Runnable::BlastDB;
use Bio::EnsEMBL::Pipeline::SeqFetcher::FetchFromBlastDB;
use Bio::SeqIO;

ok(1); #1

my $tmppepfile = 't/data/testpep.fa';
my $tmpestfile = 't/data/testest.fa';

system("cp $tmppepfile $tmppepfile.safe");
system("cp $tmpestfile $tmpestfile.safe");

my $pepfile = $tmppepfile .".safe";
my $estfile = $tmpestfile .".safe";

open(IN1,"<$pepfile");
open(IN2,"<$estfile");

my $seqio1 = new Bio::SeqIO(-fh => \*IN1,-format => 'fasta');
my $seqio2 = new Bio::SeqIO(-fh => \*IN2,-format => 'fasta');

my @sequences1;
my @sequences2;

while (my $seq = $seqio1->next_seq) {
  push(@sequences1,$seq);
}

while (my $seq = $seqio2->next_seq) {
  push(@sequences2,$seq);
}

close(IN1);
close(IN2);

ok(my $db1 = new Bio::EnsEMBL::Pipeline::Runnable::BlastDB(-sequences => \@sequences1,
							   -type      => 'PROTEIN')); #2

ok($db1->type eq 'PROTEIN'); #3
ok(scalar($db1->sequences) == 2); #4

ok($db1->run); #5

ok(my $prot_fetcher = Bio::EnsEMBL::Pipeline::SeqFetcher::FetchFromBlastDB->new(-db => $db1));

ok(my $prot_seq = $prot_fetcher->fetch($sequences1[0]->display_id)); #7

ok($prot_seq && ($prot_seq->isa("Bio::Seq")) && 
   ($prot_seq->display_id eq $sequences1[0]->display_id));

ok($db1->remove_index_files);

ok(my $db2 = new Bio::EnsEMBL::Pipeline::Runnable::BlastDB(-sequences => \@sequences2,
							   -type      => 'DNA'));

ok($db2->type eq 'DNA');
ok(scalar($db2->sequences) == 2);

ok($db2->run);

ok(my $dna_fetcher = Bio::EnsEMBL::Pipeline::SeqFetcher::FetchFromBlastDB->new(-db => $db2));

ok(my $dna_seq = $dna_fetcher->fetch($sequences2[0]->display_id)); #7

ok($dna_seq && ($dna_seq->isa("Bio::Seq")) && 
   ($dna_seq->display_id eq $sequences2[0]->display_id));


ok($db2->remove_index_files);

