use lib 't';
use strict;
use Test;

BEGIN { $| = 1; plan test => 27;}

use Bio::EnsEMBL::Pipeline::Runnable::BlastDB;
use Bio::SeqIO;

ok(1);

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
							   -type      => 'PROTEIN'));

ok($db1->type eq 'PROTEIN');
ok(scalar($db1->sequences) == 2);

ok($db1->run);

ok(my $dbfile1 = $db1->dbfile);
ok(my $dbname1 = $db1->dbname);

ok($db1->remove_index_files);


ok(my $db2 = new Bio::EnsEMBL::Pipeline::Runnable::BlastDB(-sequences => \@sequences2,
							   -type      => 'DNA'));

ok($db2->type eq 'DNA');
ok(scalar($db2->sequences) == 2);

ok($db2->run);

ok(my $dbfile2 = $db2->dbfile);
ok(my $dbname2 = $db2->dbname);

ok($db2->remove_index_files);

ok(my $db3 = new Bio::EnsEMBL::Pipeline::Runnable::BlastDB(-dbfile    => $pepfile,
							   -type      => 'PROTEIN'));

ok($db3->type eq 'PROTEIN');

ok($db3->run);

ok(my $dbfile3 = $db3->dbfile);
ok($db3->dbname eq 'testpep.fa.safe');

ok($db3->remove_index_files);

ok(my $db4 = new Bio::EnsEMBL::Pipeline::Runnable::BlastDB(-dbfile    => $estfile,
							   -type      => 'DNA'));

ok($db4->type eq 'DNA');

ok($db4->run);

ok(my $dbfile4 = $db4->dbfile);
ok($db4->dbname eq 'testest.fa.safe');

ok($db4->remove_index_files);



