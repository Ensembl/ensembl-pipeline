## Bioperl Test Harness Script for Modules
##
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'
#-----------------------------------------------------------------------
## perl test harness expects the following output syntax only!
## 1..3
## ok 1  [not ok 1 (if test fails)]
## 2..3
## ok 2  [not ok 2 (if test fails)]
## 3..3
## ok 3  [not ok 3 (if test fails)]
##
## etc. etc. etc. (continue on for each tested function in the .t file)
#-----------------------------------------------------------------------


## We start with some black magic to print on failure.
BEGIN { $| = 1; print "1..4\n"; 
	use vars qw($loaded); }

END { print "not ok 1\n" unless $loaded; }


use Bio::EnsEMBL::Pipeline::Runnable::CrossMatch;
use Bio::PrimarySeq;
use Bio::Seq;
use Bio::SeqIO;

$loaded = 1;
print "ok 1\n";    # 1st test passed.

$str1 = &set_seq;
$str2 = &set_seq2;

$seq1 = Bio::PrimarySeq->new(-id => 'seq1',-seq => $str1);
$seq2 = Bio::PrimarySeq->new(-id => 'seq2',-seq => $str2);


$runnable = Bio::EnsEMBL::Pipeline::Runnable::CrossMatch->new( -seq1 => $seq1, -seq2 => $seq2, -score => 50);


$runnable->run();

print "ok 2\n";

@fp = $runnable->output();


print "ok 3\n";



if( scalar(@fp) != 3 ) {	
  print "not ok 4\n";
  print STDERR "got ",scalar(@fp)," hits\n";
} else {
  print "ok 4\n";
}

#foreach $s ( @fp ) {
#	print STDERR "got $s ",$s->start,":",$s->end,"\n";
#	}



sub set_seq {
#embedded sequence! Because I can't create Bio::PrimarySeqs from files
my $seq = 
'cctgggctgcctggggaagcacccagggccagggagtgtgaccctgcaggctccacacaggactgccagaggcacac'.
'acctgctctgtctacccgagggcaccagagggcacgagaaggctggctccctggcgctgacacgtcaggcaactgag'.
'gcacaaggctggcatttctgaaccttgcccctctgcaaacacaagggggcgatggtggcactccaagcaaaggggcg'.
'tgtgggtgctgcaggaggagcacagagcactggcgcccctcccctcccgccctgcagatgccggaggccccgcctct'.
'gctgttggcagctgtgttgctgggcctggtgctgctggtggtgctgctgctgcttctgaggcactggggctggggcc'.
'tgtgccttatcggctggaacgagttcatcctgcagcccatccacaacctgctcatgggtgacaccaaggagcagcgc'.
'atcctgaaccatgtgctgcagcatgcggagcccgggaacgcacagagcgtgctggaggccattgacacctactgcga'.
'gcagaaggagtgggccatgaacgtgggcgacaagaaaggtggggtccgggccagcaggtgctcagctctgggacagg'.
'gacccaggaccaggcatcaaagcccttacaggagaagctgttatcaccccatttccagggggctgggaaccctggga'.
'tatgcccagatagggctggggggctcctctggagtcccagggtgccagggtccctgatgacccctgcaggccctgct'.
'gcctgctgccccaggacaacaggcccccacactcacagggtctgacggtggtgcagttccccttgaactctgttctg';

return $seq;
}



sub set_seq2 {
#embedded sequence! Because I can't create Bio::PrimarySeqs from files
my $seq = 
'cctgggctgcctggggaagcacccagggccagggagtgtgaccctgcaggctccacacaggactgccagaggcacac'.
'acctgctctgtctacccgagggcaccagagggcacgagaaggctggctccctggcgctgacacgtcaggcaactgag'.
'gcacaaggctggcatttctgaaccttgcccctctgcaaacacaagggggcgatggtggcactccaagcaaaggggcg'.
'tgtgggtgctgcaggaggagcacagagcactggcgcccctcccctcccgccctgcagatgccggaggccccgcctct'.
'tgtgccttatcggctggaacgagttcatcctgcagcccatccacaacctgctcatgggtgacaccaaggagcagcgc'.
'atcctgaaccatgtgctgcagcatgcggagcccgggaacgcacagagcgtgctggaggccattgacacctactgcga'.
'gcagaaggagtgggccatgaacgtgggcgacaagaaaggtggggtccgggccagcaggtgctcagctctgggacagg'.
'gacccaggaccaggcatcaaagcccttaggagaagctgttatcaccccatttccagggggctgggaaccctggga'.
'tatgcccagatagggctggggggctcctctggagtcccagggtgccagggtccctgatgacccctgcaggccctgct';


return $seq;
}



