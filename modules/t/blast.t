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
BEGIN { $| = 1; print "1..5\n"; 
	use vars qw($loaded); }

END { print "not ok 1\n" unless $loaded; }


use Bio::EnsEMBL::Pipeline::Runnable::Blast;
use Bio::PrimarySeq;
use Bio::Seq;
use Bio::SeqIO;
$| = 1;
$loaded = 1;
print "ok 1\n";    # 1st test passed.
my ($seq) =  set_seq();


my $clone =  Bio::PrimarySeq->new(  -seq         => $seq,
                                    -id          => 'AI937824',
                                    -accession   => 'AI937824',
                                    -moltype     => 'dna');
unless ($clone) 
{ print "not ok 2\n"; }
else
{ print "ok 2\n"; }

#create blast object    
my $pwd = `pwd`;
chomp($pwd);

my $blast = Bio::EnsEMBL::Pipeline::Runnable::Blast->new (   -CLONE => $clone,
                                                             -BLAST => 'blastn',
                                                             -DB    => "$pwd/t/data/AL035422.fa");
 
unless ($blast)
{ print "not ok 3\n"; }
else
{ print "ok 3\n"; }

#run Inverted                                                
$blast->run();
print "ok 4\n"; # 4th test passed

#get and store the output
my @results = $blast->output();

unless (@results) 
{ print "not ok 5\n"; }
else
{ print "ok 5\n"; }

my @methods = qw( seqname start end strand hseqname hstart hend hstrand
                  score percent_id p_value);
#Display output
foreach my $obj (@results)
{
    print "\n";
    foreach my $method_name (@methods) {
        my $value = $obj->$method_name();
        printf ("%10s = $value\n", $method_name);
    }
}

sub set_seq {

my $seq = 
"ttttttttttttttttagtggcagaaagcaaaaagggcctttgctgctttaatttttaaa".
"ttttcttacaaaaatttaggtgtttaccaatagtcttattttggcttatttttaatgctt".
"tttctcagtgtttttcttctgtttctgagtcacgaacagcaggcactgaaagcagtcccc".
"cagccactgccgaaggtcagtcccggaggtgctgcccaggcttggctggggcgaggcgct".
"gccaacccctgccgccagggggctccaagctccacggcacgatctgctcagggtggccct".
"tcttccacgatccaagccctaagaacaagaggctgggcctgggcccaggcatgcaagcta".
"tacccaagaatcaacgggctgaggcttagcgtcccctaccgcgtccaccagcctgaccgc".
"gggcctgctgggcccggggggaggggccttcctgctggggttcgagtttcaaggtgttcc".
"ttctcagaggcctttcctgcagagcttgcgcctcgtggtctacgctccctctggtccacc".
"gtgggatattccgtctgggcccagtgcatccttgttgatcacatctccgctgcc";
return $seq;
}
