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
BEGIN { $| = 1; print "1..6\n"; 
	use vars qw($loaded); }

END { print "not ok 1\n" unless $loaded; }


use Bio::EnsEMBL::Pipeline::Runnable::Blast;
use Bio::EnsEMBL::Pipeline::Runnable::SearchFilter;
use Bio::PrimarySeq;
use Bio::Seq;
use Bio::SeqIO;

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
my $blast = Bio::EnsEMBL::Pipeline::Runnable::Blast->new (   -query => $clone,
                                                             -program => 'wublastn',
                                                             -database    => 'dbSTS-1',
                                                             -threshold => 1,
                                                          );
 
unless ($blast)
{ print "not ok 3\n"; }
else
{ print "ok 3\n"; }


my $search = Bio::EnsEMBL::Pipeline::Runnable::SearchFilter->new( -runnable => $blast,
							         -coverage => 1, 
								 -minscore => 150);

unless ($search)
{ print "not ok 4\n"; }
else
{ print "ok 4\n"; }


$search->run('/tmp/');
print "ok 5\n"; # 4th test passed

#get and store the output
my @results = $search->output();
display (@results);

unless (@results) 
{ print "not ok 6\n"; }
else
{ print "ok 6\n"; }


#Display output
sub display {
    my @results = @_;
    #Display output
    foreach my $obj (@results)
    {
       print STDERR ($obj->gffstring."\n");
       if ($obj->sub_SeqFeature)
       {
            foreach my $exon ($obj->sub_SeqFeature)
            {
                print STDERR "Sub: ".$exon->gffstring."\n";
            }
       }
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
