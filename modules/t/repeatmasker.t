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


use Bio::EnsEMBL::Pipeline::Runnable::RepeatMasker;
use Bio::PrimarySeq;
use Bio::Seq;
use Bio::SeqIO;

$loaded = 1;
print "ok 1\n";    # 1st test passed.
my ($seq) =  set_seq();


my $clone =  Bio::PrimarySeq->new(  -seq         => $seq,
                                    -id          => 'HSAC74',
                                    -accession   => 'ACOOOO74',
                                    -moltype     => 'dna');


unless ($clone) 
{ print "not ok 2\n"; }
else
{ print "ok 2\n"; }

#create RepeatMasker object    
my $repmask = Bio::EnsEMBL::Pipeline::Runnable::RepeatMasker->new (-CLONE => $clone);
 
unless ($repmask)
{ print "not ok 3\n"; }
else
{ print "ok 3\n"; }

#run RepeatMasker                                                
$repmask->run();
print "ok 4\n"; # 4th test passed

#get and store the output
my @results = $repmask->output();

unless (@results) 
{ print "not ok 5\n"; }
else
{ print "ok 5\n"; }

my @methods = qw( seqname start end strand hseqname hstart hend hstrand );
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
#embedded sequence! Because I can't create Bio::PrimarySeqs from files
open (SEQ, "<genomic.seq");
my $seq = <SEQ>;
return $seq;
}
