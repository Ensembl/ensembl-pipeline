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


use Bio::EnsEMBL::Pipeline::Runnable::Est2Genome;
use Bio::PrimarySeq;
use Bio::Seq;
use Bio::SeqIO;

$loaded = 1;
print "ok 1\n";    # 1st test passed.
my ($genomic_seq, $est_seq) =  set_seq();


my $gen=    Bio::PrimarySeq->new(  -seq         => $genomic_seq,
                                   -id          => 'HSAC74',
                                   -accession   => 'ACOOOO74',
                                   -moltype     => 'dna');

my $est=    Bio::PrimarySeq->new(  -seq         => $est_seq,
                                   -id          => 'AI937824',
                                   -accession   => 'AI937824',
                                   -moltype     => 'dna');
unless ($gen && $est) 
{ print "not ok 2\n"; }
else
{ print "ok 2\n"; }

#create Est2Genome object    
my $e2g = Bio::EnsEMBL::Pipeline::Runnable::Est2Genome->new( -genomic => $gen,
                                                             -est => $est );
unless ($e2g) 
{ print "not ok 3\n"; }
else
{ print "ok 3\n"; }

#run Est2Genome                                                
$e2g->run();
print "ok 4\n"; # 4th test passed

#get and store the output
my @results = $e2g->output();

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
#The files genomic.seq and est.seq MUST be present for testing
open (GENOMIC, "<genomic.seq");
my $genomic = <GENOMIC>;
open (EST, "<est.seq");
my $est = <EST>;
return ($genomic, $est);
}
