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
BEGIN { $| = 1; print "1..8\n"; 
	use vars qw($loaded); }

END { print "not ok 1\n" unless $loaded; }


use Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise;
use Bio::PrimarySeq;
use Bio::Seq;
use Bio::SeqIO;

$loaded = 1;

print "ok 1\n";    # 1st test passed.

my $in  = new Bio::SeqIO(-file => 't/data/AC004663.fa',
                          -format => 'fasta');

my $gen = $in->next_seq;

my $protid = "SW:NTC3_MOUSE";

my @feature = ($protid);

my $blast = new Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise(-genomic => $gen,
    -ids => \@feature);


$blast->run;

unless ($gen) 
{ print "not ok 2\n"; }
else
{ print "ok 2\n"; }

#create minigenewise object
my $mgw = Bio::EnsEMBL::Pipeline::Runnable::MiniGenewise->new(-genomic  => $gen,
							      -features => \@features);

unless ($mgw) 
{ print "not ok 3\n"; }
else
{ print "ok 3\n"; }

#run genewise
$mgw->minirun();
print "ok 4\n"; # 4th test passed

#get and store the output
my @results = $mgw->output();

unless (@results) 
{ print "not ok 5\n"; }
else
{ print "ok 5\n"; }

#Display output
foreach my $obj (@results)
{
    print "$obj\n";

}


my @revfeatures;

push(@revfeatures,$fp5);
push(@revfeatures,$fp6);
push(@revfeatures,$fp7);


#create minigenewise object
my $mgwrev = Bio::EnsEMBL::Pipeline::Runnable::MiniGenewise->new(-genomic  => $gen,
								 -features => \@revfeatures);

unless ($mgwrev) 
{ print "not ok 6\n"; }
else
{ print "ok 6\n"; }

#run genewise
$mgwrev->minirun();
print "ok 7\n"; # 4th test passed

#get and store the output
@results = $mgwrev->output();

unless (@results) 
{ print "not ok 8\n"; }
else
{ print "ok 8\n"; }

#Display output
foreach my $obj (@results)
{
    print "$obj\n";

}















