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
#BEGIN { $| = 1; print "1..8\n"; 
BEGIN { $| = 1; print "1..3\n"; 
	use vars qw($loaded); }

END { print "not ok 1\n" unless $loaded; }

use lib 't';
use EnsTestDB;
use Bio::EnsEMBL::Pipeline::RunnableDB::Riken_BlastMiniGenewise;
use Bio::PrimarySeq;
use Bio::Seq;
use Bio::SeqIO;

$loaded = 1;
print "ok 1\n";    # 1st test passed.

my $ens_test = EnsTestDB->new();
# Load some data into the db
$ens_test->do_sql_file("t/BlastMiniGenewise.dump");

# Get an EnsEMBL db object for the test db
my $db = $ens_test->get_DBSQL_Obj;
print "ok 2\n";    

#my $id=('chr1.75000-100000');
my $id=('chr1.961228-999570');

my $runnable = 'Bio::EnsEMBL::Pipeline::RunnableDB::Riken_BlastMiniGenewise';
my $rbmg = "$runnable"->new(-dbobj    => $db,
                           -input_id => $id);	

unless ($rbmg)
{ print "not ok 3\n"; }
else
{ print "ok 3\n"; }

$rbmg->fetch_input();
$rbmg->run();

my @genes = $rbmg->output();
print STDERR "found " . scalar(@genes) . " genes\n";

# check those genes
my $count=0;
foreach my $gene(@genes){
   $count++;
   my $ecount = 0;
   print STDERR "gene $count\n";
     foreach my $exon ( $gene->each_unique_Exon() ) {
       $ecount++;
       print STDERR 	"exon $ecount\t" .
			$exon->contig_id . "\t" . 
			$exon->start     . "\t" . 
			$exon->end       . "\t" . 
			$exon->strand    . "\n";
       foreach my $sf($exon->each_Supporting_Feature) {
	  print STDERR "\tsupporting_features: " . 
                       $sf->seqname     . "\t" .
		       $sf->start       . "\t" .
		       $sf->end         . "\t" .
		       $sf->strand      . "\t" .
		       $sf->hseqname    . "\t" .
		       $sf->hstart      . "\t" .
		       $sf->hend        . "\n";
	}
     }
}

# not working yet ...
$rbmg->write_output();













