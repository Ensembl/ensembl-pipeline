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
BEGIN { $| = 1; print "1..3\n"; 
	use vars qw($loaded); }

END { print "not ok 1\n" unless $loaded; }


use Bio::EnsEMBL::Pipeline::Runnable::Genomewise;
use Bio::PrimarySeq;
use Bio::Seq;
use Bio::SeqIO;

$loaded = 1;
print "ok 1\n";

my $seqin = Bio::SeqIO->new( -file => "t/human.genomic");

my $gen = $seqin->next_seq();

my $run = Bio::EnsEMBL::Pipeline::Runnable::Genomewise->new;

$run->seq($gen);

#exon 1794 1935
#exon 2084 2180

my $t = Bio::EnsEMBL::Transcript->new();
my $exon = Bio::EnsEMBL::Exon->new;

$exon->start(1794);
$exon->end  (1935);
$exon->strand (1);

$t->add_Exon($exon);

my $exon = Bio::EnsEMBL::Exon->new;

$exon->start(2084);
$exon->end  (2180);
$exon->strand (1);

$t->add_Exon($exon);


$run->add_Transcript($t);

$run->run;

print "ok 2\n";

$seen = 0;
$error = 0;
foreach $t ( $run->output ) {
   if( !$t->isa('Bio::EnsEMBL::Transcript') ) {
      $error = 1;
   } else {
     foreach $e ( $t->each_Exon ) {
        #print "Exon ",$e->start," ",$e->end,"\n";
     }
   }

   $seen = 1;
}

if( $seen == 1 ) {
    print "ok 3\n";
} else {
    print "not ok 3\n";
}

