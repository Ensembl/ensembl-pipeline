
## We start with some black magic to print on failure.
BEGIN { $| = 1; print "1..3\n"; 
	use vars qw($loaded); }

END { print "not ok 1\n" unless $loaded; }


use Bio::EnsEMBL::Pipeline::Runnable::Hmmpfam;

use Bio::PrimarySeq;
use Bio::Seq;
use Bio::SeqIO;

$loaded = 1;
print "ok 1\n";

my $seqin = Bio::SeqIO->new( -file => "t/road.pep");

my $pep = $seqin->next_seq();

my $run = Bio::EnsEMBL::Pipeline::Runnable::Hmmpfam->new( -peptide => $pep, -database => "./t/rrm.HMM");

$run->run;

print "ok 2\n";

$count = 0;
foreach $sf ( $run->output ) {
	$count++
	#print $sf->start,":",$sf->end,"\n";
}

if( $count == 2 ) { 
	print "ok 3\n";
} else {
	print "not ok 3\n";
}
