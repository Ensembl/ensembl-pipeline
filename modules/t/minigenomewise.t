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


use Bio::EnsEMBL::Pipeline::Runnable::MiniGenomewise;
use Bio::PrimarySeq;
use Bio::Seq;
use Bio::SeqIO;

$loaded = 1;
print "ok 1\n";

# Obtain an input genomic sequence
my $seqin = Bio::SeqIO->new( -file => "t/data/human.genomic");
my $genomic_seq = $seqin->next_seq();

# Create a transcript to pass to the MiniGenomewise object.
my $transcript = Bio::EnsEMBL::Transcript->new();

  # Create a couple of exons to add to the transcript object
  my $exon = Bio::EnsEMBL::Exon->new;

  $exon->start(1787);
  $exon->end  (1935);
  $exon->strand (1);

  $transcript->add_Exon($exon);

  $exon = Bio::EnsEMBL::Exon->new;

  $exon->start(2084);
  $exon->end  (2180);
  $exon->strand (1);

  $transcript->add_Exon($exon);

# Create new MiniGenomewise object
my @transcripts = ($transcript);
my $minigenomewise = Bio::EnsEMBL::Pipeline::Runnable::MiniGenomewise->new(
									-genomic     => $genomic_seq,
									-transcripts => \@transcripts
									);

unless ($minigenomewise){
  print "not ok 3\n";
} else {
  print "ok 3\n";
}

# Run the runnable.
$minigenomewise->run;

print "ok 4\n";

# Taste the output.

my @results = $minigenomewise->output();

unless (@results){
  print "not ok 5\n";
} else {
  print "ok 5\n";
}

my $error = 0;
foreach $trans (@results ) {
   if( !$trans->isa('Bio::EnsEMBL::Transcript') ) {
       $error = 1;
   }
}

if ($error){
  print "not ok 6\n";
} else {
  print "ok 6\n";
}
