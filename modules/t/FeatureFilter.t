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


use Bio::EnsEMBL::Pipeline::Runnable::FeatureFilter;
use Bio::EnsEMBL::Pipeline::Runnable::Blast;
use Bio::Seq;
use Bio::SeqIO;
use Bio::PrimarySeq;
use strict;

$loaded = 1;
print "ok 1\n";    # 1st test passed.

# Now that were all loaded, the first thing to grab are some Feature
# objects to filter.  Using blastn for this.

# We'll need a Seq object to run our blast.

my $pwd = `pwd`;
chomp($pwd);

my $dnaseqio = new Bio::SeqIO('-file' => "$pwd/t/data/AP000074.fa", 
			      '-format' => 'fasta');
my $dnaseq   = $dnaseqio->next_seq->seq;
my $clone    =  Bio::PrimarySeq->new(  '-seq'         => $dnaseq,
				       '-id'          => 'AP000074',
				       '-accession'   => 'AP000074',
				       '-moltype'     => 'protein');
unless ($clone) 
{ print "not ok 2\n"; }
else
{ print "ok 2\n"; } # Weve got our Seq

# Run this Seq against a nucleotide db from the t/data directory.

my $blastn = Bio::EnsEMBL::Pipeline::Runnable::Blast->new
    ('-query'    => $clone,
     '-program'  => 'wublastn',
     '-database' => "$pwd/t/data/AI053588.fa",
     '-threshold' => 1e-6,
     );

unless ($blastn)
{ print "not ok 3\n"; }
else
{ print "ok 3\n"; }  # Now we're all set.

# Run our blast.
                                                
eval {
    $blastn->run;
};

if ($@) {
    print $@;
    print ("not ok 4\n");
} else {
    print "ok 4\n"; # Good, the blast worked and we have some Feature objects now.
}


# Start playing with FeatureFilter...

# Get a Bio::EnsEMBL::Pipeline::Runnable::FeatureFilter object.

my $feat_filter = Bio::EnsEMBL::Pipeline::Runnable::FeatureFilter->new(
							-coverage  => 5,
							-minscore  => 100,
							-maxevalue => 0.001,
							-prune     => 1
				                                 );

if ($feat_filter) {
  print "ok 5\n";     # Now we're set to run FeatureFilter.
}else {
  print "not ok 5\n";
}

# Try running a FeatureFilter job.

my @filtered_features = $feat_filter->run($blastn->output);

if (@filtered_features) {
  print "ok 6\n";         # It worked, we're done.
}else {
  print "not ok 6\n";
}
