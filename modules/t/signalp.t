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
BEGIN { $| = 1; print "1..4\n"; 
	use vars qw($loaded); }

END { print "not ok 1\n" unless $loaded; }


use Bio::EnsEMBL::Pipeline::Runnable::Protein::Signalp;
use Bio::PrimarySeq;
use Bio::Seq;
use Bio::SeqIO;

$loaded = 1;
print "ok 1\n";    # 1st test passed.

my $pwd = `pwd`;
chomp($pwd);

my $pepseqio = new Bio::SeqIO('-file' => "$pwd/t/data/seq_test.pep",  
			      '-format' => 'fasta');
my $pepseq   = $pepseqio->next_seq->seq;

my $peptide  =  Bio::PrimarySeq->new(  '-seq'         => $pepseq,
				       '-id'          => 'ENSP1',
				       '-accession'   => 'ENSP1',
				       '-moltype'     => 'protein');

unless ($peptide) 
{ print "not ok 2\n"; }
else
{ print "ok 2\n"; }

my $blastn = Bio::EnsEMBL::Pipeline::Runnable::Protein::Signalp->new('-clone'    => $peptide,
							   '-program'  => '/usr/local/ensembl/bin/signalp',
							   '-workdir' => "$pwd/t/data");

unless ($blastn)
{ print "not ok 3\n"; }
else
{ print "ok 3\n"; }


eval {
    $blastn->run;
};

if ($@) {
    print $@;
    print ("not ok 4\n");
} else {
    print "ok 4\n"; 
}

foreach my $out ($blastn->output) {
    print STDERR $out->gffstring . "\n";
}

