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
BEGIN { $| = 1; print "1..12\n"; 
	use vars qw($loaded); }

END { print "not ok 1\n" unless $loaded; }


use Bio::EnsEMBL::Pipeline::Runnable::Blast;
use Bio::PrimarySeq;
use Bio::Seq;
use Bio::SeqIO;

$loaded = 1;
print "ok 1\n";    # 1st test passed.

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
{ print "ok 2\n"; }

my $pepseqio = new Bio::SeqIO('-file' => "$pwd/t/data/AP000074.pep",  
			      '-format' => 'fasta');
my $pepseq   = $pepseqio->next_seq->seq;
my $peptide  =  Bio::PrimarySeq->new(  '-seq'         => $pepseq,
				       '-id'          => 'AP000074',
				       '-accession'   => 'AP000074',
				       '-moltype'     => 'protein');

unless ($peptide) 
{ print "not ok 3\n"; }
else
{ print "ok 3\n"; }

# First lets test the dna-dna blast
my $blastn = Bio::EnsEMBL::Pipeline::Runnable::Blast->new
    ('-query'    => $clone,
     '-program'  => 'wublastn',
     '-database' => "$pwd/t/data/AI053588.fa",
     '-threshold' => 1e-6,
     );

unless ($blastn)
{ print "not ok 4\n"; }
else
{ print "ok 4\n"; }

#run dna-dna                                                
eval {
    $blastn->run;
};

if ($@) {
    print $@;
    print ("not ok 5\n");
} else {
    print "ok 5\n"; 
}

foreach my $out ($blastn->output) {
    print $out->gffstring . "\n";
}

print "ok 6\n";

# Now the dna-pep blast
my $blastx = Bio::EnsEMBL::Pipeline::Runnable::Blast->new 
    ('-query'    => $clone,
     '-program'  => 'wublastx',
     '-database' => "$pwd/t/data/AP000074.pep",
     '-threshold' => 1e-6,
     );

unless ($blastx)
{ print "not ok 7\n"; }
else
{ print "ok 7\n"; }

#run dna-pep                                               
eval {
    $blastx->run;
};

if ($@) {
    print $@;
    print ("not ok 8\n");
} else {
    print "ok 8\n"; 
}

foreach my $out ($blastx->output) {
    print $out->gffstring . "\n";
}

print "ok 9\n";

# Now the pep-dna blast

my $tblastn = Bio::EnsEMBL::Pipeline::Runnable::Blast->new 
    ('-query'    => $peptide,
     '-program'  => 'wutblastn',
     '-database' => "$pwd/t/data/AI053588.fa",
     '-threshold' => 1e-6,
     );

unless ($tblastn)
{ print "not ok 10\n"; }
else
{ print "ok 10\n"; }

#run dna-pep                                               
eval {
    $tblastn->run;
};

if ($@) {
    print $@;
    print ("not ok 11\n");
} else {
    print "ok 11\n"; 
}

foreach my $out ($tblastn->output) {
    print $out->gffstring . "\n";
}

print "ok 12\n";


