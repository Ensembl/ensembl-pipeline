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


use Bio::EnsEMBL::Pipeline::Runnable::Blast;
use Bio::PrimarySeq;
use Bio::Seq;
use Bio::SeqIO;

$loaded = 1;
print "ok 1\n";    # 1st test passed.

my $seqio = new Bio::SeqIO(-file => '/nfs/croc/michele/trunk/ensembl-pipeline/scripts/AP000074.pep',
			   -format => 'fasta');

my $seq   = $seqio->next_seq->seq;

my $clone =  Bio::PrimarySeq->new(  -seq         => $seq,
                                    -id          => 'AP000074',
                                    -accession   => 'AP000074',
                                    -moltype     => 'protein');
unless ($clone) 
{ print "not ok 2\n"; }
else
{ print "ok 2\n"; }

#create blast object    
my $blast = Bio::EnsEMBL::Pipeline::Runnable::Blast->new (   -query    => $clone,
                                                             -program  => 'wutblastn',
                                                             -database => '/nfs/croc/michele/trunk/ensembl-pipeline/scripts/AI053588.fa',
                                                             -threshold => 1,
                                                             );
 
unless ($blast)
{ print "not ok 3\n"; }
else
{ print "ok 3\n"; }

#run Inverted                                                
$blast->run('/tmp/');
print "ok 4\n"; # 4th test passed

#get and store the output
my @results = $blast->output();
display (@results);

unless (@results) 
{ print "not ok 5\n"; }
else
{ print "ok 5\n"; }

my @methods = qw( seqname start end strand hseqname hstart hend hstrand
                  percent_id p_value hpercent_id hp_value score hscore);
#Display output
sub display {
    my @results = @_;
    #Display output
    foreach my $obj (@results)
    {
       print STDERR ($obj->gffstring."\n");
       if ($obj->sub_SeqFeature)
       {
            foreach my $exon ($obj->sub_SeqFeature)
            {
                print STDERR "Sub: ".$exon->gffstring."\n";
            }
       }
    }
}

