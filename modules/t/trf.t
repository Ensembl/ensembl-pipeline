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


use Bio::EnsEMBL::Pipeline::Runnable::TRF;
use Bio::PrimarySeq;
use Bio::Seq;
use Bio::SeqIO;

$loaded = 1;
print "ok 1\n";    # 1st test passed.
my ($sequence) =  set_seq();


my $seq	   =    Bio::PrimarySeq->new(	-seq         => $sequence,
					-id          => 'HS97D16',
					-accession   => 'AL009179',
					-moltype     => 'dna');


unless ($seq) 
{ print "not ok 2\n"; }
else
{ print "ok 2\n"; }


#create TRF object    
my $trf = Bio::EnsEMBL::Pipeline::Runnable::TRF->new (-QUERY => $seq);
 
unless ($trf)
{ print "not ok 3\n"; }
else
{ print "ok 3\n"; }

#run TRF                                      
$trf->run();
print "ok 4\n"; # 4th test passed

#get and store the output
my @results = $trf->output();
display(@results);

unless (@results) 
{ print "not ok 5\n"; }
else
{ print "ok 5\n"; }

sub display {
  my @results = @_;
  my @methods = qw( start end score strand );

  foreach my $obj (@results)
    {
      printf STDERR "\n";
      foreach my $method_name (@methods) {
        my $value = $obj->$method_name();
        printf STDERR ("%10s = $value\n", $method_name);
      }
    }
}

sub set_seq {
  #embedded sequence! Because I can't create Bio::PrimarySeqs from files
my $seq = 
'gtaaatttttataaaataataaattacccctacacaacacaggtgaggcttcgtggctta'.
'gctggttaaagcgcctgtctagtaaacaggagatcctgggttcgaatcccagcgaggcct'.
'ctttatttcttcccctaaacttaggtaaattcttgtcactagttaaacctgacttagatg'.
'tatacctaataaagcagcaaaccttgaatctcgatgtctgaggcactttgtagaggaacg'.
'tggtgatgggttggaaatatttacgatgaattagctaaatctatgccagagaccaaccgt'.
'ctagaaccatctagaattgagctgcaattttctacgtttaaatttttgtagctacttttc'.
'atatttctaaggaacttacagcttcacgtaacatttttgcttttaatacaattcaaaagg'.
'aaggcaaggaaggagatgagacatacatatggtataaattgagtgttgtgataacttcct'.
'ccttgcatgattctgaactctagaatttcactgaatatattgaaaagtggagtaacatga'.
'caacgtcttgagaggataagaagaatgattagtgaaacaccaaaagaaaacacaaacaaa'.
'aacaaaaaacctaacaagaggctggaagcaaatggggcctgaagttgtggatcatttatt'.
'agacagataaaaataaaatgatggaaaaaattaatcaaaccctggatacagtgctatgtt'.
'ttagctgctctttaaagaaaaaaaaattatcagttatatcttcaaagtatcacattctac'.
'aaatgccttccctcataacattccttgttgcagtattgcgctcagctttttgcttctccc'.
'tgtaaggaatgattacattggatgattactcctttgtgagcatttaggagaaaaaccagt'.
'tacgtaattgatgagtgcagtgcaaaatggaaatgttccttgtttaaggaaattaggaat'.
'ttcaaaccagcaatagcagagcattaaatcaagtgcctggccttgtgcaaatgcgtaggt'.
'cgcatgcccatgaagacagtcgtgactgggaagacggtcctgattgggaatgtggttatt'.
'tgtcaaaagactggattttagagacaccttgaaatggagtgaaatctaaagtgatagtca'.
'atgatatcaagtacagctagaaggcttgagcacttaaaagtactaatcaggggtgaggcg'. 
'cggtggctcacgcctgtaatcccagcacttcgggaggctgaggtgggcggatcacaaggt'.
'caggagtgtgagagcagcctggccaatatagtgaaaccccgtctccactaaaaatacaaa'.
'aagtagccaggcgtggtggcgcacgcctgtagtctcagctacttgggaggctgaggcagg'.
'agaatcgcttgaacccaggaggcggaggttacagccactgcaccccagcctgggcaacag'.
'agcaagactctatctcaaaaaaaaaaaaaaaagtactaagcagggaaatgttccagtata'.
'aaaagaattctgaaagtcatgcttttattagtcttctctaacaaagtaacatttcttata'.
'caattcaagttactcttttctcaggcaactccacctccttccaaatcccatcaaatggca'.
'ccagcctgggatctttcaacatagatatatgctccttcttcctttagcctcctgttgtct'.
'ttgggctcccacggggctccctctttaaatctggtgtaggaagagagcaagtcgtattcg'.
'tagtactgagttacacagcgtaaaccacttacaagcgtatagtcatttattcgaatggac'.
'ggggttctaagagcgtgagtcctaaaaacggtcgctgtaaattttaggaagacaagcggt'.
'cgctataagttttaggaagacaaaactgaggacaatacaatgttcccgcctgatttcgaa'.
'ccggggacttttcgcgtgtgaggcgaacataataaccactacactacggaaacccctaac'.
'tcatattatagaaacccctaagagaaggaactatcctgctggacctaatcaggaattccc'.
'tattaatcatcccgccttcctagaaatttgttttctttctattcaagtcagtgtctacaa'.
'tcgcattttaaaactggtacaaaacaagcaccacatcgagattaggaaatatccaaaacc'.
'cagaaattgctaagtatttgagttctatgagccttatacta'
;

return $seq;
}
