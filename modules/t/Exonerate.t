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

use Bio::EnsEMBL::Pipeline::Runnable::Exonerate;
use Bio::PrimarySeq;
use Bio::Seq;
use Bio::SeqIO;

$loaded = 1;
print "ok 1\n";    # 1st test passed.

my ($est1) =  set_est1();
my ($est2) =  set_est2();
my ($genomic) =  set_genomic();

my $estseq1 =  Bio::PrimarySeq->new(  -seq         => $est1,
           	                      -id          => 'XX12345',
                	              -accession   => 'XX12345',
                        	      -moltype     => 'dna');
my $estseq2 =  Bio::PrimarySeq->new(  -seq         => $est2,
           	                      -id          => 'M93650',
                	              -accession   => 'M93650',
                        	      -moltype     => 'dna');

my $genseq =  Bio::PrimarySeq->new(  -seq         => $genomic,
                                     -id          => 'Z83307',
                                     -accession   => 'Z83307',
                                     -moltype     => 'dna');

my @ests;
my @genomics;

push(@ests,$estseq1);
push(@ests,$estseq2);
push(@genomics,$genseq);

foreach my $es(@ests) {
	
$es->isa("Bio::PrimarySeqI") || die("argh!");
}

unless (scalar(@ests) && scalar(@genomics)) 
{ print "not ok 2\n"; }
else
{ print "ok 2\n"; }

#create Exonerate object    
my $exonerate = Bio::EnsEMBL::Pipeline::Runnable::Exonerate->new (-EST       => \@ests, 
								  -GENOMIC   => \@genomics);
 
unless ($exonerate)
{ print "not ok 3\n"; }
else
{ print "ok 3\n"; }

#run exonerate
$exonerate->run();
print "ok 4\n"; # 4th test passed

#get and store the output
my @results = $exonerate->output();
display(@results);

unless (@results) 
{ print "not ok 5\n"; }
else
{ print "ok 5\n"; }

sub display {
  my @results = @_;

  foreach my $obj (@results)
    {
      print STDERR $obj->gffstring() . "\n";
    }
}

sub set_est1 {
  #embedded sequence! Because I can't create Bio::PrimarySeqs from files
my $seq = 
'cagaggtcaggcttcgctaatgggccagtgaggagcggtggaggcgaggccggcgccgca'.
'cacacacattaacacacttgagccatcaccaatcagcataggaatctgagaattgctctc'.
'acacaccaacccagcaacatccgtggagaaaactctcaccagcaactcctttaaaacacc'.
'gtcatttcaaaccattgtggtcttcaagcaacaacagcagcacaaaaaaccccaaccaaa'.
'caaaactcttgacagaagctgtgacaaccagaaaggatgcctcataaagggggaagactt'.
'taactaggggcgcgcagatgtgtgaggccttttattgtgagagtggacagacatccgaga'.
'tttcagagccccatattcgagccccgtggaatcccgcggcccccagccagagccagcatg'.
'cagaacagtcacagcggagtgaatcagctcggtggtgtctttgtcaacgggcggccactg'.
'ccggactccacccggcagaagattgtagagctacctcacagcggggcccggccgtgcgac'.
'atttcccgaattctgcaggtgtccaacggatgtgtgagtaaaattctgggcaggtattac'.
'gagactggctccatcagacccagggcaatcggtggtagtaaaccgagagtagcgactcca'.
'gaagttgtaagcaaaatagcccagtataagcgggagtgcccgtccatctttgcttgggaa'.
'atccgagacagattactgtccgagggggtctgtaccaacgataacataccaagcgtgtca'.
'tcaataaacagagttcttcgcaacctggctagcgaaaagcaacagatgggcgcagacggc'.
'atgtatgataaactaaggatgttgaacgggcagaccggaagctggggcacccgccctggt'.
'tggtatccggggacttcggtgccagggcaacctacgcaagatggctgccagcaacaggaa'.
'ggagggggagagaataccaactccatcagttccaacggagaagattcagatgaggctcaa'.
'atgcgacttcagctgaagcggaagctgcaaagaaatagaacatcctttacccaagagcaa'.
'attgaggccctggagaaagagtttgagagaacccattatccagatgtgtttgcccgagaa'.
'agactagcagccaaaatagatctacctgaagcaagaatacaggtatggttttctaatcga'.
'agggccaaatggagaagagaagaaaaactgaggaatcagagaagacaggccagcaacaca'.
'cctagtcatattcctatcagcagtagtttcagcaccagtgtctaccaaccaattccacaa'.
'cccaccacaccggtttcctccttcacatctggctccatgttgggccgaacagacacagcc'.
'ctcacaaacacctacagcgctctgccgcctatgcccagcttcaccatggcaaataacctg'.
'cctatgcaacccccagtccccagccagacctcctcatactcctgcatgctgcccaccagc'.
'ccttcggtgaatgggcggagttatgatacctacacccccccacatatgcagacacacatg'.
'aacagtcagccaatgggcacctcgggcaccacttcaacaggactcatttcccctggtgtg'.
'tcagttccagttcaagttcccggaagtgaacctgatatgtctcaatactggccaagatta'.
'cagtaaaaaaaaaaaaaa';

return $seq;
}

sub set_genomic {
#embedded sequence! Because I can't create Bio::PrimarySeqs from files
my $seq 	= 
'gatccggagcgacttccgcctatttccagaaattaagctcaaacttgacgtgcagctagt'.
'tttattttaaagacaaatgtcagagaggctcatcatattttcccccctcttctatatttg'.
'gagcttatttattgctaagaagctcaggctcctggcgtcaatttatcagtaggctccaag'.
'gagaagagaggagaggagaggagagctgaacagggagccacgtcttttcctgggagggct'.
'gctatctaagtcggggctgcaggttggagatttttaaggaagtggaaattggcaattggc'.
'tttgtgtgtctgtggtttttggggagggggactacaaagaggggctaactccctctccct'.
'attctctaaggttggaccacagggatgaggttgtgagatacaaagataaaggagggatgg'.
'ggaacactatgatgtggtatttttcttttctgttttttctttttgataatatctatcctt'.
'ggctaaggggagggcggagttcagcggcggaaataaagcgagcagtggctggtgcgaacc'.
'gactcccggctctaagccttacttgcggtggccggacttgtctgtggctgaagccgagcc'.
'cgggctctgactctcacgtctgcactggaggtgcggaaacctggactggggttcaccagc'.
'cacatactggctgctctggttgttctctcctcttccttcttctattctaaccaaacaaaa'.
'cccaactcgatgcttgtgccaggccttggccagtttggacggtgatggaaacatttctgg'.
'attttagcgtctaagggagcactcaagggctgtaaggtttgctatcactgtcccaacact'.
'gcagagaccttgaaggttcagtgtgggatctgtagaacccatggtttgcggtgacactca'.
'gaggggatctttgggagatttatagagccggttcttagcgctttgtgggttcagaggctg'.
'gctgcagtgtttatgaagaggggcagtgggctggggacactgctggggttatggctgtag'.
'tgaggtccatgtgtacttgttccttggcatgtgtctgactctgtgttgctgctgcagtac'.
'agtgggcaggggcacggttgcttggactgggctgcccgatgtctggcatggctggtggtc'.
'ctgttgtcctttatttgatcgatagcagggaactgaccgccgaggttggcacaggttggc'.
'agggggatgaggatgcattgtggttgtctcctcctcctctccttcttcctcttccccttc'.
'ctcctctcctttctgttcttgttctccctcatcttcctcttccttcttctccctcttctt'.
'cctcttcactctgctctcttctcttcttttcccctttcctctcctctccctttcctcagg'.
'tcacagcggagtgaatcagctcggtggtgtctttgtcaacgggcggccactgccggactc'.
'cacccggcagaagattgtagagctagctcacagcggggcccggccgtgcgacatttcccg'.
'aattctgcaggtgatcctcccggcgccgccccactcgccgcccccgcggccctccactct'.
'caacgccctctcttcatttcttactgtaaacgatgctaattatggacccccccaccctca'.
'cccaccccagtccccagtccccccacccactcccctgccttcccactttcccctctcctt'.
'ccaccatccttccctatctctttcaacctgggtcagttacataagataactcatctgctt'.
'cagaaaaatattgtgtgcggatttttttaaaaaaatctttctctttttatttggtaaaga'.
'cattcattgtgggaacagtttatgaaccgtaataagctgagttatataaaccggcataat'.
'ttcattctgctctccccagcccatgcttacctcagcttttgaggtatttgtttattcttc'.
'atgtttatgaataatatatattgattttaaaaggcaaatgctatttactcctccctaact'.
'cacatttactcaattgagtattttttaaagaagaagaagtaaaaaaatcagtttattttt'.
'tacctttgttctttaaaacataacgccactttaagcaaggtcagcacaaaaataaattta'.
'tctacttcgttttgatgcatcttcaggcagtgtttaagaaaagttttttttttaaaaaaa'.
'gcttttaaattcgttttttagttcaaattgtttgaaagtatcatcatatttgtagttttt'.
'agggctacaaatgtaattttaagaaaaaaagctctctacagtaagttctcataccattga'.
'aggtatatttttgtgttatagacccatgcagatgcaaaagtccaagtgctggacaatcaa'.
'aacgtaagcttgtcattgtttaatgcatacttaaacaattttatttttgtcttgaaatta'.
'ttaataatgtggttttctgtccacttcccctatgcaggtgtccaacggatgtgtgagtaa'.
'aattctgggcaggtattacgagactggctccatcagacccagggcaatcggtggtagtaa'.
'accgagagtagcgactccagaagttgtaagcaaaatagcccagtataagcgggagtgccc'.
'gtccatctttgcttgggaaatccgagacagattactgtccgagggggtctgtaccaacga'.
'taacataccaagcgtaagttcattgagaacatctgccctccctgccctaagcccaatgct'.
'ctctcctctttacctcctcccaccctctctctccactgtctcttagtctgtgtcctctcc'.
'ctgctccacatttgtctcctttgtacctgggggaacagagaggaatgccctgacttttct'.
'ttgactgtctggaaaatgggagtcaagtgtggggagtcattcacttcatttgcatgctgc'.
'aaaacagagggcggaggcaccagggaaaggcacttgaatgaagaaggaaaatgagaacca'.
'gattgtaacttcgtcctaatccacctgccagaactttccttcaggtgtcacacatccatt'.
'tccatcctaatattaaacaatataatgaaagaaagctttacaacaagtctaaatagtttt'.
'tattttcgggcagtcttttaacaaagcaaagcaaccgtgtgagttaggtcaccagagaca'.
'ccaagcaatggtgaaggaccccctccgcccaattctctatccaactaaatttccatgccc'.
'aaagtgatagctatcattttttccacggtgtatctgcaaatccacccactgtcccggggt'.
'ggctgggagctttttaacgggttgagagttgctttttaaggttgtgggtgagctgagatg'.
'ggtgactgtgtcttcaggagacactaccatttggtttgattttggtttgatttgcaggtg'.
'tcatcaataaacagagttcttcgcaacctggctagcgaaaagcaacagatgggcgcagac'.
'ggcatgtatgataaactaaggatgttgaacgggcagaccggaagctggggcacccgccct'.
'ggttggtatccggggacttcggtgccagggcaacctacgcaaggtaaaacccaagcagcc'.
'atccacgcagctctccatatgtgcatccctttgcctgtcccctactctcccaaccatttc'.
'ctctcagggcttccatgcttggggaatgtcatgggtgagactgcatttgaaggcctggga'.
'cacatcgaccacattgtatgtaggtggtgtttgttgaaccattttgttggcctgggaatg'.
'tgagggttgtgtatgtcaggagcggtatgaaggtggctggtgggtgcgtgtgtgtgtttg'.
'tgggcagcactgtgtgcataagaacgagctgcgtggtgcactattcaaatttgacattag'.
'aatgcagggagagcccccttgtagaagtgtgacaagataacgcactgacagactgcgagg'.
'ttaacataattgtatcctgttgctttttgctgtaattgacaaattatgtgactatgagta'.
'gggtgcacatgaaaaggtgagcagggagatgcagggcgagtcggggatggggcattggaa'.
'gtggtggtctcccaactaaaaccaccccactatccttgcccctcccccgttcgcccctct'.
'gggtcccactttctcactctctgtgggtggcgactatgtggaagtggtgaggtgaagata'.
'ctgtgagcaggacggcgggatgccgggagccgctctggtctgggaggagagatgggacct'.
'gaaacaccccaggttctgtgagaagggctggaggaagaagagcgggggatggggggggtg'.
'gggcggaggggagggtcagttggctcccctccgctcccctccacaccccgaggtagcttg'.
'gtcaagaataggaaacagtttactccgaagattaatcggtatttgttgtcctcgtgacaa'.
'ggagataatatgcgggataatgggacgtggcgtctcactcccaggagcacaacggagccg'.
'aggggctcggggccaggggcggcgcgcgggactgcgaggaggctgccggcgcatccgtga'.
'gctctcggaccggctccggatgcccggcgcctccatggagtcgggcgagggaagggaaaa'.
'tcggtaacagtatgcgaaatatctggtacaaaggatgcttctgtcactcgcaagaatttg'.
'ttatgggaacaatcctgtcgctctgtaattgctcattagagaacgttttgtttctcatta'.
'aggcgacatgaataagcgtctagttgaaggagacagctgtagaaatgttccacaagagac'.
'ctctggacacgattcggcaccaattgccaattagacatgtcagttttaagagcagaaaac'.
'aatataggacggacactgggaagatagggagcaaagcgctcgcttcttgtttcccagcgc'.
'ggccgctcggccggcggaggcctcctgacgcgcggggcccgggctcccctgggcgacccc'.
'gcccgcacgggccccgaacccgccgcctccgcgagctcagagggcccccaggcccgctcc'.
'cgcctctggtgagactagcgattgaccccaccgagtgtaggcaaccccatccctgcctga'.
'cttttgaaaccgtaggattcccacttgtagacactggacacacccaagtgacacccgaac'.
'ttcgcactcacgagcgtccacccctccgcccccagcaatctgaggtcctggtgatccctc'.
'ccctgacaagagtctggcccacttagtccgagaagtcttggggagggactcagggcacgg'.
'tggcctaggcctgggagggctgagcccagagcctccactccgggctctccccggtgcgtg'.
'ggctcccatgagcttggggcgcaggatacaccccccttccgcttcttggagttgggcaga'.
'agccgaggacagcagcttccaccagaggtggcgcccgattcgggcggattttcgtcgccg'.
'gcagtctgagggcgggaagtgagagtcaaatcagctagggtccgcctgtgatggttgggg'.
'ccggaggtggagctagtccaggggcttttcactttaacgcgcctagttgctagaaaattc'.
'tgggtcaactggcaatcaccacccccccacccgcaagtgaggtttttggggggtgattgt'.
'cacctcaggcacaggaccgtcgcccctagttgcgtggggaggacacgtgggccaaacgaa'.
'gcacagagttgccggcgtgggggtagcgggcaccggtggtgagcttcgcttcacgttccc'.
'ggggaatcctcagactttttggtgctaagagggtccccctccctccctagccccctttcc'.
'cagcccaggaagcctttcttgggttttaaagagtccatactccaaaagcataaggcgggg'.
'gttaggtctcagagggagacaaatatggtttccatctgatcaacgccatgcggcggtgaa'.
'taaaatcgcgacgacgtgggctgacaacaacaagagacaactctatccagccccaattct'.
'ccggcgatttggtggttttagggatcacggagacactttttcttatctcctcccctcacc'.
'ccacccccaccccccataatacccaagccctctataggatttagtcagcctctctttcaa'.
'cagatgctacaatgttttgatccgcgctccagagctgaaaaggtgccgcaaccgggttgc'.
'tatcttttcttgcttgttttccgcctcactgattaagacactaagaaactaaggaatgat'.
'tttatttcccagcaactctctctctgtctctctctctctctctctctctcacctccttcc'.
'aggaaacaaatcctattcccatcgtcagatggtcctggggctgggataatgggaggcgct'.
'ttcactcgctttctcacggattaggctggaaggtgaacggcacccactgaaggcccctgc'.
'tttgcgcggctccgagggggcgccttttgaagaagatatcttaattgtccaccacttagg'.
'tttatcgtgggggtgggggagtcgagatccaacttctagttttattttgttaaaagcttt'.
'aagagtggaaggcaaatccactaaagtggggggaaatgtacccattttatgtaagagcag'.
'ccgctaggtcaccgcctggctgcaaacagttgccacttctaaagtaatgaactgtctccg'.
'tgccgctctccccgcggggtagagaagggatcctgcagcgacaggtttctttttgctgtg'.
'gaattccaatccccgcttcccatcttttttttttttttttttttttccagatttaaaatt'.
'gagctctactgtccctgctgacttttctctcttaagtgtcagttctggaggcagcacagg'.
'gcctggcggcgatcgcgtgctgcctgtgtaaacccgtgtgtatgtttgtgtgagacagaa'.
'catggataagaatgtgaatctccatgttttaaaattaagaaaatgtcaaagaggcaaatg'.
'agagcagagcatgctatcggcggggttcttcggaggcagattcgggcaactttgtttaat'.
'tggcgagatcgaatgtgccagaaaagggaccgagatggggcccctaggaaggcatgggca'.
'tgagtggccatggggagatgggcattgtgttctcctcagtgtcccccaccccataccaga'.
'tgtgcacagatttttgtctgttttccaaaagcaataaccctccggcctttctagttggcc'.
'aaaaggagctatttccagtcttccttcttcaccaaacccgaaaataaagcgaggggtggg'.
'ggtggtcgcttcatttcttcccttgaaaaacgctttgaaaagtccccggaaacttggggc'.
'agtaaattagcacagacgcttgtgcccgcacacgtgagcggagcgcgtctccctcagcta'.
'agtagccctcctcgacccccacctttctgatggagtgcaaagacccagagtccagatggg'.
'gctcagttactaatattctctggccctcatcctccttttcctctctctcccactcttgtt'.
'cggccctcgtgtctctttttcctctcagcattcttttctgtccttttttaaaaatacctt'.
'tctagcctctatgttttttgagccttcttttcccctccatcgccctctttttccctccat'.
'cgccctctgtctttcccctaacccaggccccccaatccttcccccttgcttccaccccca'.
'ccgctggtttccatcctgcctcccctccatctcctctatgtaggaacccagccctggtcc'.
'ccgtccactctccctgcctccccccatgctgggccctggcccccccatccaccctccatc'.
'cccccagccaagcgccctaaatagcacggaggcgcccgctcttcggacagtgattaatga'.
'tagcagagcagaggggttaacacacttcactgaaaagtctgttgactgggcttcttgtaa'.
'cacaatgtggcccgctgcacgcctcaagagaatccttttgttgtccgcgctcattgtagc'.
'ctcaaaattctgcccacgaaagtttgccaacgctcctgccccaggagtttaatagtttcc'.
'cttactcgcggggcattgtgcagcgctgaaaagcagcccctcgctattcaagtgttggtg'.
'gtcatctcaatagatctccaagggcccatatggtggccagtgccgatgaatccgcctgtt'.
'taaatgggggagaaagttggggttttaaaacatttcaaagttcctgaaaagatcccacta'.
'gatcctgtcacaattccctgaactctttgaaggcgcagcctattgtctcctggttataaa'.
'taatattcctggccaagtcgattcccaccaaggcgctcctaggggcccctccaccccgcc'.
'tctggccaagttttgaggatagggaggtgggtagcccatgctggtatcccttgggggtca'.
'tttcaggaccccagacccaggacccagcctgtagtggccctggaggcccagtcagagttt'.
'aggcaatcctgcttcccttcactgtggccttacaggcaaggtgcaagccgggagccagct'.
'agcccaagtgcacatcttgccctccaggcagcaagggaaaaccataaaactttccccctg'.
'tataatgatagaagttacattcaaagtggtggaaatgccctaattaaaaatgtaccactt'.
'aaatgatatgccaaagattcagctgcagcatagaatatttagctagttagtgaactctct'.
'actatttcttttttaaaattacacgttaaaaatttaaaagaaatcttacttttctctggt'.
'gcaacacattacaaagaatggacagtccttttatctaaatataaaattccattttcagca'.
'attatagcctgttcttggtgatgatataattacagctgtgctcgtaataggtgtcaaggc'.
'gaagcgcactcatacaatagttgaataaaactgcagaacaatgcgagcacttctaaattc'.
'acaacctttaaatgaaattgactgtgcagaattcaaataaaaactagataaatatttttt'.
'aaaaaagattacaagcttggtttggtttcttaatagcattttctacaaacctatcaacat'.
'ctaattgattagataggaattactgattatagatcaactgaaatcttgaacatcactctt'.
'ctattttcccaacacagccataggtaaagagattctgtagggagtggagtgaggcatttg'.
'gggagttggggttacttataaaagatcattttaattggaaacttcagtgcttattttttc'.
'attcaagaaatgccacttagtgtgtgtattataaagtctcatcttgataaaaacaaagaa'.
'aatggtggtcaggtaactaacatcgcaaatgtatttttttaaaagaaggctgacagttac'.
'cttgggaatgttttggtgaggctgtcgggatataatgctcttggagtttaagactacacc'.
'aggccccttttggaggctccaagttaatccaaatttctcttaccatcctattctttttgt'.
'tccagatggctgccagcaacaggaaggagggggagagaataccaactccatcagttccaa'.
'cggagaagattcagatgaggctcaaatgcgacttcagctgaagcggaagctgcaaagaaa'.
'tagaacatcctttacccaagagcaaattgaggccctggagaaaggtgatagagtttttca'.
'aagtagagaagcagtaaatcaaagtaaatgccacatcttcagtacaaagagctaaattta'.
'gccagggccctttgcatagaagaatgaaaagatttccttttttctgtctttttatttctc'.
'tgggcatcttttcagtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgt'.
'gtgtgtgtgtgtggtgtgtgtgtgtgtgtttcttcttttcatctaccagtaattcaaaga'.
'ctaaatgtctgacttataaggaaaaatgatgatttggctatttcaggccacagaaaggtc'.
'actgaatgccattccaaagaaaatttaacttggttctggtgggaaagttcttccaagtac'.
'agtcaacactagaagcattttaaagggaattggttggaggtaatgggagtggggaggtgg'.
'gaaccagtttgatgcacagtttggtcaacatattttgtgtagttctggcacaatatggaa'.
'aatcaacttactctttcagagtttgagagaacccattatccagatgtgtttgcccgagaa'.
'agactagcagccaaaatagatctacctgaagcaagaatacaggtaccgagagactgtgca'.
'gtttcacactttgtgattcataccatttgtctttcctagagacagaggtgcttgtacaga'.
'gtactatttatttataggactaatataataaaaaggttcagtctgctaaatgctctgctg'.
'ccatgggcgtggggagggcagcagtggaggtgccaaggtggggctgggctcgacgtagac'.
'acagtgctaacctgtcccacctgatttccaggtatggttttctaatcgaagggccaaatg'.
'gagaagagaagaaaaactgaggaatcagagaagacaggccagcaacacacctagtcatat'.
'tcctatcagcagtagtttcagcaccagtgtctaccaaccaattccacaacccaccacacc'.
'gggtaatttgaaatactaatactacgaatcaatgtctttaaacctgtttgctccgggctc'.
'tgactctcactctgactactgtcatttctcttgccctcagtttcctccttcacatctggc'.
'tccatgttgggccgaacagacacagccctcacaaacacctacagcgctctgccgcctatg'.
'cccagcttcaccatggcaaataacctgcctatgcaagtaagtgcggctggtggtggcctg'.
'cataacccaggccccagagaagtgaggagtggctcagggcctgcggacctcattggctgt'.
'gtctgcacccttgagagcttttcgcactacagtgattggcttgaccagtcaagtcggaga'.
'cagtcaatcccatcacttttaagtgattgactcattaattcatgccctaaaaaaatgagt'.
'aataaaaatctgtccagttttgtcaggttgatctgccttttattatactgttaccttgat'.
'aatgttgggtggtggtggggcatgttttgggtaccagggagccttgcaccagaaagtgga'.
'aataatgctggcacattatcagatagcagattagtagtttaaaatttgggtttatattaa'.
'tgtgtttgtatgctaaatatagaatctgtgcacgcatttggggcattactttgggtatat'.
'gtgataaactagtgagaaagaaaaaaggatcagaaatgggattcatatttacatggtgag'.
'atatacataatataatgagaatgctagttttctgtctgtatctacaataagaaaaggcat'.
'agcaggtatttgctggaaatttagtgtgtctttgctgtgaatggtgtgacgagtttgtgg'.
'ccctcctagctgcctgggaagcttgatgctattcacttggtatgacagcctgcctctcct'.
'cttagttctgtcccaaatatctattagcccctacatttagaggtcctgcactaggttcac'.
'cctttatgatgtaagttggataaggcagatggtttgtactagacctttgttgctggatgg'.
'attcttgataggaaaaatgtctgtcttctggtaggcctttcccagtggttttcctagaac'.
'tcctgtttgtgcaacaattagagatattagatggtacgatattggccagcatgagcctct'.
'gctggaaacagttctggggctacactgattgtttattctccattgaacattttttgctgg'.
'atttcaaatccaaataacagcaaaataaatgtttcacagtcttcagactaatataggagc'.
'agctagataagcaacttcagaggaattattcacatgtttatttttattgcatctggatat'.
'tgttggccatagtgcaattgatgtaaattaagggattaacagcccattagttggtgttgc'.
'tataactgcgttggaattttccagaagtcaggttgcctagaggaactcattgcaggaatt'.
'agaaacaaatgcaagctgaaattctggcaggccctcaaggcctgtttgcctctttgaact'.
'tgatgttagcattcatcatctgactttaataaggccacagagtgtctcgttcagtttcat'.
'gtattataacaacatccaggtttctgtgaagatagcaaaatgtgtatgtgagaaaataat'.
'aagacgaacaagtagatgctgcaattatattagggctgtttccacatacagagcggtatg'.
'gggaatcaatctatttcaaacactgacttttaaaaattacataatttgtataattcaaaa'.
'agtacatgcattcattaagcataaagaaaatcaaaatgatcaataactttactcttcaga'.
'gaaaactactctttgaattttggagtgggtgagcccgattgtctatccatttattcattc'.
'tcttgaccaaaatcatcattttacaaatggtgggtgtcttcctctggaccttctctgata'.
'tgcttccctctctgagcatatggattttggaaattaaacatttctccagtttgcagagaa'.
'gagaaatgatagtatactgtactatcatctgacctgcttgttcccgacacacttctttta'.
'attcatcccaagtttgctgccatggcatagtaccatgagagctcatgaggcgttttgtta'.
'ggagaaaggtttacctcctcagtgctgccattccgaacagttcgagggtagcaacagttt'.
'tctagttataagtcagtctgggcttttgggtctgttaggaagtcatggttaattctcagg'.
'atggaggtcgttatatttactaacacttctgatcactttaacttggggtcattacagatc'.
'tgcttcttcaaagaatctttaatcccatcagtgaaaggttcccaaggtccctgaaccata'.
'ctttactgagagctcctcccactgctctctgtcccaaagcttgagatatgggctctgggg'.
'gtcatagggttcccaaataaatccagatttgcagggagaggggatgtgttttgatgaagg'.
'tcctcatgctataggttcttcaaagatgctagcagaggttagagacaaaaatcctattaa'.
'tttatggatagtggcaaccatctagctggacagttgtcagaacctaaagtgctttagagg'.
'cttgatacataggcagctttcttctagctgtggccagtggaaggactagctcgaggccca'.
'atcttagatttatcatatggaattccagtacttcacgtgaaggcatctttaatgatcaga'.
'cttgttggcagagttcctcgggaggagggagcctggggctgtggctgtgtgatgtgttcc'.
'tcagtaaccacaggtttgcctctctcctcacagcccccagtccccagccagacctcctca'.
'tactcctgcatgctgcccaccagcccttcggtgaatgggcggagttatgatacctacacc'.
'cccccacatatgcagacacacatgaacagtcagccaatgggcacctcgggcaccacttca'.
'acaggtgagccactgctttctgcaggctgcacagaggcgatctctcttcactagaagttt'.
'acccaaacagaatctcctggtcttatgggagggcgtgtttaactccttgctttccttgtc'.
'cctgggggatggggattgaaaagggaaattcagttaagctaattagtaactttacaccat'.
'atagacaaaaactaaaattgtttttcctgaatttggtcacaaaagttgtgtatgaagaca'.
'aggcctgagactgcaagttttctgaggacagattattagacgaagctcagtagggggccc'.
'actgagctgtaggtgcgtgcttgttgaaatgcttcttgccctcatagctcctctagacct'.
'tttgctggaaataaaaagtgacacattggttttccagagacagctttattgtaaaagttc'.
'caaacatgcaaacaaacagaggattttttttttcttttcctttggattggggtggggggt'.
'acttgggatccaataggtatatatacatatattgtctagtttctgaaggtgctactttta'.
'tttgtaacaattgaagtgattttaatacagtaaaaaatgttagaaagtattagttttttt'.
'ttttttttttttttttgtaaacctataaatttgtattccatgtctgtttctcaaagggaa'.
'tatctacatggctatttctttcatccacttctaggactcatttcccctggtgtgtcagtt'.
'ccagttcaagttcccggaagtgaacctgatatgtctcaatactggccaagattacagtaa'.
'aaaaaaaaaaaaaaaaaaaaaggaaaggaaatattgtgttaattcagtcagtgactatgg'.
'ggacacaacagttgagctttcaggaaagaaagaaaaatggctgttagagccgcttcagtt'.
'ctacaattgtgtcctgtattgtaccactggggaaggaatggacttgaaacaaggaccttt'.
'gtatacagaaggcacgatatcagttggaacaaatcttcattttggtatccaaacttttat'.
'tcattttggtgtattatttgtaaatgggcatttgtatgttataatgaaaaaaagaacaat'.
'gtagactggatggatgtttgatctgtgttggtcatgaagttgttttttttttttttaaaa'.
'agaaaaccatgatcaacaagctttgccacgaatttaagagttttatcaagatatatcgaa'.
'tacttctacccatctgttcatagtttatggactgatgttccaagtttgtatcattccttt'.
'gcatataattaaacctggaacaacatgcactagatttatgtcagaaatatctgttggttt'.
'tccaaaggttgttaacagatgaagtttatgtgcaaaaaagggtaagatataaattcaagg'.
'aagaaaaaaagttgatagctaaaaggtagagtgtgtcttcgatataatccaatttgtttt'.
'atgtcaaaatgtaagtatttgtcttccctagaaatcctcagaatgatttctataataaag'.
'ttaatttcatttatatttgacaagaatatagatgttttatacacattttcatgcaatcat'.
'acgtttcttttttggccagcaaaagttaattgttcttagatatagttgtattactgttca'.
'cggtccaatcattttgtgcatctagagttcattcctaatcaattaaaagtgcttgcaaga'.
'gttttaaacttaagtgttttgaagttgttcacaactacatatcaaaattaaccattgttg'.
'attgtaaaaaaccatgccaaagcctttgtatttcctttattatacagttttctttttaac'.
'cttatagtgtggtgttacaaattttatttccatgttagatcaacattctaaaccaatggt'.
'tactttcacacacactctgttttacatcctgatgatccttaaaaaataatccttatagat'.
'accataaatcaaaaacgtgttagaaaaaaattccacttacagcagggtgtagatctgtgc'.
'ccatttatacccacaacatatatacaaaatggtaacatttcccagttagccatttaattc'.
'taaagctcaaagtctagaaataatttaaaaatgcaacaagcgattagctaggaattgttt'.
'tttgaattaggactggcattttcaatctgggcagatttccattgtcagcctatttcaaca'.
'atgatttcactgaagtatattcaaaagtagatttcttaaaggagactttctgaaagctgt'.
'tgcctttttcaaataggccctctcccttttctgtctccctcccctttgcacaagaggcat'.
'catttcccattgaaccactacagctgttcccatttgaatcttgctttctgtgcggttgtg'.
'gatggttggagggtggaggggggatgttgcatgtcaaggaataatgagcacagacacatc'.
'aacagacaacaacaaagcagactgtgactggccggtgggaattaaaggccttcagtcatt'.
'ggcagcttaagccaaacattcccaaatctatgaagcagggcccattgttggtcagttgtt'.
'atttgcaatgaagcacagttctgatcatgtttaaagtggaggcacgcagggcaggagtgc'.
'ttgagcccaagcaaaggatggaaaaaaataagcctttgttgggtaaaaaaggactgtctg'.
'agactttcatttgttctgtgcaacatataagtcaatacagataagtcttcctctgcaaac'.
'ttcactaaaaagcctgggggttctggcagtctagattaaaatgcttgcacatgcagaaac'.
'ctctggggacaaagacacacttccactgaattatactctgctttaaaaaaatccccaaaa'.
'gcaaatgatcagaaatgtagaaattaatggaaggatttaaacatgaccttctcgttcaat'.
'atctactgttttttagttaaggaattacttgtgaacagataattgagattcattgctccg'.
'gcatgaaatatactaataattttattccaccagagttgctgcacatttggagacaccttc'.
'ctaagttgcagtttttgtatgtgtgcatgtagttttgttcagtgtcagcctgcactgcac'.
'agcagcacatttctgcaggggagtgagcacacatacgcactgttggtacaattgccggtg'.
'cagacatttctacctcctgacattttgcagcctacattccctgagggctgtgtgctgagg'.
'gaactgtcagagaagggctatgtgggagtgcatgccacagctgctggctggcttacttct'.
'tccttctcgctggctgtaatttccaccacggtcaggcagccagttccggcccacggttct'.
'gttgtgtagacagcagagactttggagacccggatgtcgcacgccaggtgcaagaggtgg'.
'gaatgggagaaaaggagtgacgtgggagcggagggtctgtatgtgtgcacttgggcacgt'.
'atatgtgtgctctgaaggtcaggattgccagggcaaagtagcacagtctggtatagtctg'.
'aagaagcggctgctcagctgcagaagccctctggtccggcaggatgggaacggctgcctt'.
'gccttctgcccacaccctagggacatgagctgtccttccaaacagagctccaggcactct'.
'cttggggacagcatggcaggctctgtgtggtagcagtgcctgggagttggccttttactc'.
'attgttgaaataatttttgtttattatttatttaacgatacatatatttatatatttatc'.
'aatggggtatctgcagggatgttttgacaccatcttccaggatggagattatttgtgaag'.
'acttcagtagaatcccaggactaaacgtctaaattttttctccaaacttgactgacttgg'.
'gaaaaccaggtgaatagaataagagctgaatgttttaagtaataaacgttcaaactgctc'.
'taagtaaaaaaatgcattttactgcaatgaatttctagaatatttttcccccaaagctat'.
'gcctcctaacccttaaatggtgaacaactggtttcttgctacagctcactgccatttctt'.
'cttactatcatcactaggtttcctaagattcactcatacagtattatttgaagattcagc'.
'tttgttctgtgaatgtcatcttaggattgtgtctatattcttttgcttatttctttttac'.
'tctgggcctctcatactagtaagattttaaaaagccttttcttctctgtatgtttggctc'.
'accaaggcgaaatatatattcttctctttttcatttctcaagaataaacctcatctgctt'.
'ttttgtttttctgtgttttggcttggtactgaatgactcaactgctcggttttaaagttc'.
'aaagtgtaagtacttagggttagtactgcttatttcaataatgttgacggtgactatctt'.
'tggaaagcagtaacatgctgtcttagaaatgacattaataatgggcttaaacaaatgaat'.
'aggggggtccccccactctccttttgtatgcctatgtgtgtctgatttgttaaaagatgg'.
'acagggaattgattgcagagtgtcgcttccttctaaagtagttttattttgtctactgtt'.
'agtatttaaagatcctggaggtggacataaggaataaatggaagagaaaagtagatattg'.
'tatggtggctactaaaaggaaattcaaaaagtcttagaacccgagcacctgagcaaactg'.
'cagtagtcaaaatatttatctcatgttaaagaaaggcaaatctagtgtaagaaatgagta'.
'ccatatagggttttgaagttcatatactagaaacacttaaaagatatcatttcagatatt'.
'acgtttggcattgttcttaagtatttatatctttgagtcaagctgataattaaaaaaaat'.
'ctgttaatggagtgtatatttcataatgtatcaaaatggtgtctatacctaaggtagcat'.
'tattgaagagagatatgtttatgtagtaagttattaacataatgagtaacaaataatgtt'.
'tccagaagaaaggaaaacacattttcagagtgcgtttttatcagaggaagacaaaaatac'.
'acacccctctccagtagcttatttttacaaagccggcccagtgaattagaaaaacaaagc'.
'acttggatatgatttttggaaagcccaggtacacttattattcaaaatgcacttttactg'.
'agtttgaaaagtttcttttatatttaaaataagggttcaaatatgcatattcaattttta'.
'tagtagttatctatttgcaaagcatatattaactagtaattggctgttaattttatagac'.
'atggtagccagggaagtatatcaatgacctattaagtattttgacaagcaatttacatat'.
'ctgatgacctcgtatctctttttcagcaagtcaaatgctatgtaattgttccattgtgtg'.
'ttgtataaaatgaatcaacacggtaagaaaaaggttagagttattaaaataataaactga'.
'ctaaaatactcatttgaatttattcagaatgttcataatgctttcaaaggacatagcaga'.
'gcttttgtggagtatccgcacaacattatttattatctatggactaaatcaattttttga'.
'agttgctttaaaatttaaaagcacctttgcttaatataaagccctttaattttaactgac'.
'agatcaattctgaaactttattttgaaaagaaaatggggaagaatctgtgtctttagaat'.
'taaaagaaatgaaaaaaataaacccgacattctaaaaaaatagaataagaaacctgattt'.
'ttagtactaatgaaatagcgggtgacaaaatagttgtctttttgattttgatcacaaaaa'.
'ataaactggtagtgacaggatatgatggagagatttgacatcctggcaaatcactgtcat'.
'tgattcaattattctaattctgaataaaagctgtatacagtatgtgtttatgctacagtg'.
'ggtttttttaagtgactgacattcatcatattggttagacagttttaaaaacctagtctt'.
'tgttttctacaatgttgtcaagaacaaataaagtataatttgtggtatattaaaagcaag'.
'ctgcttaaaaatgctagttataactgcttcaaaatattaaatacatttgaaaatgtattt'.
'tggaaattcagtttctcccaatggatgttaaacaccttttaaaaaatcaagactcttaaa'.
'tatgcaagttacttttcaactttcattttttatccattatctgtcaaaagctttgagaat'.
'atagaattgtaattaaaacacctaatgtttcatatgaaaactggttcgatatgctatacc'.
'ccccaaagacatgtttacatggttagaaatttacacactcaaccaatgttcattattaac'.
'aaaatatttatgtgatgcaatggaaatctgagttagtttcattttttcactttgtctcag'.
'acctatgttctgtaaaatctctcaacaagcttcagtatatgtttttctcctggtatccta'.
'gcattaatgtagctgtcagtttgaatttcaatcctgaatgatgcctgttgacataatttc'.
'atagtactgacttcaattgtttaggattactaaactaaagaccatccacttcttttttct'.
'atgtttcctatcactgaagttctttatatttaataactactatggtagactaattatctc'.
'tagtggaatctaaatcccaaaatattactttgtctttatcaatttaataactggatacat'.
'aaatgcaaatttattttatatattagtaaacacattttcatttaaaaaacagcatttttg'.
'gtggcaagtgactgtttctgtaataaagagaaagtatgttttcacatacaatcttgaatt'.
'tcatcatgttttgtatgctacctgctgtgacatgcactttaccatgacactaaatgcctt'.
'agtacacacacgaatgaccactgtcaggaatcttgagattttgacattttttggcgttta'.
'tttaaaaataaaattattttctctagagcatgaaggaaaaatttaaaagtatttctgctg'.
'tgctagttctgtaaaaagttaaaaaatgagaaatggcggcaccattttataatcctgtgt'.
'ttttcaaggagtggagggtcaaattgttaagaaaaatatttacataagctattagcaatc'.
'ataattagagtgtcatagaattctgcatgcagcgactaaggaggaatccctagaagtcca'.
'ggtgcttcttgcctccggccatcatgccacagcctgggcccagccgcttggcggattctg'.
'ccagatccattttgcttgagcggctcactgtgtctgacaagtctggaggcaaatgcagtc'.
'gctatggcagaagatttcaggagaagaaaaaaagaagaaagaagaaaaaaaattaatttt'.
'gctttggctagcatacttcagattacaatattatgctctgttgatcagcatttagacttg'.
'catataaaattcctgctaaaacaattcttcagactgaattccagtataatcacctcctgt'.
'tcctactggattttacaaattacagccatcacagaaaccgtaccaggaatcttattggaa'.
'cattttcagtgtccaaaccaaaatgtattaagatcttttctcattaggtagaagtgaaca'.
'catgtaaacacataaataacagcaacgttggaaatatcgcataaataacggcagatgtct'.
'tttttctttttggctttgttttattttgccagagctccaacatgcaatttttaaggtcaa'.
'atgtagcttttcagcaattatagattcccatgaatacttgtttgtctttgtctattagtt'.
'gacctgacgattctttttcaaggttatactttctctaggaacttaaaattttcacagaaa'.
'tgttgacaacataaaaaaatgataccaaccggagggtcaattttttcatatgttcaggta'.
'agtcagccacagcatctattaccggtcttggtaagtttgtttctatgcataatgttgggt'.
'catgttttcaaaagtatagtgttgctttgaacaccctcccaacccccgcccatcccaact'.
'gtttttctgcagatacaggctgatttaaactgaccccagcaactgaccacaattacagaa'.
'gttactgcaattagggagaaacatccatatttcaaaataacatttttctgttttcaaaag'.
'aatatcaaattcatttaactctccgtgctcccagctcgcaaaatttatttcataaaaatc'.
'caaacttaaaggaacttatctgttgtgtgagacacaaagtggtgtgggaggctaaagata'.
'aggcagcatagggctccccactgatgactacacacagccttctagaggagaactgacctg'.
'gagacaacctagctgaaacccttcatttggaaaattctcttcaatatgggggggaaaaac'.
'agatgaaaaaggggagaaccatattattttggtcaaaatattgtggtccacaagcatatg'.
'ctccagttagtttctttcttgaataaaggctttttattgtcatgtaaacacaagctgtgt'.
'gcacatgatcaaaatattttaaaactaaaaataatttatgaaaaaatattcttccttgat'.
'ttcaacctgcctgtacttatttttaatacaaatatatctaggataaaagatactattata'.
'caaatgcatgatcaaggaagatgtcagaaaggttaacggggtcaagaaaagctgtaacac'.
'tcatagagtaatatccatacagaactattccttagtatccatgggacccagcc';

return $seq;
}

sub set_est2 {
  #embedded sequence! Because I can't create Bio::PrimarySeqs from files
my $seq =
'CAGAGGTCAGGCTTCGCTAATGGGCCAGTGAGGAGCGGTGGAGGCGAGGCCGGCGCCGCA'.
'CACACACATTAACACACTTGAGCCATCACCAATCAGCATAGGAATCTGAGAATTGCTCTC'.
'ACACACCAACCCAGCAACATCCGTGGAGAAAACTCTCACCAGCAACTCCTTTAAAACACC'.
'GTCATTTCAAACCATTGTGGTCTTCAAGCAACAACAGCAGCACAAAAAACCCCAACCAAA'.
'CAAAACTCTTGACAGAAGCTGTGACAACCAGAAAGGATGCCTCATAAAGGGGGAAGACTT'.
'TAACTAGGGGCGCGCAGATGTGTGAGGCCTTTTATTGTGAGAGTGGACAGACATCCGAGA'.
'TTTCAGAGCCCCATATTCGAGCCCCGTGGAATCCCGCGGCCCCCAGCCAGAGCCAGCATG'.
'CAGAACAGTCACAGCGGAGTGAATCAGCTCGGTGGTGTCTTTGTCAACGGGCGGCCACTG'.
'CCGGACTCCACCCGGCAGAAGATTGTAGAGCTAGCTCACAGCGGGGCCCGGCCGTGCGAC'.
'ATTTCCCGAATTCTGCAGGTGTCCAACGGATGTGTGAGTAAAATTCTGGGCAGGTATTAC'.
'GAGACTGGCTCCATCAGACCCAGGGCAATCGGTGGTAGTAAACCGAGAGTAGCGACTCCA'.
'GAAGTTGTAAGCAAAATAGCCCAGTATAAGCGGGAGTGCCCGTCCATCTTTGCTTGGGAA'.
'ATCCGAGACAGATTACTGTCCGAGGGGGTCTGTACCAACGATAACATACCAAGCGTGTCA'.
'TCAATAAACAGAGTTCTTCGCAACCTGGCTAGCGAAAAGCAACAGATGGGCGCAGACGGC'.
'ATGTATGATAAACTAAGGATGTTGAACGGGCAGACCGGAAGCTGGGGCACCCGCCCTGGT'.
'TGGTATCCGGGGACTTCGGTGCCAGGGCAACCTACGCAAGATGGCTGCCAGCAACAGGAA'.
'GGAGGGGGAGAGAATACCAACTCCATCAGTTCCAACGGAGAAGATTCAGATGAGGCTCAA'.
'ATGCGACTTCAGCTGAAGCGGAAGCTGCAAAGAAATAGAACATCCTTTACCCAAGAGCAA'.
'ATTGAGGCCCTGGAGAAAGAGTTTGAGAGAACCCATTATCCAGATGTGTTTGCCCGAGAA'.
'AGACTAGCAGCCAAAATAGATCTACCTGAAGCAAGAATACAGGTATGGTTTTCTAATCGA'.
'AGGGCCAAATGGAGAAGAGAAGAAAAACTGAGGAATCAGAGAAGACAGGCCAGCAACACA'.
'CCTAGTCATATTCCTATCAGCAGTAGTTTCAGCACCAGTGTCTACCAACCAATTCCACAA'.
'CCCACCACACCGGTTTCCTCCTTCACATCTGGCTCCATGTTGGGCCGAACAGACACAGCC'.
'CTCACAAACACCTACAGCGCTCTGCCGCCTATGCCCAGCTTCACCATGGCAAATAACCTG'.
'CCTATGCAACCCCCAGTCCCCAGCCAGACCTCCTCATACTCCTGCATGCTGCCCACCAGC'.
'CCTTCGGTGAATGGGCGGAGTTATGATACCTACACCCCCCCACATATGCAGACACACATG'.
'AACAGTCAGCCAATGGGCACCTCGGGCACCACTTCAACAGGACTCATTTCCCCTGGTGTG'.
'TCAGTTCCAGTTCAAGTTCCCGGAAGTGAACCTGATATGTCTCAATACTGGCCAAGATTA'.
'CAGTAAAAAAAAAAAAAA';

return $seq;
}
