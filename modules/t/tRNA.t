use lib 't';
use strict;
use Test;

BEGIN { $| = 1; plan test => 5;}

use Bio::EnsEMBL::Pipeline::Runnable::tRNAscan_SE;
use Bio::PrimarySeq;
use Bio::Seq;
use Bio::SeqIO;

ok(1);

ok(my $seqstr =  set_seq);

ok(my $seq =  Bio::PrimarySeq->new(-seq         => $seqstr,
				   -id          => 'HS97D16',
				   -accession   => 'AL009179',
				   -moltype     => 'dna'));



ok(my $tRNA = Bio::EnsEMBL::Pipeline::Runnable::tRNAscan_SE->new (-QUERY => $seq));
 
ok($tRNA->run);

ok(my @results = $tRNA->output);

ok(display(@results));


sub display {
  my @results = @_;

  my @methods = qw( seqname start end score strand );

  foreach my $obj (@results) {

      print "\n";
      foreach my $method_name (@methods) {
        my $value = $obj->$method_name();
        printf ("%10s = $value\n", $method_name);
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
