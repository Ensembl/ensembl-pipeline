use lib 't';
use Test;
use strict;

BEGIN { $| = 1; plan test => 7;}

use Bio::EnsEMBL::Pipeline::Runnable::TRF;
use Bio::PrimarySeq;
use Bio::Seq;
use Bio::SeqIO;

ok(1);

ok(my $sequence =  set_seq());

ok(my $seq =    Bio::PrimarySeq->new(	-seq         => $sequence,
					-id          => 'HS97D16',
					-accession   => 'AL009179',
					-moltype     => 'dna'));


ok(my $trf = Bio::EnsEMBL::Pipeline::Runnable::TRF->new (-QUERY => $seq));
 
ok($trf->run);

ok(my @results = $trf->output);

ok(display(@results));

sub display {
  my @results = @_;
  my @methods = qw( start end score strand );

  foreach my $obj (@results) {
    printf "\n";
    foreach my $method_name (@methods) {
      my $value = $obj->$method_name();
      printf ("%10s = $value\n", $method_name);
    }
  }
  1;
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
