use lib 't';
use strict;
use Test;

BEGIN { $| = 1; plan test => 7;}

use Bio::EnsEMBL::Pipeline::Runnable::RepeatMasker;
use Bio::PrimarySeq;
use Bio::Seq;
use Bio::SeqIO;

ok(1);
ok(my ($seq) =  set_seq());


ok(my $clone =  Bio::PrimarySeq->new(-seq         => $seq,
				     -id          => 'HSAC74',
				     -accession   => 'ACOOOO74',
				     -moltype     => 'dna'));



ok(my $repmask = Bio::EnsEMBL::Pipeline::Runnable::RepeatMasker->new (-QUERY => $clone));

ok($repmask->run());

ok(my @results = $repmask->output());

my @methods = qw( seqname start end strand repeat_consensus hstart hend hstrand );

foreach my $obj (@results) {
    foreach my $method_name (@methods) {
      my $value = "";
      
      if (defined($obj->$method_name)) {
	$value = $obj->$method_name();
      }
      printf ("%10s = $value\n", $method_name);
    }
  }

ok(1);

sub set_seq {
  #embedded sequence! Because I can't create Bio::PrimarySeqs from files
	  my $seq = 
	  'cctgggctgcctggggaagcacccagggccagggagtgtgaccctgcaggctccacacaggactgccagaggcacac'.
	  'acctgctctgtctacccgagggcaccagagggcacgagaaggctggctccctggcgctgacacgtcaggcaactgag'.
  'ggcaccaccaagggctcaaaggctggggagctccagggagagggctgcagaggctgccagtgcagtgtcttctcccc'.
  'tgcacactcagcccacccactccgcccccaatgatgagcatgaagtttgtaggagtgtggagaggtgtagtgcccta'.
  'ctaagtgggccggacattggtggcggtggctggcccattgcctccacggcttaccgggacgggcacccaccccgcaa'.
  'ggtgccgtgtggcccatgccctgccttagggagggcaggcgcagtggaggcggaggcggggtgtaaaaggcccagac'.
  'atggagctgagccccatcctctgtactctggggtgggcacagggaccaaagcatccctgaggagcgaactgacaggg'.
  'ttcacagctaggctgggatcagaggcagacagccacagggcctgggagccctgccctaccatgccctacggccccac'.
  'tgaggcactgtaggagacagctgtccaccaggctgaggggcctcagcgggcaggtgcagcagtcaccacgctggcca'.
  'ccagcctgagtgcacagggccgtccgcctgcaacgccctttcccccctcattggaagacggggttagaccccttctt'.
  'cccagcccctctatcaccatggaaacacccaggcttctcctgagggacagagcctctttgtggtgacacaaagcccc'.
  'tttgaggggaggtggggtgcagcccctgggagcctgtggccaaagggcccctggggggtaaagcagggcccgccggc'.
  'ccagcccccctcccaggtttcggagccaaggaagctgaggtaggagagccttcaggttaggccagcccgccccagct'.
  'gcccacccagaagggtctccggatggatgagccggacagcccaggaggccaggctgaggctgggcagagacagcacc'.
  'gagaagcctctctgcaggagctggggggttcctgctttcctctctgttgttttgacaaagaacacactctgcttaga'.
  'taagtggaaaacagcactggcaggggagctcccgccagtcaggaaatgcccaggatggctgcaggaacagtgcctgg'.
  'tgggagatggggacaggaagtcctccgtgggagggcaggcaggggccgtgagacggcaggtcagacagtgggaggtc'.
  'tggcaccaggtgtgccgagaggactggccctggcatctgggccccactcagtttactgtcctctcacttctcaagcc'.
  'tcagcctccgtcttcacctatccgcccagggggaaatgtcactggccttgctgcagccacaactctgggacaggcca'.
  'gagttggataaaacggaacctcaccatggcttggcttattccagggctggggaaatggaaccctgaggatgatgtgg'.
  'gtccctggcacagtgacagcatggctgtcccatccatgctacagccccacagccagggaggcagaggaggccaggcc'.
  'acatttccccagactcccagccctggagcccctcctgtactggcctgagagtccctgcccccacccactgactgctg'.
  'ctgaccaggcacctggcagctccctgcatgtgggcgccctgcctagaggccctcccgagttgacctcactctgaccc'.
  'acagcagctttgctaactaagcccctggggtggggttggggtggggatcactggccccctgcggcagcacgaggagg'.
  'cactggcccagggagagatccatccttttgtaggcaaatctacctgactactagggagtggggccaagggctgcccc'.
  'cacagggcagtgcctgggtcagcaagcacttgctgcacccccaatggaggccctcaccacacccctcccaacctgcc'.
  'ttccccactggcaccccactctctgctgacacctgcaggcacaggggtgggggttgggaccctccagctgccgtggc'.
  'agggagtggccatagggtggttgctgtggtggacgtaggccatatggttcagaatggatcagaggatggtgtgggtg'.
  'ggagaattaagtctgggctgaccgcaagggatcaggcctgggtggctggagggaggggccgggaccctgagtgggag'.
  'ggtgcaaagggagtctctgaggctggagaggctgagggaccactagatgcatgaagttctggaagaaaggggttcag'.
  'gctggagacagaaccaccggcaccttcagcctgtgaagccgagagataagacaggtgcaagccccagtgcctcctcc'.
  'aaaggaccctcccctgaccttctagaacaaggggacaaagtcactttcccagccatcatcttccctgcactggaagg'.
  'ggcagaacgaggacccagcactgagaggtgagcctcaggaaggaggcagtgctgccagcccttggggacaacagcct'.
  'gtcccctgaaagtcccaaaacccttgggaaaggtgctcctcactgtacccttctcccacccagcacatcagaggagc'.
  'cagtggtgtgagatcttgacccagcatgaagtagagcccaccctgtcatatgcagggaggttccatgtcccctttgg'.
  'tctgcccttgccccagggtcctgtgtggcctctggccagcgctggctgctcccagggtgcctggggcaccgctcatc'.
  'ccttccagagcctgagctggatgctggccttgccgctcccaaccccctcctctgctgggagccagctccactccgga'.
  'gtgcccatggcaacaggcctagacaccaacctccagcaagtgcacaaagtcctgcctttcatatgcaagcaacacag'.
  'agggccaggcagggggctgcaccaccccatccagggcagattccaggccctgcagcctcccatctggctcccaccgc'.
  'cacataccaggagagcctcgggacttgactgtgtgtccagtccccactacaaagaggtgggagaagtgctgacagtt'.
  'tggggccctcctcaaaggcacctactggcacatgtgtgcagccacccacactgctacggtgtgaatgctagtgttcc'.
  'ccaaattcacattaaaacctaatctccaatgtgacagtattaagaggtagggcctttaagaggtgattagtgccctt'.
  'atatgatcagtgtccttatcaaagcggctcttgagggctggtttaccccttccaccatgagaggacacacagcgcca'.
  'tctatcaatagccctcaccagacacggaatctgctggcgcctcgatcttggacttcccagcctccagcactgtgaac'.
  'agtatgtttctgttgtttataaaggactccatctaaggtatttggtggcagagcgtgaatgaactaagacacacaca'.
  'gacacacacacaagcacacccacgcacagggacccacagctgccccacaggcacacgcaccagcaggcacactcgtc'.
  'cacacccatatggctgctgcttatcccacctagaagtcccgacaggccctttggggggtgggggacatcttggctcc'.
  'tggtgaagatcagcatgctaggcagaccccattcctgtcttggacaccgggcaaggcagaccccctctgagacccct'.
  'gaggccgggtcaagcctgcccctgacaccaagatctgcaggaccccctgcctctgatcccacatgaggttcctctga'.
  'gacaggtgggagtggctgggcctgcagaccccgagcttgcacagacccacctgcagggacaagcacaacttacagca'.
  'gcctagtgcccagggaagacccttccaccatgggaggactcaggaagccgcttcccgctgagtcacatttagcgcat'.
  'ttcttaggaaacagctggtatgtgccctagaaggttaaagacggtgttcaactccgctggccacagaccacaatgaa'.
  'aaggagctgctggctgggcacgatggctcatacctataatcccagcactttggaaggccaaggtgggaggatcacct'.
  'gaggttgggagttcgagacctgcctgaccaacatggagaaaccccatctctactaaagatacaaacaaaattagctg'.
  'ggtgtggtggtgcatgcctgtaatcccagctactcaggactgaggctgaggcaggaggcagaggttgtggtgagcca'.
  'ggatcacctcactgtactccagcctgggcaacaaaagcaaaacttcatctcaaaaaaagaaaagaaagaaagaaaag'.
  'gagccgctataaataatgaggatcccatggagatttaaaggaaatcccagcagaggcagtaaacagcagactgtaga'.
  'aaacacgcagtgtgaggctggagaaagtacactaagacggtcctcccaaccccaacagcaaacaggagagatgacaa'.
  'gaggatggtaggagccacagaagacagagcccagggccaaacttccacatgataggagtttccagagcgagaacgac'.
  'agacagtgctcttccaaccacagaaaaacggccagaaaccacctcaacaggagtacgtccaacggggctatcagcgt'.
  'cagggctggaatctgagttagggtgtgagccttctgtccagaaaagaagcccaggaagtcacggatgtgtctacctc'.
  'tggcttttacaattcagtctcatttcctctggctcctaaatgcgctgtgtccacttcgaggtgtagggaggagggtg'.
  'cccagctctgacctgcgggctgtgctgaagcctgggcctctctctataccaggccgtctgcaggcctgcaggggctg'.
  'catactgtccctctgtaaatgtcagcctcaacacgctgctgcctgctgcctggctcccgtctctcacagctcagaag'.
  'ctgttttatgataaaaaagagctgacagctttgttaaaacaaacaaacgcctctcctctctacttgagatcttttac'.
  'aaagttaacatcgctctgaaaccccataaggaaaaattatggtgaggtttctggattttcagaggcagcaaactcat'.
  'cttacctagagccctggactacccagagaagtagttttaagaaccaaaagaagagaagaattcaggagtaaatgagg'.
  'cctggataacccaagggatcctgtgcttgcactgctgggattcaaatttgggggctgccaccaccacaagtttgtga'.
  'gaagcccacagcctgcgcatgaagattgtgtcaagggaatagcagctgctgatgtgtggaattgattttttcagcaa'.
  'gggcctgcgattctgaatgcattaaggcagctccagttgtctgtagcagtcctctgtcagctgcccattggaagcat'.
  'ttttactggctcaggaaattcatgggtttgtggactgatgactccaacagattaggtacaaatccaagaaccacaaa'.
  'acaaaatctaaccaagccaaacaacaatcaaacaataaaacagaacaaaaaaaaacctcaatttgttttccttttct'.
  'ttgtataacattaaaataaaaagatggtaacagactaccataaagaattatacaccaacaaattagataacccagat'.
  'gaaatggaaaaatttcctagaaacacacaacctaccaagactgaatcatgaagaaatagaaaacattgattagattt'.
  'acaactagcaacaaaactgaatcagtaatcaaaaacctccagcaaagatggtttcactggtgtattctaccaaacat'.
  'tcaagtaattttaacaccaatcgtcctcaaacttttcccaaaaataaagaggaaggaacacttcctaattcatttta'.
  'tgaggtcagccttaccctaataccaaagccaggaaaagacactacaagaaaagaaaaccacaggccaataaccctta'.
  'tgaatactgaaacaaaaattctcaaaacaaacaaacaaacaaaaaaacactagcaaaccaaattcagctgcatttta'.
  'caagaattataaacctaaccaagtgggattttactcctggaatgcaaggatggttccatgtacaaaaatcaatcaat'.
  'gtaatacgccacagtaatagaataaaggacaaaaatcacaatcatgtcaaatgatgcaggaaaaacatttaacaaaa'.
  'ttcaatagtctgtcatgataaaaacactcaaactaggaaaagaaggaaactttctcaacatgatgaaagctatatgt'.
  'gaaaagcccataacatcatactcaatggtgaaagactgaaagcttttctcctaagatcataaacaagacaaagatgg'.
  'ctgttttgaccacttccattcaacatagtactagaagttttagacagagtaattaggcaagaaaaagaaataaaaga'.
  'catccaaattagaaagaaagaaaaaccatctctgttcatagacgacataaccttatataaaaccgtaaagattctac'.
  'atgcattttttaaaaactctattagatgtaaaaaatgaaatcagcaaaggtgcaggataaaaaatgaatacccagaa'.
  'atcagttgcatttctatatatgaacaatgaacaacctgaaagggaaatttaaaaaaactattccatttacaatacca'.
  'ttaaaaagaatagaatacataatactaaacttaagaaacttaagatctgtacactgaaaactataaaacattaccaa'.
  'gaaaaatttaaaaagacctaaatcagtggaaagatatcccatgttcatggctgggaagatttaccactgtcacaatg'.
  'acgatactactgaaagggatctacagattcagtgcaattcctattaaaatcccaatgacagccaggcatagtagcat'.
  'gtatctgtaatcctagctagtcagaggctgaggctaaagaactgcttgaacccaggggttcaaggccgcagtaagct'.
  'aggatcacaccattgcattccagcctagatgacagaacaagacctgtcaaaaaaagaaaagaaaagaaaaaaaaaat'.
  'cccaatgacattgtttgcagaaacagagaaatctatcataagatttatatcgatctttaaggagcccaaagagccaa'.
  'aacaatctggaaaaagaacaacatttaaggactcacacttcctgatttcaaaacttactacaaagctacagtaatca'.
  'aaacagtgtggtattagcataaggacagacacaaagaccaatggaatagaatacagagcccagaaacaaaccctcat'.
  'gaatacggtcaaatgaccttcaataaaggtgccaagaccattcattagggaaatgatggtattttcaacaaatgatg'.
  'gaactggatacacacacgtaaatgaatgaactcacacccttaccttattacatacacaaaaattaactcaaaataga'.
  'tgaaagatctaaacataagacctcaaactataaaactcttagaagaaaatataggggaaagtcttcatgacactgga'.
  'tttggcaatgatatcttagataagataccaaaaacacaggcaacagaaggaaaaataaattgtactacatccaaatt'.
  'taaaatgtatgttcatcaaaggacacaatcaatatagtgaaaaagcaaccaacagaatgggagaaaagatttgtaaa'.
  'tcatatatctgaaaaggaattcatatccagaatatgtaaagaactcccaaaacttaacaacagcaaaaaaacacaaa'.
  'taacccaatttaaaaatgagtggagttgaactgacatttctccaaagaagatatataaacggccaacaagcacgtta'.
  'aaaaatttgctcaacatttctatcattagggaagagcaaatcaaggctgcaatgaaatatcaccttacacctattag'.
  'gatggctagtctaaaaaaaaaaaaaaagttttggcagagatacagacaaattaagacccttgtgcactgctggtggg'.
  'aatgtaagatcagtatggcaattcctcaaaaaattaaaaatagaattaccatatggtccagcaatcagaaaattctg'.
  'cattatgtatctaaaagaattgaggccaggcgtggtggctcacgcctgtaatcccaacactttgggaggccgaggcg'.
  'ggtggatcacaaggtcaggagatcgagaccatcctggctagcacggtgaaaccccgtctctactaaaaacacaaaaa'.
  'aattagccgggcgtggtgacaggcgcctgtagtcccagctacttgtgaacttgggaggcggagcttgcagtgagcca'.
  'agatcgtgccactgcactccagcctgggcggcagagcgagactccgtctaaaaaaaaaaaaaaaaagaattgaaagc'.
  'agagtctccaagagatatgtgcacaatcatattcacagcagcattattcacaatagccaaaaaggggaagcaaccca'.
  'agtgtccaccagtggaaaaatggataaacaaaatgtggtatagatatacaatggaatatcaggcatcaaattctgac'.
  'atatgatattctcatgatactaagtgaaataagccagtcacaaaaagacagatactgtgagattctatttatgtaag'.
  'gcatctaaagtagtcaaacaataggaaaagtagaagcacttcaagaaaggggcaaacggggagttttctaatgggta'.
  'cagagtttctattttgcaagaataaaaaaaagttctgaaattggttgcacagcaatgtgaacaaacaacactactga'.
  'attgtacacttgaaaatggttaagatagtaaattttgttatgtgtattctaccataattaaaaattaaaattaaaaa'.
  'agaggctccctggccaggcgtggtggctcacgcctgtaatctcagcactttgggaggccaaggcgggcggatcacga'.
  'ggtcaggagatcgagaccatcctggctaacatggtgaaaccccgtctctactaaaagtacaaaaaaaattagccggg'.
  'tgtggcggtgagcgcctgtagttccagctacttgggaggctgaggcaggagaatggtgtgaacccaggaggcggagc'.
  'ttgcagtgagccgagatcgcgccactgcactccagcctgggcggtagagcgtgactctgtctcaaaaaaaaaaagag'.
  'gctccctcacataataagtaagaaattgcataatacccttttatttatttattttttgagacggagtcttgccctgt'.
  'cgccaggctagagtgcagtgg';

return $seq;
}
