use lib 't';
use Test;
use strict;

BEGIN { $| = 1; plan test => 5;}

use Bio::EnsEMBL::Pipeline::Runnable::NewExonerate;
use Bio::SeqIO;

ok(1);

ok(my $estseqio = new Bio::SeqIO(-file   => 't/data/testest.fa', 
                                 -format => 'fasta'));

my @estseqs;

while (my $seq = $estseqio->next_seq) {
  push(@estseqs,$seq);
}

ok(my $exonerate = Bio::EnsEMBL::Pipeline::Runnable::NewExonerate->new(
   -query_seqs  => \@estseqs,
   -query_type  => 'DNA',
   -database    => 't/data/AC099340.fa.masked', #$genseq,
   -target_type => 'DNA',
   -exonerate   => '/usr/local/ensembl/bin/exonerate-0.6.7',
#   -exonerate   => '/usr/local/ensembl/bin/exonerate-0.8.2',
#   -options     => ' --dnahspthreshold 30 --forcegtag TRUE --maxintron 500000'
   ));

$exonerate->_verbose(0);

ok($exonerate->run);

ok(my @transcripts = $exonerate->output());

ok(display(@transcripts));

sub display {
  my @transcripts = @_;

  foreach my $transcript (@transcripts) {
    my $exons = $transcript->get_all_Exons;
    my $dafs = $exons->[0]->get_all_supporting_features;

    print "Match : (hseqname) " . 
    $dafs->[0]->hseqname . "\t" .
    "(length) " . $transcript->length . "\t" . 
    "(exons) " . scalar @$exons . "\t" . 
    "(strand) " . $transcript->strand . "\t" . 
     "\n";
  }

  return 1
}

