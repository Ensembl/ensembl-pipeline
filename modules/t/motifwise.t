use strict;
use lib 't';
use Test;

BEGIN { $| = 1; plan test => 12;}

use Bio::EnsEMBL::Pipeline::Runnable::Motifwise;
use Bio::SeqIO;

ok(1);

ok(my $seqio = Bio::SeqIO->new( -file => "t/data/PKD_testdata.fa"));
ok(my $genomic_seq = $seqio->next_seq());


ok(my $mw = Bio::EnsEMBL::Pipeline::Runnable::Motifwise->new());

ok($mw->query_seq($genomic_seq));
ok($mw->parameters('-tfm_cutoff 11.0'));
ok($mw->workdir('/tmp'));
ok($mw->executable('motifwise'));
ok($mw->motif_file('/ecs2/work1/dta/motifwise/motif_dropout.lr'));

ok($mw->run);

ok(my @output = $mw->output);

ok(display(@output));

sub display {
    my @results = @_;
    #Display output
    foreach my $obj (@results) {
      print ($obj->gffstring."\n");
    }
    return 1
}
