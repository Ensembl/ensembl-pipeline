use lib 't';
use strict;
use Test;

BEGIN { $| = 1; plan test => 16; }

use EnsTestDB;
use Bio::EnsEMBL::Pipeline::RunnableDB::Genscan;
use Bio::EnsEMBL::Analysis;

ok(1);


ok(my $ens_test = EnsTestDB->new);

ok($ens_test->do_sql_file("t/dumps/runnabledb.dump"));
    
ok(my $db = $ens_test->get_DBSQL_Obj);

ok(my $runnable = 'Bio::EnsEMBL::Pipeline::RunnableDB::Genscan');

ok(my $ana_adaptor = $db->get_AnalysisAdaptor);

ok(my $ana = $ana_adaptor->fetch_by_logic_name('genscan'));

my $id ='AB015752.00001';

ok($ana_adaptor->exists( $ana ));

ok(my $runobj = "$runnable"->new(-db      => $db,
				 -input_id   => $id,
				 -analysis   => $ana ));


ok($runobj->fetch_input);

ok($runobj->run);

ok(my @out = $runobj->output);

ok(display(@out));

ok($runobj->write_output());

ok(my $contig = $db->get_RawContigAdaptor()->fetch_by_name($id));

ok(my $prediction_transcripts = $contig->get_all_PredictionTranscripts());

foreach my $pred_trans (@{$prediction_transcripts}) {
  print STDERR "Transcript  " . $pred_trans->translate . "\n";
}    

sub display {
  my @results = @_;

  foreach my $obj (@results) {
    if ($obj->get_all_Exons) {
      foreach my $exon (@{$obj->get_all_Exons}) {
	print STDERR "Trancr/Exon: ".$exon->gffstring."\n";
      }
    }
  }
  return 1;
}
