use lib 't';
use Test;
use strict;

BEGIN { $| = 1; plan test => 19;
      }

use EnsTestDB;
use Bio::EnsEMBL::Pipeline::RunnableDB::EPCR;
use Bio::EnsEMBL::Analysis;

ok(1);

ok(my $ens_test = EnsTestDB->new);

ok($ens_test->do_sql_file("t/dumps/runnabledb.dump"));
    
ok(my $db = $ens_test->get_DBSQL_Obj);

ok(my $runnable = 'Bio::EnsEMBL::Pipeline::RunnableDB::EPCR');

ok(my $sts_db   = $::pipeConf{'datadir'} . "/mapprimer_maps_110");
  
ok(my $ana_adaptor = $db->get_AnalysisAdaptor);

ok(my $ana = $ana_adaptor->fetch_by_logic_name('e-PCR'));

ok($ana->db_file($sts_db));

my $id = 'AB015752.00001'; 

ok(my $runobj = "$runnable"->new(  -db       => $db,
				   -input_id => $id,
				   -analysis => $ana ));

ok($runobj->fetch_input);

ok($runobj->run);

ok(my @out = $runobj->output);

ok(display(@out));

ok($runobj->write_output);

ok(my $contig =  $db->get_RawContigAdaptor->fetch_by_name($id));

ok(my $mfa = $db->get_MarkerFeatureAdaptor);

ok(my @features = @{$mfa->fetch_all_by_RawContig_and_priority($contig)});

ok(display(@features));

sub display {
  my @results = @_;
 
  foreach my $obj (@results) {
    print ($obj->gffstring."\n");

    if ($obj->sub_SeqFeature) {

      foreach my $exon ($obj->sub_SeqFeature) {
	print "Sub: ".$exon->gffstring."\n";
      }

    }
  }
  return 1;
}









