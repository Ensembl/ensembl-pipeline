use lib 't';
use Test;
use strict;

BEGIN {$| = 1; plan test => 15;}

use EnsTestDB;
use Bio::EnsEMBL::Pipeline::RunnableDB::tRNAscan_SE;
use Bio::EnsEMBL::Analysis;

ok(1);
    
ok(my $ens_test = EnsTestDB->new);

ok($ens_test->do_sql_file("t/dumps/runnabledb.dump"));
    
ok(my $db = $ens_test->get_DBSQL_Obj);

ok(my $runnable = 'Bio::EnsEMBL::Pipeline::RunnableDB::tRNAscan_SE');

ok(my $ana_adaptor = $db->get_AnalysisAdaptor);

ok(my $ana = $ana_adaptor->fetch_by_logic_name('tRNAscan'));

ok($ana_adaptor->exists( $ana ));

my $id = 'AL009179.00001';

ok(my $runobj = "$runnable"->new(  -db         => $db,
				   -input_id   => $id,
				   -analysis   => $ana ));

ok($runobj->fetch_input());

ok($runobj->run());

ok(my @out = $runobj->output());

foreach my $obj (@out) {
  print ($obj->gffstring."\n");
}

ok($runobj->write_output());

ok(my $contig =  $db->get_RawContigAdaptor->fetch_by_name($id));

ok(my @features = @{$db->get_SimpleFeatureAdaptor->fetch_all_by_RawContig($contig)});

foreach my $obj (@features) {
  print ($obj->gffstring."\n");
}
