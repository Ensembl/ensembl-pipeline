use lib 't';
use strict;
use Test;

BEGIN { $| = 1; plan test => 6; }

use EnsTestDB;
use Bio::EnsEMBL::Pipeline::RunnableDB::EponineTSS;
use Bio::EnsEMBL::Analysis;

ok(1);

ok(my $ens_test = EnsTestDB->new());
ok($ens_test->do_sql_file("t/runnabledb.dump"));
ok(my $db = $ens_test->get_DBSQL_Obj());

ok(my $runnable = 'Bio::EnsEMBL::Pipeline::RunnableDB::EponineTSS');
ok(my $ana_adaptor = $db->get_AnalysisAdaptor());
ok(my $ana = $ana_adaptor->fetch_by_logic_name('Eponine'));

my $id =  'AL009179.00001';  

$ana_adaptor->exists( $ana );

my $runobj = "$runnable"->new(  -db         => $db,
																-input_id   => $id,
                                -analysis   => $ana );
ok($runobj);

$runobj->fetch_input();

ok(1);

$runobj->run();

ok(1);

ok(my @out = $runobj->output());

$runobj->write_output();

ok(1);

ok(my $contig =  $db->get_RawContigAdaptor->fetch_by_name($id));
ok(my @features = @{$db->get_SimpleFeatureAdaptor->fetch_all_by_Contig($contig, 'Eponine')});

foreach my $obj (@features) {
	print ($obj->gffstring."\n");
}

ok(1);
