use lib 't';
use strict;
use Test;

BEGIN { $| = 1; plan test => 16;}

use EnsTestDB;
use Bio::EnsEMBL::Pipeline::RunnableDB::CPG;
use Bio::EnsEMBL::Analysis;

ok(1);

ok(my $ens_test = EnsTestDB->new());
ok($ens_test->do_sql_file("t/runnabledb.dump"));
    
ok(my $db = $ens_test->get_DBSQL_Obj());

ok(my $runnable = 'Bio::EnsEMBL::Pipeline::RunnableDB::CPG');
ok(my $ana_adaptor = $db->get_AnalysisAdaptor());
ok(my $ana = $ana_adaptor->fetch_by_logic_name('cpg'));

my $id ='AB016897.00001'; 

ok($ana_adaptor->exists( $ana ));

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

ok(my $contig   =  $db->get_RawContigAdaptor->fetch_by_name($id));
ok(my @features = @{$db->get_SimpleFeatureAdaptor->fetch_all_by_RawContig($contig, 'cpg')});

display(@features);

ok(1);

sub display {
    my @results = @_;
    #Display output
    foreach my $obj (@results)
    {
       print ($obj->gffstring."\n");
    }
}
