use lib 't';
use Test;
use strict;


BEGIN { $| = 1; plan test => 16;}

use EnsTestDB;
use Bio::EnsEMBL::Pipeline::RunnableDB::TRF;
use Bio::EnsEMBL::Analysis;

ok(1);

ok(my $ens_test = EnsTestDB->new);

ok($ens_test->do_sql_file("t/dumps/runnabledb.dump"));
    
ok(my $db          = $ens_test->get_DBSQL_Obj);

ok(my $runnable    = 'Bio::EnsEMBL::Pipeline::RunnableDB::TRF');

ok(my $ana_adaptor = $db->get_AnalysisAdaptor);

ok(my $ana         = $ana_adaptor->fetch_by_logic_name('trf'));

my $id =  'AL009179.00001';  

ok($ana_adaptor->exists( $ana ));

ok(my $runobj = "$runnable"->new(  -db      => $db,
				   -input_id   => $id,
				   -analysis   => $ana ));

ok($runobj->fetch_input);

ok($runobj->run);

ok(my @out = $runobj->output);

ok($runobj->write_output);

ok(my $contig =  $db->get_RawContigAdaptor->fetch_by_name($id));

ok(my @features =@{$contig->get_all_RepeatFeatures('TRF')});

ok(display(@features));

sub display {
    my @results = @_;

    foreach my $obj (@results) {
       print ($obj->gffstring."\n");
    }
    1;
}
