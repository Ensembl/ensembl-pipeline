use lib 't';
use Test;
use strict;

BEGIN { $| = 1; plan test => 7;}

use Bio::EnsEMBL::Pipeline::RunnableDB::GC;
use EnsTestDB;
use Bio::EnsEMBL::Pipeline::RunnableDB::Genscan;
use Bio::EnsEMBL::Analysis;

ok(1);
ok(my $ens_test = EnsTestDB->new);

ok($ens_test->do_sql_file("t/dumps/analysis.dump"));
ok($ens_test->do_sql_file("t/dumps/dna.dump"));
    
ok(my $db = $ens_test->get_DBSQL_Obj);

ok(my $runnable = 'Bio::EnsEMBL::Pipeline::RunnableDB::GC');

ok(my $ana_adaptor = $db->get_AnalysisAdaptor);

ok(my $ana = $ana_adaptor->fetch_by_logic_name('gc'));

my $id ='AC099340';

ok($ana_adaptor->exists( $ana ));

ok(my $runobj = "$runnable"->new(-db         => $db,
				 -input_id   => $id,
				 -analysis   => $ana ));


ok($runobj->fetch_input);

ok($runobj->run);

ok(my @out = $runobj->output);

ok(display(@out));

ok($runobj->write_output());

ok(my $contig = $db->get_RawContigAdaptor()->fetch_by_name($id));

ok (my @sf    = @{$contig->get_all_SimpleFeatures});

ok(display(@sf));


ok(1);

sub display {

  my (@f) = @_;

  foreach my $f (@f) {
    print $f->gffstring . "\n";
  }
  return 1;
}

