use lib 't';
use strict;
use Test;

BEGIN { $| = 1; plan test => 21;}

use EnsTestDB;

use Bio::EnsEMBL::Pipeline::RunnableDB::Motifwise;
use Bio::EnsEMBL::Analysis;

ok(1);

### Set up a test database that will be deleted.

ok(my $ens_test = EnsTestDB->new);

ok($ens_test->do_sql_file("t/dumps/chr12_snippet.dump"));

ok(my $db = $ens_test->get_DBSQL_Obj);


### Make runnable object and other objects to be passed to runnable.


ok(my $runnable = 'Bio::EnsEMBL::Pipeline::RunnableDB::Motifwise');

ok(my $ana_adaptor = $db->get_AnalysisAdaptor);

ok(my $ana_tss   = $ana_adaptor->fetch_by_logic_name('motifwise_tss'));
ok(my $ana_motif = $ana_adaptor->fetch_by_logic_name('motifwise_motif'));

my $id = '12.63894423-63914407';
my $chr = '12';
my $start = 63894423;
my $end = 63914407;

ok($ana_adaptor->exists( $ana_tss ));
ok($ana_adaptor->exists( $ana_motif ));

ok(my $runobj = "$runnable"->new(-db             => $db,
				 -input_id       => $id,
				 -analysis_tss   => $ana_tss,
				 -analysis_motif => $ana_motif));

### Use the runnable.

ok($runobj->fetch_input);

ok($runobj->run);

ok(my @out = $runobj->output());

ok((scalar @out) > 0);

ok($runobj->write_output);

### Check that the runnable is writing the correct things to the database.

ok(my $slice   =  $db->get_SliceAdaptor->fetch_by_chr_start_end($chr, $start, $end));

ok(my @tss_features = @{$db->get_SimpleFeatureAdaptor->fetch_all_by_Slice($slice, 'motifwise_tss')});
ok(my @motif_features = @{$db->get_SimpleFeatureAdaptor->fetch_all_by_Slice($slice, 'motifwise_motif')});

ok(display(@tss_features));
ok(display(@motif_features));

sub display {
    my @results = @_;
    #Display output
    foreach my $obj (@results) {
      print ($obj->gffstring."\n");
    }
    return 1;
}
