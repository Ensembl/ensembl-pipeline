use lib 't';
use strict;
use Test;

BEGIN { $| = 1; plan test => 9;} 

use EnsTestDB;
use Bio::EnsEMBL::Pipeline::RunnableDB::Slice_EPCR;

ok(1);

my $ens_test = EnsTestDB->new();
# Load some data into the db
$ens_test->do_sql_file("t/dumps/chr12_snippet.dump");
$ens_test->do_sql_file("t/dumps/slice_epcr.dump");
    
# Get an EnsEMBL db object for the test db
my $db = $ens_test->get_DBSQL_Obj;
ok(1);

my $ana_adaptor = $db->get_AnalysisAdaptor;
my $ana = $ana_adaptor->fetch_by_logic_name('EPCR');

ok($ana);

my($chr, $chr_start, $chr_end) = ('12', 63894423, 63914407);

ok(my $epcr = Bio::EnsEMBL::Pipeline::RunnableDB::Slice_EPCR->new(
    -db       => $db,
    -analysis => $ana,
    -input_id => "$chr.$chr_start.$chr_end",
));

$epcr->fetch_input;
ok(1);

$epcr->run;
ok(1);

$epcr->write_output();
ok(1);

my $sfa   = $db->get_SimpleFeatureAdaptor;
my $slice = $db->get_SliceAdaptor->fetch_by_chr_start_end($chr, $chr_start, $chr_end);
my @features = @{$sfa->fetch_all_by_Slice($slice, 'EPCR')};
ok(1);

my @methods = qw(seqname start end strand score);

foreach my $window (@features) {
    print STDERR "\n";
    foreach my $method_name (@methods) {
        my $value = $window->$method_name();
        printf STDERR "%10s = $value\n", $method_name;
    }
}
ok(1);

