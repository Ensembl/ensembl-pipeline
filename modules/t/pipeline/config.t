use lib 't/pipeline';
use strict;

BEGIN { $| = 1;
	use Test ;
	plan tests => 3;
}

use MultiTestDB;
use TestUtils qw(debug test_getter_setter);

use Bio::EnsEMBL::Pipeline::Config;

my $multi_test_db = new MultiTestDB();

my $db_adaptor = $multi_test_db->get_DBAdaptor('pipeline');

ok($db_adaptor);

my $config = new Bio::EnsEMBL::Pipeline::Config(-DB => $db_adaptor);

ok($config);

ok($config->get_parameter("DEFAULT", "default_key1") eq "value1");

