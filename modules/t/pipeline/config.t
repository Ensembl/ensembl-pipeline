use lib 't/pipeline';
use strict;

BEGIN { $| = 1;
	use Test ;
	plan tests => 10;
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

eval { $config->get_parameter("abcdef", "0193734") };
ok($@);   # previous eval should have thrown exception, will be in $@

my @headers = $config->get_headers();
ok(@headers > 0);

my @keys = $config->get_keys("REPEATMASK_TASK");
ok(@keys > 0);

my $host = $config->get_parameter("PIPELINE_DATABASE", "host");
ok($host eq "127.0.0.1");

eval { my @badkeys = $config->get_keys("askjfwhjgfwruihtrueh") };
ok($@);

eval { my $config_broken = new Bio::EnsEMBL::Pipeline::Config(-FILES => "duplicate.conf") };
ok($@);

eval { new Bio::EnsEMBL::Pipeline::Config(-FILES => "test.conf", -DB    =>"some adaptor") };
ok($@);


