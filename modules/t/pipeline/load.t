use lib 't/pipeline';
use strict;

use Bio::EnsEMBL::Pipeline::Load;


BEGIN { $| = 1;
	use Test ;
	plan tests => 2;
}


my $load = new Bio::EnsEMBL::Pipeline::Load();
ok($load);

my @uptimes = $load->get();



