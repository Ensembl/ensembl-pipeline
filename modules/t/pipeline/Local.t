# Test the local submission system

use lib 't/pipeline';
use strict;

BEGIN { $| = 1;
	use Test ;
	plan tests => 1;
}

use MultiTestDB;
use TestUtils qw(debug test_getter_setter);
use Bio::EnsEMBL::Pipeline::Config;
use Bio::EnsEMBL::Pipeline::PipelineManager;

my $config = new Bio::EnsEMBL::Pipeline::Config(-FILE => "local_ss.conf");
ok($config);

my $pm = new Bio::EnsEMBL::Pipeline::PipelineManager($config);
ok($pm);

my $lss = Bio::EnsEMBL::Pipeline::SubmissionSystem::Local->new($pm);
ok($lss);

