use lib 't/pipeline';
use strict;
use warnings;

BEGIN { $| = 1;
	use Test ;
	plan tests => 5;
}

use TestUtils qw(debug test_getter_setter);
use MultiTestDB;
use Bio::EnsEMBL::Pipeline::Job;


my $multi_test_db = new MultiTestDB();

my $db_adaptor = $multi_test_db->get_DBAdaptor('pipeline');

ok($db_adaptor);

my $job_adaptor = $db_adaptor->get_JobAdaptor;

ok($job_adaptor);

my $job = Bio::EnsEMBL::Pipeline::Job->new(
					   -input_id => 'c011601284.1.6205',
					   -taskname => 'RepeatMasker',
					   -module => 'Bio::EnsEMBL::Pipeline::RunnableDB::RepeatMasker',
					  );




ok($job);


ok($job_adaptor->store($job));


ok($job->set_current_status('SUBMITTED'));

#need a test for run but need more framework in place before its possible

