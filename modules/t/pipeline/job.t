use lib 't/pipeline';
use strict;
use warnings;

BEGIN { $| = 1;
	use Test ;
	plan tests => 17;
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

ok(test_getter_setter($job, 'dbID', 1234));
ok(test_getter_setter($job, 'adaptor', $job_adaptor));
ok(test_getter_setter($job, 'taskname', 'test_task'));
ok(test_getter_setter($job, 'input_id', 'AL504321.1'));
ok(test_getter_setter($job, 'submission_id', 21421));
ok(test_getter_setter($job, 'job_name', 323));
ok(test_getter_setter($job, 'array_index', 1234));
ok(test_getter_setter($job, 'parameters', 'ensdb:3307:ecs1c'));
ok(test_getter_setter($job, 'module', 'Bio::EnsEMBL::Pipeline::RunnableDB::Test'));
ok(test_getter_setter($job, 'stderr_file', '/acari/nfs/tmp/12.out'));
ok(test_getter_setter($job, 'stdout_file', '/acari/nfs/tmp/1234.err'));
ok(test_getter_setter($job, 'status', 'WRITING'));
ok(test_getter_setter($job, 'retry_count', 1));


$multi_test_db->save('pipeline', 'job_status', 'job');

ok($job_adaptor->store($job));
ok($job->set_current_status('SUBMITTED'));

$multi_test_db->restore('pipeline');

#need a test for run but need more framework in place before its possible

