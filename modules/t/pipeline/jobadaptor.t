use lib 't/pipeline';
use strict;
use warnings;

BEGIN { $| = 1;
	use Test ;
	plan tests => 13;
}

use TestUtils qw(debug test_getter_setter);
use MultiTestDB;

our $verbose = 0;

my $multi_test_db = new MultiTestDB();

my $db_adaptor = $multi_test_db->get_DBAdaptor('pipeline');

ok($db_adaptor);

my $job_adaptor = $db_adaptor->get_JobAdaptor;

ok($job_adaptor);

my $job = $job_adaptor->fetch_by_dbID(1);

ok($job);

my $status = $job_adaptor->fetch_current_status($job);
ok($status eq 'FAILED');

my @jobs = @{$job_adaptor->fetch_all_by_taskname('RepeatMasker')};

ok(@jobs == 2);


$job = Bio::EnsEMBL::Pipeline::Job->new(
					-input_id => 'c011601284.1.6205',
					-taskname => 'RepeatMasker',
					-module => 'Bio::EnsEMBL::Pipeline::RunnableDB::RepeatMasker',
				       );




ok($job);

$multi_test_db->save('pipeline', 'job', 'job_status');

ok($job_adaptor->store($job));

$status = $job_adaptor->fetch_current_status($job);

ok($status eq 'CREATED');


ok($job_adaptor->update_status($job, 'SUBMITTED'));

$status = $job_adaptor->fetch_current_status($job);

ok($status eq 'SUBMITTED');

my @lists = @{$job_adaptor->list_current_status};

ok(@lists);

if($verbose) {
  foreach my $l(@lists){
    print STDERR join ', ', @$l;
    print STDERR "\n";
  }
}

#
# Test that the update method works
#
$job->job_name('new_test');
$job->submission_id('test');

$job_adaptor->update($job);

$job = $job_adaptor->fetch_by_dbID($job->dbID);
ok($job->job_name eq 'new_test');
ok($job->submission_id eq 'test');


$multi_test_db->restore('pipeline');


