use lib 't/pipeline';
use strict;
use warnings;

BEGIN { $| = 1;
	use Test ;
	plan tests => 11;
}

use TestUtils qw(debug test_getter_setter);
use MultiTestDB;



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


ok($job_adaptor->store($job));

$status = $job_adaptor->fetch_current_status($job);

ok($status eq 'CREATED');


ok($job_adaptor->update_status($job, 'SUBMITTED'));

$status = $job_adaptor->fetch_current_status($job);

ok($status eq 'SUBMITTED');

my @lists = @{$job_adaptor->list_current_status};

ok(@lists);

#foreach my $l(@lists){
#  print STDERR join ', ', @$l;
#  print STDERR "\n";
#}
