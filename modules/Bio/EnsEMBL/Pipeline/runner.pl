#!/usr/local/ensembl/bin/perl

# default pipeline runner script
# this script is passed as part of the batch submission request
# and is therefore the first script to be executed on the
# remote side

# note that you may need to alter the #! line above to suit your
# local perl installation

use strict;

use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Sys::Hostname;
use Getopt::Long;


#parameters for Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor

my $host;
my $dbname;
my $dbuser;

my $port            = '3306';
my $pass            = undef;
my $jobname;
my $index;
my $check;
my $output_dir;

GetOptions(
    'host=s'       => \$host,
    'port=n'       => \$port,
    'dbname=s'     => \$dbname,
    'dbuser=s'     => \$dbuser,
    'pass=s'       => \$pass,
    'index'        => \$index,
    'jobname=s'      => \$jobname,
    'check!'       => \$check,
    'output_dir=s' => \$output_dir
)
or die ("Couldn't get options");

if( $check ) {
  my $host = hostname();
  if ( ! -e $output_dir ) {
    die "output dir $output_dir doesn't exist according to host $host";
  }

  exit 0;

}

print STDERR "have jobname ".$jobname."\n";
$jobname or die("jobname not defined");

my $db = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new(
    -host   => $host,
    -user   => $dbuser,
    -dbname => $dbname,
    -pass   => $pass,
    -port   => $port);


my $job_adaptor = $db->get_JobAdaptor();

$index ||= '';
print STDERR "Fetching job [$jobname] index [$index]\n";
my $jobs;

if($index) {
  my $array_index = $ENV{LSB_JOBINDEX};
  $jobs = $job_adaptor->fetch_all_by_job_name($jobname, $array_index);
} else {
  $jobs = $job_adaptor->fetch_all_by_job_name($jobname);
}

if(@$jobs > 1) {
  die "expected single job with [$jobname] index [$index] but got [" .
    scalar(@$jobs) . "] instead";
}

my ($job) = @$jobs;

die "Couldn't recreate job [$jobname] index [$index]\n" if(!$job);

print STDERR "Running job [$jobname] index [$index]\n";
print STDERR "Module is [".$job->module."]\n";
print STDERR "Input id is [".$job->input_id."]\n";
print STDERR "Files are [".$job->stdout_file."] [".$job->stderr_file."]\n";

eval {
  $job->run;
};

if ($@) {
  print STDERR "Job [$jobname] index [$index] failed: [$@]\n";
}

print STDERR "Finished job [$jobname] index [$index]\n";
