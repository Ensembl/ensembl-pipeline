#!/usr/local/ensembl/bin/perl

# default pipeline runner script
# this script is passed as part of the batch submission request
# and is therefore the first script to be executed on the
# remote side

# note that you may need to alter the #! line above to suit your
# local perl installation


use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Config::General;
use Bio::EnsEMBL::Pipeline::Config::BatchQueue;
use Sys::Hostname;
use Getopt::Long;


#parameters for Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor

my $host;
my $dbname;
my $dbuser;

my $port            = '3306';
my $pass            = undef;
my $job_id;
my $queue_manager;
my $cleanup = 0;

GetOptions(
           'dbhost=s'        => \$host,
           'dbport=n'        => \$port,
           'dbname=s'        => \$dbname,
           'dbuser=s'        => \$dbuser,
           'dbpass=s'        => \$pass,
           'check!'          => \$check,
           'output_dir=s'    => \$output_dir,
           'queue_manager=s' => \$queue_manger,
           'cleanup!'       => \$cleanup,
)
or die ("Couldn't get options");

if( $check ) {
  my $host = hostname();
  if ( ! -e $output_dir ) {
    die "output dir $output_dir doesn't exist according to host $host";
  }
  my $deadhostfile = "$output_dir/deadhosts";
  open( FILE, $deadhostfile ) or exit 0;
  while( <FILE> ) {
    chomp;
    if( $host eq $_ ) {
      die "Cant use this host";
    }
  }
  exit 0;

}
my $queue_manager = $QUEUE_MANAGER unless($queue_manager);
my $batch_q_module = 
  "Bio::EnsEMBL::Pipeline::BatchSubmission::$queue_manager";




my $file = "$batch_q_module.pm";
$file =~ s{::}{/}g;
eval {
    require "$file";
};
if ($@) {
    print STDERR "Error trying to load $batch_q_module;\ncan't find $file\n";
    exit 1;
}

my $batch_q_object = $batch_q_module->new();
my $db = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new(
    -host   => $host,
    -user   => $dbuser,
    -dbname => $dbname,
    -pass   => $pass,
    -port   => $port
)
or die ("Failed to create Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor to db $dbname \n");

my $job_adaptor = $db->get_JobAdaptor();

# Print a list of jobs in this 'batch' - this is perhaps useful later
# to parse these output files for CPU/Mem usage etc.
print STDOUT "LSF Batch summary\n";
print STDOUT "Time ", scalar localtime time, " (", time, ")\n";
print STDOUT "Job ID\tinput ID\tanalysis ID\n";
my @jobs;
my $submission_id_not = 0;
foreach my $id (@ARGV) {
    my $job = $job_adaptor->fetch_by_dbID($id);
    print STDERR "lost job ".$id." while trying to run\n" if(!$job);
    exit(0) if(!$job);
    print STDOUT join ("\t", $id, $job->input_id, $job->analysis->logic_name), "\n";
    if(!$job->submission_id){
      $submission_id_not = 1;
      $job->submission_id($ENV{'LSB_JOBID'});
    }
    $job->set_status('WAITING');
    push(@jobs, $job);
}
if($submission_id_not){
  $job_adaptor->update(@jobs);
}

my $hostname = [ split(/\./, hostname()) ];
my $host = shift(@$hostname);
JOB:foreach my $job(@jobs) {

  #print STDERR "Fetching job " . $job_id . "\n";
  
  #my $job = $job_adaptor->fetch_by_dbID($job_id);
  $job->execution_host($host);
  eval{
    $job_adaptor->update($job);
  };

  if($@){
    print STDERR "Job $job_id failed: [$@]";
    if($batch_q_object->can('copy_output')){
      $batch_q_object->copy_output($job->stderr_file, $job->stdout_file);
    }else{
      print STDERR "the batch_q_object ".$queue_manger." needs to implement ".
        " the copy_output method\n";
    }
    #if($batch_q_object->can('delete_output')){
    #  $batch_q_object->delete_output() unless(!$cleanup);
    #}
  }
  if( !$job) {
    print STDERR ( "Couldnt recreate job $job_id\n" );
    next;
  }
  print STDERR "Running job $job_id\n";
  print STDERR "Module is " . $job->analysis->module . "\n";
  print STDERR "Input id is " . $job->input_id . "\n";
  print STDERR "Analysis is ".$job->analysis->logic_name."\n";
  print STDERR "Files are " . $job->stdout_file . " " . $job->stderr_file . "\n";

  eval {
    $job->run_module;
  };
  $pants = $@;
  print STDERR "have batchqueue module ".$batch_q_object."\n";
  if ($pants) {
    print STDERR "Job $job_id failed: [$pants]";
    if($batch_q_object->can('copy_output')){
      $batch_q_object->copy_output($job->stderr_file, $job->stdout_file);
    }else{
      print STDERR "the batch_q_object ".$queue_manger." needs to implement ".
        " the copy_output method\n";
    }
    next JOB;
  }
  if(!$cleanup){
    if($batch_q_object->can('copy_output')){
      $batch_q_object->copy_output($job->stderr_file, $job->stdout_file);
    }else{
      print STDERR "the batch_q_object ".$queue_manger." needs to implement ".
        " the copy_output method\n";
    }
  }
  if ($job->current_status->status eq "SUCCESSFUL"){
    $job->adaptor->remove( $job );
  }
  print STDERR "Finished job $job_id\n";
}

#if($batch_q_object->can('delete_output')){
#  $batch_q_object->delete_output();
#}


