#!/usr/local/ensembl/bin/perl 

# default pipeline runner script
# this script is passed as part of the batch submission request
# and is therefore the first script to be executed on the
# remote side

# note that you may need to alter the #! line above to suit your
# local perl installation

use strict;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Config::General;
use Bio::EnsEMBL::Pipeline::Config::BatchQueue;
use Sys::Hostname;
use Getopt::Long;

#parameters for Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor

my $dbname;
my $host;
my $dbuser;
my $port            = '3306';
my $pass            = undef;
my $job_id;
my $queue_manager;
my $cleanup = 0;
my $check;
my $output_dir;
my $pants;

GetOptions(
           'dbhost=s'        => \$host,
           'dbport=n'        => \$port,
           'dbname=s'        => \$dbname,
           'dbuser=s'        => \$dbuser,
           'dbpass=s'        => \$pass,
           'check!'          => \$check,
           'output_dir=s'    => \$output_dir,
           'queue_manager=s' => \$queue_manager,
           'cleanup!'        => \$cleanup,
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

$queue_manager = $QUEUE_MANAGER unless($queue_manager);

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
    #print "SETTING STATUS to WAITING\n";
    $job->set_status('WAITING');
    push(@jobs, $job);
}
if($submission_id_not){
  $job_adaptor->update(@jobs);
}

my $hostname = [ split(/\./, hostname()) ];
$host = shift(@$hostname);


#if($batch_q_object->can('delete_output')){
#  $batch_q_object->delete_output();
#}


if($cleanup){
  &run_jobs_with_selfcopy(\@jobs, $host);
}else{
  &run_jobs_with_lsfcopy(\@jobs, $host);
}

sub run_jobs_with_selfcopy{
  my ($jobs, $host) = @_;
 
  # fix needed here; when running a batch of jobs, a jobs failure causes
  # all output so far to be copied to the output directory. Upon completion
  # of the whole batch, LSF then appends output for the whole batch  to the 
  # end of the earlier-created file in the output directory. This causes 
  # duplication of all output up to the point of the job failure. If you 
  # intend to use for output for anything (e.g. data mining) then this will
  # obviously affect the results!

 JOB:foreach my $job(@$jobs) {
    
    my $job_id = $job->dbID;
    
    $job->execution_host($host);
    
    eval{
      $job_adaptor->update($job);
    };
    
    if($@){
      print STDERR "Job $job_id failed: [$@]";
      if($batch_q_object->can('copy_output')){
        $batch_q_object->copy_output($job->stderr_file, $job->stdout_file);
      }else{
        print STDERR "the batch_q_object ".$queue_manager.
          " needs to implement  the copy_output method\n";
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
        print STDERR "the batch_q_object ".$queue_manager." needs to implement ".
          " the copy_output method\n";
      }
      next JOB;
    }
    if ($job->current_status->status eq "SUCCESSFUL"){
      $job->adaptor->remove( $job );
    }
    print STDERR "Finished job $job_id\n";
  }
}

sub run_jobs_with_lsfcopy{
  my ($jobs, $host) = @_;

 JOB:foreach my $job(@$jobs) {
    
    my $job_id = $job->dbID;
    
    $job->execution_host($host);
    
    eval{
      $job_adaptor->update($job);
    };
    if($@){
      print STDERR "Job update ".$job->dbID." failed: [$@]";  
    }
    print STDERR "Running job $job_id\n";
    print STDERR "Module is " . $job->analysis->module . "\n";
    print STDERR "Input id is " . $job->input_id . "\n";
    print STDERR "Analysis is ".$job->analysis->logic_name."\n";
    print STDERR "Files are " . $job->stdout_file . " " .
      $job->stderr_file . "\n";
    
    eval {
      $job->run_module;
    };
    $pants = $@;
    
    if ($pants) {
      print STDERR "Job $job_id failed: [$pants]";
    }
    
    print STDERR "Finished job $job_id\n";
    if ($job->current_status->status eq "SUCCESSFUL"){
      $job->adaptor->remove( $job );
    }
  }
  
}
