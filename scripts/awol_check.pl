#!/usr/local/ensembl/bin/perl -w

use strict;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Config::BatchQueue;

my $dbhost    = $ENV{'ENS_DBHOST'};
my $dbname    = $ENV{'ENS_DBNAME'};
my $dbuser    = $ENV{'ENS_DBUSER'};
my $dbpass    = $ENV{'ENS_DBPASS'};
my $dbport    = $ENV{'ENS_DBPORT'} || 3306;
my $ignore_lock;
my $queue_manager;
my $verbose;
my $help;
my @command_args = @ARGV;
GetOptions(
           'dbhost=s'      => \$dbhost,
           'dbname=s'      => \$dbname,
           'dbuser=s'      => \$dbuser,
           'dbpass=s'      => \$dbpass,
           'dbport=s'      => \$dbport,
           'verbose!' => \$verbose,
           'ignore_lock!' => \$ignore_lock,
           'queue_manager=s' => \$queue_manager,
           'help!' => \$help,
          ) or useage(\@command_args);


if(!$dbhost || !$dbname || !$dbuser){
  $help = 1;
}

if($help){
  useage(\@command_args);
}
my $db = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new(
                                                       -host   => $dbhost,
                                                       -dbname => $dbname,
                                                       -user   => $dbuser,
                                                       -pass   => $dbpass,
                                                       -port   => $dbport,
                                                      );


if(!$ignore_lock){
  if (my $lock_str = $db->pipeline_lock) {
    my($user, $host, $pid, $started) = 
      $lock_str =~ /(\w+)@(\w+):(\d+):(\d+)/;
    $started = scalar localtime $started;
    
    print STDERR ("Error: There is a RuleManager already running\n. ".
                  "You probably don't want to run this script then\n\n".
                  "      db       $dbname\@$dbhost\n".
                  "      pid      $pid on host $host\n".
                  "      started  $started\n\n".
                  "If the process doesn't exist the lock can be deleted\n".
                  "using this sql\n\n".
                  "      delete from meta where meta_key = ".
                  "'pipeline.lock';\n\n".
                  "If the process does exist are you want to run this ".
                  "script anyway use the -ignore_lock commandline ".
                  "option\n\n".
                  "Thank you\n");
    exit 1;
  }
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
    print STDERR "$@\nError trying to load $batch_q_module;\ncan't find $file\n";
    exit 1;
}

if($batch_q_module->can('check_existance')){
&job_existance($batch_q_module, $verbose, $db->get_JobAdaptor);
}else{
  print STDERR "Your batch_q_module $batch_q_module need to have the ".
    "method check_existance implemented for this script to work";
  exit(1);
}

sub job_existance{
  my ($batch_q_module, $verbose, $job_adaptor) = @_;
  my @jobs = $job_adaptor->fetch_all;
  my @valid_status = qw(SUBMITTED RUNNING READING WRITING WAITING);
  my %valid_status = map { $_ => 1 } @valid_status;
  my %job_submission_ids;
  JOB:foreach my $job(@jobs){
    my $status = $job->current_status->status;
    if($valid_status{$status}){
      if(!$job->submission_id){
        if($status eq 'SUBMITTED'){
          next JOB;
        }else{
          print STDERR ("Job ".$job->dbID." doesn't have a submission id ".
                        "We can't check if it still exists and the ".
                        "status is aparently ".$status."\n");
        }
      }
      if(!$job_submission_ids{$job->submission_id}){
        $job_submission_ids{$job->submission_id} = [];
      }
      push(@{$job_submission_ids{$job->submission_id}}, $job);
    }
  }
  
  my @awol_jobs = @{$batch_q_module->check_existance
                 (\%job_submission_ids, $verbose)};
  
  foreach my $job(@awol_jobs){
    if($valid_status{$job->current_status->status}){
      $job->set_status('AWOL');
    }
  }
}

sub useage{
  my ($command_args) = @_;
  print "Your commandline was :\n".
    "awol_check.pl ".join("\t", @$command_args), "\n\n";
	exec('perldoc', $0);
	exit;
}

=pod

=head1 NAME

awol_check.pl

=head1 SYNOPSIS

awol_check.pl, a script to check if jobs have gone missing from
your batch_submission system

=head1 OPTIONS

DB Connection Details


   -dbhost     The host where the pipeline database is.
   -dbport     The port.
   -dbuser     The user to connect as.
   -dbpass     The password to use.
   -dbname     The database name.


Other Useful Options

  -ignore_lock, will ignore any lock strings already present in the
  pipeline database
  -verbose will print out information about the checks being carried out
  -queue_manager specifies which Bio::EnsEMBL::Pipeline::BatchSubmission
   module to use


=head1 EXAMPLES

./awol_check -dbhost ecs2b -dbuser ensadmin -dbpass ****
 -dbname pipeline_db


this will get all the jobs out of the database. Any jobs which have
valid statuses, SUBMITTED, WAITING, READING, WRITING, RUNNING will have
their submission id checked to make sure they still exist

those which no longer exist will be marked as awol
