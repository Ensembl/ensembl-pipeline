package Bio::EnsEMBL::Pipeline::RuleManager;


=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RuleManager

=head1 SYNOPSIS


=head1 DESCRIPTION

The RuleManager object is to provide functionailty for creating Jobs, checking if they can run and checking their status and existence

=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


use vars qw(@ISA);
use strict;
use warnings;
use File::Copy;
use Sys::Hostname;
use Socket;

use Bio::EnsEMBL::Pipeline::Config::General;
use Bio::EnsEMBL::Pipeline::Config::BatchQueue;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Pipeline::Job;

@ISA = qw();


=head2 new

  Arg [1]   : Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
  Function  : create a Bio::EnsEMBL::Pipeline::RuleManager object 
  Returntype: Bio::EnsEMBL::Pipeline::RuleManager
  Exceptions: throws if not passed in a DBAdaptor
  Example   : my $rulemanager = Bio::EnsEMBL::Pipeline::RuleManager 
              new->( -DB => $db);

=cut



sub new{
  my ($class,@args) = @_;

  my $self = bless {},$class;

  &verbose('WARNING');

  my ($db, $input_ids, $rules, $queue_manager,
      $awol_mark, $rename, $max_sleep,
      $min_sleep, $base_sleep, $verbose,
      $runner, $output_dir, $job_limit,$delete_lock,$number_output_dirs) = rearrange (
         ['DB', 
          'INPUT_IDS', 
          'RULES', 
          'QUEUE_MANAGER', 
          'MARK_AWOL',
          'RENAME_ON_RETRY',
          'MAX_JOB_SLEEP',
          'MIN_JOB_SLEEP',
          'SLEEP_PER_JOB',
          'VERBOSE',
          'RUNNER',
          'OUTPUT_DIR',
          'JOB_LIMIT',
          'UNLOCK',
          'NUMBER_OUTPUT_DIRS',
         ],@args);

  if(!$db){
    throw("Can't run the RuleManager without a dbadaptor");
  }

  $self->db($db);
  $self->delete_lock if $delete_lock ; 
  $self->is_locked;
  $self->create_lock;
  $self->number_output_dirs($number_output_dirs) ;
  if(!$queue_manager){
    $queue_manager = $QUEUE_MANAGER; #found in BatchQueue.pm
  }

  my $batch_q_module = "Bio::EnsEMBL::Pipeline::BatchSubmission::$queue_manager";

  my $file = "$batch_q_module.pm";
  $file =~ s{::}{/}g;

  eval {
    require "$file";
  };

  if ($@) {
    throw("Can't find $file [$@]");
  }

  if (!defined($awol_mark)) {
    $awol_mark = $MARK_AWOL_JOBS; #found in BatchQueue.pm;
  }

  if (!defined($rename)) {
    $rename = $RENAME_ON_RETRY; #found in General.pm;
  }

  if (!defined($max_sleep)) {
    $max_sleep = $MAX_JOB_SLEEP; #found in BatchQueue.pm
  }

  if (!defined($min_sleep)) {
    $min_sleep = $MIN_JOB_SLEEP;
  }

  if (!defined($base_sleep)) {
    $base_sleep = $SLEEP_PER_JOB;
  }
  if (!defined($job_limit)) {
    $job_limit = $JOB_LIMIT;
  }

  if (!$runner) {
    $runner = $DEFAULT_RUNNER;
  }

  $self->batch_q_module($batch_q_module);
  $self->input_ids($input_ids) if($input_ids);
  $self->rules($rules) if($rules);
  $self->mark_awol_jobs($awol_mark);
  $self->rename_on_retry($rename);
  $self->max_job_sleep($max_sleep);
  $self->min_job_sleep($min_sleep);
  $self->sleep_per_job($base_sleep);
  $self->be_verbose($verbose);
  $self->runner($runner);
  $self->output_dir($output_dir);
  $self->job_limit($job_limit);
  return $self;
}


###################
#container methods#
###################


=head2 db

  Arg [1]   : Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
  Function  : stores the DBadaptor for the object
  Returntype: Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
  Exceptions: throws if argument passed in isn't a 
              Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor'
  Example   : my $rules_adaptor = $self->db->get_RulesAdaptor;

=cut

sub db{
  my ($self, $db) = @_;

  if ($db) {
    if (!$db->isa('Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor')) {
      throw("Can't run the RuleManager with $db you need a ".
            "Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor");
    }
    $self->{'dbadaptor'} = $db;
  }

  return $self->{'dbadaptor'};
}


=head2 adaptors

  Arg [1]   : Bio::EnsEMBL::Pipeline::DBSQL::Adaptormodule
  Function  : store and/or create appropriate adaptors
  Returntype: the requested adaptor
  Exceptions: none
  Example   : my @rules = $self->rule_adaptor->fetch_all;

=cut

sub job_adaptor{
  my ($self, $adaptor) = @_;

  if ($adaptor) {
    $self->{'job_adaptor'} = $adaptor;
  }

  if (!$self->{'job_adaptor'}) {
    my $job_adaptor = $self->db->get_JobAdaptor;
    $self->{'job_adaptor'} = $job_adaptor;
  }
  return $self->{'job_adaptor'};
}


sub analysis_adaptor{
  my ($self, $adaptor) = @_;

  $self->{'analysis_adaptor'} = $adaptor;

  if (!$self->{'analysis_adaptor'}) {
    $self->{'analysis_adaptor'} = $self->db->get_AnalysisAdaptor;
  }
  return $self->{'analysis_adaptor'};
}

sub rule_adaptor{
  my ($self, $adaptor) = @_;

  if ($adaptor) {
    $self->{'rule_adaptor'} = $adaptor;
  }
  if (!$self->{'rule_adaptor'}) {
    my $rule_adaptor = $self->db->get_RuleAdaptor;
    $self->{'rule_adaptor'} = $rule_adaptor;
  }
  return $self->{'rule_adaptor'};
}


sub stateinfocontainer{
  my ($self, $adaptor) = @_;

  if ($adaptor) {
    $self->{'stateinfocontainer'} = $adaptor;
  }
  if (!$self->{'stateinfocontainer'}) {
    my $stateinfocontainer = $self->db->get_StateInfoContainer;
    $self->{'stateinfocontainer'} = $stateinfocontainer;
  }
  return $self->{'stateinfocontainer'};
}

=head2 input_ids

  Arg [1]   : Hash references to hash of input_ids
  Function  : container for hash of input ids keyed on input_id_type
  Returntype: Hash ref
  Exceptions: throws if not passed a hash ref
  Example   : my %input_ids = %{$self->input_ids};

=cut

sub input_ids{
  my ($self, $input_ids) = @_;

  if ($input_ids) {
    throw("Must have a hash ref of input_ids not a $input_ids ") 
      unless(ref($input_ids) eq 'HASH');
    $self->{'input_ids'} = $input_ids;
  }

  if (!$self->{'input_ids'}) {
    $self->{'input_ids'} = $self->stateinfocontainer->get_all_input_id_analysis_sets;
  }
  return $self->{'input_ids'};
}


sub empty_input_ids{
  my ($self) = @_;
  $self->{'input_ids'} = undef;
}

=head2 rules

  Arg [1]   : array ref to array of Bio::EnsEMBL::Pipeline::Rule
  Function  : container for array of rules
  Returntype: arrayref
  Exceptions: throws if not passed an array ref;
  Example   : my @rules = @{$self->rules};

=cut



sub rules{
  my ($self, $rules) = @_;
  if ($rules) {
    throw("Must has an array ref of rules not a $rules ") 
      unless(ref($rules) eq 'ARRAY');
    $self->{'rules'} = $rules;
  }
  if (!$self->{'rules'}) {
    my @rules = $self->rule_adaptor->fetch_all;
    $self->{'rules'} = \@rules;
  }
  return $self->{'rules'};
}

sub empty_rules{
  my ($self) = @_;
  $self->{'rules'} = undef;
}

=head2 boolean toggles 

  Arg [1]   : int
  Function  : boolean toggle as to whether perform a certain action
  Returntype: int
  Exceptions: none
  Example   : if($self->rename_on_retry)

=cut


#whether to rename job's output files when retried
sub rename_on_retry{
  my $self = shift;

  $self->{'rename_on_retry'} = shift if (@_);
  return $self->{'rename_on_retry'};
}

#whether to mark jobs which have disappeared from the submission system
sub mark_awol_jobs{
  my $self = shift;

  $self->{'mark_awol_jobs'} = shift if (@_);
  return $self->{'mark_awol_jobs'};
}


=head2 batch_q_module

  Arg [1]   : string, perl path to batch_q_module
  Function  : container for name of batch_q_module
  Returntype: string
  Exceptions: none
  Example   : $self->batch_q_module->job_stats;

=cut

sub batch_q_module{
  my $self = shift;

  $self->{'batch_q_module'} = shift if (@_);
  return $self->{'batch_q_module'};
}


=head2 sleep_values

  Arg [1]   : int
  Function  : containers for variables concerned with how long 
  the RuleManager sleeps for, or the maximum number of jobs in job_limit's
  case'
  Returntype: int
  Exceptions: none
  Example   : 

=cut


sub max_job_sleep{
  my $self = shift;
  $self->{'max_job_sleep'} = shift if (@_);
  return $self->{'max_job_sleep'};
}

sub min_job_sleep{
  my $self = shift;

  $self->{'min_job_sleep'} = shift if (@_);
  return $self->{'min_job_sleep'};
}

sub sleep_per_job{
  my $self = shift;

  $self->{'sleep_per_job'} = shift if (@_);
  return $self->{'sleep_per_job'};
}

sub job_limit{
  my $self = shift;

  $self->{'job_limit'} = shift if (@_);
  return $self->{'job_limit'};
}

=head2 runner/output_dir

  Arg [1]   : string, path to file
  Function  : container for path to output_dir/runner script
  Returntype: string
  Exceptions: none
  Example   : 
  Notes     : Do not set output_dir if you want the value to be take
  from the BatchQueue.pm configuration
=cut

sub runner{
  my $self = shift;

  $self->{'runner'} = shift if (@_);
  return $self->{'runner'};
}

sub output_dir{
  my $self = shift;

  $self->{'output_dir'} = shift if (@_);
  return $self->{'output_dir'};
}




=head2 be_verbose

  Arg [1]   : int
  Function  : flag for prints
  Returntype: int
  Exceptions: none
  Example   : 

=cut

#called be_verbose so not to interfere with
#Bio::EnsEMBL::Utils::Exception::verbose
sub be_verbose{
  my $self = shift;

  $self->{'verbose'} = shift if (@_);
  return $self->{'verbose'};
} 

#################
#Utility methods#
#################


=head2 read_id_file

  Arg [1]   : filename
  Function  : reads a file in the format input_id input_id_type and
  places entires in an 2D hash ref intially keyed on input_id_type
  then on input_id
  Returntype: hash ref
  Exceptions: throws if can't open file'
  Example   : my %ids_not_to_run = $rulemanager->read_id_file($skip_ids);

=cut

sub read_id_file{
  my ($self, $file) = @_;

  my %ids;
  open IDS, "< $file" or throw("Can't open $file");
  while (<IDS>) {
    chomp;
    my ($id, $type) = split;
    if (!$type) {
      throw(" id ".$id ." type ".$type."\n need to know what type the ".
            "input ids in the file are  format should be input_id type\n");
    }
    $ids{$type}{$id} = 1;
  }
  return \%ids
}


=head2 process_types

  Arg [1]   : arrayref to array of strings
  Function  : processes a list of types into a hash
  Returntype: hashref
  Exceptions: throws if not passed an array ref
  Example   : %types_to_run = %{$rulemanager->process_types(\@types)};

=cut

sub process_types{
  my ($self, $types) = @_;

  throw("Can't process $types if not array ref") unless(ref($types) eq 'ARRAY');

  my %types;
  foreach my $type (@$types){
    $types{$type} = 1;
  }
  return \%types;
}



=head2 starts_from_input_ids

  Arg [1]   : arrayref to array of Bio::EnsEMBL::Pipeline::Analysis objects
  Function  : uses these analysis objects to generate hash of input_ids
  Returntype: 2d hashref
  Exceptions: throws if not passed a array ref
  Example   : my $input_id_hash = 
  $rulemanager->starts_from_input_ids(\@analyses);

=cut


sub starts_from_input_ids{
  my ($self, $analyses) = @_;

  throw("Can't process $analyses if not array ref") unless(ref($analyses) eq 'ARRAY');

  my %ids;
  foreach my $analysis (@$analyses) {
    my @ids = @{$self->stateinfocontainer->list_input_ids_by_analysis($analysis->dbID)};
    print "analysis ".$analysis->logic_name." has got ".@ids." ids\n";
    foreach my $id (@ids){
      $ids{$analysis->input_id_type}{$id} = 1;
    }
  }
  return \%ids;
}

=head2 create_and_store_job

  Arg [1]   : string, input_id
  Arg [2]   : Bio::EnsEMBL::Pipeline::Analysis
  Arg [3]   : string, directory path  (optional)
  Arg [4]   : string, runner script path (optional)
  Arg [5]   : int, for a boolean flag to mark verbosity (optional)
  Function  : Create a job based on the input_id and analysis passed in.
  Store it in the database and submit to be run
  Returntype: Bio::EnsEMBL::Pipeline::Job
  Exceptions: throws if not passed an input_id or an analysis object and
  if fails to store the job 
  Example   : my $job = $rulemanager->create_and_store_job('filename', 
                                                           $analysis
                                                           'path/to/dir', 
                                                           'runnerscript',
                                                           1);

=cut

sub create_and_store_job{
  my ($self, $input_id, $analysis) = @_;

  if(!$input_id || !$analysis){
    throw("Can't create job without an input_id $input_id or analysis ".
          "$analysis");
  }  

  my $job = Bio::EnsEMBL::Pipeline::Job->new
    (
     -input_id => $input_id,
     -analysis => $analysis,
     -output_dir => $self->output_dir,
     -runner => $self->runner, 
     -number_output_dirs => $self->number_output_dirs, 
    );

  eval{
    $self->job_adaptor->store($job);
  };

  if ($@) {
    throw("Failed to store job ".$job->input_id." ".
          $job->analysis->logic_name." ".$@);
  } else {
    print "Stored ".$job->dbID." ".$job->input_id." ".
      $job->analysis->logic_name."\n" if ($self->be_verbose);
  }
  return $job;
}

=head2 can_run_job

  Arg [1]   : string, input_id
  Arg [2]   : Bio::EnsEMBL::Pipeline::Analysis
  Arg [3]   : string, directory path  (optional)
  Arg [4]   : string, runner script path (optional)
  Arg [5]   : int, for a boolean flag to mark verbosity (optional)
  Function  : Check if a job can be created for an input_id and analysis
              If a job already exists check if it needs to be retried
  Returntype: int
  Exceptions: throws if not passed an input_id or an analysis object and
              if fails to submit the job 
  Example   : $rulemanager->can_run_job('filename', $analysis
                                        'path/to/dir', 
                                        'path/to/runner', 1);

=cut

sub can_job_run{
  my ($self, $input_id, $analysis, $current_jobs) = @_;

  if (!$input_id || !$analysis) {
    throw("Can't create job without an input_id $input_id or analysis ".
          "$analysis");
  }

  my $job;

  if ($current_jobs->{$analysis->dbID}) {
    my $cj = $current_jobs->{$analysis->dbID};
    #my $status = $cj->current_status->status;

    # This is a hack to get the status at the same time as the job was retrieved
    # _status should be private but current_status is not useful to us here as 
    # it will do a new query. Here we don't need the absolutely up-to-the-second
    # status, we can always get the status on the next round of checks. This
    # will also significantly reduce the number of queries of the job table.
    #
    # One problem is jobs which are in the middle of updating their current status.
    # These can be in a state with no current status, which will mean they don't get 
    # listed in the current_jobs hash.
    my $status = $cj->{_status}->status;

    if (($status eq 'FAILED' || $status eq 'AWOL') && $cj->can_retry) {
      print "\nRetrying job with status $status!!!!\n" if $self->be_verbose;

      if ($self->rename_on_retry) {
        $self->rename_files($cj);
      }
      $cj->set_status('CREATED');
      $job = $cj;
    }
  } else {
    $job = $self->create_and_store_job($input_id, $analysis);
  }

  if ($job) {
    eval {
      print "\tBatch running job\n" if $self->be_verbose;
      $job->batch_runRemote;
    };
    if ($@) {
      throw("ERROR running job " . $job->dbID . " ".
            $job->analysis->logic_name." " . $job->stderr_file . " [$@]");
    }
    return 1;
  }

  return 0;
}


=head2 rename_files

  Arg [1]   : Bio::EnsEMBL::Pipeline::Job
  Function  : rename the stderr and stdout files of a job when it is 
  retried
  Returntype: int
  Exceptions: throws if not passed a Job object or if rename fails
  Example   : $self->rename_files($cj);

=cut


sub rename_files{
  my ($self, $job) = @_;
  
  if (!$job || !$job->isa("Bio::EnsEMBL::Pipeline::Job")) {
    throw("Need a job object " . $job . " to rename its files\n");
  }
  eval{
    move($job->stdout_file, $job->stdout_file.'.retry'.$job->retry_count);
    move($job->stderr_file, $job->stderr_file.'.retry'.$job->retry_count);
  };

  if ($@) {
    throw("Couldn't rename job ".$job->id."'s output files $@ ");
  }
  return 1;
}


=head2 job_stats/job_limit_check

  Arg [1]   : int, max number of jobs (optional)
  Arg [2]   : array_ref to array of Bio::EnsEMBL::Pipeline::Job objects
  (optional)
  Function  : gets statistics from BatchSubmission module about what
  jobs are running and what their status is then take action on this 
  information. job_limit_check checks how many jobs of a defined status
  there are and sleeps if appropriate, job_stats will mark awol jobs
  as well
  Returntype: int
  Exceptions: throws if batch_submission module can't do the method
  job stats'
  Example   : $rulemanager->job_stats;('1000', \@jobs);

=cut


sub job_stats{
  my ($self, $job_limit, $jobs) = @_;

  if (!$job_limit) {
    $job_limit = $self->job_limit;
  }


  # Do job_stats call before getting jobs
  if (!$self->batch_q_module->can('job_stats')) {
    throw($self->batch_q_module." doesn't have the job_stats method");
  }
  my %statuses_to_count = map{$_, 1} @{$JOB_STATUSES_TO_COUNT}; #found in
                                                                #BatchQueue.pm
  my %job_stats = %{$self->batch_q_module->job_stats};

  my @jobs;
  if (!$jobs) {
    @jobs = $self->job_adaptor->fetch_all;
  } else {
    @jobs = @$jobs;
  }

  my @awol_jobs;
  my $job_count = 0;
  JOB:foreach my $job (@jobs) {
    if (!$job_stats{$job->submission_id}) {
      push(@awol_jobs, $job);
      next JOB;
    }
    if ($statuses_to_count{$job_stats{$job->submission_id}}) {
      $job_count++;
    }
  }
  if ($self->mark_awol_jobs) {
    foreach my $awol (@awol_jobs){
      if ($self->valid_statuses_for_awol->{$awol->current_status->status}) {
        $awol->set_status('AWOL');
      }
    }
  }
  if ($job_count >= $job_limit) {
    $self->sleep($job_count, $job_limit);
  }
  return 1;
}



sub job_limit_check{
  my ($self, $job_limit, $jobs) = @_;
  
  if (!$job_limit) {
    $job_limit = $self->job_limit;
  }

  my @jobs;
  if (!$jobs) {
    @jobs = $self->job_adaptor->fetch_all;
  } else {
    @jobs = @$jobs;
  }

  if (!$self->batch_q_module->can('job_stats')) {
    throw($self->batch_q_module." doesn't have the job_stats method");
  }
  my %statuses_to_count = map{$_, 1} @{$JOB_STATUSES_TO_COUNT}; #found in
                                                                #BatchQueue.pm
  my %job_stats = %{$self->batch_q_module->job_stats};
  my $job_count = 0;

  JOB:foreach my $job (@jobs) {
    if ($statuses_to_count{$job_stats{$job->submission_id}}) {
      $job_count++;
    }
  }

  if ($job_count >= $job_limit) {
    $self->sleep($job_count, $job_limit);
  }
  return 1;
}

=head2 valid_statuses_for_awol/kill

  Arg [1]   : none
  Function  : produces a hash of valid statuses to check when either
  marking jobs as awol or killing jobs for running for too long
  Returntype: hashref
  Exceptions: none
  Example   : if($self->
                 valid_statuses_for_awol->{$job->current_status->status})

=cut

sub valid_statuses_for_awol {
  my ($self)  = @_;
  if(!$self->{'status_for_awol'}) {
    my %statuses = map{$_, 1} ('SUBMITTED', 'RUNNING', 'READING',
                               'WRITING', 'WAITING');
      $self->{'status_for_awol'} = \%statuses;
  }
  return $self->{'status_for_awol'};
}

sub valid_statuses_for_kill{
  my ($self)  = @_;

  if (!$self->{'status_for_kill'}) {
    my %statuses = map{$_, 1} ('RUNNING', 'WAITING');
      $self->{'status_for_kill'} = \%statuses;
  }
  return $self->{'status_for_kill'};
}



=head2 sleep

  Arg [1]   : int, number of jobs
  Arg [2]   : int, job_limit
  Function  : sleep for an appropriate amount of time given the number
  of jobs over the maximum there is
  Returntype: int the amount of time sleeped for
  Exceptions: none
  Example   : $self->sleep($job_count, $job_limit);

=cut

sub sleep{
  my ($self, $job_number, $job_limit) = @_;
  
  my $extra_jobs = $job_number - $job_limit;
  my $sleep      = $extra_jobs * $self->sleep_per_job;

  $sleep = $self->max_job_sleep if ($sleep > $self->max_job_sleep);
  $sleep = $self->min_job_sleep if ($sleep < $self->min_job_sleep);
 
  sleep($sleep);

  return $sleep;
}




=head2 create_lock

  Arg [1]   : none
  Function  : create a string to lock the pipeline database with
  Returntype: string
  Exceptions: none
  Example   : $lock_str = $self->create_lock;

=cut


sub create_lock{
  my ($self) = @_;

  my $host = $self->qualify_hostname(hostname());
  my $user = scalar getpwuid($<);
  my $lock_str = join ":", "$user\@$host", $$, time();

  $self->db->pipeline_lock($lock_str);

  return $lock_str;
}


=head2 delete_lock 

  Arg [1]   : none
  Function  : removes entry from meta-table of database 'pipeline.lock'  
  Returntype: string
  Exceptions: none
  Example   : $lock_str = $self->delete_lock;

=cut


sub delete_lock {
  my ($self) = @_;

  my $host = $self->qualify_hostname(hostname());
  my $user = scalar getpwuid($<);
  my $lock_str = join ":", "$user\@$host", $$, time();

  my $arg = $self->db->pipeline_unlock();
  return $arg;
}

 
=head2 qualify_hostname

  Arg [1]   : string, output of Sys::Hostname::hostname
  Function  : produces fully qualified host name
  Returntype: string
  Exceptions: none
  Example   : my $host = $self->qualify_hostname(hostname());

=cut


sub qualify_hostname{
  my ($self, $hostname) = @_;

  my $addr = gethostbyname($hostname);
  my $host = gethostbyaddr($addr, AF_INET);
  return $host;
}



=head2 is_locked

  Arg [1]   : none
  Function  : to establish if the pipeline is locked, ie if there is  lock
  string already present in the database
  Returntype: int
  Exceptions: throws if a lock str is present
  Example   : $self->is_locked

=cut

sub is_locked{
  my ($self) = @_;

  if (my $lock_str = $self->db->pipeline_lock) {
    my($user, $host, $pid, $started) = $lock_str =~ /(\w+)@(\w+):(\d+):(\d+)/;

    $started = scalar localtime $started;

    my $dbname = $self->db->dbname;
    my $dbhost = $self->db->host;

    my $error_str = ("Error: this pipeline appears to be running!\n\n".
                     "\tdb       $dbname\@$dbhost\n".
                     "\tpid      $pid on ".
                     "host $host\n\tstarted  $started\n\n".
                     "The process above must be terminated before this ".
                     "script can be run.\nIf the process does not exist, ".
                     "remove the lock by removing the lock from the ".
                     "pipeline database:\n\ndelete from meta where ".
                     "meta_key = 'pipeline.lock';\n\n\n\n" . 
                     "\tYou could also use the -unlock option to remove the lock" . 
                     "\n\n\n Thank you !\n\n");

    print STDERR $error_str;
    throw("Can't run RuleManager there may be another rulemanager ".
          "running look in ".$self->db->dbname." meta table ");
  }
  return 1;
}

=head2 fetch_complete_accumulators

  Arg [1]   : none
  Function  : fetches the analysis of accumulators which have already run
  Returntype: hash ref
  Exceptions: warns if an analysis with input_id ACCUMULATOR doesn't have 
  the type accumulator'
  Example   : my %complete_accumulators = %{$self->
                                              fetch_complete_accumulators};

=cut

sub fetch_complete_accumulators{
  my ($self) = @_;
  
  $self->{'complete_accumulators'} = {};

  my @accumulators = @{$self->stateinfocontainer->fetch_analysis_by_input_id('ACCUMULATOR')};

  foreach my $analysis (@accumulators) {
    if ($analysis->input_id_type eq 'ACCUMULATOR') {
      $self->{'complete_accumulators'}->{$analysis->logic_name} = 1;
    } else {
      warn(" analysis " . $analysis->logic_name . " must have input id " .
           "type ACCUMULATOR");
    }
  }
  return $self->{'complete_accumulators'};
}


=head2 cleanup_waiting_jobs

  Arg [1]   : none
  Function  : ensures any jobs which are sat waiting to be run at the end of
              a pipeline get left unrun
  Returntype: 
  Exceptions: 
  Example   : 

=cut

sub cleanup_waiting_jobs{
  my ($self) = @_;

  my ($a_job) = $self->job_adaptor->fetch_by_Status("CREATED");

  if ($a_job) {
    $a_job->flush_runs($self->job_adaptor);
  } else {
    print STDERR "have no jobs to clean up\n" if ($self->be_verbose);
  } 
}


=head2 add_created_jobs_back

  Arg [1]   : none
  Function  : ensures at created by unsubmitted jobs get readded to the
  queuing system
  Returntype: int
  Exceptions: none
  Example   : 

=cut

sub add_created_jobs_back{
  my ($self) = @_;
  my @created_jobs = $self->job_adaptor->fetch_by_Status("CREATED");

  foreach my $j (@created_jobs) {
    $j->batch_runRemote;
  }
  return 1;
}

=head2 rules_setup

  Arg [1]   : hashref, keyed on analysis_id of analyses which are to run
  Arg [2]   : hashref keyed on analysis_id of analyses to skip
  Arg [3]   : arrayref of all rules
  Arg [4]   : hashref of all accumulator analyses
  Arg [5]   : hashref of incomplete accumulators
  Function  : to setup the rules array on the basis of specified analyses
              to either run or skip
  Returntype: arrayref
  Exceptions: throws if no rules are produces at the end
  Example   : 

=cut

sub rules_setup{
  my ($self, $analyses_to_run, $analyses_to_skip, $all_rules,
     $accumulator_analyses, $incomplete_accumulators) =@_;
  my @rules;
  #print "SETTING UP RULES\n";
  if (keys(%$analyses_to_run)) {
    foreach my $rule (@$all_rules) {
      #print "CHECKING RULE ".$rule->goalAnalysis->logic_name."\n";
      if (exists($analyses_to_run->{$rule->goalAnalysis->dbID})) {
        #print "ADDING RULE\n";
        push (@rules, $rule);
      } elsif ($accumulator_analyses->{$rule->goalAnalysis->logic_name}) {
        $incomplete_accumulators->{$rule->goalAnalysis->logic_name} = 1;
      }
    }
  } elsif (keys(%$analyses_to_skip)) {
    foreach my $rule (@$all_rules) {
      if (!exists($analyses_to_skip->{$rule->goalAnalysis->dbID})) {
        push (@rules, $rule);
      }
    }
  } else {
    @rules = @$all_rules;
  }
  if (scalar(@rules) == 0) {
    throw("Something is wrong with the code or your commandline setup ".
          "rules_setup has returned no rules");
  }
  $self->rules(\@rules);
  return \@rules;
}



=head2 logic_name2dbID

  Arg [1]   : arrayref of Bio::EnsEMBL::Pipeline::Analysis
  Function  : produce a hash keyed on analysis dbID based on array passed 
  in
  Returntype: hashref 
  Exceptions: none
  Example   : 

=cut

sub logic_name2dbID {
  my ($self, $analyses) = @_;
  my %analyses;
  
  foreach my $ana (@$analyses) {
    if ($ana =~ /^\d+$/) {
      $analyses{$ana} = 1;
    } else {
      my $id ; 
      eval {  $id = $self->analysis_adaptor->fetch_by_logic_name($ana)->dbID; }  ; 
      if ($@){
      	print STDERR "\n\nERROR: The analysis (logic_name) you've specified does not exist in the db\n"."\n"x10; 
      	throw($@) ; 
      }
      if ($id) {
        $analyses{$id} = 1;
      } else {
        print STDERR "Could not find analysis $ana\n";
      }
    }
  }
  return \%analyses;
}


=head2 input_ids_setup

  Arg [1]   : string, filename pointing file of input_ids to run
  Arg [2]   : string, filename pointing file of input_ids to skip
  Arg [3]   : array ref, array of input_id_types to run
  Arg [4]   : array ref, array of input_id_types to skip
  Function  : prepare the hash of input_ids to use when running 
  Returntype: hashref
  Exceptions: throws if the id_hash is empty
  Example   : 
  Notes     :both files are expected in the format input_id input_id_type
=cut


sub input_id_setup{
  my ($self, $ids_to_run, $ids_to_skip, 
      $types_to_run, $types_to_skip, $starts_from) = @_;
  my $id_hash;

  $self->empty_input_ids;
  
  if ($ids_to_run) {
    
    $id_hash = $self->read_id_file($ids_to_run);
  } elsif (@$starts_from) {
    
    my @analyses;

    foreach my $logic_name(@$starts_from){
      my $analysis = $self->analysis_adaptor->
        fetch_by_logic_name($logic_name);
      push(@analyses, $analysis);
    }
    $id_hash = $self->starts_from_input_ids(\@analyses);
  }elsif(scalar(@$types_to_run)){
    my %id_type_hash;
    foreach my $type(@$types_to_run){
      my $ids = $self->stateinfocontainer
        ->list_input_ids_by_type($type);
      foreach my $id  (@$ids) {
        $id_type_hash{$type}{$id} = 1;
      }
    }
    $id_hash = \%id_type_hash;
  } else {
    print "Getting standard set\n";
    $id_hash = $self->input_ids;
  }

  if ($ids_to_skip) {
    my $skip_id_hash = $self->read_id_file($ids_to_skip);
    foreach my $type (keys(%$skip_id_hash)){
      foreach my $id (keys(%{$skip_id_hash->{$type}})){
        delete($skip_id_hash->{$type}->{$id});
      }
    }
  }

  if(@$types_to_skip){
    my %types_to_skip = map{$_, 1} @$types_to_skip;
    foreach my $type (keys(%$id_hash)){
      if($types_to_skip{$type}){
        delete($id_hash->{$type});
      }
    }
  }

  if(keys(%$id_hash) == 0){
    throw("Something is wrong with the code or the commandline ".
          "input_ids_setup has produced no ids\n");
  }

  $self->input_ids($id_hash);

  return $id_hash;
}

=head2 check_if_done

  Arg [1]   : none
  Function  : check if the pipeline is finished running
  Returntype: int
  Exceptions: none
  Example   : 

=cut

sub check_if_done {
  my ($self) = @_;
  my @jobs = $self->job_adaptor->fetch_all;
  my $continue;

 JOB: 
  foreach my $job (@jobs) {
    my $status = $job->current_status->status;

    if ($status eq 'KILLED' || $status eq 'SUCCESSFUL') {
      next JOB;
    } elsif ($status eq 'FAILED' || $status eq 'AWOL') {
      if (!$job->can_retry) {
        next JOB;
      } else {
        return 1;
      }
    } else {
      return 1;
    }
  }

  return 0;
}

sub number_output_dirs { 
  my ($self,$arg) = @_ ; 
  $self->{number_output_dirs} = $arg if $arg ;  
  return $self->{number_output_dirs} ; 
} 

1;
