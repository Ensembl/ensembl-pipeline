
# Object for storing details of an analysis job
#
# Cared for by Michele Clamp  <michele@sanger.ac.uk>
#
# Copyright Michele Clamp
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Job

=head1 SYNOPSIS

=head1 DESCRIPTION

Stores run and status details of an analysis job

=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::Job;


use vars qw(@ISA $SAVE_RUNTIME_INFO);
use strict;
use Bio::EnsEMBL::Pipeline::Config::BatchQueue;
use Bio::EnsEMBL::Pipeline::Config::General;

@ISA = qw(Bio::EnsEMBL::Root);


# dynamically load appropriate queue manager (e.g. LSF)

my $batch_q_module = "Bio::EnsEMBL::Pipeline::BatchSubmission::$QUEUE_MANAGER";

my $file = "$batch_q_module.pm";
$file =~ s{::}{/}g;
require "$file";

# BATCH_QUEUES is a package variable which stores available
# 'queues'. this allows different analysis types to be sent
# to different nodes etc. it is keyed on analysis logic_name
# and has different parameters such as resource (e.g. node set),
# number of jobs to be batched together etc.
#
# this may be better encapsulated as a separate class, but
# i can't think at the moment how best to do this.

my %BATCH_QUEUES = &set_up_queues;

sub new {
    my ($class, @args) = @_;
    my $self = bless {},$class;
    my ($p, $f, $l) = caller;
    my ($adaptor,$dbID,$submission_id,$input_id,$analysis,$stdout,$stderr,$retry_count, $output_dir, $runner) 
	= $self->_rearrange([qw(ADAPTOR
				ID
				SUBMISSION_ID
				INPUT_ID
				ANALYSIS
				STDOUT
				STDERR
				RETRY_COUNT
				OUTPUT_DIR
				RUNNER
				)],@args);
    $dbID    = -1 unless defined($dbID);
    $submission_id   = -1 unless defined($submission_id);
    $input_id   || $self->throw("Can't create a job object without an input_id");
    $analysis   || $self->throw("Can't create a job object without an analysis object");

    $analysis->isa("Bio::EnsEMBL::Analysis") ||
	$self->throw("Analysis object [$analysis] is not a Bio::EnsEMBL::Analysis");
    
    $self->dbID($dbID);
    $self->adaptor($adaptor);
    $self->input_id($input_id);
    $self->analysis($analysis);
    $self->stdout_file($stdout);
    $self->stderr_file($stderr);
    $self->retry_count($retry_count);
    $self->submission_id($submission_id);
    $self->output_dir($output_dir);
    if($self->output_dir){
      $self->make_filenames;
    }else{
      my $dir;
      if(!$BATCH_QUEUES{$analysis->logic_name}){
        $dir =  $DEFAULT_OUTPUT_DIR;
      }else{
        $dir = $BATCH_QUEUES{$analysis->logic_name}{output_dir}  || $DEFAULT_OUTPUT_DIR;
      }
      $self->throw("need an output directory passed in from RuleManager or from Config/BatchQueue $!") unless($dir);
      $self->output_dir($dir);
      $self->make_filenames;
    }
    $self->runner($runner);
    if(!$self->runner){
      $self->runner($PIPELINE_RUNNER_SCRIPT);
    }
    return $self;
}


=head2 create_by_analysis_input_id

  Title   : create_by_analysis_input_id
  Usage   : $class->create_by.....
  Function: Creates a job given an analysis object and an input_id
            Recommended way of creating job objects!
  Returns : a job object, not connected to db
  Args    : 

=cut

sub create_by_analysis_input_id {
  my ($dummy, $analysis, $inputId, $output_dir, $auto_update) = @_;
  
  $dummy->warn("Bio::EnsEMBL::Pipeline::Job->create_by_analysis_input_id is deprecated ".(caller)." should now call the constructor directly");
  my $job = Bio::EnsEMBL::Pipeline::Job->new(
					     -input_id    => $inputId,
					     -analysis    => $analysis,
					     -output_dir => $output_dir,
					     -auto_update => $auto_update,
					     -retry_count => 0
					    );
  #$job->make_filenames;
  return $job;
}


=head2 dbID

  Title   : dbID
  Usage   : $self->dbID($id)
  Function: get set the dbID for this object, only used by Adaptor
  Returns : int
  Args    : int

=cut

sub dbID {
    my ($self, $arg) = @_;

    if ($arg) {
	$self->{'_dbID'} = $arg;
    }
    return $self->{'_dbID'};

}

sub runner {
    my ($self, $arg) = @_;

    if ($arg) {
	$self->{'_runner'} = $arg;
    }
    return $self->{'_runner'};

}


=head2 adaptor

  Title   : adaptor
  Usage   : $self->adaptor
  Function: get database adaptor, set only for constructor and adaptor usage. 
  Returns : 
  Args    : 

=cut

sub adaptor {
    my ($self, $arg) = @_;

    if ($arg) {
	$self->{'_adaptor'} = $arg;
    }
    return $self->{'_adaptor'};
}


=head2 input_id

  Title   : input_id
  Usage   : $self->input_id($id)
  Function: Get/set method for the id of the input to the job
  Returns : string
  Args    : string

=cut

sub input_id {
    my ($self, $arg) = @_;

    if (defined($arg)) {
	$self->{'_input_id'} = $arg;
    }
    return $self->{'_input_id'};
}


=head2 analysis

  Title   : analysis
  Usage   : $self->analysis($anal);
  Function: Get/set method for the analysis object of the job
  Returns : Bio::EnsEMBL::Analysis
  Args    : Bio::EnsEMBL::Analysis

=cut

sub analysis {
    my ($self, $arg) = @_;

    if (defined($arg)) {
	$self->throw("[$arg] is not a Bio::EnsEMBL::Analysis object" ) 
            unless $arg->isa("Bio::EnsEMBL::Analysis");

	$self->{'_analysis'} = $arg;
    }
    return $self->{'_analysis'};
}


=head2 flush_runs

  Title   : flush_runs
  Usage   : $job->flush_runs( jobadaptor, [queue] );
  Function: Issue all jobs in the queue and empty the queue.
    Set LSF id in all jobs. Uses the given adaptor for connecting to
    db. Uses last job in queue for stdout/stderr. 
  Returns : 
  Args    : 

=cut

sub flush_runs {
  my ($self, $adaptor, $queue) = @_;

  # flush_runs is optionally sent a queue to deal with
  # @analyses is a list of logic_names (strings)
  my @analyses = ($queue) || (keys %BATCH_QUEUES);
  
  if( !defined $adaptor ) {
    $self->throw( "Cannot run remote without db connection" );
  }
  
  local *FILE;

  my $db       = $adaptor->db;
  my $host     = $db->host;
  my $username = $db->username;
  my $dbname   = $db->dbname;
  my $pass     = $db->password;
  my $port     = $db->port;

  # runner.pl: first look at value set in RuleManager ($RUNNER_SCRIPT)
  # then in same directory as Job.pm,
  # and fail if not found

  my $runner = $self->runner;
 
  unless (-x $runner) {
    $runner = __FILE__;
    $runner =~ s:/[^/]*$:/runner.pl:;
    my $caller = caller(0);
    $self->throw("runner ".$runner." not found - needs to be set in $caller\n") unless -x $runner;
  }

  # $anal is a logic_name
 
 ANAL:for my $anal (@analyses) {
   
    my $queue = $BATCH_QUEUES{$anal};
    
    my @job_ids = @{$queue->{'jobs'}};
    if(!@job_ids){
      next ANAL;
    }
    my $this_runner = $queue->{'runner'};
    $this_runner = (-x $this_runner) ? $this_runner : $runner;
   
    my $lastjob = $adaptor->fetch_by_dbID($job_ids[-1]);
    
    if( ! $lastjob ) {
      $self->throw( "Last batch job not in db" );
    }
    
    my $pre_exec = $this_runner." -check -output_dir ".$self->output_dir;
   

    my $batch_job = $batch_q_module->new
      (
       -STDOUT     => $lastjob->stdout_file,
       -PARAMETERS => $queue->{'sub_args'},
       -PRE_EXEC   => $pre_exec,
       -QUEUE      => $queue->{'queue'},
       -JOBNAME    => $dbname . ':' . $anal,
       -NODES      => $queue->{'nodes'},
       -RESOURCE   => $queue->{'resource'}
      );
    my $cmd;
    if(!$self->cleanup){
      $batch_job->stderr_file($lastjob->stderr_file);
    }
    

    # check if the password has been defined, and write the
    # "connect" command line accordingly otherwise -pass gets the
    # first job id as password, instead of remaining undef
    if ($pass) {
      $cmd = $runner." -dbhost $host -dbuser $username -dbname $dbname -dbpass $pass -dbport $port";
    }
    else {
      $cmd = $runner." -dbhost $host -dbuser $username -dbname $dbname -dbport $port";
    }
    $cmd .= " -output_dir ".$self->output_dir;
    $cmd .= " -queue_manager $QUEUE_MANAGER  " ;
    if($self->cleanup){
      $cmd .= " -cleanup "
    }
    $cmd .= " @job_ids";
    
    $batch_job->construct_command_line($cmd);
    eval {
	#print STDERR "Submitting: ", $batch_job->bsub, "\n";
	$batch_job->open_command_line();
    };
    if ($@) {
	print STDERR "Couldnt batch submit @job_ids \n[$@]\n";
	foreach my $job_id (@job_ids) {
	    my $job = $adaptor->fetch_by_dbID($job_id);
	    $job->set_status( "FAILED" );
	}
    } 
    else {
	#print STDERR "have submitted ".@job_ids." jobs with ".$batch_job->id
	#."\n";
	my @jobs = $adaptor->fetch_by_dbID_list(@job_ids);
	foreach my $job (@jobs) {
	    if( $job->retry_count > 0 ) {
		for ( $job->stdout_file, $job->stderr_file ) {
		    open( FILE, ">".$_ ); close( FILE );
		}
	    }
	    #print STDERR "altering stderr file to ".$lastjob->stderr_file."\n";
	    if ($batch_job->id) {
		$job->submission_id( $batch_job->id );
	    } else {
		# submission seems to have succeeded, but we didnt
		# get a job ID. Safest NOT to raise an error here,
		# (a warning would have already issued) but flag
		print STDERR "Job: Null submission ID for the following, but continuing: @job_ids\n";
		$job->submission_id( 0 );		
	    }
	    
	    $job->retry_count( $job->retry_count + 1 );
	    $job->set_status( "SUBMITTED" );
	    $job->stdout_file($lastjob->stdout_file);
	    $job->stderr_file($lastjob->stderr_file);
	}
	$adaptor->update(@jobs);
    }
    $queue->{'jobs'} = [];
    $queue->{'last_flushed'} = time;
  }
}


=head2 batch_runRemote

  Title   : batch_runRemote
  Usage   : $job->batch_runRemote
  Function: Issue more than one pipeline Job in one batch job because 
    job submission is very slow
  Returns : 
  Args    : Is static, private function, dont call with arrow notation.

=cut

sub batch_runRemote {
  my ($self) = @_;

  my $queue;
  
  if (!$BATCH_QUEUES{$self->analysis->logic_name}) {
      $queue = 'default';
  }
  else {
      $queue = $self->analysis->logic_name;
  }
  
  push @{$BATCH_QUEUES{$queue}{'jobs'}}, $self->dbID;
 
  if (scalar(@{$BATCH_QUEUES{$queue}{'jobs'}}) >=
               $BATCH_QUEUES{$queue}{'batch_size'}) {

    $self->flush_runs($self->adaptor, $queue);
  }
}


=head2 runLocally

  Title   : running
  Usage   : $self->run...;
  Function: runLocally doesnt submit to LSF
            run_module is like runLocally, but doesnt redirect STDOUT and 
            STDERR. 
            runRemote submits to LSF via the runner.pl script.
  Returns : 
  Args    : 

=cut

sub runLocally {
  my $self = shift;
 
  print STDERR "Running locally " . $self->stdout_file . " " . $self->stderr_file . "\n"; 

  local *STDOUT;
  local *STDERR;

  if( ! open ( STDOUT, ">".$self->stdout_file )) {
    $self->set_status( "FAILED" );
    return;
  }
        
  if( ! open ( STDERR, ">".$self->stderr_file )) {
    $self->set_status( "FAILED" );
    return;
  }
       print STDERR "Running inLSF\n"; 
  $self->run_module();
}


# question, when to submit the success report to the db?
# we have to parse the output of LSF anyway....
sub run_module {
  my $self = shift;
  my $module = $self->analysis->module;
  print STDERR "Running ".$module." with ".$self."\n";
  my $rdb;
  my ($err, $res);
  my $autoupdate = $AUTO_JOB_UPDATE;
  my $hash_key = $self->analysis->logic_name;
  if(!$BATCH_QUEUES{$hash_key}){
    $hash_key = 'default';
  }
  my $runnable_db_path = 
    $BATCH_QUEUES{$hash_key}{runnabledb_path};
  my $perl_path;
  print STDERR "Getting ".$hash_key." batchqueue value\n";
  if($module =~ /::/){
    print STDERR "Module contains path info already\n";
    $module =~ s/::/\//g;
    $perl_path = $module;
  }else{
    $perl_path = $runnable_db_path."/".$module;
  }
  print STDERR "have perlpath ".$perl_path."\n";
 STATUS: 
  { 
    eval {
      require $perl_path.".pm";
      $perl_path =~ s/\//::/g;
      $rdb = $perl_path->new( -analysis => $self->analysis,
                              -input_id => $self->input_id,
                              -db => $self->adaptor->db );
    };
      
    if ($err = $@) {
      print (STDERR "CREATE: Lost the will to live Error\n");
      $self->set_status( "FAILED" );
      $self->throw( "Problems creating runnable $module for " . 
                    $self->input_id . " [$err]\n");
    }
    
    # "READING"
    eval {   
      $self->set_status( "READING" );
      $res = $rdb->fetch_input;
    };
    if ($err = $@) {
      $self->set_status( "FAILED" );
      print (STDERR "READING: Lost the will to live Error\n");
      $self->throw( "Problems with $module fetching input for " . 
                    $self->input_id . " [$err]\n");
    }
    
    if ($rdb->input_is_void) {
      $self->set_status( "VOID" );
    }
    else {
      # "RUNNING"
      
      eval {
        $self->set_status( "RUNNING" );
        $rdb->db->disconnect_when_inactive(1); 
        $rdb->run;
        $rdb->db->disconnect_when_inactive(0); 
      };
      
      if ($err = $@) {
        if(my $err_state = $rdb->failing_job_status){
          $self->set_status( $err_state );
        }else{
          $self->set_status( "FAILED" ); # default to just failed these jobs get retried
        }
        print (STDERR "RUNNING: Lost the will to live Error\n");
        $self->throw("Problems running $module for " . $self->input_id . " [$err]\n");
      }
      
      # "WRITING"
      eval {
        
        $self->set_status( "WRITING" );
        $rdb->write_output;
	      if($rdb->can('db_version_searched')){
          my $new_db_version = $rdb->db_version_searched();
          my $analysis = $self->analysis();
          my $old_db_version = $analysis->db_version();
          $analysis->db_version($new_db_version);
          # where is the analysisAdaptor??
          # $self->adaptor->get_AnalysisAdaptor->store($analysis);
          # if $new_db_version gt $old_db_version;
	      } else {
          $SAVE_RUNTIME_INFO = 0;
	      }
	      $self->set_status( "SUCCESSFUL" );
      }; 
      if ($err = $@) {
	      $self->set_status( "FAILED" );
	      print (STDERR "WRITING: Lost the will to live Error\n");
	      $self->throw("Problems for $module writing output for " . $self->input_id . " [$err]" );
      }
    }
    
  }    
  
  # update job in StateInfoContainer
  if ($autoupdate) {
    eval {
	    my $sic = $self->adaptor->db->get_StateInfoContainer;
	    # -------------------------------------------------------------
	    $sic->store_input_id_analysis(
                                    $self->input_id,
                                    $self->analysis,
                                    $self->execution_host,
                                    $SAVE_RUNTIME_INFO
                                   );
	    # -------------------------------------------------------------
    };
    if ($err = $@) {
	    print STDERR "Error updating successful job ".$self->dbID .
        "[$err]\n";
	    $self->throw("Problems for updating sucessful job for " . 
                   $self->input_id . " [$err]" );
    } else {
      print STDERR "Updated successful job ".$self->dbID."\n";
    }
  }
}


=head2 set_status

  Title   : set_status
  Usage   : my $status = $job->set_status
  Function: Sets the job status
  Returns : nothing
  Args    : Bio::EnsEMBL::Pipeline::Status

=cut

sub set_status {
  my ($self, $arg) = @_;
  
  $self->throw("No status input" ) unless defined($arg);
  
  
  if (!$self->adaptor) {
    $self->warn("No database connection.  Can't set status to $arg");
    return;
  }
  
  return $self->adaptor->set_status( $self, $arg );
}


=head2 current_status

  Title   : current_status
  Usage   : my $status = $job->current_status
  Function: Get/set method for the current status
  Returns : Bio::EnsEMBL::Pipeline::Status
  Args    : Bio::EnsEMBL::Pipeline::Status

=cut

sub current_status {
  my ($self, $arg) = @_;
  
  if( ! defined( $self->adaptor )) {
    return undef;
  }
  my $status;
  eval{
    $status = $self->adaptor->current_status( $self, $arg );
  };
  if($@){
    $self->throw("Failed to get status for ".$self->dbID." ".$self->input_id.
                 " ".$self->analysis->logic_name." error $@");
  }
  return $status;
}


=head2 get_all_status

  Title   : get_all_status
  Usage   : my @status = $job->get_all_status
  Function: Get all status objects associated with this job
  Returns : @Bio::EnsEMBL::Pipeline::Status
  Args    : @Bio::EnsEMBL::Pipeline::Status

=cut

sub get_all_status {
  my ($self) = @_;
  
  if( $self->adaptor ) {
    return $self->adaptor->get_all_status( $self );
  } else {
    return undef;
  }
}


=head2 get_last_status

  Title   : get_last_status
  Usage   : my @status = $job->get_all_status ($status)
  Function: Get latest status object associated with this job
  Returns : Bio::EnsEMBL::Pipeline::Status
  Args    : status string

=cut

sub get_last_status {
  my ($self) = @_;
  
  if( $self->adaptor ) {
    return $self->adaptor->get_last_status( $self );
  } else {
    return undef;
  }
}


=head2 make_filenames

=cut

sub make_filenames {
  my ($self) = @_;
  
  my $num = int(rand(10));
  
  my $dir = $self->output_dir . "/$num/";
  if( ! -e $dir ) {
    system( "mkdir $dir" );
  }

  my $stub = $self->input_id.".";
  $stub .= $self->analysis->logic_name.".";
  $stub .= int(rand(1000));

 
  $self->stdout_file($dir.$stub.".out") unless($self->stdout_file);
  $self->stderr_file($dir.$stub.".err") unless($self->stderr_file);
}


=head2 stdout_file

  Title   : stdout_file
  Usage   : my $file = $self->stdout_file
  Function: Get/set method for stdout.
  Returns : string
  Args    : string

=cut

sub stdout_file {
    my ($self, $arg) = @_;

    if (defined($arg)) {
	$self->{'_stdout_file'} = $arg;
    }
    return $self->{'_stdout_file'};
}


=head2 stderr_file

  Title   : stderr_file
  Usage   : my $file = $self->stderr_file
  Function: Get/set method for stderr.
  Returns : string
  Args    : string

=cut

sub stderr_file {
    my ($self, $arg) = @_;

    if ($arg) {
      if($arg !~ /err/){
        my ($p, $f, $l) = caller;
        $self->throw("You can't set stderr file to ".$arg." $f:$l\n");
      }
      $self->{'_stderr_file'} = $arg;
    }
    return $self->{'_stderr_file'};
}


=head2 submission_id

  Title   : submission_id
  Usage   : 
  Function: Get/set method for the submission ID
  Returns : 
  Args    : 

=cut

sub submission_id {
  my ($self, $arg) = @_;

  if(defined $arg){
    $self->{'_submission_id'} = $arg ;
  }
  return $self->{'_submission_id'};
}


sub output_dir{
 my ($self, $arg) = @_;

 if($arg){
   
   $self->{'_output_dir'} = $arg;
 }

 return $self->{'_output_dir'};
}



=head2 retry_count

  Title   : retry_count
  Usage   : 
  Function: Get/set method for the retry_count
  Returns : 
  Args    : 

=cut

sub retry_count {
  my ($self, $arg) = @_;
  if($arg) {
    $self->{'_retry_count'} = $arg; 
   }
  $self->{'_retry_count'};
}

sub can_retry{
  my ($self, $logic_name) = @_;

  $logic_name = $self->analysis->logic_name if(!$logic_name);
  if(!$BATCH_QUEUES{$logic_name}){
     $logic_name = 'default';
  }
  my $max_retry = $BATCH_QUEUES{$logic_name}{'retries'};
  if($self->retry_count <= $max_retry){
     return 1;
  }else{
    return 0;
  }
}

sub cleanup{
  my ($self, $logic_name) = @_;

  $logic_name = $self->analysis->logic_name if(!$logic_name);
  if(!$BATCH_QUEUES{$logic_name}){
     $logic_name = 'default';
  }
  if($BATCH_QUEUES{$logic_name}{'cleanup'} eq 'yes'){
    return 1;
  }else{
    return 0;
  }
}

=head2 remove

  Arg [1]   : STRING, analysis logic_name
  Function  : remove job and delete output
  Returntype: none
  Exceptions: none
  Caller    : $self
  Example   : $self->remove($self->analysis->logic_name);

=cut



sub remove {
  my $self = shift;
  my $logic_name = shift;

  if(!$BATCH_QUEUES{$logic_name}){
     $logic_name = 'default';
  } 
  if($BATCH_QUEUES{$logic_name}{'cleanup'} eq 'yes'){
    if( -e $self->stdout_file ) { unlink( $self->stdout_file ) };
    if( -e $self->stderr_file ) { unlink( $self->stderr_file ) };
  }

   if( defined $self->adaptor ) {
     $self->adaptor->remove( $self );
   }
}


=head2 set_up_queues

=cut

# scp: set up package variable for queue config
# i'm not sure this is the best way of doing this
# should ideally have this stuff in object(s)

sub set_up_queues {
    my %q;

    foreach my $queue (@$QUEUE_CONFIG) {
      my $ln = $queue->{'logic_name'};
      next unless $ln;
      delete $queue->{'logic_name'};
      while (my($k, $v) = each %$queue) {
        $q{$ln}{$k} = $v;
        $q{$ln}{'jobs'} = [];
        $q{$ln}{'last_flushed'} = undef;
        $q{$ln}{'batch_size'} ||= $DEFAULT_BATCH_SIZE;
        $q{$ln}{'queue'} ||= $DEFAULT_BATCH_QUEUE;
        $q{$ln}{'retries'} ||= $DEFAULT_RETRIES;
        $q{$ln}{'cleanup'} ||= $DEFAULT_CLEANUP;
        $q{$ln}{'runnabledb_path'} ||= $DEFAULT_RUNNABLEDB_PATH;
      }

	# a default queue for everything else
	unless (defined $q{'default'}) {
	    $q{'default'}{'batch_size'} = $DEFAULT_BATCH_SIZE;
            $q{'default'}{'retries'} ||= $DEFAULT_RETRIES;
	    $q{'default'}{'last_flushed'} = undef;
	    $q{'default'}{'queue'} = $DEFAULT_BATCH_QUEUE;
            $q{'default'}{'jobs'} = [];
            $q{'default'}{'cleanup'} = $DEFAULT_CLEANUP;
       $q{'default'}{'runnabledb_path'} ||= $DEFAULT_RUNNABLEDB_PATH;
	}
    }
  
    return %q;
}

sub execution_host {
  my ($self, $arg) = @_;
  
  if ($arg) {
    $self->{'_execution_host'} = $arg;
  }
  return $self->{'_execution_host'} || '';
}


sub temp_dir {
  my ($self, $arg) = @_;

  if($arg){
    $self->{'_temp_dir'} = $arg;
  }
  if(!$self->{'_temp_dir'}){
    $self->{'_temp_dir'} = $self->batch_q_object->temp_filename;
  }

  return $self->{'_temp_dir'} || '';
}


sub batch_q_object {
    my ($self, $arg) = @_;

    if ($arg) {
      $self->{'_batch_q_object'} = $arg;
    }
    if(!$self->{'_batch_q_object'}){
      my $object = $batch_q_module->new();
      $self->{'_batch_q_object'} = $object;
    }


    return $self->{'_batch_q_object'};
}



1;
