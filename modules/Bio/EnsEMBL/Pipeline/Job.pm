=head1 LICENSE


# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Pipeline::Job - 

=head1 SYNOPSIS


=head1 DESCRIPTION

  Stores run and status details of an analysis job

=head1 METHODS

=cut

=head1 APPENDIX

  The rest of the documentation details each of the object methods.
  Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::Job;


use warnings ;
use vars qw(@ISA $SAVE_RUNTIME_INFO);
use strict;
use File::Copy;
use Bio::EnsEMBL::Pipeline::Config::BatchQueue;
use Bio::EnsEMBL::Pipeline::Config::General;
use Bio::EnsEMBL::Analysis::Tools::Logger;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use List::Util qw( min );    # Used in flush_runs()

@ISA = qw();


# dynamically load appropriate queue manager (e.g. LSF)

my $batch_q_module = "Bio::EnsEMBL::Pipeline::BatchSubmission::$QUEUE_MANAGER";

my $file = "$batch_q_module.pm";
$file =~ s{::}{/}g;
require "$file";

# This is to catch in which state our job died. It will print the status of the running job
# and its input_id
my $SIGINT_INPUT_ID = 'NULL';
my $SIGINT_STATUS = 'NULL';
$SIG{INT} = \&_inthandler;


# BATCH_QUEUES is a package variable which stores available 'queues'.
# This allows different analysis types to be sent to different nodes
# etc.  It is keyed on analysis logic_name and has different parameters
# such as resource (e.g. node set), number of jobs to be batched
# together etc.
#
# This may be better encapsulated as a separate class, but I can't think
# at the moment how best to do this.

my %BATCH_QUEUES = set_up_queues();

sub new {
  my ( $class, @args ) = @_;

  my $self = bless {}, $class;
  my ( $p, $f, $l ) = caller;

  my ( $adaptor,    $dbID,     $submission_id,
       $input_id,   $analysis, $stdout,
       $stderr,     $rename,   $retry_count,
       $output_dir, $runner,   $number_output_dirs )
    = rearrange( [ 'ADAPTOR',       'ID',
                   'SUBMISSION_ID', 'INPUT_ID',
                   'ANALYSIS',      'STDOUT',
                   'STDERR',        'RENAME_ON_RETRY',
                   'RETRY_COUNT',   'OUTPUT_DIR',
                   'RUNNER',        'NUMBER_OUTPUT_DIRS',
                 ],
                 @args );

  if ( !defined($dbID) )          { $dbID          = -1 }
  if ( !defined($submission_id) ) { $submission_id = -1 }

  if ( !defined($input_id) ) {
    throw("Can't create a job object without an input_id");
  }
  elsif ( !defined($analysis) ) {
    throw("Can't create a job object without an analysis object");
  }

  if ( !$analysis->isa("Bio::EnsEMBL::Analysis") ) {
    throw(
        sprintf( "Analysis object [%s] is not a Bio::EnsEMBL::Analysis",
                 $analysis ) );
  }

  if (!defined($rename)) {
    $rename = $RENAME_ON_RETRY; #found in General.pm;
  }

  $self->dbID($dbID);
  $self->adaptor($adaptor);
  $self->input_id($input_id);
  $self->analysis($analysis);
  $self->stdout_file($stdout);
  $self->stderr_file($stderr);
  $self->retry_count($retry_count);
  $self->submission_id($submission_id);
  $self->number_output_dirs($number_output_dirs || 10);
  $self->rename_on_retry($rename);

  if ( !defined( $output_dir ) ) {
    if ( exists $BATCH_QUEUES{ $analysis->logic_name() } and
         exists $BATCH_QUEUES{ $analysis->logic_name() }{'output_dir'} and
         defined $BATCH_QUEUES{ $analysis->logic_name()}{'output_dir'} )
    {
      $output_dir = $BATCH_QUEUES{ $analysis->logic_name() }{'output_dir'};
    }
    else {
      $output_dir = $BATCH_QUEUES{'default'}{'output_dir'}.'/'.$analysis->logic_name();
    }
  }

  $self->output_dir($output_dir);
  $self->make_filenames() if (!$stdout and !$stderr);

    if (not defined $runner) {
    if ( exists $BATCH_QUEUES{ $analysis->logic_name() } and
         exists $BATCH_QUEUES{ $analysis->logic_name() }{'runner'} and
         defined $BATCH_QUEUES{ $analysis->logic_name()}{'runner'} )
    {
      $runner = $BATCH_QUEUES{ $analysis->logic_name() }{'runner'};
    }
    else {
      $runner = $BATCH_QUEUES{'default'}{'runner'};
    }
  }

  $self->runner($runner);

  return $self;
} ## end sub new


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
  
  warning("Bio::EnsEMBL::Pipeline::Job->create_by_analysis_input_id is deprecated ". (caller) .
               " should now call the constructor directly");
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
    throw("[$arg] is not a Bio::EnsEMBL::Analysis object" ) 
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
  
  if ( !defined $adaptor) {
    throw( "Cannot run remote without db connection" );
  }
  
  local *FILE;

  my $dbc      = $adaptor->db->dbc;
  my $host     = $dbc->host;
  my $username = $dbc->username;
  my $dbname   = $dbc->dbname;
  my $pass     = $dbc->password;
  my $port     = $dbc->port;

  # runner.pl: first look at value set in RuleManager ($RUNNER_SCRIPT)
  # then in same directory as Job.pm,
  # and fail if not found

  my $runner = $self->runner;

  if (!$runner || ! -x $runner) {
    $runner = __FILE__;
    $runner =~ s:/[^/]*$:/runner.pl:;
    my $caller = caller(0);
    throw("runner ".$runner." not found - needs to be set in ".
                 "$caller\n") unless -x $runner;
  }

  # $anal is a logic_name
 
 ANAL:
  for my $anal (@analyses) {

    my $queue = $BATCH_QUEUES{$anal};
   
    
    my @job_ids = @{$queue->{'jobs'}};
    if (!@job_ids) {
      next ANAL;
    }

    my $this_runner = $queue->{'runner'};
    #print "This runner is ".$this_runner."\n";
    $this_runner = (-x $this_runner) ? $this_runner : $runner;
    
    my $lastjob = $adaptor->fetch_by_dbID($job_ids[-1]);
    
    if ( ! $lastjob) {
      throw( "Last batch job not in db" );
    }
    $lastjob->make_filenames;

    my $run_index = $self->retry_count();


    # jhv LD_LIBRARY_PATH lsf-pre-exec hack  
    # We're running the runner.pl with /usr/local/ensembl32/bin/perl - this is a bit un-usual , but it's 
    # currently ( July 2010 ) needed to overcome some lsf bug / lsf security feature. 
    # The probelm is that if you're using the 64bit perl( to use more than 4 gigs of memory ), the job
    # will fail in the pre-exec command stage; ( you can see this with "bhist -l <jobnr>" 
    # This is because the bsub pre-exec command runs in a different environment than your normal one, 
    # and the  LD_LIBRARY_PATH is unset. so you can't use the libraries you want, ie if you you're using 
    #  /software/intel_cce_80/bin/iccvars.csh
  
    my $pre_exec_perl = "perl"; 
    if ( defined $queue->{lsf_pre_exec_perl} ) {
      $pre_exec_perl = $queue->{lsf_pre_exec_perl}; 
       
      if (! -x $pre_exec_perl ) {  
         throw("Your LSF_PRE_EXEC_PERL or DEFAULT_LSF_PRE_EXEC_PERL [$pre_exec_perl] is not executable\n");
      } 
    } 
    #my $pre_exec = '/usr/local/ensembl32/bin/perl ' . $this_runner." -check -output_dir ".$self->output_dir;
    my $pre_exec = $pre_exec_perl." ". $this_runner." -check -output_dir ".$self->output_dir;

      
    my $system_queue = $queue->{'queue'};
    my $parameters   = $queue->{'sub_args'};
    my $resources    = $queue->{'resource'};
    my $memory       = $queue->{'memory'};

    # The $system_queue, $parameters, $resources, and $memory settings
    # may be arrays where the first entry corresponds to the first run
    # of the job and the following entries are used for consequent
    # retries.
    #
    # The arrays may be shorter than the number of required retries, in
    # which case the last entry should be re-used.
    # 
    # If the setting is not an array, the code reverts to using the old
    # '_retry' setting, if present.

    if ( ref($system_queue) eq 'ARRAY' ) {
      # The 'queue' setting is an array, pick out the one corresponding
      # to the current run, or the last item if we have retried more
      # times than the number of items in the list.

      $system_queue =
        $system_queue->[ min( $#{$system_queue}, $run_index ) ];
    }
    elsif ( $run_index > 0 && defined( $queue->{'retry_queue'} ) ) {
      # We're retrying, and there's a 'retry_queue' specified, use it.
      $system_queue = $queue->{'retry_queue'};
    }
    ## else {
    ##   # Nothing to do, just use $system_queue as it is.
    ## }

    # Now do the same for $parameters, $resources, and $memory:

    # $parameters:
    if ( ref($parameters) eq 'ARRAY' ) {
      $parameters = $parameters->[ min( $#{$parameters}, $run_index ) ];
    }
    elsif ( $run_index > 0 && defined( $queue->{'retry_sub_args'} ) ) {
      $parameters = $queue->{'retry_sub_args'};
    }

    # $resources:
    if ( ref($resources) eq 'ARRAY' ) {
      $resources = $resources->[ min( $#{$resources}, $run_index ) ];
    }
    elsif ( $run_index > 0 && defined( $queue->{'retry_resource'} ) ) {
      $resources = $queue->{'retry_resource'};
    }

    # $memory:
    if ( ref($memory) eq 'ARRAY' ) {
      $memory = $memory->[ min( $#{$memory}, $run_index ) ];
    }
    elsif ( $run_index > 0 && defined( $queue->{'retry_memory'} ) ) {
      $memory = $queue->{'retry_memory'};
    }

    my $batch_job =
      $batch_q_module->new( -STDOUT     => $lastjob->stdout_file(),
                            -PARAMETERS => $parameters,
                            -PRE_EXEC   => $pre_exec,
                            -QUEUE      => $system_queue,
                            -JOBNAME    => $dbname . ':' . $anal,
                            -NODES      => $queue->{'nodes'},
                            -RESOURCE   => $resources,
                            -MEMORY     => $memory, );

    my $cmd;

    if (!$self->cleanup) {
      $batch_job->stderr_file($lastjob->stderr_file);
    }
    

    # check if the password has been defined, and write the
    # "connect" command line accordingly otherwise -pass gets the
    # first job id as password, instead of remaining undef 
   
    my $exec_perl ="perl" ; 
 
    if ( defined $queue->{lsf_perl} ) {
      $exec_perl = $queue->{lsf_perl};  
      if (! -x $exec_perl ) {  
         warning("using lsf_perl : $exec_perl\n");
         throw("Your LSF_PERL or DEFAULT_LSF_PERL [$exec_perl] is not executable\n");
      } 
    } 

    if ($pass) {  
      $cmd = "$exec_perl $runner  -dbhost $host -dbuser $username -dbname $dbname -dbpass $pass -dbport $port";
      #$cmd = "/usr/local/ensembl64/bin/perl $runner  -dbhost $host -dbuser $username -dbname $dbname -dbpass $pass -dbport $port";
    } else {
      $cmd = "$exec_perl $runner  -dbhost $host -dbuser $username -dbname $dbname -dbport $port";
    }
  
 
    $cmd .= " -output_dir ".$self->output_dir;
    $cmd .= " -queue_manager $QUEUE_MANAGER  " ;
    if ($self->cleanup) {
      $cmd .= " -cleanup "
    }
    $cmd .= " @job_ids";
    #print "Job.pm-cmd : $cmd\n";   
    $batch_job->construct_command_line($cmd);
    return "DB overloaded" if ($main::CHECK_DBLOAD and $batch_job->is_db_overloaded($queue->{load_pending_cost} || 10));

    eval {
      # SMJS LSF Specific for debugging
      # print STDERR "Submitting: ", $batch_job->command, "\n"; # only if you run jobs with LSF /batch submission; does not work if you run jobs locally
      $batch_job->open_command_line();
    };

    if ($@) {
      print STDERR "Couldnt batch submit @job_ids \n[$@]\n";
      print STDERR "Using ".$batch_job->command."\n";
      foreach my $job_id (@job_ids) {
        my $job = $adaptor->fetch_by_dbID($job_id);
        $job->set_status( "FAILED" );
      }
    } else {
      #print STDERR "have submitted ".@job_ids." jobs with ".$batch_job->id
      #."\n";
      my @jobs = $adaptor->fetch_by_dbID_list(@job_ids);
      foreach my $job (@jobs) {

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
  
  if (!exists($BATCH_QUEUES{$self->analysis->logic_name})) {
    $queue = 'default';
  } else {
    $queue = $self->analysis->logic_name;
  }
  
  push @{$BATCH_QUEUES{$queue}{'jobs'}}, $self->dbID;
 
  if (scalar(@{$BATCH_QUEUES{$queue}{'jobs'}}) >=
               $BATCH_QUEUES{$queue}{'batch_size'}) {

    $self->flush_runs($self->adaptor, $queue);
  } else {
    # print STDERR "Not running ".$self->dbID." yet as batch size is ".
    #              $BATCH_QUEUES{$queue}{'batch_size'} . " and there are ".
    #               scalar(@{$BATCH_QUEUES{$queue}{'jobs'}})." jobs\n";
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
 
  #print STDERR "Running locally " . $self->stdout_file . " " 
  #. $self->stderr_file . "\n"; 

  local *STDOUT;
  local *STDERR;

  $self->make_filenames;

  if ( ! open(STDOUT, ">".$self->stdout_file)) {
    $self->set_status( "FAILED" );
    return;
  }
        
  if( ! open(STDERR, ">".$self->stderr_file)) {
    $self->set_status( "FAILED" );
    return;
  }

  # print STDERR "Running in LSF\n"; 

  $self->run_module();
  if ($self->current_status->status eq "SUCCESSFUL"){
    $self->adaptor->remove( $self );
    if($self->cleanup){
      unlink $self->stderr_file if(-e $self->stderr_file);
      unlink $self->stdout_file if(-e $self->stdout_file);
    }
  }
}


# question, when to submit the success report to the db?
# we have to parse the output of LSF anyway....
sub run_module {
  my $self = shift;

  $SIGINT_INPUT_ID = $self->input_id;
  my $module     = $self->analysis->module;
  my $hash_key   = $self->analysis->logic_name;
  my $rdb;
  my ($err, $res);

  #print STDERR "Running ".$module." with ".$self."\n";

  if (!exists($BATCH_QUEUES{$hash_key})) {
    $hash_key = 'default';
  }

  my $runnable_db_path = $BATCH_QUEUES{$hash_key}{runnabledb_path};
  my $verbosity =  $BATCH_QUEUES{$hash_key}{verbosity};

  my $perl_path;
  #print STDERR "Getting ".$hash_key." batchqueue value\n";

  if ($module =~ /::/) {
    #print STDERR "Module contains path info already\n";
    $module =~ s/::/\//g;
    $perl_path = $module;
  } elsif($runnable_db_path) {
    $perl_path = $runnable_db_path."/".$module;
  }else{
    $perl_path = $module;
  }

  my $current_verbosity = logger_verbosity;
  my $start = 0;
  my $end = 0;
  logger_verbosity($verbosity);

  #print STDERR "have perlpath ".$perl_path."\n";
 STATUS: 
  { 
    eval {
      require $perl_path.".pm";
      $perl_path =~ s/\//::/g;
      $rdb = $perl_path->new( -analysis => $self->analysis,
                              -input_id => $self->input_id,
                              -db => $self->adaptor->db,
                              -verbosity => $verbosity,
                            );
    };
    
    if ($err = $@) {
      print (STDERR "CREATE: Lost the will to live Error\n");
      $self->set_status( "FAILED" );
      throw( "Problems creating runnable $module for " . 
                    $self->input_id . " [$err]\n");
    }
    
    # "READING"
    eval {   
      $SIGINT_STATUS = 'READING';
      $self->set_status( "READING" );
      $res = $rdb->fetch_input;
    };
    if ($err = $@) {
      $self->set_status( "FAILED" );
      print (STDERR "READING: Lost the will to live Error\n");
      throw( "Problems with $module fetching input for " . 
                    $self->input_id . " [$err]\n");
    }
    
    if ($rdb->input_is_void) {
      $SIGINT_STATUS = 'VOID';
      $self->set_status( "VOID" );
    } else {
      # "RUNNING"
      eval {
        $SIGINT_STATUS = 'RUNNING';
        $self->set_status( "RUNNING" );
        $start = $self->get_last_status->created();
        $rdb->db->dbc->disconnect_when_inactive(1); 
        $rdb->run;
        $rdb->db->dbc->disconnect_when_inactive(0); 
      };
      if ($err = $@) {
        print STDERR $@ . "\n";
        
        if (my $err_state = $rdb->failing_job_status) {
          $self->set_status( $err_state );
        } else {
          $self->set_status( "FAILED" ); # default to just failed 
          #these jobs get retried
        }

        print (STDERR "RUNNING: Lost the will to live Error\n");
        throw("Problems running $module for " . 
                     $self->input_id . " [$err]\n");
      }
      
      # "WRITING"
      eval {
        $SIGINT_STATUS = 'WRITING';
        $self->set_status( "WRITING" );
        $rdb->write_output;

        if ($rdb->can('db_version_searched')) {
          my $new_db_version = $rdb->db_version_searched();
          my $analysis       = $self->analysis();
          my $old_db_version = $analysis->db_version();

          $analysis->db_version($new_db_version);

          # where is the analysisAdaptor??
          # $self->adaptor->get_AnalysisAdaptor->store($analysis);
          # if $new_db_version gt $old_db_version;
        } else {
          $SAVE_RUNTIME_INFO = 0;
        }

        $SIGINT_STATUS =  'SUCCESSFUL';
        $self->set_status("SUCCESSFUL");
        $end = $self->get_last_status->created();
      }; 
      if ($err = $@) {
        $self->set_status( "FAILED" );
        print (STDERR "WRITING: Lost the will to live Error\n");
        throw("Problems for $module writing output for " . 
        $self->input_id . " [$err]" );
      }
    }
  }    

  logger_verbosity($current_verbosity);
  
  # update job in StateInfoContainer
  eval {
    my $sic = $self->adaptor->db->get_StateInfoContainer;
    $SAVE_RUNTIME_INFO = $end - $start;
    $sic->store_input_id_analysis(
                                  $self->input_id,
                                  $self->analysis,
                                  $self->execution_host,
                                  $SAVE_RUNTIME_INFO
                                 );
  };
  if ($err = $@) {
    my $error_msg = "Job finished successfully, but could not be ".
      "recorded as finished.  Job : [" . $self->input_id . "]\n[$err]";
    eval {
      $self->set_status("FAIL_NO_RETRY");
    };
    $error_msg .= ("(And furthermore) Encountered an error in updating the job to status failed_no_retry.\n[$@]") if $@;
    throw($error_msg);
  } else {
    print STDERR "Updated successful job ".$self->dbID."\n";
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
  
  throw("No status input" ) unless defined($arg);
  
  if (!$self->adaptor) {
    warning("No database connection.  Can't set status to $arg");
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
  
  if ( ! defined( $self->adaptor)) {
    return undef;
  }

  my $status;
  eval{
    $status = $self->adaptor->current_status( $self, $arg );
      
  };
  if ($@) {
    throw("Failed to get status for ".$self->dbID." ".$self->input_id.
                 " ".$self->analysis->logic_name." error $@");
  }

  #if($status->status eq 'SUCCESSFUL'){
  #  my ($p, $f, $l) = caller;
  #  print STDERR $f.":".$l."\n";
  #}
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
  
  if ($self->adaptor) {
    return $self->adaptor->get_all_status( $self );
  } else {
    return undef;
  }
}


=head2 get_last_status

  Title   : get_last_status
  Usage   : my @status = $job->get_last_status ($status)
  Function: Get latest status object associated with this job
  Returns : Bio::EnsEMBL::Pipeline::Status
  Args    : status string

=cut

sub get_last_status {
  my ($self) = @_;
  
  if ($self->adaptor) {
    return $self->adaptor->get_last_status( $self );
  } else {
    return undef;
  }
}


=head2 make_filenames

=cut 

sub make_filenames {
  my ($self) = @_;
 
  my $output_dir_number ; 
  if ( $self->number_output_dirs ) {   
   $output_dir_number =$self->number_output_dirs ; ;  
  }  else { 
   $output_dir_number =10 ;  
  } 
  my $num = int(rand($output_dir_number));
  
  my $dir = $self->output_dir . "/$num/"; 
  #print "MY dir ".$dir."\n"; 
  
  if ( ! -e $self->output_dir) {
    warning("Your output-directory " . $self->output_dir . " does not exist - i'll create it now.\n") ;  
    eval{
      system("mkdir -p " . $self->output_dir ) ;  
    };
    if($@){
      throw("Failed to make dir " . $self->output_dir . "$@" );
    }
  }  

  if ( ! -e $dir) {
    eval{
      my $command  = "mkdir $dir";
      system( $command );
    };
    if($@){
      throw("Failed to make dir ".$dir." $@");
    }
  }

  my $fname = $self->input_id;
  $fname =~ s/\/+$//;
  if ($fname =~ /\//) {
  my ($suffix) = $fname =~ /\/([^\/]+)$/;
  $fname = $suffix;
  }
  $fname .= "." . $self->analysis->logic_name;
  $fname .= "." . int(rand(1000));
  $self->stdout_file($dir.$fname.'.out') unless $self->stdout_file;
  $self->stderr_file($dir.$fname.'.err') unless $self->stderr_file;
  
  if ($self->retry_count > 0 and $self->can_retry($self->analysis->logic_name) and $self->rename_on_retry) {
      if (-e $self->stdout_file) {
          move($self->stdout_file, $self->stdout_file.'.retry'.$self->retry_count) || throw('Failed to copy log file '.$self->stdout_file);
      }
      if (-e $self->stderr_file) {
          move($self->stderr_file, $self->stderr_file.'.retry'.$self->retry_count) || throw('Failed to copy log file '.$self->stderr_file);
      }
    $self->stdout_file($dir.$fname.$self->retry_count.'.out');
    $self->stderr_file($dir.$fname.$self->retry_count.'.err');
  }
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

     # modification for find_false_introns to shorten long file-names to deal with long input_ids 
     #my @tmp = split /\:/,$arg ; 
     #if ( @tmp > 7 ) {
     #   @tmp  = splice(@tmp, 0, 7 );
     #   $self->{'_stderr_file'} = join(":",@tmp).".out" ;
     #} 
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
    $self->{'_stderr_file'} = $arg;

     # modification for find_false_introns to shorten long file-names to deal with long input_ids 
     #my @tmp = split /\:/,$arg ; 
     #if ( @tmp > 7 ) {
     #   @tmp  = splice(@tmp, 0, 7 );
     #   $self->{'_stderr_file'} = join(":",@tmp).".err" ;
     #} 

    if ($arg !~ /err/){
      my ($p, $f, $l) = caller;
      throw("You can't set stderr file to ".$arg." $f:$l\n");
    }
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

  if (defined $arg) {
    $self->{'_submission_id'} = $arg ;
  }
  return $self->{'_submission_id'};
}


sub output_dir{
  my ($self, $arg) = @_;

  if ($arg) {
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
  if ($arg) {
    $self->{'_retry_count'} = $arg; 
  }
  $self->{'_retry_count'} || 0;
}

sub can_retry{
  my ($self, $logic_name) = @_;

  $logic_name = $self->analysis->logic_name if(!$logic_name);
  
  if (!exists($BATCH_QUEUES{$logic_name})) {
    $logic_name = 'default';
  }
  my $max_retry = $BATCH_QUEUES{$logic_name}{'retries'};
  if (!$self->retry_count) {
    $BATCH_QUEUES{$logic_name}{'batch_size'} = $BATCH_QUEUES{$logic_name}{'retry_batch_size'} if (exists $BATCH_QUEUES{$logic_name}{'retry_batch_size'});
    return 1;
  }
  if ($self->retry_count && $self->retry_count <= $max_retry) {
    $BATCH_QUEUES{$logic_name}{'batch_size'} = $BATCH_QUEUES{$logic_name}{'retry_batch_size'} if (exists $BATCH_QUEUES{$logic_name}{'retry_batch_size'});
    return 1;
  } else {
    return 0;
  }
}


sub cleanup{
  my ($self, $logic_name) = @_;

  $logic_name = $self->analysis->logic_name if(!$logic_name);

  if (!exists($BATCH_QUEUES{$logic_name})) {
    $logic_name = 'default';
  }

  if ($BATCH_QUEUES{$logic_name}{'cleanup'} eq 'yes') {
    return 1;
  } else {
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

  if (!exists($BATCH_QUEUES{$logic_name})) {
     $logic_name = 'default';
  } 
  if ($BATCH_QUEUES{$logic_name}{'cleanup'} eq 'yes') {
    if ( -e $self->stdout_file) { unlink( $self->stdout_file ) };
    if ( -e $self->stderr_file) { unlink( $self->stderr_file ) };
  }

   if (defined $self->adaptor) {
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
    my $ln = $queue->{logic_name};
    #print "LOOKING AT ".$ln."\n";
    next unless $ln;

    while (my($k, $v) = each %$queue) {
      next if $k eq 'logic_name';
      #print "PREPARING ".$k." ".$v." entry \n";
      $q{$ln}{$k}     = $v;
    }
    $q{$ln}{jobs} = [];
    $q{$ln}{last_flushed} = undef;
    $q{$ln}{batch_size} = $DEFAULT_BATCH_SIZE if !defined($q{$ln}{batch_size});
    $q{$ln}{queue}     ||= $DEFAULT_BATCH_QUEUE;
    $q{$ln}{resource}  ||= $DEFAULT_RESOURCE;
                #Allow 0 retries
    $q{$ln}{retries} = $DEFAULT_RETRIES if !defined($q{$ln}{retries});
    $q{$ln}{cleanup}                            ||= $DEFAULT_CLEANUP;
    $q{$ln}{runnabledb_path}||= $DEFAULT_RUNNABLEDB_PATH;
    #$q{$ln}{output_dir}                ||= $DEFAULT_OUTPUT_DIR;
    $q{$ln}{runner}                             ||= $DEFAULT_RUNNER;
    $q{$ln}{verbosity}                  ||= $DEFAULT_VERBOSITY;
    $q{$ln}{sub_args}                   ||= $DEFAULT_SUB_ARGS;
    $q{$ln}{retry_queue}                ||= $DEFAULT_RETRY_QUEUE;
    $q{$ln}{retry_resource} ||= $DEFAULT_RETRY_RESOURCE;
    $q{$ln}{retry_sub_args} ||= $DEFAULT_RETRY_SUB_ARGS; 

    $q{$ln}{lsf_pre_exec_perl} ||= $DEFAULT_LSF_PRE_EXEC_PERL;
    $q{$ln}{lsf_perl} ||= $DEFAULT_LSF_PERL; 
    #print "HAVE SUB ARGS ".$q{$ln}{sub_args}."\n";
  }

  # a default queue for everything else
  if ( ! exists($q{default})) {
    $q{default}{jobs} = [];
    $q{default}{last_flushed} = undef;
  }
  # Need these set, do the ||= thing (but check for numbers that are 0) 
  $q{default}{batch_size} = $DEFAULT_BATCH_SIZE if !defined($q{default}{batch_size});
  $q{default}{queue}     ||= $DEFAULT_BATCH_QUEUE;
  $q{default}{resource}  ||= $DEFAULT_RESOURCE;
  $q{default}{retries} = $DEFAULT_RETRIES if !defined($q{default}{retries});
  $q{default}{cleanup}    ||= $DEFAULT_CLEANUP;
  $q{default}{runnabledb_path}||= $DEFAULT_RUNNABLEDB_PATH;
  $q{default}{output_dir} ||= $DEFAULT_OUTPUT_DIR;
  $q{default}{runner}     ||= $DEFAULT_RUNNER;
  $q{default}{verbosity}  ||= $DEFAULT_VERBOSITY;
  $q{default}{sub_args}   ||= $DEFAULT_SUB_ARGS;
  $q{default}{retry_queue}||= $DEFAULT_RETRY_QUEUE;
  $q{default}{retry_resource} ||= $DEFAULT_RETRY_RESOURCE;
  $q{default}{retry_sub_args} ||= $DEFAULT_RETRY_SUB_ARGS; 
  $q{default}{lsf_pre_exec_perl} ||= $DEFAULT_LSF_PRE_EXEC_PERL;
  $q{default}{lsf_perl} ||= $DEFAULT_LSF_PERL; 
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

  if ($arg) {
    $self->{'_temp_dir'} = $arg;
  }
  if (!$self->{'_temp_dir'}) {
    $self->{'_temp_dir'} = $self->batch_q_object->temp_filename;
  }

  return $self->{'_temp_dir'} || '';
}


sub batch_q_object {
  my ($self, $arg) = @_;

  if ($arg) {
    $self->{'_batch_q_object'} = $arg;
  }

  if (!$self->{'_batch_q_object'}) {
    my $object = $batch_q_module->new();
    $self->{'_batch_q_object'} = $object;
  }

  return $self->{'_batch_q_object'};
}


sub number_output_dirs {  
   my ($self, $arg ) = @_ ; 
   $self->{number_output_dirs} = $arg if $arg ; 
   return $self->{number_output_dirs}; 
} 

=head2 rename_on_retry 

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

sub _inthandler {
    print STDERR "ZUT! The job died while $SIGINT_STATUS for $SIGINT_INPUT_ID\n";
}

1;
