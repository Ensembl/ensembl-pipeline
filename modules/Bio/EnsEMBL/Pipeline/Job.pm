
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


use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object;

use Bio::EnsEMBL::Pipeline::Config::BatchQueue;

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
    #print STDERR "calling Job->new from ".(caller)."\n";
    #print STDERR "@args\n";
    my ($adaptor,$dbID,$submission_id,$input_id,$analysis,$stdout,$stderr,$retry_count, $output_dir) 
	= $self->_rearrange([qw(ADAPTOR
				ID
				SUBMISSION_ID
				INPUT_ID
				ANALYSIS
				STDOUT
				STDERR
				RETRY_COUNT
				OUTPUT_DIR
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
      my $dir = $BATCH_QUEUES{$analysis->logic_name}->{'output_dir'} || $DEFAULT_OUTPUT_DIR;
      $self->throw("need an output directory passed in from RuleManager or from Config/BatchQueue $!") unless($dir);
      $self->output_dir($dir);
      $self->make_filenames;
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

    if (defined($arg)) {
	$self->{'_dbID'} = $arg;
    }
    return $self->{'_dbID'};

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

    if (defined($arg)) {
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

  local *SUB;
  local *FILE;

  my $db       = $adaptor->db;
  my $host     = $db->host;
  my $username = $db->username;
  my $dbname   = $db->dbname;
  my $pass     = $db->password;

  # runner.pl: first look at value set in RuleManager ($RUNNER_SCRIPT)
  # then in same directory as Job.pm,
  # and fail if not found

  my $runner = $::RUNNER_SCRIPT || undef;

  unless (-x $runner) {
    $runner = __FILE__;
    $runner =~ s:/[^/]*$:/runner.pl:;
    my $caller = caller(0);
    $self->throw("runner not found - needs to be set in $caller\n") unless -x $runner;
  }

  # $anal is a logic_name
  for my $anal (@analyses) {

    my $queue = $BATCH_QUEUES{$anal};
   
    my @job_ids = @{$queue->{'jobs'}};
    next unless @job_ids;

    my $this_runner = $queue->{'runner'};
    $this_runner = (-x $this_runner) ? $this_runner : $runner;
 
    my $lastjob = $adaptor->fetch_by_dbID($job_ids[-1]);

    if( ! defined $lastjob ) {
      $self->throw( "Last batch job not in db" );
    }
  
    my $pre_exec = $this_runner." -check -output_dir ".$self->output_dir;
    my $batch_job = $batch_q_module->new(
	-STDOUT     => $lastjob->stdout_file,
	-STDERR     => $lastjob->stderr_file,
	-PARAMETERS => $queue->{'sub_args'},
	-PRE_EXEC   => $pre_exec,
	-QUEUE      => $queue->{'queue'},
	-JOBNAME    => $dbname . '/' . $anal,
	-NODES      => $queue->{'nodes'},
	-RESOURCE   => $queue->{'resource'}
    );
    my $cmd;
  
    

    # check if the password has been defined, and write the
    # "connect" command line accordingly (otherwise -pass gets the
    # first job id as password, instead of remaining undef)
    if ($pass) {
      $cmd = $runner." -host $host -dbuser $username -dbname $dbname -pass $pass";
    }
    else {
      $cmd = $runner." -host $host -dbuser $username -dbname $dbname";
    }
    $cmd .= " -output_dir ".$self->output_dir;
    $cmd .= " @job_ids";

    $batch_job->construct_command_line($cmd);
    $batch_job->open_command_line();
    if( ! defined $batch_job->id ) {
      print STDERR ( "Couldnt batch submit @job_ids" );
      foreach my $job_id (@job_ids) {
        my $job = $adaptor->fetch_by_dbID($job_id);
        $job->set_status( "FAILED" );
      }
    } else {
    
      my @jobs = $adaptor->fetch_by_dbID_list(@job_ids);
      foreach my $job (@jobs) {
        if( $job->retry_count > 0 ) {
          for ( $job->stdout_file, $job->stderr_file ) {
            open( FILE, ">".$_ ); close( FILE );
          }
        }
	$job->submission_id( $batch_job->id );
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

  if ($BATCH_QUEUES{$self->analysis->logic_name}) {
      $queue = $self->analysis->logic_name;
  }
  else {
      $queue = 'default';
  }
  
  push @{$BATCH_QUEUES{$queue}{'jobs'}}, $self->dbID;

  if (scalar(@{$BATCH_QUEUES{$queue}{'jobs'}}) >=
               $BATCH_QUEUES{$queue}{'batch_size'}) {

    $self->flush_runs($self->adaptor, $queue);
  }
}


=head2 runLocally
=head2 runRemote( boolean withDB, queue )
=head2 run_module

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
  my $rdb;
  my ($err, $res);
  my $autoupdate = $AUTO_JOB_UPDATE;

  STATUS:
  { 
    # "SUBMITTED"
    eval {
      if( $module =~ /::/ ) {
          $module =~ s/::/\//g;
          require "${module}.pm";
          $module =~ s/\//::/g;

          $rdb = "${module}"->new
              ( -analysis => $self->analysis,
                -input_id => $self->input_id,
                -db => $self->adaptor->db );
      } else {
          require "Bio/EnsEMBL/Pipeline/RunnableDB/${module}.pm";
          $module =~ s/\//::/g;
          $rdb = "Bio::EnsEMBL::Pipeline::RunnableDB::${module}"->new
              ( -analysis => $self->analysis,
                -input_id => $self->input_id,
                -db => $self->adaptor->db );
      }
    };
    
    if ($err = $@) {
      print (STDERR "CREATE: Lost the will to live Error\n");
      $self->set_status( "FAILED" );
      $self->throw( "Problems creating runnable $module for " . $self->input_id . " [$err]\n");
    }

    # "READING"
    eval {   
      $self->set_status( "READING" );
      $res = $rdb->fetch_input;
    };
    if ($err = $@) {
      $self->set_status( "FAILED" );
      print (STDERR "READING: Lost the will to live Error\n");
      die "Problems with $module fetching input for " . $self->input_id . " [$err]\n";
    }
    if ($res eq "") {
    }
    if ($rdb->input_is_void) {
      $self->set_status( "VOID" );
      return;
    }

    # "RUNNING"
    eval {
      $self->set_status( "RUNNING" );
      $rdb->run;
    };
    if ($err = $@) {
      $self->set_status( "FAILED" );
      print (STDERR "RUNNING: Lost the will to live Error\n");
      die "Problems running $module for " . $self->input_id . " [$err]\n";
    }

    # "WRITING"
    eval {
      $self->set_status( "WRITING" );
      $rdb->write_output;
      $self->set_status( "SUCCESSFUL" );
    }; 
    if ($err = $@) {
      $self->set_status( "FAILED" );
      print (STDERR "WRITING: Lost the will to live Error\n");
      die "Problems for $module writing output for " . $self->input_id . " [$err]" ;
    }
  } # STATUS

  # update job in StateInfoContainer
  if ($autoupdate) {
    eval {
      my $sic = $self->adaptor->db->get_StateInfoContainer;
      $sic->store_input_id_analysis(
        $self->input_id,
        $self->analysis
      );
    };
    if ($err = $@) {
      print STDERR "Error updating successful job ".$self->dbID ."[$err]\n";
    }
    else {
      print STDERR "Updated successful job ".$self->dbID."\n";
      eval {
        $self->remove;
      };
      if ($err = $@) {
         print STDERR "Error deleting job ".$self->dbID." [$err]\n";
      }
      else {
         print STDERR "Deleted job ".$self->dbID."\n";
      }
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
  
  
  if (!(defined($self->adaptor))) {
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
 
  return $self->adaptor->current_status( $self, $arg );
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


=head 2 make_filenames

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

    if (defined($arg)) {
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

 my ($p, $f, $l) = caller;
 if($arg){
   #print STDERR $f." ".$l." is calling output_dir with ".$arg."\n";
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
  if(defined $arg) {
    #print STDERR "setting retry count to ".$arg."\n"; 
    $self->{'_retry_count'} = $arg; 
   }
  #print STDERR "retry count is = ".$self->{'_retry_count'}."\n";;
  $self->{'_retry_count'};
}


=head2 remove
 
=cut

sub remove {
  my $self = shift;
  
  if( -e $self->stdout_file ) { unlink( $self->stdout_file ) };
  if( -e $self->stderr_file ) { unlink( $self->stderr_file ) };
  

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
	}

	# a default queue for everything else
	unless (defined $q{'default'}) {
	    $q{'default'}{'batch_size'} = $DEFAULT_BATCH_SIZE;
            $q{'default'}{'retries'} ||= $DEFAULT_RETRIES;
	    $q{'default'}{'last_flushed'} = undef;
	    $q{'default'}{'queue'} ||= $DEFAULT_BATCH_QUEUE;
            $q{'default'}{'jobs'} = [];
	}
    }
    return %q;
}

1;
