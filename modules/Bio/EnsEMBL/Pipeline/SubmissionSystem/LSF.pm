use strict;
use warnings;

package Bio::EnsEMBL::Pipeline::SubmissionSystem::LSF;

use vars qw(@ISA);

use Bio::EnsEMBL::Pipeline::SubmissionSystem;
use Bio::EnsEMBL::Pipeline::Job;

@ISA = qw(Bio::EnsEMBL::Pipeline::SubmissionSystem);

# LSF_MAX_BATCH_SIZE could go as high as 65536, but we need a way to seperate
# job_array output files into seperate dirs to ease filesystem burden
use constant LSF_MAX_BATCH_SIZE => 1000;
use constant MAX_INT => 2147483647; # (2^31)-1 max 32 bit signed int

=head2 kill

  Arg [1]    : Bio::EnsEMBL::Pipeline::Job
  Example    : $lsf_sub_system->kill($job);
  Description: kills a job that has been submitted already
  Returntype : none
  Exceptions : none
  Caller     : general

=cut

sub kill {
  my $self = shift;
  my $job  = shift;

	my $job_id    = $job->dbID;
  my $array_idx = $job->array_index;
	my $taskname  = $job->taskname;
	my $sub_id    = $job->submission_id;
	
	if(!$sub_id) {
		$self->warn("Cannot kill job [$job_id] - does not have submission id");
		return;
	}

	my $arg;
	if($array_idx) {
		$arg = '"' . $job_id . "[$array_idx]" . '"';
	} else {
		$arg = $job_id;
	}

	my $path = $self->get_Config->get_parameter('LSF', 'path');

	my $rc = system($path.'bkill', $arg);

	if($rc & 0xffff) {
		$self->warn("bkill returned non-zero exit status $!");
		return;
	}

	$job->update_status('KILLED');
}



=head2 flush

  Arg [1]    : (optional) $taskname
               Gives a hint as to which queues should be flushed.  If no
               taskname is supplied then it is implied that all internal
               queues should be flushed.
  Example    :
  Description: Flushes this submission systems queues
  Returntype : none
  Exceptions : none
  Caller     : Bio::EnsEMBL::Pipeline::PipelineManager

=cut

sub flush {
	my $self = shift;
	my $taskname = shift;

  #print STDERR "FLUSHING taskqueue $taskname\n";

	#keep track of the number of submissions to guarentee unique names
  $self->{'sub_count'} = 0 if(!defined($self->{'sub_count'}));
	$self->{'sub_count'}++;

	#don't want to exceed max int size
	$self->{'sub_count'} = 1 if($self->{'sub_count'} == MAX_INT);

	my @tasknames;

	#if no task name specified, flush all task queues
	if($taskname) {
		@tasknames = ($taskname);
	} else {
		@tasknames = keys %{$self->_task_queue()};
	}

	my $config = $self->get_Config();
	my $job_adaptor = $config->get_DBAdaptor->get_JobAdaptor();

  #construct command line db options for runner script
  my $db = $config->get_DBAdaptor();
  $self->throw('could not obtain pipeline database connection') if(!$db);
  my $dbargs = join(' ',
                    '-host',   $db->host,
                    '-port',   $db->port,
                    '-dbname', $db->dbname,
                    '-dbuser', $db->username,
                    '-pass',   $db->password);

  my $runner = $config->get_parameter('LSF', 'runner');

	foreach $taskname (@tasknames) {

		#extract the queue name and free form parameters from the 'where'
		# the first part of the 'where' string is LSF which we already know
		my $where = $config->get_parameter($taskname, 'where');
		my ($queue, $other_parms);
		(undef, $queue, $other_parms) = split(':', $where, 3);

    $other_parms ||= '';
		
		$queue ||
			$self->throw("Could not determine LSF queue for task [$taskname]");

		my @jobs = @{$self->_task_queue->{$taskname} || []};

    #print STDERR "GOT " . scalar(@jobs) . " jobs from taskqueue $taskname\n";

		next if(@jobs == 0);

    #empty the queue
    $self->_task_queue->{$taskname} = [];

		#submission needs to be uniquely named
		my $lsf_job_name = join('_', $taskname, time(),$self->{'sub_count'});

		my $dir_prefix = $self->_dir_prefix($taskname);

		#add the preexec to the arguments
		my @args = ('-E', "\"$runner -check -output_dir $dir_prefix\"");
		#add the queue name
		push @args, ('-q', $queue);
		

    my $command;

		#
		# If there is only a single job submit it normally
		# otherwise submit it as a job array
		#
		if(@jobs == 1) {
			my ($job) = @jobs;

			$job->job_name($lsf_job_name);
      my $stdout = ($dir_prefix) ?
        "$dir_prefix$lsf_job_name.out" : '/dev/null';
      my $stderr = ($dir_prefix) ?
        "$dir_prefix$lsf_job_name.err" : '/dev/null';
			$job->stdout_file($stdout);
			$job->stderr_file($stderr);
			$job->array_index(undef);

			$job_adaptor->update($job);
			$job->set_current_status('SUBMITTED');

			#add the output dirs to the stdout list
			push @args, ('-o', $stdout);
			push @args, ('-e', $stderr);
			push @args, ('-J', $lsf_job_name);

      $command = "$runner -jobname $lsf_job_name $dbargs";
		} else {
			my $array_index = 0;

			foreach my $job (@jobs) {
        $array_index++;

        my $job_stdout = ($dir_prefix) ? 
          "$dir_prefix${lsf_job_name}_${array_index}.out" : '/dev/null';
        my $job_stderr = ($dir_prefix) ?
          "$dir_prefix${lsf_job_name}_${array_index}.err" : '/dev/null';

        $job->stdout_file($job_stdout);
        $job->stderr_file($job_stderr);
				$job->job_name($lsf_job_name);
				$job->array_index($array_index);

				# write entry into job table
				$job_adaptor->update($job);

				# set each job status to submitted
				$job->set_current_status('SUBMITTED');
			}	
			#add the output dirs to the stdout list
      my $stdout = ($dir_prefix) ? 
        "$dir_prefix$lsf_job_name".'_%I.out' : '/dev/null';
      my $stderr = ($dir_prefix) ?
        "$dir_prefix$lsf_job_name".'_%I.err' : '/dev/null';

			push @args, ('-o', $stdout);
			push @args, ('-e', $stderr);
			
			#add the job name and array index
			push @args, ('-J', '"'.$lsf_job_name.'[1-'.scalar(@jobs).']"');
      $command = "$runner -jobname $lsf_job_name -index %I $dbargs";
		}

		#execute the bsub to submit the job or job_array
		#would rather use system() since it doesn't require a fork, but
		#need to get job_id out of stdout :(
		my $bsub = 'bsub ' . join(' ', @args, $other_parms, $command);

    print STDERR "LSF: EXECUTING COMMAND:\n$bsub\n";
 #   my $sub_id = int(rand(100000));

		open(SUB, $bsub." 2>&1 |") or
			$self->throw("could not execute command [$bsub]");
		my $sub_id;
    my $err;
		while(<SUB>) {
			if (/Job <(\d+)>/) {
				$sub_id = $1;
			} else {
        $err .= $_;
      } 
		}
		close(SUB);

    #if the bsub failed, set the status of all the jobs to failed
    #otherwise update them with submission ids
    if(!$sub_id) {
      $self->warn("bsub command [$bsub] failed: [$err]");
      foreach my $job (@jobs) {
        $job->set_current_status('FAILED');
      }
    } else {

      #update the job table so that the submission id is also stored
      foreach my $job(@jobs) {
        $job->submission_id($sub_id);
        $job_adaptor->update($job);
      }
    }
  }
	
	return;
}


=head2 submit

  Arg [1]    : Bio::EnsEMBL::Pipeline::Job::LSFJob $job
	Arg [2]    :
  Example    : $lsf->submit($job);
  Description: Submits an LSF job.
  Returntype : none
  Exceptions : thrown if
  Caller     : PipelineManager

=cut

sub submit {
	my $self = shift;
	my $job  = shift;

	unless(ref($job) && $job->isa('Bio::EnsEMBL::Pipeline::Job')) {
		$self->throw('expected Bio::EnsEMBL::Pipeline::Job argument');
	}

	# retrieve batch size from config
	my $config = $self->get_Config();
	my $taskname = $job->taskname();
	my $batch_size = $config->get_parameter($taskname, 'batchsize');

	$batch_size || $self->throw('Could not determine batch size for task ' .
															"[$taskname] from config");


	if($batch_size > LSF_MAX_BATCH_SIZE) {
		$self->warn("Maximum batch size [" . LSF_MAX_BATCH_SIZE .
								"] exceeded for LSF task [$taskname]");
		$batch_size = LSF_MAX_BATCH_SIZE;
	}

	#place the job into a queue specific to each task and module
	$self->_task_queue()->{$taskname} ||= [];
  #print STDERR "SUBMITTED job to taskqueue $taskname\n";
	push(@{$self->_task_queue()->{$taskname}}, $job);


	#
	# if the queue is full submit all of the jobs in it
	#
	if(@{$self->_task_queue()->{$taskname}} >= $batch_size) {
		$self->flush($taskname);
	}
}




=head2 create_Job

  Arg [1]    : string $taskname
  Arg [2]    : string $module
  Arg [3]    : string $input_id
  Arg [4]    : string $parameter_string
  Example    :
  Description: Factory method.  Creates a new LSF job.
  Returntype : Bio::EnsEMBL::Job::LSFJob;
  Exceptions : none
  Caller     : PipelineManager

=cut

sub create_Job {
	my ($self, $taskname, $module, $input_id, $parameter_string) = @_;

	my $config = $self->get_Config();
	my $job_adaptor = $config->get_DBAdaptor->get_JobAdaptor();

	my $job = Bio::EnsEMBL::Pipeline::Job->new
    (-TASKNAME   => $taskname,
     -MODULE     => $module,
     -INPUT_ID   => $input_id,
     -PARAMETERS => $parameter_string);


	#store the job and set its status to created
	$job_adaptor->store($job);

	return $job;
}


sub _task_queue {
	my $self = shift;

	$self->{'_task_queue'} ||= {};

	return $self->{'_task_queue'};
}



#
# creates temp directories like:
# /tmp/repeat_masker_task/1/
# /tmp/genscan_task/4/
# etc..
#
sub _dir_prefix {
	my $self = shift;
	my $taskname  = shift;

  #distribute temp files evenly into 10 different dirs so that we don't
  #get too many files in the same dir
  $self->{'dir_num'} = 0 if(!defined($self->{'dir_num'}));
	$self->{'dir_num'} = ($self->{'dir_num'}+1) % 10;
	
	#
	# get temp dir from config
	#
	my $config = $self->get_Config();
	my $temp_dir = $config->get_parameter('LSF', 'output_dir');

	if(!$temp_dir) {
    $self->warn("could not determine output dir for task [$taskname]" .
               " : using /dev/null");
    return undef;
  }
	
	my $task_dir = "$temp_dir/$taskname";

  if(! -e $task_dir) {
    if(!mkdir($task_dir)) {
      $self->warn("could not create output dir [$task_dir] : using /dev/null");
      return undef;
    }
  }
	
	$temp_dir = "$task_dir/".$self->{'dir_num'};

	#create the dir if it doesn't exist
  if(! -e $temp_dir) {
    if(!mkdir($temp_dir)) {
      $self->warn("could not create output dir [$temp_dir] : using /dev/null");
      return undef;
    }
  }

  return "$temp_dir/";
}


1;
