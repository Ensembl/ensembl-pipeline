use strict;
use warnings;

package Bio::EnsEMBL::Pipeline::SubmissionSystem::LSF;

use constant LSF_MAX_BATCH_SIZE => 65536;


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

	my @tasknames;

	if($taskname) {
		@tasknames = ($taskname);
	} else {
		@tasknames = keys %{$self->_task_queue()};
	}

	my $config = $self->get_Config();
	my $lsf_job_adaptor = $config->db()->get_LSFJobAdaptor();

	foreach $taskname (@tasknames) {

		#extract the queue name and free form parameters from the 'where'
		# the first part of the 'where' string is LSF which we already know
		my $where = $config->get_parameter($taskname, 'where');
		my ($queue, $other_parms);
		(undef, $queue, $other_parms) = split(':', $where, 3);
		
		$queue ||
			$self->throw("Could not determine LSF queue for task [$taskname]");

		my @jobs = $self->_task_queue->{$taskname};

		next if(@jobs == 0);

		#
		# If there is only a single job submit it normally
		# otherwise submit it as a job array
		#
		#array index needs to be uniquely named
		my $lsf_job_name = $job->taskname() . '_' . time();

		my $bsub = 'bsub';
		my $pre_exec = "runLSF.pl -pre_exec ";
		$bsub .= "-E $pre_exec";
		
		
		
		if(@jobs == 1) {
			my $job = $jobs[0];
			my $file_prefix = $self->_generate_file_prefix($job);

			$job->lsf_job_id($lsf_job_name);
			$job->stdout("${file_prefix}.out");
			$job->stderr("${file_prefix}.err");
			$job->array_index(undef);

			$lsf_job_adaptor->store($job);
			$job->update_status('SUBMITTED');

		} else {
			my $array_index = 1;

			foreach my $job (@jobs) {
				my $file_prefix = $self->_generate_file_prefix($job);

				$job->stdout("${file_prefix}.out");
				$job->stderr("${file_prefix}.err");
				$job->lsf_job_id($lsf_job_name);
				$job->array_index($array_index++);

				# write entry into LSF job table
				$lsf_job_adaptor->store($job);

				# set each job status to submitted
				$job->update_status($submitted);
			}	
			
			$bsub .= " -o " . $job->stdout . " -e " . $job->stderr .
				" -E $pre_exec -q $queue -J${lsf_job_name}[".scalar(@jobs).']';
		}

		#execute the bsub to submit the job or job_array
			
			#update the LSF job table to include the lsf job id
			$
		}
	}


		



	return;
}


=head2 submit

  Arg [1]    : Bio::EnsEMBL::Pipeline::Job::LSFJob $job
	Arg [2]    :
  Example    : $lsf->submit($job);
  Description: Submits an LSF job.
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub submit {
	my $self = shift;
	my $job  = shift;

	unless(ref($job) && $job->isa('Bio::EnsEMBL::Job::LSFJob')) {
		$self->throw('expected Bio::EnsEMBL::Job::LSFJob argument');
	}

	# retrieve batch size from config
	my $config = $self->get_Config();
	my $taskname = $job->taskname();
	my $batch_size = $config->get_parameter($taskname, 'batch_size');

	$batch_size || $self->throw('Could not determine batch size for task ' .
															"[$taskname] from config");


	if($batch_size > LSF_MAX_BATCH_SIZE) {
		$self->warn("Maximum batch size [" . LSF_MAX_BATCH_SIZE .
								"] exceeded for LSF task [$taskname]");
		$batch_size = LSF_MAX_BATCH_SIZE;
	}

	#place the job into a queue specific to each task and module
	$self->_task_queue()->{$taskname} ||= [];
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
	my ($self, $taskname, $module, $input_id, $parameter_string);

	my $config = $self->get_Config();
	my $job_adaptor = $config->db()->get_JobAdaptor();

	my $job = new Bio::EnsEMBL::Job::LSFJob->new($taskname, $module,
																							 $input_id, $parameter_string);


	$job_adaptor->store($job);
	$job->update_status('CREATED');
}


sub _task_queue {
	my $self = shift;

	$self->{'_task_queue'} ||= {};

	return $self->{'_task_queue'};
}


sub _generate_filename_prefix {
	my $self = shift;
	my $job  = shift;

	#
	# get temp dir from config
	#
	my $taskname = 
	my $config = $self->get_Config();
	my $temp_dir = $config->get_parameter($job->taskname(), 'temp_dir');

	$temp_dir || $self->throw('Could not determine temp dir for task ['.
														$task->taskname() . ']');
	
	my $time = localtime(time());
	$time =~ tr/ :/_./;

	return "$temp_dir/" . $job->taskname . "_job" . $job->dbID() . "$time";
}


1;
