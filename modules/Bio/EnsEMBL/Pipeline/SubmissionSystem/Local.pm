package Bio::EnsEMBL::Pipeline::SubmissionSystem::Local;
use vars qw(@ISA);
use strict;
use warnings;

use Bio::EnsEMBL::Pipeline::SubmissionSystem;
use Bio::EnsEMBL::Pipeline::Job;

@ISA = qw(Bio::EnsEMBL::Pipeline::SubmissionSystem);

my $children = 0;
my @_queue;

=head2 new

  Arg [1]    : Bio::EnsEMBL::Pipeline::PipelineManager
  Example    : 
  Description: Constructor.  Creates a new local submission system.  This class
               is a singleton, and only one instance is ever in existance.
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;

  my $self = bless {}, $class;

  $SIG{CHLD} = \&sig_chld;

  return $self->SUPER::new(@_);
}


=head2 submit

  Arg [1]    : Bio::EnsEMBL::Pipeline::Job $job
  Example    : 
  Description: This is used to submit the job.  For the Local submission
               system this simply means forking and running the job locally.
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub submit{

  my $self = shift;
  my $job  = shift;
  print STDERR "have ".$self." and job ".$job."\n";
  # Check params
  unless(ref($job) && $job->isa('Bio::EnsEMBL::Pipeline::Job')) {
    $self->throw('expected Bio::EnsEMBL::Pipeline::Job argument');
  }


  unshift @_queue, $job;

  while (@_queue > 0) {
    if ($children < 1) {
      my $job = shift @_queue;

      if (my $pid = fork) { # fork returns PID of child to parent, 0 to child
        # PARENT
        $children++;
        print "Size of children is now " . $children . "; " . 
          scalar(@_queue) . " jobs left in queue\n";
	
      } else {
        #CHILD

	my $config = $self->get_Config();
	my $temp_dir = $config->get_parameter('LOCAL', 'tmpdir');
        my $file_prefix = $job->_generate_filename_prefix($job, $temp_dir);

	# redirect stdout/stderr
        $job->stdout_file("${file_prefix}.out");
        $job->stderr_file("${file_prefix}.err");
	close(STDERR);
	close(STDOUT);
	open(STDERR, ">" . $job->stderr_file()) || warn "Error redirecting STDERR to " .  $job->stderr_file();
        open(STDOUT, ">" . $job->stdout_file()) || warn "Error redirecting STDOUT to " .  $job->stdout_file();

        #print "Executing $job with PID $$\n";
        $job->submission_id($$);
        $job->adaptor->update($job);
        $job->set_current_status('SUBMITTED');

        $job->run();
        exit(0); # child process is finished now!
      }

    } else {

      sleep(2);

    }

  }

  # wait for any children to complete
  while ($children > 0) {
    print "$children children still running, waiting ...\n";
    sleep(2);
  }


}


sub sig_chld {
  $children--;
  print "Child died; children now $children\n";
}

=head2 create_job

  Arg [1]    : string $taskname
  Arg [2]    : string $module
  Arg [3]    : string $input_id
  Arg [4]    : string $parameter_string
  Example    : 
  Description: Factory method.  Creates a job.
  Returntype :
  Exceptions : 
  Caller     : 

=cut

sub create_Job {

  my ($self, $taskname, $module, $input_id, $parameter_string) = @_;

  my $config = $self->get_Config();
  my $job_adaptor = $config->get_DBAdaptor()->get_JobAdaptor();
  
  my $job = Bio::EnsEMBL::Pipeline::Job->new(
					     -TASKNAME => $taskname, 
					     -MODULE => $module, 
					     -INPUT_ID => $input_id, 
					     -PARAMETERS => $parameter_string);

  $job_adaptor->store($job);

  return $job;

}


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

  my $job_id = $job->dbID;

  my $rc = system('kill', '-3', $job_id);     # is -3 enough?

  if($rc & 0xffff) {
    $self->warn("kill of job $job_id returned non-zero exit status $!");
    return;
  }

  $job->update_status('KILLED');
}


=head2 flush

  Arg [1]    : 
  Example    : 
  Description: Present so this is polymorphic with all submission systems
               flush() does nothing for the local submission system.
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub flush {

  my $self = shift;

  return;

}

sub _generate_filename_prefix {

  my $self = shift;
  my $job = shift;

  # get temp dir from config
  my $config = $self->get_Config();
  my $temp_dir = $config->get_parameter('LOCAL', 'output_dir');
  $temp_dir || $self->throw('Could not determine output dir for job ' . $job->taskname() . ' ID ' . $job->dbID() . '\n');

  # have a subdirectory for each type of task
  $temp_dir .= "/" .$self->taskname;

  #distribute temp files evenly into 10 different dirs so that we don't
  #get too many files in the same dir
  $self->{'dir_num'} = 0 if(!defined($self->{'dir_num'}));
  $self->{'dir_num'} = $self->{'dir_num'} +1 % 10;
	
  $temp_dir .= "/" . $self->{'dir_num'};

  #create the dir if it doesn't exist
  mkdir($temp_dir) if(! -e $temp_dir);

  my $time = localtime(time());
  $time =~ tr/ :/_./;

  return "$temp_dir/" . "_job_" . $self->dbID() . "$time";

}

1;
