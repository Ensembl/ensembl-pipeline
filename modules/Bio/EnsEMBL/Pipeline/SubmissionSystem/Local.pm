package Bio::EnsEMBL::Pipeline::SubmissionSystem::Local;
use vars qw(@ISA);
use strict;
use warnings;

use Bio::EnsEMBL::Pipeline::SubmissionSystem;
use Bio::EnsEMBL::Pipeline::Job;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use POSIX;

@ISA = qw(Bio::EnsEMBL::Pipeline::SubmissionSystem);

my $children = 0;

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

  $self->{'_queue'} = [];

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

sub submit {

  my $self = shift;
  my $job  = shift;

  # Check params
  unless(ref($job) && $job->isa('Bio::EnsEMBL::Pipeline::Job')) {
    $self->throw('expected Bio::EnsEMBL::Pipeline::Job argument');
  }

  my $config = $self->get_Config();
  my $max_jobs = $config->get_parameter("LOCAL", "maxjobs"); # default???

  #put the just received job onto the end of the queue
  push @{$self->{'_queue'}}, $job;

  # if there aren't too many jobs executing take one off the start of the
  # queue and run it
  if ($children < $max_jobs) {

    $self->_start_job(shift @{$self->{'_queue'}});

  }

}


sub sig_chld {
  $children--;
  print "Child finished; children now $children\n";
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
  Description: If there are jobs on the queue and < max number of jobs 
               are currently running, start a new one.
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub flush {

  my $self = shift;
  my $config = $self->get_Config();
  my $max_jobs = $config->get_parameter("LOCAL", "maxjobs"); # default???

  if (scalar(@{$self->{'_queue'}}) > 0) {

    if ($children < $max_jobs) {
      my $job = shift @{$self->{'_queue'}};
      $self->_start_job($job);
    }

  }

  return;

}

sub _generate_filename_prefix {
  my $self = shift;
  my $job = shift;

  # get temp dir from config
  my $config = $self->get_Config();
  my $temp_dir = $config->get_parameter('LOCAL', 'output_dir');
  $temp_dir || $self->throw('Could not determine output dir for job ' . 
			    $job->taskname() . ' ID ' . $job->dbID() . '\n');

  # have a subdirectory for each type of task
  $temp_dir .= "/" .$job->taskname;

  if(! -e $temp_dir) {
    mkdir($temp_dir);
  }


  #distribute temp files evenly into 10 different dirs so that we don't
  #get too many files in the same dir
  $self->{'dir_num'} = 0 if(!defined($self->{'dir_num'}));
  $self->{'dir_num'} = $self->{'dir_num'} +1 % 10;
	
  $temp_dir .= "/" . $self->{'dir_num'};

  if(! -e $temp_dir) {
    mkdir($temp_dir);
  }


  my $time = localtime(time());
  $time =~ tr/ :/_./;

  return "$temp_dir/" . "job_" . $job->dbID() . "$time";
 
}


# Actually start running a job
sub _start_job {

  my $self = shift;
  my $job = shift;

  my $db = $job->adaptor()->db();
  my $dbname = $db->dbname();
  my $host   = $db->host();
  my $pass   = $db->password();
  my $user   = $db->username();
  my $port   = $db->port();

  if (my $pid = fork) {	 # fork returns PID of child to parent, 0 to child
    # PARENT
    $children++;
    print "Size of children is now " . $children . "; " . 
      scalar(@{$self->{'_queue'}}) . " jobs left in queue\n";
	
  } else {
    #CHILD

    #The child needs its own connection to the database. We don't want the
    #parent's database handle cleaned up when the child exits and vice-versa

    $db = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new
      (-dbname   => $dbname,
       -host     => $host,
       -user     => $user,
       -pass     => $pass,
       -port     => $port);

    #re-retrive the job so it is using the new dbhandle

    $job = $db->get_JobAdaptor->fetch_by_dbID($job->dbID);

    my $file_prefix = $self->_generate_filename_prefix($job);

    $job->stdout_file("${file_prefix}.out");
    $job->stderr_file("${file_prefix}.err");

    POSIX::setsid();  # make session leader, and effectively a daemon

    # redirect stdout/stderr to files
    # read stdin from /dev/null
    close(STDERR);
    close(STDOUT);
    close(STDIN);
    open(STDERR, "+>" . $job->stderr_file()) 
      || warn "Error redirecting STDERR to " .  $job->stderr_file();
    open(STDOUT, "+>" . $job->stdout_file()) 
      || warn "Error redirecting STDOUT to " .  $job->stdout_file();
    open(STDIN,  "+>/dev/null");

    #print "Executing $job with PID $$\n";
    $job->submission_id($$);
    $job->adaptor->update($job);
    $job->set_current_status('SUBMITTED');

    $job->run();
    exit(0);			# child process is finished now!
  }

  return;

}

1;
