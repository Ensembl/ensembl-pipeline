use strict;
use warnings;

package Bio::EnsEMBL::Pipeline::SubmissionSystem;


=head2 new

  Arg [1]    : Bio::EnsEMBL::Pipeline::Config $config
  Example    : $ss = Bio::EnsEMBL::Pipeline::SubmissionSystem::LSF->new($pm);
  Description: Constructor.  This class is abstract and only implementing
               subclasses should actually be instantiated.
               Creates a new submission system.
  Returntype : Bio::EnsEMBL::Pipeline::SubmissionSystem
  Exceptions : none
  Caller     : general

=cut

sub new {
  my $caller = shift;
  my $config = shift;

	my $class = ref($caller) || $caller;

	$self->throw('config argument required') if(!$config);
	$self = bless( {'config' => $config} $class);
	
	return $self;
}


=head2 get_Config

  Arg [1]    : none
  Example    : $config = $sub_system->get_Config();
  Description: Getter for the config controlling this
               submission system.
  Returntype : Bio::EnsEMBL::Pipeline::Config
  Exceptions : none
  Caller     : implementing subclasses

=cut

sub get_Config {
	my $self = shift;

	return $self->{'config'};
}



=head2 flush

  Arg [1]    : (optional) $taskname
               Gives a hint as to which queues should be flushed.  If no
               taskname is supplied then it is implied that all internal
               queues should be flushed.
  Example    :
  Description: Does nothing unless overridden by subclass.  Should be used
               to define internal queue flushing behaviour (such as for batches
               of jobs) for submitted jobs
  Returntype : none
  Exceptions : none
  Caller     : Bio::EnsEMBL::Pipeline::PipelineManager

=cut

sub flush {
	my $self = shift;

	return;
}


=head2 submit

  Arg [1]    : Bio::EnsEMBL::Pipeline::Job $job
  Example    :
  Description: Abstract.  Submits a job. Implementation must be defined
               by subclass.
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub submit {
	my $self = shift;
	my $job  = shift;

	$self->throw('Abstract method should have been implemented by subclass');
}



=head2 create_Job

  Arg [1]    : string $taskname
  Arg [2]    : string $module
  Arg [3]    : string $input_id
  Arg [4]    : string $parameter_string
  Example    :
  Description: Abstract. Factory method.  Creates a job of an appropriate
               type.  Must be implemented by subclasses. An LSF implementation
               should create an LSF job.
  Returntype :
  Exceptions :
  Caller     :

=cut

sub create_Job {
	my $self = shift;

	$self->throw('Abstract method should have been implemented by subclass');
}


=head2 kill

  Arg [1]    : job
  Example    : $submission_system->kill($job);
  Description: Attempts to kill a job that has been submitted to this
               submission system. This base class implementation does nothing,
               and the method should be overridden by all submission
               system subclasses that have the ability to kill a submitted job.
  Returntype : none
  Exceptions : none
  Caller     : PipelineManager

=cut

sub kill {
	my $self = shift;
	my $job = shift;
}


1;
