#
# Manual.pm - A checkpoint or manual submission system interface
#
# 
# You may distribute this module under the same terms as perl itself
#

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::SubmissionSystem::Manual - 
A checkpoint / manual submission system

=head1 SYNOPSIS


=head1 DESCRIPTION

This is an implementation of the common submission system interface for
checkpoint or manual intervention jobs. It does very little except to prompt
the user of the pipelinemanager that a job needs to be completed manually.


=head1 CONTACT

ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


use strict;
use warnings;

package Bio::EnsEMBL::Pipeline::SubmissionSystem::Manual;

use vars qw(@ISA);

use Bio::EnsEMBL::Pipeline::SubmissionSystem;

@ISA = qw(Bio::EnsEMBL::Pipeline::SubmissionSystem);


=head2 submit

  Arg [1]    : Bio::EnsEMBL::Pipeline::Job $job
  Example    :
  Description: When a job is submitted to the manual submission system
               it is set to a submitted state and the value of its parameters
               are presented to the pipelinemanager user as a prompt.  The
               job will never become successful without manual intervention,
               Though it may timeout and be resubmnitted, thereby reprompting
               the user.
  Returntype : none
  Exceptions : none
  Caller     : PipelineManager::run

=cut

sub submit {
  my $self = shift;
  my $job  = shift;

  my $config = $self->get_Config();
  my $job_adaptor = $config->get_DBAdaptor()->get_JobAdaptor();

  $job->set_current_status('SUBMITTED');

  #ensure all of the job details are updated in db including retry count
  $job_adaptor->update($job);

  my $line = '=' x 80;

  print STDERR "$line\nMANUAL JOB SUBMITTED - must be manually completed:\n" . 
    $job->parameters . "\n$line\n";



  return;
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
  my ($self, $taskname, $module, $input_id, $parameter_string) = @_;

  my $config = $self->get_Config();
  my $job_adaptor = $config->get_DBAdaptor()->get_JobAdaptor();
  
  my $job = Bio::EnsEMBL::Pipeline::Job->new(
					     -TASKNAME => $taskname, 
					     -MODULE => 'MANUAL JOB', 
					     -INPUT_ID => $input_id, 
					     -PARAMETERS => $parameter_string);
  
  $job_adaptor->store($job);

  return $job;
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

  $job->set_current_status('KILLED');
}


1;
