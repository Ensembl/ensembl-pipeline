use strict;
use warnings;


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
  my $pm = shift;

  my $class = ref($caller) || $caller;

  # if the singleton has already been created, return it,
  # otherwise create the singleton first
  my $singleton;
  unless (defined $singleton){
    $singleton = bless {}, $class;
  }

  return $singleton;
	
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

sub create_job {
	
  my $self = shift;

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


1;
