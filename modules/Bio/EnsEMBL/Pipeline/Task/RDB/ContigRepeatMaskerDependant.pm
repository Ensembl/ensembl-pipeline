use strict;
use warnings;
package Bio::EnsEMBL::Pipeline::Task::RDB::ContigRepeatMaskerDependant;

use vars qw(@ISA);

use Bio::EnsEMBL::Pipeline::Task::RDB;
use Bio::EnsEMBL::Pipeline::IDSet;


@ISA = ('Bio::EnsEMBL::Pipeline::Task::RDB');


#the new is used from the base class as this constructor wouldn't need
#to do any additional work


=head2 can_start

  Arg [1]   : none
  Function  : figures out what contig based input ids have sucessfully 
  finished the repeatmasker_task and as such can start thet next job
  Returntype: none
  Exceptions: none
  Caller    : 
  Example   : 

=cut



sub can_start{
  my $self = shift;
 
  my $input_ids = $self->get_input_ids;
  my $repeatmask_success = $self->get_TaskStatus('repeatmasker_task')->get_successful;
  my $can_start = $input_ids->and($repeatmask_success);
  $self->input_ids_to_start($can_start);
}



=head2 input_id_to_start

  Arg [1]   : Bio::EnsEMBL::Pipeline::IDSet
  Function  : compares the passed in IDSet to the existing IDSet to produce
  a list which contains the union of both lists
  Returntype: Bio::EnsEMBL::Pipeline::IDSet
  Exceptions: none
  Caller    : 
  Example   : 

=cut



sub input_ids_to_start{
  my $self = shift;

  if(@_){
    my $can_start = shift;
    if(!$self->{'can_start'}){
      $self->{'can_start'} = Bio::EnsEMBL::Pipeline::IDSet->new();
    }
    my $unique = $self->{'can_start'}->or($can_start);
    $self->{'can_start'} = $unique;
  }
  return $self->{'can_start'};
}




1;
