use strict;
use warnings;
package Bio::EnsEMBL::Pipeline::Task::RDB::NoDependancies;

use vars qw(@ISA);

use Bio::EnsEMBL::Pipeline::Task::RDB;



@ISA = ('Bio::EnsEMBL::Pipeline::Task::RDB');


#the new is used from the base class as this constructor wouldn't need
#to do any additional work



=head2 can_start

  Arg [1]   : none
  Function  : to inform if the task can start, it this case it always 
  returns 1 as it can always start
  Returntype: 1
  Exceptions: none
  Caller    : 
  Example   : 

=cut


sub can_start{
  my $self = shift;
  return 1;
}


=head2 input_ids_to_start

  Arg [1]   : none
  Function  : to return an IDSet of ids which can be started in this case
  because there are no dependancies it can always start
  Returntype: Bio::EnsEMBL::Pipeline::IDSet
  Exceptions: none
  Caller    : 
  Example   : 

=cut


sub input_ids_to_start{
  my $self = shift;
  return $self->get_input_ids;
}



1;
