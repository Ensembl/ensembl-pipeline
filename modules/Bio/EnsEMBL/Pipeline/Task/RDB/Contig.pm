use strict;
use warnings;
package Bio::EnsEMBL::Pipeline::Task::RDB::Contig;

use vars qw(@ISA);

use Bio::EnsEMBL::Pipeline::Task::RDB;
use Bio::EnsEMBL::Pipeline::IDSet;


@ISA = ('Bio::EnsEMBL::Pipeline::Task::RDB');


#the new is used from the base class as this constructor wouldn't need
#to do any additional work


=head2 generate_input_ids

  Arg [1]   : none
  Function  : fetchs a list of contig names and creates an IDSet using
  them
  Returntype: Bio::EnsEMBL::Pipeline::IDSet 
  Exceptions: none
  Caller    : 
  Example   : my $names = $self->generate_input_idsl

=cut



sub generate_input_ids{
  my ($self) = @_;

  my $rawcontig_adaptor = $self->db->get_RawContigAdaptor;

  my $names = $rawcontig_adaptor->fetch_all_names;

  my $idset = Bio::EnsEMBL::Pipeline::IDSet(
					    -id_list => $names,
					   );

  return $idset;

}


=head2 get_input_ids

  Arg [1]   : none
  Function  : returns an IDSet of contig names
  Returntype: Bio::EnsEMBL::Pipeline::IDSet
  Exceptions: none
  Caller    : 
  Example   : my input_ids = $self->get_input_ids

=cut



sub get_input_ids{
  my ($self) = @_;

  if(!$self->{'input_ids'}){
    my $idset = $self->generate_input_ids;
    $self->{'input_ids'} = $idset;
  }
  
  return $self->{'input_ids'};
}
