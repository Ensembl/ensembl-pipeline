package Bio::EnsEMBL::Pipeline::Analysis;

use  vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Analysis;


@ISA = qw(Bio::EnsEMBL::Analysis);





sub new {
  my($class,@args) = @_;
  
  my $self = $class->SUPER::new(@args);
  
  my ($type) = $self->_rearrange([qw(TYPE)], @args);

  $self->input_id_type($type);

  return $self;

}


sub input_id_type{
  my ($self, $type) = @_;
  if($type){
    $self->{'type'} = $type;
  }
  
  return $self->{'type'};
}
