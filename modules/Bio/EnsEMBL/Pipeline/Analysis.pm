package Bio::EnsEMBL::Pipeline::Analysis;

use  vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( deprecate warning throw );

@ISA = qw(Bio::EnsEMBL::Analysis);





sub new {
  my($class,@args) = @_;
  
  my $self = $class->SUPER::new(@args);
  
  my ($type) = rearrange([qw(INPUT_ID_TYPE)], @args);

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
