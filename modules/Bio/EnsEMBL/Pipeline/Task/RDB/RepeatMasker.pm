use strict;
use warnings;
package Bio::EnsEMBL::Pipeline::Task::RDB::RepeatMasker;

use vars qw(@ISA);

use Bio::EnsEMBL::Pipeline::Task::RDB::Contig;


@ISA = ('Bio::EnsEMBL::Pipeline::Task::RDB::Contig');



#this should beable to use the new method from RDB


sub can_start{
  my $self = shift;
  return 1;
}


sub run{
  my $self = shift;

  my $parameters = $self->parameter_string;
  my $potential = $self->get_input_ids;
  my $existing = $self->get_TaskStatus->get_existing;
  my $id_set = $potential->not($existing);
  $self->create_Jobs('Bio::EnsEMBL::Pipeline::RunnableDB::RepeatMasker', $id_set, $parameters);
}

sub is_finished{
  my $self = shift;

  my $potential = $self->get_input_ids;
  my $successful = $self->get_TaskStatus->get_successful;

  if(!$potential || !$successful){
    return undef;
  }elsif($potential->count == $successful->count){
    return 1;
  }else{
    return undef;
  }
}

sub name{
  my $self = shift;
  return 'repeatmasker_task';
}

sub logic_name{
  my $self = shift;
  return 'RepeatMask';
}
