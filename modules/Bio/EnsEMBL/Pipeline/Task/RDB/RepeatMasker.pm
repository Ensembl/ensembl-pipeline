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
  my $module = $self->module;
  my $potential = $self->get_input_ids;
  my $existing = $self->get_TaskStatus->get_existing;
  my $id_set = $potential->not($existing)->subset($self->max_create);
  $self->create_Jobs($module, 
		     $id_set, $parameters);

  return 2; #TASK_DONE
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
    return 0;
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

sub module{
  my $self = shift;
  return 'Bio::EnsEMBL::Pipeline::RunnableDB::RepeatMasker';
}

sub description{
  print STDERR "RepeatMasker runs the runnable RepeatMasker and has no dependancies\n";
}
