use strict;
use warnings;
package Bio::EnsEMBL::Pipeline::Task::RDB::BestPmatch;

use vars qw(@ISA);

use Bio::EnsEMBL::Pipeline::Task::RDB::NoDependancies;

@ISA = ('Bio::EnsEMBL::Pipeline::Task::RDB::NoDependancies');



sub can_start{
  my $self = shift;
  return $self->get_TaskStatus('pmatch_task')->is_finished;
}

sub name{
  my $self = shift;
  return 'bestpmatch_task';
}

sub module{
  my $self = shift;
  return 'Bio::EnsEMBL::Pipeline::RunnableDB::BestPmatch';
}

sub description{
  my $self = @_;
  return $self->name." ".$self->logic_name." runs the runnableDB BestPmatch and is dependant on all Pmatch analysis being finished\n";
}


sub get_input_ids{
  my ($self) = @_;
  my @input_ids = ('genome');
  my $idset = Bio::EnsEMBL::Pipeline::IDSet->new(
						 -id_list => \@input_ids,
						);
  return $idset;
}
