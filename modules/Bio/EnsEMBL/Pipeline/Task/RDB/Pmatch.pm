use strict;
use warnings;
package Bio::EnsEMBL::Pipeline::Task::RDB::Pmatch;

use vars qw(@ISA);

use Bio::EnsEMBL::Pipeline::Task::RDB::RepeatMaskerDependant;

@ISA = ('Bio::EnsEMBL::Pipeline::Task::RDB::RepeatMaskerDependant');



#this should beable to use the new method from RDB

sub name{
  my $self = shift;
  return 'pmatch_task';
}

sub module{
  my $self = shift;
  return 'Bio::EnsEMBL::Pipeline::RunnableDB::Pmatch';
}

sub description{
  my $self = @_;
  return $self->name." ".$self->logic_name." runs the runnabledb Pmatch and is dependant on RepeatMasker being finished\n";
}
