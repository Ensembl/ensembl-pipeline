use strict;
use warnings;
package Bio::EnsEMBL::Pipeline::Task::RDB::Genscan;

use vars qw(@ISA);

use Bio::EnsEMBL::Pipeline::Task::RDB::RepeatMaskerDependant;

@ISA = ('Bio::EnsEMBL::Pipeline::Task::RDB::RepeatMaskerDependant');



#this should beable to use the new method from RDB

sub name{
  my $self = shift;
  return 'genscan_task';
}

sub module{
  my $self = shift;
  return 'Bio::EnsEMBL::Pipeline::RunnableDB::Genscan';
}

sub description{
  my $self = @_;
  return $self->name." ".$self->logic_name." runs the runnabledb Genscan and is dependant on RepeatMasker being finished\n";
}
