use strict;
use warnings;
package Bio::EnsEMBL::Pipeline::Task::RDB::Fgenesh;

use vars qw(@ISA);

use Bio::EnsEMBL::Pipeline::Task::RDB::ContigRepeatMaskerDependant;

@ISA = ('Bio::EnsEMBL::Pipeline::Task::RDB::ContigRepeatMaskerDependant');


#this should beable to use the new method from RDB


sub name{
  my $self = shift;
  return 'fgenesh_task';
}

sub logic_name{
  my $self = shift;
  return 'Fgenesh';
}

sub module{
  my $self = shift;
  return 'Bio::EnsEMBL::Pipeline::RunnableDB::Fgenesh';
}

sub description{
  my $self = @_;
  print STDERR $self->name." ".$self->logic_name." runs the runnabledb Fgenesh and is dependant on RepeatMasker being finished\n";
}
