use strict;
use warnings;
package Bio::EnsEMBL::Pipeline::Task::RDB::Blast;

use vars qw(@ISA);

use Bio::EnsEMBL::Pipeline::Task::RDB::RepeatMaskerDependant;

@ISA = ('Bio::EnsEMBL::Pipeline::Task::RDB::RepeatMaskerDependant');



#this should beable to use the new method from RDB


sub module{
  my $self = shift;
  return 'Bio::EnsEMBL::Pipeline::RunnableDB::Blast';
}

sub description{
  my $self = @_;
  print STDERR $self->name." ".$self->logic_name." runs the runnabledb Blast and is dependant on RepeatMasker being finished\n";
}
