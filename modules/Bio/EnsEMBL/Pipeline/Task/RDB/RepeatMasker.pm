use strict;
use warnings;
package Bio::EnsEMBL::Pipeline::Task::RDB::RepeatMasker;

use vars qw(@ISA);

use Bio::EnsEMBL::Pipeline::Task::RDB::NoDependancies;


@ISA = ('Bio::EnsEMBL::Pipeline::Task::RDB::NoDependancies');



#this should beable to use the new method from RDB

sub name{
  my $self = shift;
  return 'repeatmasker_task';
}

sub module{
  my $self = shift;
  return 'Bio::EnsEMBL::Pipeline::RunnableDB::RepeatMasker';
}

sub description{
  return "RepeatMasker runs the runnabledb RepeatMasker and has no dependancies\n";
}
