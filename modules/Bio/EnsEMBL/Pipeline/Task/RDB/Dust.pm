use strict;
use warnings;
package Bio::EnsEMBL::Pipeline::Task::RDB::Dust;

use vars qw(@ISA);

use Bio::EnsEMBL::Pipeline::Task::RDB::NoDependancies;


@ISA = ('Bio::EnsEMBL::Pipeline::Task::RDB::NoDependancies');

#this should beable to use the new method from RDB

sub name{
  my $self = shift;
  return 'dust_task';
}

sub logic_name{
  my $self = shift;
  return 'Dust';
}

sub module{
  my $self = shift;
  return 'Bio::EnsEMBL::Pipeline::RunnableDB::Dust';
}

sub description{
  print STDERR "Dust runs the runnabledb Dust and has no dependancies\n";
}
