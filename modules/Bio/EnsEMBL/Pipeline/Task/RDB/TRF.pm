use strict;
use warnings;
package Bio::EnsEMBL::Pipeline::Task::RDB::TRF;

use vars qw(@ISA);

use Bio::EnsEMBL::Pipeline::Task::RDB::NoDependancies;


@ISA = ('Bio::EnsEMBL::Pipeline::Task::RDB::NoDependancies');

#this should beable to use the new method from RDB

sub name{
  my $self = shift;
  return 'trf_task';
}

sub module{
  my $self = shift;
  return 'Bio::EnsEMBL::Pipeline::RunnableDB::TRF';
}

sub description{
  print STDERR "TRF runs the runnabledb TRF and has no dependancies\n";
}
