use strict;
use warnings;
package Bio::EnsEMBL::Pipeline::Task::RDB::Slice_Dust;

use vars qw(@ISA);

use Bio::EnsEMBL::Pipeline::Task::RDB::NoDependancies;


@ISA = ('Bio::EnsEMBL::Pipeline::Task::RDB::NoDependancies');

#this should beable to use the new method from RDB

sub name{
  my $self = shift;
  return 'slice_dust_task';
}

sub module{
  my $self = shift;
  return 'Bio::EnsEMBL::Pipeline::RunnableDB::Slice_Dust';
}

sub description{
  return "Slice_dust runs the runnable Slice_Dust and has no dependancies\n";
}
