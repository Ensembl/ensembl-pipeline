use strict;
use warnings;
package Bio::EnsEMBL::Pipeline::Task::RDB::Slice_CPG;

use vars qw(@ISA);

use Bio::EnsEMBL::Pipeline::Task::RDB::NoDependancies;


@ISA = ('Bio::EnsEMBL::Pipeline::Task::RDB::NoDependancies');

#this should beable to use the new method from RDB

sub name{
  my $self = shift;
  return 'slice_cpg_task';
}

sub logic_name{
  my $self = shift;
  return 'CpG';
}

sub module{
  my $self = shift;
  return 'Bio::EnsEMBL::Pipeline::RunnableDB::Slice_CPG';
}

sub description{
  print STDERR "SliceCPG runs the runnable Slice_CPG and has no dependancies\n";
}
