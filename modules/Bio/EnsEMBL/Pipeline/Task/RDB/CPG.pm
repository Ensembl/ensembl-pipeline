use strict;
use warnings;
package Bio::EnsEMBL::Pipeline::Task::RDB::CPG;

use vars qw(@ISA);

use Bio::EnsEMBL::Pipeline::Task::RDB::NoDependancies;


@ISA = ('Bio::EnsEMBL::Pipeline::Task::RDB::NoDependancies');

#this should beable to use the new method from RDB

sub name{
  my $self = shift;
  return 'cpg_task';
}

sub module{
  my $self = shift;
  return 'Bio::EnsEMBL::Pipeline::RunnableDB::CPG';
}

sub description{
  print STDERR "CPG runs the runnable CPG and has no dependancies\n";
}
