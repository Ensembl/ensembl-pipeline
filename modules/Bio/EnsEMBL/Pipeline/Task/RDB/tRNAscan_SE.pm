use strict;
use warnings;
package Bio::EnsEMBL::Pipeline::Task::RDB::tRNAscan_SE;

use vars qw(@ISA);

use Bio::EnsEMBL::Pipeline::Task::RDB::NoDependancies;


@ISA = ('Bio::EnsEMBL::Pipeline::Task::RDB::NoDependancies');

#this should beable to use the new method from RDB

sub name{
  my $self = shift;
  return 'tRNAscan_task';
}

sub logic_name{
  my $self = shift;
  return 'tRNAscan';
}

sub module{
  my $self = shift;
  return 'Bio::EnsEMBL::Pipeline::RunnableDB::tRNAscan_SE';
}

sub description{
  print STDERR "tRNAscan runs the runnable tRNAscan_SE and has no dependancies\n";
}
