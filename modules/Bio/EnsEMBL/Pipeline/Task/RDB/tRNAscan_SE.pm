use strict;
use warnings;
package Bio::EnsEMBL::Pipeline::Task::RDB::tRNAscan_SE;

use vars qw(@ISA);

use Bio::EnsEMBL::Pipeline::Task::RDB::NoDependancies;


@ISA = ('Bio::EnsEMBL::Pipeline::Task::RDB::NoDependancies');

#this should beable to use the new method from RDB

sub name{
  my $self = shift;
  return 'trnascan_task';
}

sub module{
  my $self = shift;
  return 'Bio::EnsEMBL::Pipeline::RunnableDB::tRNAscan_SE';
}

sub description{
  return "tRNAscan runs the runnable tRNAscan_SE and has no dependancies\n";
}
