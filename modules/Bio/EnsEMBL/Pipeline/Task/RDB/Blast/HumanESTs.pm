use strict;
use warnings;
package Bio::EnsEMBL::Pipeline::Task::RDB::Blast::HumanESTs;

use vars qw(@ISA);

use Bio::EnsEMBL::Pipeline::Task::RDB::Blast;

@ISA = ('Bio::EnsEMBL::Pipeline::Task::RDB::Blast');



sub name {
  my $self = shift;
  return 'humanests_task';
}
