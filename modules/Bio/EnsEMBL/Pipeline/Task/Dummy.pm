use strict;
use warnings;
package Bio::EnsEMBL::Pipeline::Task::Dummy;

use vars qw(@ISA);

use Bio::EnsEMBL::Pipeline::Task;


@ISA = ('Bio::EnsEMBL::Pipeline::Task');


sub can_start{
  my $self = shift;
  return 1;
}


sub run{
  my $self = shift;
  $self->create_Job('Bio::EnsEMBL::Pipeline::DummyModule',
                    1,
                    '');

  return 'TASK_OK';
}

sub is_finished{
  my $self = shift;

  return 0;
}

sub name{
  my $self = shift;
  return 'dummy_task';
}

1;
