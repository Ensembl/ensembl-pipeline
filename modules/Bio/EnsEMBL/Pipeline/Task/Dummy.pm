use strict;
use warnings;
package Bio::EnsEMBL::Pipeline::Task::Dummy;

use vars qw(@ISA);

use Bio::EnsEMBL::Pipeline::Task;


@ISA = ('Bio::EnsEMBL::Pipeline::Task');


sub can_start{
  my $self = shift;

  $self->{'remaining'} = 20;
  return 1;
}


sub run{
  my $self = shift;

  if($self->{'remaining'}) {

    $self->create_Job('Bio::EnsEMBL::Pipeline::DummyModule',  1, '');
    $self->{'remaining'}--;
    return 'TASK_OK';
  }

  return 'TASK_DONE';
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
