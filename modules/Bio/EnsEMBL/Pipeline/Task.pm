# Bio::EnsEMBL::Pipeline::Task
#
# Cared for by EnsEMBL (ensembl-dev@ebi.ac.uk)
#
# Copyright (c) EnsEMBL
#
# You may distribute this module under the same terms as perl itself
#

=pod

=head1 NAME

Bio::EnsEMBL::Pipeline::Task

=head1 SYNOPSIS

# Create some pending tasks

$pending{'RepeatMasker'} =
	new Bio::EnsEMBL::Pipeline::Task::RepeatMasker('RepeatMasker',
																								 $pipeline_manager);
$pending{'Genscan'} =
	new Bio::EnsEMBL::Pipeline::Task::Genscan('Genscan',
																						$pipeline_manager);

...

while(1) {
	# move ready tasks to the running queue
	foreach $ptname (keys %pending) {
    my $pt = $pending{$name};
    if($pt->can_start()) {
      %running{$ptname} = $pt;
			delete $pending{$ptname};
		}
	}

	# move finished tasks to the finished queue
	# give a time share for creating jobs to tasks that are still running
	foreach my $rtname (keys %running) {
		my $rt = $running{$name};
		if($rt->is_finished()) {
			$finished{$rtname} = $rt;
			delete $running{$rtname};
		} else {
			$rt->run();
		}
	}
}

=head1 DESCRIPTION

This is the base class for all pipeline Tasks.  Tasks are responsible for
implementing control flow logic that defines when a particular peice of
work needs to be run by the pipeline system.  Work is broken into smaller
chunks (or 1 or more larger chunks as applicable) known as Jobs.

The three responsibilities of all tasks are:

  (1) Deciding when to start creating jobs for a particular peice of work
  (2) Creating the jobs
  (3) Deciding when the work is completed


=head1 CONTACT

Post general questions to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods have an underscore prefix.

=cut

use strict;
package Bio::EnsEMBL::Pipeline::Task;

use vars qw(@ISA @EXPORT_OK);

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::TaskStatus;

@ISA = ('Bio::EnsEMBL::Root');


=head2 new

  Arg [1]    : string $name - The name of the task being created.
	             This is read in from the config file and should match the
               name returned by the name method.
  Arg [2]    : Bio::EnsEMBL::Pipeline::PipelineManager $plm
               The pipelinemanager which instantiated this task
  Example    : $task = Bio::EnsEMBL::Pipeline::Task->new($name, $pm);
  Description: Creates a new Task object.  This is an abstract superclass
               and should not be instantiated directly.
  Returntype : Bio::EnsEMBL::Pipeline::Task
  Exceptions : thrown if the name passed to the constructor does not
               match the name defined by the method 'name'
  Caller     : general

=cut

sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;
  
  my $self = bless({}, $class);

  my ($taskname, $plm)=$self->_rearrange([qw(TASKNAME PIPELINE_MANAGER)], @_);
	
  #
  # Verify that the configuration is using the correct name
  #
  if($taskname ne $self->name) {
    $self->throw("Task name in config [$taskname] != [" . $self->name() .
		 "] as specified by the module " . ref($self));
  }

  #
  # verify the pipeline manager argument
  #
  unless(ref($plm) && $plm->isa('Bio::EnsEMBL::Pipeline::PipelineManager')) {
    $self->throw('PipelineManager argument is required');
  }
  
  $self->{'pipeline_manager'} = $plm;

  #create a new status object for this task
  $self->{'TaskStatus'} = Bio::EnsEMBL::Pipeline::TaskStatus->new();
  
  return $self;
}



=head2 name

  Arg [1]    : none
  Example    : $name = $task->name();
  Description: Abstract. Must be implemented by subclass.  This should return
               the name of this task and must match the name of the task in
               the configuration file.
  Returntype : string
  Exceptions : thrown if not implemented by subclass
  Caller     : pipelineManager, constructor

=cut

sub name {
  my $self = shift;

  $self->throw('abstract method name not implemented by subclass');
}


=head2 run

  Arg [1]    : none
  Example    : $retval = $task->run();
  Description: 
  Returntype : int
  Exceptions : 
  Caller     : 

=cut

sub run {
  my $self = shift;
  
  $self->throw('abstract method run not implmented by subclass');
}




sub can_start {
  my $self = shift;

  $self->throw('abstract method can_start not implmented by subclass');
}


sub is_finished {
  my ($self) = @_;
  $self->throw('abstract method is_finished not implmented by subclass');
}


sub get_TaskStatus {
  my $self = shift;
  my $taskname = shift;

  if($taskname && $taskname ne $self->name()) {
    return $self->get_PipelineManager($taskname);
  }
  
  return $self->{'TaskStatus'};
}

sub get_PipelineManager{
  my ($self) = @_;

  return $self->{'pipeline_manager'};
}



sub create_Job {
  my $self = shift;
  my $module = shift;
  my $id = shift;
  my $parameters = shift;

  return $self->get_PipelineManager->create_Job($self->name(), $module, 
						$id, $parameters); 
}


sub create_Jobs {
  my $self = shift;
  my $module = shift;
  my $id_set = shift;
  my $parameters = shift;

  return $self->get_PipelineManager->create_Jobs($self->name(), $module, 
						 $id_set, $parameters);
}

1;
