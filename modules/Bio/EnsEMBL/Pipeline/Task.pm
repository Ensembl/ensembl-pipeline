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
    new Bio::EnsEMBL::Pipeline::Task::RepeatMasker(
                                       -TASKNAME => 'RepeatMasker',
                                       -PIPELINE_MANAGER => $pipeline_manager);
  $pending{'Genscan'} =
    new Bio::EnsEMBL::Pipeline::Task::Genscan(
                                       -TASKNAME => 'Genscan',
				       -PIPELINE_MANAGER => $pipeline_manager);

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


This base class is a general implementation that is quite flexible and reads
control flow information from configuration. If requirements are more
specific this class may be extended and the can_start(), run() and
is_finished() can be overridden.


=head1 CONTACT

Post general questions to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods have an underscore prefix.

=cut

use strict;
package Bio::EnsEMBL::Pipeline::Task;

use vars qw(@ISA);

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::TaskStatus;

@ISA = ('Bio::EnsEMBL::Root');


=head2 new

  Arg [1]    : string $name - The name of the task being created.
               This is read in from the config file and should match the
               name returned by the name method (which it always will unless
               the name method has been overridden to ensure strict naming).
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
	
  $self->{'name'} = $taskname;

  #
  # Verify that the configuration is using the correct name.  The name method
  # may be overridden to ensure that a certain name is always used
  #
  if($taskname ne $self->name) {
    $self->throw("Task name in config [$taskname] != [" . $self->name() .
		 "] as specified by the module " . ref($self));
  }

  # verify the pipeline manager argument
  unless(ref($plm) && $plm->isa('Bio::EnsEMBL::Pipeline::PipelineManager')) {
    $self->throw('PipelineManager argument is required');
  }

  $self->{'pipeline_manager'} = $plm;

  #create a new status object for this task
  $self->{'TaskStatus'} = Bio::EnsEMBL::Pipeline::TaskStatus->new();

  $self->_done_creating(0);

  return $self;
}



=head2 name

  Arg [1]    : none
  Example    : $name = $task->name();
  Description: Returns the name of the task specified in the configuration file
               but may be overridden by subclass.  This may be useful to ensure
               that a specific name always be used for example.  The Task
               constructor validates that the name returned by this method
               corresponds to the name in the config file.  This will always
               be the case if this method is not overridden.
  Returntype : string
  Exceptions : none
  Caller     : pipelineManager, constructor

=cut

sub name {
  my $self = shift;

  return $self->{'name'};
}


=head2 run

  Arg [1]    : none
  Example    : $retval = $task->run();
  Description: Abstract method must be implemented by subclass.
               This method is responsible for creating jobs to be run by 
               the pipeline manager.  Jobs may be created in batches of 
               arbitrary size through calls to the create_Job() and 
               create_Jobs() methods.  This method will be repeatedly called 
               until the the is_finished method returns true.

               This method must return one of the following strings:
               'TASK_OK'     - normal operation, still creating more jobs
               'TASK_DONE'   - normal operation, all finished creating jobs 
               'TASK_FAILED' - something is wrong, jobs cannot be created

               'TASK_OK' should be returned until all jobs have been created
               after which time 'TASK_DONE' should be returned.  The 
               'TASK_DONE' job is used to indicate to the pipeline manager that
               no more jobs are to be created by this job and internal 
               submission system queues should be flushed for this job.  If
               'TASK_DONE' is not returned unnecessary delays may occur in the
               submission of jobs.

  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub run {
  my $self = shift;

  #save some work if we already know we are done
  return 'TASK_DONE' if($self->_done_creating());

  my $config = $self->get_PipelineManager->get_Config();
  my $dtype = $config->get_parameter($self->name(), 'dependency_type');

  my $idset;
  #use a conservative batchsize of 50 if none was specified
  my $batchsize = $config->get_parameter($self->name(), 'batchsize') || 50;

  my $failed = $self->{'_create_failed_idset'};

  my $failed_count = ($failed) ? $failed->count : 0;

  my $needed = $batchsize - $failed_count;

  if($dtype && uc($dtype) eq 'PARTIAL') {
    $idset = $self->_get_partial_id_set($config, $needed);
  } else {
    $idset = $self->_get_full_id_set($config, $needed);
  }

  #retry jobs which were not created
  if($idset) {
    $idset = $failed->or($idset) if($failed);
  } else {
    $idset = $failed;
  }

  #if undef was returned no more ids are left
  if(!$idset) {
    $self->_done_creating(1);
    return 'TASK_DONE';
  }

  #create a batch of jobs
  if(!$self->create_Jobs($self->module(), $idset, $self->parameters())) {

    #if there was partial/full failure of job submission collect the list of
    #jobs that were not created
    $idset = $idset->not($self->get_TaskStatus()->get_existing());

    $self->{'_create_failed_idset'} = $idset;
  }

  return 'TASK_OK';
}


#
# Gets ids to create new jobs with that have been completed by tasks
# that this task is dependent on
#

sub _get_partial_id_set {
  my $self = shift;
  my $config = shift;
  my $needed = shift;

  my $name = $self->name();

  #get list of modules we depend on from config
  my $dstr = $config->get_parameter($self->name(), 'dependency_list') || '';
  my @dlist = split(';', $dstr);

  if(!@dlist) {
    $self->throw("Task [$name] config must define value for key " .
                 "'dependency_list' or override run() method");
  }

  my $idf_name = $config->get_parameter($self->name(), 'id_factory');
  my $idset;

  my $all_finished = 1;

  #get the union of completed jobs of the tasks this task depends on
  foreach my $dname (@dlist) {

    #make sure that jobs that this one depends on all use same id factory
    if($idf_name ne $config->get_parameter($dname, 'id_factory')) {
      $self->throw("Task [$name] has partial dependency on task [$dname] " .
                   "but they use different id factories");
    }

    my $ts = $self->get_TaskStatus($dname);
    if(!$idset) {
      $idset = $ts->get_successful();
    } else {
      $idset = $idset->and($ts->get_successful());
    }

    $all_finished = 0 if(!$ts->is_finished());
  }


  my $existing = $self->get_TaskStatus()->get_existing();

  $idset = $idset->not($existing)->subset($needed);

  #return undefined if all possible ids have been created already
  if($all_finished && ($needed > 0) && ($idset->count() == 0)) {
    return undef;
  }

  return $idset;
}


#
# Gets ids to create new jobs with
#
sub _get_full_id_set {
  my $self = shift;
  my $config = shift;
  my $needed = shift;

  my $id_factory = $self->get_IDFactory();

  #
  #retrieve uncreated ids from the id_factory, getting the minimum
  #necessary each time and discarding already ones from already existing jobs
  #

  my $idset = Bio::EnsEMBL::Pipeline::IDSet->new(-ID_LIST => []);
  my $next;
  my $existing = $self->get_TaskStatus->get_existing();
  while(($idset->count() < $needed) &&
        ($next = $id_factory->next($needed - $idset->count))) {
    $idset = $idset->or($next->not($existing));
  }

  #no more id can be created
  return undef if(!defined($next) && $idset->count == 0);

  return $idset;
}



=head2 can_start

  Arg [1]    : none
  Example    : run_task($task) if($task->can_start());
  Description: This method returns true when this task is able to begin
               creating jobs to be executed.  This method will be repeatedly
               called by the manager until such time as it returns true,
               when the task will be moved from the pending task list to the
               running task list.

               When an active pipeline is restarted this method will be
               called at least once in order to determine the tasks state.
               A task will not progress to a running or finished state without
               first returning true to this method.

               The default implementation of this method assumes the config
               for this task specifies a dependency_type and a list of
               dependencies.  This method may be overridden to provide
               alternative behavior.

               The dependency type must be one of: FULL, NONE, or PARTIAL:

               FULL - all of the tasks in the dependency list must completely 
               finish before this task can start.

               PARTIAL - this task can start working on identifiers which
               have been completed by all of the tasks listed in the 
               dependencies. Note that this task must create jobs with
               coincident ids to the tasks on which it depends.

               NONE - this task can start immediatly since it has no
               dependencies.

               The following is an example of how dependencies can be
               encoded in a configuration file:

               #----------------------
               [REPEAT_MASKER_TASK]
               dependency_type = NONE
               id_factory = contig_id_factory

               [GENSCAN_TASK]
               dependency_type = PARTIAL
               dependency_list = repeat_masker_task
               id_factory = contig_id_factory

               [BLAST_TASK]
               dependency_type = PARTIAL
               dependency_list = genscan_task
               id_factory = contig_id_factory

               [GENEWISE_TASK]
               dependency_type = FULL
               dependency_list = genscan_task;blast_task
               id_factory = slice_id_factory
               #-----------------------

  Returntype : boolean
  Exceptions : none
  Caller     : general

=cut

sub can_start {
  my $self = shift;

  my $name = $self->name();

  my $config = $self->get_PipelineManager->get_Config();

  my $dtype = $config->get_parameter($name, 'dependency_type');

  if(!$dtype) {
    $self->throw("[$name] configuration must specify value for key " .
                 "'dependency' or override method can_start.");
  }

  $dtype = uc($dtype);

  if($dtype eq 'FULL') {
    return $self->_full_can_start($name, $config);
  } elsif($dtype eq 'PARTIAL') {
    return $self->_partial_can_start($name, $config);
  } elsif($dtype eq 'NONE') {
    print STDERR "[$name] can start (no dependency)\n";
    return 1; #can always start tasks with no dependencies
  } else {
    $self->throw("Unknown dependency type [$dtype]");
  }
}


#
# helper method for can_start (FULL dependency)
#
sub _full_can_start {
  my ($self,$name,$config) = @_;

  my @dlist = split(/;/,$config->get_parameter($name, 'dependency_list')||'');

  if(!@dlist) {
    $self->throw("[$name] must specify value for key 'dependency_list'")
  }

  # if any of the dependent tasks are not complete then we are not ready
  # to start
  foreach my $dep (@dlist) {
    #make sure we know about task
    if(!$config->get_parameter('TASKS', $dep)) {
      $self->throw("[$name] configuration specifies dependency on " .
                   "unknown task '$dep'");
    }

    print STDERR "[$name] cannot start (full dependency)\n";

    return 0 if(!$self->get_TaskStatus($dep)->is_finished);
  }

  print STDERR "[$name] can start (full dependency)\n";

  return 1;
}

#
# helper method for can_start (PARTIAL dependency)
#
sub _partial_can_start {
  my ($self, $name, $config) = @_;

  my @dlist = split(/;/,$config->get_parameter($name, 'dependency_list')||'');

  if(!@dlist) {
    $self->throw("[$name] must specify value for key 'dependency_list'")
  }

  #we can start if all of the dependent tasks share a completed job
  #with the same identifier
  my $id_set;
  foreach my $dep(@dlist) {
    if(!$config->get_parameter('TASKS', $dep)) {
      $self->throw("[$name] configuration specifies dependency on " .
                   "unknown task '$dep'");
    }


    if(!$id_set) {
      $id_set = $self->get_TaskStatus($dep)->get_successful();
    } else {
      $id_set = $id_set->and($self->get_TaskStatus($dep)->get_successful());
    }

    #if the union of completed jobs is at any point emptyset we know we
    #cannot start
    if($id_set->count() == 0) {
      print STDERR "[$name] cannot start (partial dependency)\n";
      return 0;
    }
  }

  print STDERR "[$name] can start (partial dependency)\n";

  return 1;
}


=head2 is_finished

  Arg [1]    : none
  Example    : run_task($task) if($task->can_start());
  Description: This method returns true when this task is completely finished.
               This method will be repeatedly called on running tasks by
               the pipeline manager until such time as it returns true and
               the task is moved to the finished queue.

               The default implementation of this method returns true when
               all of this tasks jobs are complete.
               This method may be overridden.
  Returntype : boolean
  Exceptions : none
  Caller     : PipelineManager::run

=cut

sub is_finished {
  my ($self) = @_;

  my $ts = $self->get_TaskStatus;

  my $existing = $ts->get_existing();
  my $successful = $ts->get_successful();
  my $done = $self->_done_creating();

  #
  # We are finished when all jobs are created and successfully completed
  #

  if($done && ($existing->count() == $successful->count())) {
    print STDERR '[',$self->name, "] finished";
    print STDERR " (done=$done existing=",$existing->count();
    print STDERR " successful=",$successful->count(),")\n";
    return 1;
  }

  print STDERR '[',$self->name, "] not finished";
  print STDERR " (done=$done existing=",$existing->count();
  print STDERR " successful=",$successful->count(),")\n";
  return 0;
}



=head2 module

  Arg [1]    : none
  Example    : $module_name = $self->module();
  Description: Returns the name of the module that jobs will run.  The
               module must implement the methods fetch_input(), run(),
               write_output().  The default implementation of this method
               retrieves the name of this module from the configuration
               but this behaviour may be overridden.
  Returntype : string
  Exceptions : none
  Caller     : run() method

=cut

sub module {
  my $self = shift;

  my $module = $self->{'_module'};

  if(!$module) {
    my $config = $self->get_PipelineManager->get_Config();
    my $name = $self->name();

    $module = $config->get_parameter($name, 'module');

#    if(!$module) {
#      $self->throw("Task [$name] configuration must define value for key "
#                   "'module' or override the module() and/or run() methods");
#    }

    $self->{'_module'} = $module;
  }

  return $module;
}


=head2 parameters

  Arg [1]    : none
  Example    : $parameters = $self->parameters();
  Description: Returns any additional parameters to be passed to jobs that
               are to be run. The default implementation of this method
               retrieves the parameters from the configuration
               but this behaviour may be overridden.
  Returntype : string
  Exceptions : none
  Caller     : run() method

=cut

sub parameters {
  my $self = shift;

  my $config = $self->get_PipelineManager->get_Config();
  my $name = $self->name();

  return $config->get_parameter($name, 'parameters') || '';
}


=head2 get_IDFactory

  Arg [1]    : none
  Example    : $id_factory = $self->get_IDFactory()
  Description: Retrieves an id factory for this task.  The default
               implementation of this method reads information to create
               this id factory fromi the configuration but this
               behaviour can be overridden by subclasses.
  Returntype : Bio::EnsEMBL::Pipeline::IDFactory
  Exceptions : none
  Caller     : run() method

=cut

sub get_IDFactory {
  my $self = shift;

  my $id_factory = $self->{'_id_factory'};

  if(!$id_factory) {
    my $config = $self->get_PipelineManager->get_Config;
    my $name = $self->name();

    my $idf_header = $config->get_parameter($name, 'id_factory');

    if(!$idf_header) {
      $self->throw("Task [$name] configuration must define value for key " .
                   "'id_factory' or override the run and/or " .
                   "get_IDFactorymethod");
    }

    my $idf_module = $config->get_parameter($idf_header, 'module');

    if(!$idf_module) {
      $self->throw("IDFactory [$idf_header] configuration must define value" .
                   " for key 'module'");
    }

    eval "require $idf_module";

    if($@) {
      $self->throw("$idf_module cannot be found for id factory [$idf_header]\n"
                   . "Exception: $@\n");
    }

    $id_factory = "$idf_module"->new(-NAME => $idf_header, -CONFIG => $config);

    $self->{'_id_factory'} = $id_factory;
  }

  return $id_factory;
}



=head2 get_TaskStatus

  Arg [1]    : (optional) string $taskname
               The name of the task to obtain the TaskStatus of.  If not
               provided this tasks status is returned.
  Example    : my $status = $self->get_TaskStatus();
  Description: Retrieves the task status of this Task or of another tasks.
  Returntype : Bio::EnsEMBL::Pipeline::TaskStatus
  Exceptions : none
  Caller     : PipelineManager, internal

=cut

sub get_TaskStatus {
  my $self = shift;
  my $taskname = shift;

  if($taskname && $taskname ne $self->name()) {
    return $self->get_PipelineManager->get_TaskStatus($taskname);
  }

  return $self->{'TaskStatus'};
}


=head2 get_PipelineManager

  Arg [1]    : none
  Example    : $pm = $self->get_PipelineManager();
  Description: Getter for the pipeline manager object which is running the 
               pipeline
  Returntype : none
  Exceptions : none
  Caller     : internal

=cut

sub get_PipelineManager{
  my ($self) = @_;

  return $self->{'pipeline_manager'};
}



=head2 create_Job

  Arg [1]    : string $module
               The module which executes the jobs for this task
  Arg [2]    : string $id
               The input_id for the job to be created
  Arg [3]    : (optional) string $parameters
               Any additional parameters to be passed to this Job
  Example    : $self->create_Job(
                        'Bio::EnsEMBL::Pipeline::RunnableDB::RepeatMasker',
			'12341', $params);

  Description: Creates a single job for submission.  If the job is successfully
               created by the submission system the Job is returned.  Otherwise
               undef is returned.
  Returntype : Bio::EnsEMBL::Pipeline::Job
  Exceptions : none
  Caller     : run method

=cut

sub create_Job {
  my $self = shift;
  my $module = shift;
  my $id = shift;
  my $parameters = shift;

  return $self->get_PipelineManager->create_Job($self->name(), $module,
						$id, $parameters);
}


=head2 create_Jobs

  Arg [1]    : string $module
  Arg [2]    : Bio::EnsEMBL::Pipeline::IDSet
               The set of input ids to create new jobs from
  Arg [3]    : (optional) string $parameters
               Any additional parameters to be passed to the job
  Example    : $self->create_Jobs(
                            'Bio::EnsEMBL::Pipeline::RunnableDB::Genscan',
                            $idset, $parameters); 
  Description: Creates a set of new jobs to be submitted
  Returntype : none
  Exceptions : none
  Caller     : run method

=cut

sub create_Jobs {
  my $self = shift;
  my $module = shift;
  my $id_set = shift;
  my $parameters = shift;

  return $self->get_PipelineManager->create_Jobs($self->name(), $module,
                                                 $id_set, $parameters);
}


=head2 description

  Arg [1]    : none
  Example    : print STDERR $task->description();
  Description: Returns a string which describes this task. The description
               is taken from the description value specified in the 
               configuration file.  This method may be overridden to 
               provide a specific description.
  Returntype : string 
  Exceptions : none
  Caller     : general

=cut

sub description{
  my $self = @_;

  my $name = $self->name();
  my $desc = $self->get_Config->get_parameter($name, 'description');

  $desc ||= "no description";

  return "$name: $desc";
}



=head2 get_Config

  Arg [1]    : none
  Example    : $conf = $self->get_Config;
  Description: Obtains the config object for the currently running pipeline
  Returntype : Bio::EnsEMBL::Pipeline::Config
  Exceptions : none
  Caller     : general

=cut

sub get_Config{
  my ($self) = @_;

  if(!$self->{'config'}){
    $self->{'config'} = $self->get_PipelineManager->get_Config;
  }

  return $self->{'config'};
}



sub _done_creating{
  my $self = shift;
  $self->{'_done_creating'} = shift if(@_);
  return $self->{'_done_creating'};
}


1;
