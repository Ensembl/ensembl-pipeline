use strict;
use warnings;
package Bio::EnsEMBL::Pipeline::Task::RDB;

use vars qw(@ISA);

use Bio::EnsEMBL::Pipeline::Task;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Task::Utils::InputIDFactory;

@ISA = ('Bio::EnsEMBL::Pipeline::Task');


=head2 new

  Arg [1]   : none
  Function  : Create a new RDB object and instantiate its InputIDFactory
  Returntype: Bio::EnsEMBL::Pipeline::Task::RDB
  Exceptions: none
  Caller    : 
  Example   : 

=cut



sub new {
    my ($class,@args) = @_;

    my $self = $class->SUPER::new(@args);    

    $self->{'input_id_factory'} = undef;
    $self->{'input_ids'} = undef;
    $self->{'config'} = undef;

    my $inputidfac = Bio::EnsEMBL::Pipeline::Task::Utils::InputIDFactory->new(
  -CONFIG => $self->get_Config,
  -TASKNAME => $self->name,
  -DB => $self->db,
 );

   $self->input_id_factory($inputidfac); 

   return $self;
}



=head2 input_id_factory

  Arg [1]   : Bio::EnsEMBL::Pipeline::Task::Utils::InputIDFactory
  Function  : getter/setter
  Returntype: 
  Exceptions: none
  Caller    : 
  Example   : $self->input_id_factory($input_id_fact);

=cut



sub input_id_factory{
  my $self = shift;

  if(@_){
    $self->{'input_id_factory'} = shift;
  }

  return  $self->{'input_id_factory'};
}

=head2 abstract methods

  Arg [1]   : none
  Function  : these methods should all be implemented in the child
  classes logicname should return the analysis logic_name, module should
  return the full perl declatation for the module which the Task needs to
  instantiate
  Returntype: string 
  Exceptions: This method will throw as it is an abstract method but
  sub classes should implement this method and as such shouldn't throw'
  Caller    : 
  Example   : my $logic_name = $self->logic_name

=cut


sub logic_name{
  my ($self) = @_;

  $self->throw("logic_name should be implemented by subclass $!");
}

sub module{
  my ($self) = @_;
  
  $self->throw("module should be implemented by subclass $!");
}


sub input_ids_to_start{
 my ($self) = @_;
  
 $self->throw("input_ids_to_start should be implemented by subclass $!"); 
}



=head2 get_input_ids

  Arg [1]   : none
  Function  : calls the InputIDFactory method generate_input_ids
  Returntype: Bio::EnsEMBL::Pipeline::IDSet
  Exceptions: none
  Caller    : 
  Example   : my $ids = $self->get_input_ids;

=cut



sub get_input_ids{
 my ($self) = @_;
  
 if(!$self->{'input_ids'}){
   my $idset = $self->input_id_factory->generate_input_ids;
   $self->{'input_ids'} = $idset; 
 }
 return $self->{'input_ids'}; 
}

=head2 parameter_string

  Arg [1]   : none
  Function  : returns a string which contains information about
  database connection and analysis type required by RunnableDBs
  Returntype: string
  Exceptions: throws if PipelineManager has no config object
  Caller    : 
  Example   : my $parameters = $self->parameter_string

=cut


sub parameter_string{
  my ($self) = @_;


  if(!$self->{'parameter_string'}){
    my $config = $self->get_Config;

    if(!$config){
      $self->throw("PipelineManager ".$self->get_PipelineManager.
                   " seems to be missing its config $!");
    }

    my $dbhost = $config->get_parameter('ensembl_database', 'host');
    my $dbuser = $config->get_parameter('ensembl_database', 'user');
    my $dbpass = $config->get_parameter('ensembl_database', 'pass');
    my $dbname = $config->get_parameter('ensembl_database', 'dbname');
    my $dbport = $config->get_parameter('ensembl_database', 'port');
    my $logic_name = $self->logic_name;

    my $string = "$dbhost:$dbport:$dbuser:$dbpass:$dbname:$logic_name";
    $self->{'parameter_string'} = $string;
    return $self->{'parameter_string'};
  }

  return $self->{'parameter_string'};
}


=head2 db

  Arg [1]   : none
  Function  : instantiates a core dbadaptor and returns it
  Returntype: Bio::EnsEMBL::DBSQL::DBAdaptor
  Exceptions: throws if the PipelineManager has no config object
  Caller    : 
  Example   : my $db = $self->db;

=cut

sub db{
  my ($self) = @_;

  if(!$self->{'core_db'}){
    my $config = $self->get_Config;
    
    if(!$config){
      $self->throw("PipelineManager ".$self->get_PipelineManager.
		   " seems to be missing its config $!");
    }

    my $dbhost = $config->get_parameter('ensembl_database', 'host');
    my $dbuser = $config->get_parameter('ensembl_database', 'user');
    my $dbpass = $config->get_parameter('ensembl_database', 'pass');
    my $dbname = $config->get_parameter('ensembl_database', 'dbname');
    my $dbport = $config->get_parameter('ensembl_database', 'port');
    
    my $dbadaptor = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
							-dbname => $dbname,
							-host => $dbhost,
							-user => $dbuser,
							-pass => $dbpass,
							-port => $dbport
						       );

    $self->{'core_db'} = $dbadaptor;
  }

  return $self->{'core_db'};
}


=head2 max_create

  Arg [1]   : none
  Function  : returns a number which represents the maximum number of ids
  a Task should submit at once
  Returntype: integer
  Exceptions: none
  Caller    : 
  Example   : my $id_set = $potential->not($existing)->subset($self->max_create); 

=cut


sub max_create{
  my ($self) = @_;

  if(!$self->{'max_create'}){
    my $config = $self->get_Config;
    my $max_create = $config->get_parameter($self->name, 'max_create');
    $self->{'max_create'} = $max_create;
  }
  return $self->{'max_create'};
}



=head2 run

  Arg [1]   : none
  Function  : calls to create_Jobs
  Returntype: TASK_DONE
  Exceptions: none
  Caller    : 
  Example   : $task->run;

=cut



sub run{
  my $self = shift;

  my $parameters = $self->parameter_string;
  my $module = $self->module;
  my $potential = $self->input_ids_to_start;
  my $existing = $self->get_TaskStatus->get_existing;
  my $id_set = $potential->not($existing)->subset($self->max_create);
  eval{
    $self->create_Jobs($module, 
		       $id_set, $parameters);
  };

  if($@){
    return 'TASK_FAILED';
  }

  return 'TASK_OK'; 
}



=head2 is_finished

  Arg [1]   : none
  Function  : checks if all the input_ids that exist for this task have
  been successfully completed
  Returntype: 1/0 depending on if the Task is finished or not
  Exceptions: none
  Caller    : 
  Example   : if($task->is_finished){}

=cut



sub is_finished{
  my $self = shift;

  my $potential = $self->get_input_ids;
  my $successful = $self->get_TaskStatus->get_successful;

  if(!$potential || !$successful){
    return undef;
  }elsif($potential->count == $successful->count){
    return 1;
  }else{
    return 0;
  }
}

sub get_Config{
  my ($self) = @_;

  if(!$self->{'config'}){
    $self->{'config'} = $self->get_PipelineManager->get_Config;
  }

  return $self->{'config'};
}


1;
