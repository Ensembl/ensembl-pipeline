use strict;
use warnings;
package Bio::EnsEMBL::Pipeline::Task::RDB::RepeatMaskerDependant;

use vars qw(@ISA);

use Bio::EnsEMBL::Pipeline::Task::RDB;
use Bio::EnsEMBL::Pipeline::IDSet;


@ISA = ('Bio::EnsEMBL::Pipeline::Task::RDB');


#the new is used from the base class as this constructor wouldn't need
#to do any additional work


=head2 can_start

  Arg [1]   : none
  Function  : figures out what contig based input ids have sucessfully 
  finished the repeatmasker_task and as such can start thet next job
  Returntype: none
  Exceptions: none
  Caller    : 
  Example   : 

=cut



sub can_start{
  my $self = shift;
 
  my $input_ids = $self->get_input_ids;
  my $config = $self->get_Config;
  my $type = $self->input_id_type;
  if($type eq 'contig'){
    my $repeatmask_success = $self->get_TaskStatus('repeatmasker_task')->get_successful;
    
   
    return $repeatmask_success->count;
  }else{
    #$self->warn("this task ".$self->name." uses a input_id type 
    #which isn't "."contig but ".$type." its going to wait till 
    #all RepeatMasker jobs are finished");
    return $self->get_TaskStatus('repeatmasker_task')->is_finished;
  }
}

sub update_input_ids{
  my $self = shift;
  
  my $input_ids = $self->get_input_ids;
  my $type = $self->input_id_type;
  if($type eq 'contig'){
    my $repeatmask_success = $self->get_TaskStatus('repeatmasker_task')->get_successful;
    
    my $can_start = $input_ids->and($repeatmask_success);
    $self->input_ids_to_start($can_start);
  }else{
    if($self->get_TaskStatus('repeatmasker_task')->is_finished){
      $self->input_ids_to_start($input_ids);
    }
  }
}


=head2 input_id_to_start

  Arg [1]   : Bio::EnsEMBL::Pipeline::IDSet
  Function  : compares the passed in IDSet to the existing IDSet to produce
  a list which contains the union of both lists
  Returntype: Bio::EnsEMBL::Pipeline::IDSet
  Exceptions: none
  Caller    : 
  Example   : 

=cut



sub input_ids_to_start{
  my $self = shift;

  if(@_){
    my $can_start = shift;
    if(!$self->{'can_start'}){
      $self->{'can_start'} = Bio::EnsEMBL::Pipeline::IDSet->new();
    }
    my $unique = $self->{'can_start'}->or($can_start);
    $self->{'can_start'} = $unique;
  }
  return $self->{'can_start'};
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
  $self->update_input_ids;
  return $self->start;
}

1;
