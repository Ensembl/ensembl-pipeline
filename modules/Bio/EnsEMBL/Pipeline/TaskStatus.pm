# Object for storing an the listf of input-ids which have a particular status for a particular jon
#
# Cared for by ensembl (ensembl-dev@ebi.ac.uk)
#
# Copyright ensembl
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod

=head1 NAME

Bio::EnsEMBL::Pipeline::TaskStatus

=head1 SYNOPSIS

=head1 DESCRIPTION

stores IDSets for a Task for all of the different statuses jobs can have from that task
these status incule Created, Submitted, Reading, Running, Writing, Successful, Failed, Fatal, Exists, Killed

=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::TaskStatus;

use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::IDSet;

@ISA = qw(Bio::EnsEMBL::Root);



=head2 new

  Arg [1]   : list of args @args
  Function  : create a new Bio::EnsEMBL::Pipeline::TaskStatus
  Returntype: Bio::EnsEMBL::Pipeline::TaskStatus
  Exceptions: none
  Caller    : 
  Example   : my $taskstatus = Bio::EnsEMBL::Pipeline::TaskStatus->new();

=cut


sub new{
  my ($class, @args) = @_;
  my $self = bless {}, $class;
  my ($created, $submitted, $reading, $running, $writing,$successful, $failed, $fatal, $killed, $existing) = $self->_rearrange([qw(CREATED SUBMITTED READING RUNNING WRITING SUCCESSFUL FAILED FATAL KILLED EXISTING)], @args);

  $self->{'_created'} = undef;
  $self->{'_submitted'} = undef;
  $self->{'_reading'} = undef;
  $self->{'_running'} = undef;
  $self->{'_writing'} = undef;
  $self->{'_killed'} = undef;
  $self->{'_successful'} = undef;
  $self->{'_failed'} = undef;
  $self->{'_fatal'} = undef;
  $self->{'_existing'} = undef;

  $self->add_created($created) if($created);
  $self->add_submitted($submitted) if($submitted);
  $self->add_reading($reading) if($reading);
  $self->add_running($running) if($running);
  $self->add_writing($writing) if($writing);
  $self->add_killed($killed) if($killed);
  $self->add_succesful($successful) if($successful);
  $self->add_failed($failed) if($failed);
  $self->add_fatal($fatal) if($fatal);
  $self->add_existing($existing) if($existing);

  return $self;
}

=head2 _check

  Arg [1]   : Bio::EnsEMBL::Pipeline::IDSet or arrayref
  Function  : to check what the argument is anc if its a listref create an IDSet
 Returntype: Bio::EnsEMBL::Pipeline::IDSet
  Exceptions: if the argument isn't of the required type'
  Caller    : 
  Example   : $self->_check

=cut

sub _check{
    my ($self, $arg) = @_;
    
    if(ref($arg) eq 'ARRAY'){
      
       my $idset = Bio::EnsEMBL::Pipeline::IDSet->new(
                                                       -ID_LIST => $arg,
                                                     );
       return $idset;
    }
    
    if($arg->isa("Bio::EnsEMBL::Pipeline::IDSet")){
       return $arg;
    }
     
    $self->throw("Must pass either Bio::EnsEMBL::Pipeline::IDSet or array refs to TaskStatus add methods not $arg : $!");
    

}






=head2 add methods

  Arg [1]   : Bio::EnsEMBL::Pipeline::IDSet or list ref
  Function  : add the arg to the appropriate hash element
  Returntype: none
  Exceptions: none
  Caller    : 
  Example   : $self->add_created($idset);

=cut



sub add_created{
    my ($self, $arg) = @_;

    my $idset = $self->_check($arg);
    if(!$self->{'_created'}){
      $self->{'_created'} = $idset;
    }else{
      $self->{'_created'} = ($self->{'_created'}) ?
        $self->{'_created'}->or($idset) : $idset;
    }
}

sub add_submitted{
    my ($self, $arg) = @_;

    my $idset = $self->_check($arg);
    if(!$self->{'_submitted'}){
      $self->{'_submitted'} = $idset;
    }else{
      $self->{'_submitted'} = ($self->{'_submitted'}) ?
        $self->{'_submitted'}->or($idset) : $idset;
    }
}

sub add_reading{
    my ($self, $arg) = @_;

    my $idset = $self->_check($arg);
    if(!$self->{'_reading'}){
      $self->{'_reading'} = $idset;
    }else{
      $self->{'_reading'} = ($self->{'_reading'}) ?
        $self->{'_reading'}->or($idset) : $idset;
    }
}

sub add_running{
    my ($self, $arg) = @_;

    my $idset = $self->_check($arg);
    if(!$self->{'_running'}){
      $self->{'_running'} = $idset;
    }else{
    $self->{'_running'} = ($self->{'_running'}) ?
      $self->{'_running'}->or($idset) : $idset;
    }
}



sub add_successful{
  my ($self, $arg) = @_;

  my $idset = $self->_check($arg);
  if(!$self->{'_successful'}){
    $self->{'_successful'} = $idset;
  }else{
    $self->{'_successful'} = ($self->{'_successful'}) ?
      $self->{'_successful'}->or($idset) : $idset;
  }
}

sub add_writing{
    my ($self, $arg) = @_;

    my $idset = $self->_check($arg);
    if(!$self->{'_writing'}){
      $self->{'_writing'} = $idset;
    }else{
      $self->{'_writing'} = ($self->{'_writing'}) ?
        $self->{'_writing'}->or($idset) : $idset;
    }
}

sub add_failed{
    my ($self, $arg) = @_;

    my $idset = $self->_check($arg);
    if(!$self->{'_failed'}){
      $self->{'_failed'} = $idset;
    }else{
      $self->{'_failed'} = ($self->{'_failed'}) ?
        $self->{'_failed'}->or($idset) : $idset;
    }
}


sub add_fatal{
    my ($self, $arg) = @_;

    my $idset = $self->_check($arg);
    if(!$self->{'_fatal'}){
      $self->{'_fatal'} = $idset;
    }else{
      $self->{'_fatal'} = ($self->{'_fatal'}) ?
        $self->{'_fatal'}->or($idset) : $idset;
    }
}

sub add_killed{
    my ($self, $arg) = @_;

    my $idset = $self->_check($arg);
    if(!$self->{'_killed'}){
      $self->{'_killed'} = $idset;
    }else{
      $self->{'_killed'} = ($self->{'_killed'}) ?
        $self->{'_killed'}->or($idset) : $idset;
    }
}

sub add_existing{
    my ($self, $arg) = @_;
    my $idset = $self->_check($arg);
    if(!$self->{'_existing'}){
      $self->{'_existing'} = $idset;
    }else{
      $self->{'_existing'} = ($self->{'_existing'}) ?
        $self->{'_existing'}->or($idset) : $idset;
    }
}


=head2 get methods

  Arg [1]   : none
  Function  : return the requested IDSet
  Returntype: Bio::EnsEMBL::Pipeline::IDSet
  Exceptions: none 
  Caller    : 
  Example   : my $idset = $self->get_created

=cut



sub get_created{
    my ($self) = @_;
    return $self->{'_created'};
}

sub get_submitted{
    my ($self) = @_;

    return $self->{'_submitted'};
}

sub get_reading{
    my ($self) = @_;

    return $self->{'_reading'};
}

sub get_running{
    my ($self) = @_;

    return $self->{'_running'};
}

sub get_writing{
    my ($self) = @_;

    return $self->{'_writing'};
}

sub get_successful{
    my ($self) = @_;

    return $self->{'_successful'};
}

sub get_failed{
    my ($self) = @_;

    return $self->{'_failed'};
}

sub get_fatal{
    my ($self) = @_;

    return $self->{'_fatal'};
}

sub get_killed{
    my ($self) = @_;

    return $self->{'_killed'};
}

sub get_existing{
    my ($self) = @_;
    if(!$self->{'_existing'}){
      $self->create_existing;
    }
    return $self->{'_existing'};
}



=head2 create_existing

  Arg [1]   : none
  Function  : this takes all IDSets the object holds and performs unions
              to produce the total set of IDs
  Returntype: Bio::EnEMBL::Pipeline::IDSet
  Exceptions: none
  Caller    : 
  Example   : my $total = $self->create_exisiting

=cut



sub create_existing{
    my ($self) = @_;
    
    my $total_ids = Bio::EnsEMBL::Pipeline::IDSet->new;

    if($self->get_created){
       my $idset = $total_ids->or($self->get_created);
       $total_ids = $idset;
    }
    if($self->get_submitted){
       my $idset = $total_ids->or($self->get_submitted);
       $total_ids = $idset;
    }
    if($self->get_reading){
       my $idset = $total_ids->or($self->get_reading);
       $total_ids = $idset;
    }
    if($self->get_running){
       my $idset = $total_ids->or($self->get_running);
       $total_ids = $idset;
    }
    if($self->get_writing){
       my $idset = $total_ids->or($self->get_writing);
       $total_ids = $idset;
    }
    if($self->get_successful){
       my $idset = $total_ids->or($self->get_successful);
       $total_ids = $idset;
    }
    if($self->get_failed){
       my $idset = $total_ids->or($self->get_failed);
       $total_ids = $idset;
    }
    if($self->get_fatal){
       my $idset = $total_ids->or($self->get_fatal);
       $total_ids = $idset;
    }
    if($self->get_fatal){
       my $idset = $total_ids->or($self->get_fatal);
       $total_ids = $idset;
    }
    $self->add_existing($total_ids);

    return $total_ids;
}



=head2 status_report

  Arg [1]   : none
  Function  : prints the number of each status type the object holds
  Returntype: none
  Exceptions: none
  Caller    : 
  Example   : $self->status_report

=cut


sub status_report{
    my ($self) = @_;
    
    print STDERR "Created tasks\t".$self->get_created->count."\n" if($self->get_created);
    print STDERR "Submitted tasks\t".$self->get_submitted->count."\n" if($self->get_submitted);
    print STDERR "Reading tasks\t".$self->get_reading->count."\n" if($self->get_reading);
    print STDERR "Running tasks\t".$self->get_running->count."\n" if($self->get_running);
    print STDERR "Writing tasks\t".$self->get_writing->count."\n" if($self->get_writing);
    print STDERR "Successful tasks\t".$self->get_successful->count."\n" if($self->get_successful);
    print STDERR "Failed tasks\t".$self->get_failed->count."\n" if ($self->get_failed);
    print STDERR "Fatal tasks\t".$self->get_fatal->count."\n" if($self->get_fatal);
    print STDERR "Killed tasks\t".$self->get_killed->count."\n" if($self->get_killed);
    print STDERR "Exisiting tasks\t".$self->get_existing->count."\n" if($self->get_existing);
}



=head2 clean

  Arg [1]   : none 
  Function  : emptys all variables
  Returntype: none
  Exceptions: none
  Caller    : 
  Example   : $taskstatus->clean;

=cut



sub clean{
  my ($self) = @_;

  $self->{'_created'} = undef;
  $self->{'_submitted'} = undef;
  $self->{'_reading'} = undef;
  $self->{'_running'} = undef;
  $self->{'_writing'} = undef;
  $self->{'_killed'} = undef;
  $self->{'_successful'} = undef;
  $self->{'_failed'} = undef;
  $self->{'_fatal'} = undef;
  $self->{'_existing'} = undef;
}


=head2 clean methods

  Arg [1]   : Bio::EnsEMBL::Pipeline::IDSet or array ref (optional);
  Function  : to empty the variable and replace with an argument if its 
  passed in
  Returntype: none
  Exceptions: none
  Caller    : 
  Example   : $self->clean_created($created_ids);

=cut



sub clean_created{
  my ($self, $arg) = @_;

  if($arg){
    my $idset = $self->_check($arg);
    $self->{'_created'} = $idset;
  }else{
    $self->{'_created'} = undef;
  }
}

sub clean_submitted{
  my ($self, $arg) = @_;

  if($arg){
    my $idset = $self->_check($arg);
    $self->{'_submitted'} = $idset;
  }else{
    $self->{'_submitted'} = undef;
  }
}

sub clean_reading{
  my ($self, $arg) = @_;

  if($arg){
    my $idset = $self->_check($arg);
    $self->{'_reading'} = $idset;
  }else{
    $self->{'_reading'} = undef;
  }
}

sub clean_running{
  my ($self, $arg) = @_;

  if($arg){
    my $idset = $self->_check($arg);
    $self->{'_running'} = $idset;
  }else{
    $self->{'_running'} = undef;
  }
}

sub clean_writing{
  my ($self, $arg) = @_;

  if($arg){
    my $idset = $self->_check($arg);
    $self->{'_writing'} = $idset;
  }else{
    $self->{'_writing'} = undef;
  }
}

sub clean_successful{
  my ($self, $arg) = @_;

  if($arg){
    my $idset = $self->_check($arg);
    $self->{'_successful'} = $idset;
  }else{
    $self->{'_successful'} = undef;
  }
}

sub clean_failed{
  my ($self, $arg) = @_;

  if($arg){
    my $idset = $self->_check($arg);
    $self->{'_failed'} = $idset;
  }else{
    $self->{'_failed'} = undef;
  }
}

sub clean_fatal{
  my ($self, $arg) = @_;

  if($arg){
    my $idset = $self->_check($arg);
    $self->{'_fatal'} = $idset;
  }else{
    $self->{'_fatal'} = undef;
  }
}

sub clean_killed{
  my ($self, $arg) = @_;

  if($arg){
    my $idset = $self->_check($arg);
    $self->{'_killed'} = $idset;
  }else{
    $self->{'_killed'} = undef;
  }
}

sub clean_existing{
  my ($self, $arg) = @_;

  if($arg){
    my $idset = $self->_check($arg);
    $self->{'_existing'} = $idset;
  }else{
    $self->{'_existing'} = undef;
  }
}

1;
