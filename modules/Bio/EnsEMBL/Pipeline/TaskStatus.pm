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

Stores IDSets for a Task for all of the different statuses jobs can have from
that task. These statuses include Created, Submitted, Reading, Running,
Writing, Successful, Failed, Fatal, Exists, Killed

=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::TaskStatus;

use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::IDSet;
use warnings;

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
  my ($created, $submitted, $reading, $running, $writing, $successful,
      $failed, $fatal, $killed, $existing, $retried) =
    $self->_rearrange([qw(CREATED SUBMITTED READING RUNNING WRITING 
			  SUCCESSFUL FAILED FATAL KILLED EXISTING RETRIED)],
		      @args);

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
  $self->{'_retried'} = undef;
  $self->{'is_finished'} = undef;

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
  $self->add_retried($retried) if($retried);

  return $self;
}

=head2 _check

  Arg [1]   : Bio::EnsEMBL::Pipeline::IDSet or arrayref
  Function  : to check what the argument is anc if its a listref create an 
              IDSet
  Returntype: Bio::EnsEMBL::Pipeline::IDSet
  Exceptions: if the argument isn't of the required type'
  Caller    :
  Example   : $self->_check

=cut

sub _check{
  my ($self, $arg) = @_;

  if($arg){
    if(ref($arg) eq 'ARRAY'){
	
      my $idset = Bio::EnsEMBL::Pipeline::IDSet->new(
						     -ID_LIST => $arg,
						    );
      return $idset;
    }

    if($arg->isa("Bio::EnsEMBL::Pipeline::IDSet")){
      return $arg;
    }
  }
  $self->throw("Must pass either Bio::EnsEMBL::Pipeline::IDSet or array" .
	       " refs to TaskStatus add methods not $arg");

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
  $self->{'_created'} = ($self->{'_created'}) ?
    $self->{'_created'}->or($idset) : $idset;
}

sub add_submitted{
  my ($self, $arg) = @_;

  my $idset = $self->_check($arg);
  $self->{'_submitted'} = ($self->{'_submitted'}) ?
    $self->{'_submitted'}->or($idset) : $idset;
}

sub add_reading{
  my ($self, $arg) = @_;

  my $idset = $self->_check($arg);
  $self->{'_reading'} = ($self->{'_reading'}) ?
    $self->{'_reading'}->or($idset) : $idset;
}

sub add_running{
  my ($self, $arg) = @_;

  my $idset = $self->_check($arg);
  $self->{'_running'} = ($self->{'_running'}) ?
    $self->{'_running'}->or($idset) : $idset;
}



sub add_successful{
  my ($self, $arg) = @_;

  my $idset = $self->_check($arg);
  $self->{'_successful'} = ( $self->{'_successful'}) ?
    $self->{'_successful'}->or($idset) : $idset;
}

sub add_writing{
  my ($self, $arg) = @_;

  my $idset = $self->_check($arg);
  $self->{'_writing'} = ($self->{'_writing'}) ?
    $self->{'_writing'}->or($idset) : $idset;
}

sub add_failed{
  my ($self, $arg) = @_;

  my $idset = $self->_check($arg);
  $self->{'_failed'} = ($self->{'_failed'}) ?
    $self->{'_failed'}->or($idset) : $idset;
}


sub add_fatal{
  my ($self, $arg) = @_;

  my $idset = $self->_check($arg);
  $self->{'_fatal'} = ($self->{'_fatal'}) ?
    $self->{'_fatal'}->or($idset) : $idset;
}

sub add_killed{
  my ($self, $arg) = @_;

  my $idset = $self->_check($arg);
  $self->{'_killed'} = ($self->{'_killed'}) ?
    $self->{'_killed'}->or($idset) : $idset;
}

sub add_existing{
  my ($self, $arg) = @_;

  my $idset = $self->_check($arg);

  $self->{'_existing'} = ($self->{'_existing'}) ?
    $self->{'_existing'}->or($idset) : $idset;
}

sub add_retried{
  my ($self, $arg) = @_;
  my $idset = $self->_check($arg);

  $self->{'_retried'} = ($self->{'_retried'}) ?
    $self->{'_retried'}->or($idset) : $idset;
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
    if(!$self->{'_created'}){
      $self->{'_created'} = Bio::EnsEMBL::Pipeline::IDSet->new;
    }
    return $self->{'_created'};
}

sub get_submitted{
    my ($self) = @_;
    if(!$self->{'_submitted'}){
      $self->{'_submitted'} = Bio::EnsEMBL::Pipeline::IDSet->new;
    }
    return $self->{'_submitted'};
}

sub get_reading{
    my ($self) = @_;
    if(!$self->{'_reading'}){
      $self->{'_reading'} = Bio::EnsEMBL::Pipeline::IDSet->new;
    }
    return $self->{'_reading'};
}

sub get_running{
    my ($self) = @_;
    if(!$self->{'_running'}){
      $self->{'_running'} = Bio::EnsEMBL::Pipeline::IDSet->new;
    }
    return $self->{'_running'};
}

sub get_writing{
    my ($self) = @_;
    if(!$self->{'_writing'}){
      $self->{'_writing'} = Bio::EnsEMBL::Pipeline::IDSet->new;
    }
    return $self->{'_writing'};
}

sub get_successful{
    my ($self) = @_;
    if(!$self->{'_successful'}){
      $self->{'_successful'} = Bio::EnsEMBL::Pipeline::IDSet->new;
    }
    return $self->{'_successful'} ;
}

sub get_failed{
    my ($self) = @_;
    if(!$self->{'_failed'}){
      $self->{'_failed'} = Bio::EnsEMBL::Pipeline::IDSet->new;
    }
    return $self->{'_failed'};
}

sub get_fatal{
    my ($self) = @_;
    if(!$self->{'_fatal'}){
      $self->{'_fatal'} = Bio::EnsEMBL::Pipeline::IDSet->new;
    }
    return $self->{'_fatal'};
}

sub get_killed{
    my ($self) = @_;
    if(!$self->{'_killed'}){
      $self->{'_killed'} = Bio::EnsEMBL::Pipeline::IDSet->new;
    }
    return $self->{'_killed'};
}

sub get_existing{
    my ($self) = @_;
    if(!$self->{'_existing'}){
      $self->{'_existing'} = Bio::EnsEMBL::Pipeline::IDSet->new;
    }

    return $self->{'_existing'};
}

sub get_retried{
    my ($self) = @_;
    if(!$self->{'_retried'}){
      $self->{'_retried'} = Bio::EnsEMBL::Pipeline::IDSet->new;
    }
    return $self->{'_retried'};
}


=head2 status_report

  Arg [1]   : none
  Function  : returns a formatted string listing the number of each status
              type the object holds.  Useful for debugging or displaying
              status information in a monitor tool.
  Returntype: none
  Exceptions: none
  Caller    : general
  Example   : print STDERR $self->status_report

=cut


sub status_report{
    my ($self) = @_;

    return join("\n",
                "  Created   : ".$self->get_created->count,
                "  Submitted : ".$self->get_submitted->count,
                "  Reading   : ".$self->get_reading->count,
                "  Running   : ".$self->get_running->count,
                "  Writing   : ".$self->get_writing->count,
                "  Successful: ".$self->get_successful->count,
                "  Failed    : ".$self->get_failed->count,
                "  Fatal     : ".$self->get_fatal->count,
                "  Killed    : ".$self->get_killed->count,
                "  Retried   : ".$self->get_retried->count,
                "  Existing  : ".$self->get_existing->count);
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
  $self->{'_retried'} = undef;
  $self->{'_writing'} = undef;
  $self->{'_killed'} = undef;
  $self->{'_successful'} = undef;
  $self->{'_failed'} = undef;
  $self->{'_fatal'} = undef;
  $self->{'_existing'} = undef;
}



sub is_finished{
  my ($self) = shift;
  $self->{'is_finished'} = shift if(@_);
  return $self->{'is_finished'};
}

1;
