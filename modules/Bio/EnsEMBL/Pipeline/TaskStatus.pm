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
  my ($created, $submitted, $reading, $running, $writing,

      $successful, $failed, $fatal, $killed, $existing) = $self->_rearrange([qw(CREATED, SUBMITTED, READING, RUNNING,
                                                                               WRITING, SUCCESSFUL, FAILED, FATAL,
                                                                               KILLED, EXISTING)], @args);

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
    
    if(ref($arg) == 'ARRAY'){
       my $idset = Bio::EnsEMBL::Pipeline::IDSet->new(
                                                       -ID_LIST => $arg,
                                                     );
       return $idset;
    }elsif($arg->isa("Bio::EnsEMBL::Pipeline::IDSet")){
       return $arg;
    }else{
       $self->throw("Must pass either Bio::EnsEMBL::Pipeline::IDSets or array refs to TaskStatus add methods not $arg : $!");
    }

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
    $self->{'_created'} = $idset;
}

sub add_submitted{
    my ($self, $arg) = @_;

    my $idset = $self->_check($arg);
    $self->{'_submitted'} = $idset;
}

sub add_reading{
    my ($self, $arg) = @_;

    my $idset = $self->_check($arg);
    $self->{'_reading'} = $idset;
}


sub add_running{
    my ($self, $arg) = @_;

    my $idset = $self->_check($arg);
    $self->{'_running'} = $idset;
}

sub add_writing{
    my ($self, $arg) = @_;

    my $idset = $self->_check($arg);
    $self->{'_writing'} = $idset;
}

sub add_successful{
    my ($self, $arg) = @_;

    my $idset = $self->_check($arg);
    $self->{'_successful'} = $idset;
}

sub add_failed{
    my ($self, $arg) = @_;

    my $idset = $self->_check($arg);
    $self->{'_failed'} = $idset;
}

sub add_fatal{
    my ($self, $arg) = @_;

    my $idset = $self->_check($arg);
    $self->{'_fatal'} = $idset;
}

sub add_killed{
    my ($self, $arg) = @_;

    my $idset = $self->_check($arg);
    $self->{'_killed'} = $idset;
}

sub add_existing{
    my ($self, $arg) = @_;

    my $idset = $self->_check($arg);
    $self->{'_existing'} = $idset;
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

    return $self->{'created'};
}

sub get_submitted{
    my ($self) = @_;

    return $self->{'submitted'};
}

sub get_reading{
    my ($self) = @_;

    return $self->{'reading'};
}

sub get_running{
    my ($self) = @_;

    return $self->{'running'};
}

sub get_writing{
    my ($self) = @_;

    return $self->{'writing'};
}

sub get_successful{
    my ($self) = @_;

    return $self->{'successful'};
}

sub get_failed{
    my ($self) = @_;

    return $self->{'failed'};
}

sub get_fatal{
    my ($self) = @_;

    return $self->{'fatal'};
}

sub get_killed{
    my ($self) = @_;

    return $self->{'killed'};
}

sub get_existing{
    my ($self) = @_;

    return $self->{'existing'};
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



sub create_exisiting{
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
    
    $self->add_exisiting($total_ids);
    
    return $total_ids;
}



=head2 status_report

  Arg [1]   : none
  Function  : prints the number od each status type the object holds
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


