
#
# BioPerl module for Bio::EnsEMBL::Pipeline::DBSQL::Clone
#
# Cared for by Elia Stupka <elia@ebi.ac.uk>
#
# Copyright Elia Stupka
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::DBSQL::Clone - Object representing a Pipeline clone

=head1 SYNOPSIS

    # $db is Bio::EnsEMBL::DB::Obj 

    $clone = $db->get_Clone();

=head1 DESCRIPTION

Represents information on one Pipeline Clone and its update status

=head1 CONTACT

e-mail: elia@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Pipeline::DBSQL::Clone;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::RootI

use Bio::Root::RootI;
#use Bio::EnsEMBL::DB::CloneI;
#use Bio::EnsEMBL::DBSQL::Clone;

@ISA = qw(Bio::Root::RootI Bio::EnsEMBL::DB::CloneI);

# new() is inherited from Bio::Root::RootI

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;

  my ($dbobj,$disk_id) = $self->_rearrange([qw(DBOBJ
					  DISK_ID
					  )],@args);

  $disk_id || $self->throw("Cannot make clone db object without a disk_id");
  $dbobj   || $self->throw("Cannot make clone db object without db object");
  $dbobj->isa('Bio::EnsEMBL::Pipeline::DBSQL::Obj') || $self->throw("Cannot make contig db object with a $dbobj object");

  $self->disk_id($disk_id);
  $self->_dbobj($dbobj);

# set stuff in self from @args
  return $make; # success - we hope!
}

=head2 disk_id

 Title   : disk_id
 Usage   : $clone->disk_id()
 Function: Get/set methd for disk_id, which is the primary key for a clone 
 Example : $clone->disk_id()
 Returns : disk_id string
 Args    : none


=cut

sub disk_id{
    my ($self,$value) = @_;
    
    if( defined $value) {
      $self->{'disk_id'} = $value;
    }
    return $self->{'disk_id'};

}

=head2 clone_group

 Title   : clone_group
 Usage   : $clone->clone_group()
 Function: Gives the value of clone_group
           (Please note: this method will probably soon disappear)
 Example : $clone->clone_group()
 Returns : clone_group string
 Args    : none


=cut

sub clone_group{
   my ($self) = @_;
   return $self->_select('clone_group');
}

=head2 chromosome

 Title   : chromosome
 Usage   : $obj->chromosome()
 Function: Gives the chromosome number for this clone or 'unk', i.e.unknown
 Returns : chromosome string
 Args    : none


=cut

sub chromosome{
   my ($self) = @_;
   return $self->_select('chromosome');
}

=head2 last_check

 Title   : last_check
 Usage   : $obj->last_check()
 Function: Gives the last_check unix date
 Returns : last_check unix date
 Args    : none

=cut

sub last_check{
   my ($self) = @_;
   return $self->_select_date('last_check');
}

=head2 dna_update_state

 Title   : dna_update_state
 Usage   : $obj->dna_update_state()
 Function: Gives the current dna update state number
 Returns : dna_update_state string
 Args    : none


=cut

sub dna_update_state{
   my ($self) = @_;
   return $self->_select('dna_update_state');
}

=head2 update_state

 Title   : update_state
 Usage   : $obj->update_state()
 Function: Gives the currrent general update state number
 Returns : update_state string
 Args    : none


=cut

sub update_state{
   my ($self) = @_;
   return $self->_select('update_state');
}

=head2 internal_lock

 Title   : internal_lock
 Usage   : $obj->internal_lock()
 Function: Checks whether the clone is locked internally
 Note    : If a clone is locked internally, it must be locked externally!
 Returns : 1 if locked, 0 if unlocked
 Args    : none


=cut

sub internal_lock{
   my ($self) = @_;
   my $lock = $self->_select('internal_lock');
   if ($lock && !$self->external_lock) {
       $self->throw("A clone locked internally must have an external lock!");
   }
   return $lock;
}

=head2 external_lock

 Title   : external_lock
 Usage   : $obj->external_lock()
 Function: Checks wether the clone is locked to the outside world
 Returns : 1 if locked, 0 if unlocked
 Args    : none


=cut

sub external_lock{
   my ($self) = @_;
   return $self->_select('external_lock');
}

=head2 update_label

 Title   : update_label
 Usage   : $obj->update_label()
 Function: Gives the current update label
 Returns : update_label string
 Args    : none


=cut

sub update_label{
   my ($self) = @_;
   return $self->_select('update_label');
}

=head2 update_date

 Title   : update_date
 Usage   : $obj->update_date()
 Function: Gives the update_date unix date
 Returns : update_date unix date
 Args    : none


=cut

sub update_date{
   my ($self) = @_;
   return $self->_select_date('update_date');
}

=head2 created

 Title   : created
 Usage   : $obj->created()
 Function: Gives the created unix date
 Returns : created unix date
 Args    : none


=cut

sub created{
   my ($self) = @_;
   return $self->_select_date('created');
}

=head2 modified

 Title   : modified
 Usage   : $obj->modified()
 Function: Gives the modified unix date
 Returns : modified unix date
 Args    : none


=cut

sub modified{
   my ($self) = @_;
   return $self->_select_date('modified');
}


=head2 _select

 Title   : _select
 Usage   : $obj->_select
 Function: internal method for mysql select statement for a specific row name
 Example : $obj->_select('clone_group')
 Returns : result of mysql select statement
 Args    : row name


=cut

sub _select{
    my ($self,$row) = @_;

    my $disk_id = $self->disk_id();
    
    my $sth = $self->_dbobj->prepare("select $row from clone where disk_id = '$disk_id'");
    $sth->execute();
    my $rowhash = $sth->fetchrow_hashref();
    return $rowhash->{$row};
}

=head2 _select_date

 Title   : _select_date
 Usage   : $obj->_select_date
 Function: internal method for mysql select statement for a specific datetime row name
 Example : $obj->_select_date('clone_group')
 Returns : result of mysql select statement in unix time
 Args    : row name


=cut

sub _select_date{
    my ($self,$row) = @_;

    my $disk_id = $self->disk_id();
    
    my $sth = $self->_dbobj->prepare("select $row from clone where disk_id = '$disk_id'");
    $sth->execute();
    my $rowhash = $sth->fetchrow_hashref();
    my $datetime = $rowhash->{$row};
    $sth = $self->_dbobj->prepare("select UNIX_TIMESTAMP('".$datetime."')");
    $sth->execute();
    $rowhash = $sth->fetchrow_arrayref();
    return $rowhash->[0];
}

=head2 _dbobj

 Title   : _dbobj
 Usage   : $obj->_dbobj($newval)
 Function: 
 Example : 
 Returns : value of _dbobj
 Args    : newvalue (optional)


=cut

sub _dbobj{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'_dbobj'} = $value;
    }
    return $obj->{'_dbobj'};

}


sub add_Contig {
    my ($self,$contig) = @_;

    if (!(defined($self->{_contigs}))) {
	$self->{_contigs} = [];
    }
    push(@{$self->{_contigs}},$contig);

}

sub get_all_Contigs {
    my ($self) = @_;

    return @{$self->{_contigs}};
}
1;


