#
# Object for storing sequence analysis details
#
# Cared for by Michele Clamp  <michele@sanger.ac.uk>
#
# Copyright Michele Clamp
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Status - small object storing job status tags

=head1 SYNOPSIS

    my $obj    = new Bio::EnsEMBL::Pipeline::Status
    ('-jobid'              => $jobid,
     '-status'             => $status,
     '-created'            => $created,
     );

=head1 DESCRIPTION

Stores the status of a job at a certain time

=head1 CONTACT

ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::Status;

use vars qw(@ISA);
use strict;

use Bio::Root::RootI;

@ISA = qw(Bio::Root::RootI);


sub new {
  my($class,@args) = @_;
  
  my $self = $class->SUPER::new(@args);

  my ($jobid,$status,$created)  =
      $self->_rearrange([qw(JOBID
			    STATUS
			    CREATED
			    )],@args);

  $jobid   || $self->throw("Can't create a status object with no jobid");
  $status  || $self->throw("Can't create a status object with no status string");
  $created || $self->throw("Can't create a status object with no created time");

  $self->jobid             ($jobid);
  $self->status            ($status);
  $self->created           ($created);

  return $self;
}


=head2 jobid

  Title   : jobid
  Usage   : $self->jobid
  Function: Get/set method for the jobid
  Returns : int
  Args    : int

=cut

sub jobid {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_jobid} = $arg;
    }

    return $self->{_jobid};
}

=head2 status

  Title   : status
  Usage   : $self->status
  Function: Get/set method for the status string
  Returns : string
  Args    : string

=cut

sub status {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_status} = $arg;
    }

    return $self->{_status};
}

=head2 created

  Title   : created
  Usage   : $self->created
  Function: Get/set method for the created time
  Returns : int
  Args    : int

=cut

sub created {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_created} = $arg;
    }

    return $self->{_created};
}

1;
