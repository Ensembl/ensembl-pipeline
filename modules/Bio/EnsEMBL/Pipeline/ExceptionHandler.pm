# Cared for by Simon Potter <scp@sanger.ac.uk>
#
# Copyright GRL/EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code


=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::ExceptionHandler

=head1 SYNOPSIS

  my $exch = Bio::EnsEMBL::Pipeline::ExceptionHandler->new();
  $exch->text($str);
  $exch->parse();
  print $exch->id;

=head1 DESCRIPTION

This object is a simple parser for exceptions thrown from pipeline
Runnable's. The text of the exception is matched against a list of
known messages so that appropriate action can be taken. This could
be either writing something more informative than 'FAILED' in the job
status or used to determine determine whether to rerun the job, e.g.
if it failed as a result of a DB timeout.

=head1 AUTHOR

Simon Potter: scp@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut
#' # keep emacs happy

package Bio::EnsEMBL::Pipeline::ExceptionHandler;

use strict;
use Bio::EnsEMBL::Root;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Root);


=head2 new

    Title   : new
    Usage   : ExceptionHandler->new
    Function: Create new exception handler
    Returns : nothing
    Args    : none

=cut

# FIXME
# This method is a tad crude - much nicer to get this information
# from a table in the database than being hardwired like this. That
# would also fit nicely with that point-and-click pipeline GUI. TODO...

sub new {                   
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    $self->{'_exceptions'} = undef; 

    $self->add_exception('FAILED_Blast_Short_peptide',
     q{Need at least 3 nucleotides});

    return $self;
}


=head2 text

    Title   : text
    Usage   : ExceptionHandler->text($@)
    Function: Set text string of exception
    Returns : nothing
    Args    : string

=cut

sub text {
    my ($self, $arg) = @_;

    $self->{'_text'} = $arg if defined $arg;
    return $self->{'_text'};
}


=head2 id

    Title   : id
    Usage   : ExceptionHandler->id
    Function: Return exception id, e.g. 'FAILED_xxx'
    Returns : string
    Args    : nothing

=cut

sub id {
    my ($self, $arg) = @_;

    $self->{'_status'} = $arg if defined $arg;
    return $self->{'_status'};
}


=head2 parse

    Title   : parse
    Usage   : ExceptionHandler->parse
    Function: Parse the message string
    Returns : nothing
    Args    : nothing

=cut

sub parse {
    my ($self) = @_;

    $self->warn("Need to set text before using parse")
     unless $self->text;

    foreach my $exc (@{$self->{'_exceptions'}}) {
	my $str = $exc->{'str'};
	if ($self->text =~ /$str/) {
	    $self->id($exc->{'id'});
	    return;
	}
    }

    return;
}


=head2 add_exception

    Title   : add_exception
    Usage   : ExceptionHandler->add_exception
    Function: Add an exception
    Returns : nothing
    Args    : string, string

=cut

sub add_exception {
    my ($self, $id, $str) = @_;

    $self->throw("Need to specify exception id and string")
     unless (defined $id and defined $str);

    push @{$self->{'_exceptions'}}, {
	'id'  => $id,
	'str' => $str
    };
}

1;
