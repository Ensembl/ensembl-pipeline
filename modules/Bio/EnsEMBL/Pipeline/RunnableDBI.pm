#
# Interface for running programs
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

Bio::EnsEMBL::Pipeline::RunnableDBI

=head1 SYNOPSIS

=head1 DESCRIPTION

Interface for running external programs.  This
inherits from
C<Bio::EnsEMBL::Pipeline::RunnableI>, and, unlike
a pure C<Bio::EnsEMBL::Pipeline::RunnableI>
module, defines methods for reading and writing
input and output to and from the database which
stores the analysis.  This separation is so that

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the
object methods. Internal methods are usually
preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::RunnableDBI;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::Root::Object;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI Bio::Root::Object);

=head1 CONCRETE METHODS

=head2 dbobj

    Title   :   dbobj
    Usage   :   $self->dbobj($obj);
    Function:   Gets or sets the value of dbobj
    Returns :   A Bio::EnsEMBL::Pipeline::DB::ObjI compliant object
                (which extends Bio::EnsEMBL::DB::ObjI)
    Args    :   A Bio::EnsEMBL::Pipeline::DB::ObjI compliant object

=cut

sub dbobj {
    my( $self, $value ) = @_;
    
    if ($value) {
        $value->isa("Bio::EnsEMBL::Pipeline::DB::ObjI")
            || $self->throw("Input [$value] isn't a Bio::EnsEMBL::Pipeline::DB::ObjI");
        $self->{'_dbobj'} = $value;
    }
    return $self->{'_dbobj'};
}

=head1 ABSTRACT METHODS

These methods need to be defined in any module
implementing C<Bio::EnsEMBL::Pipeline::RunnableI>.

=head2 fetch_input

    $self->fetch_input($id);

Fetches the input (selected via C<$id>) for a job
from the database (accessed via C<dbobj>), and
stores it in the object.

=head2 fetch_output

    $self->fetch_output($file);



=head2 write_output

=cut

sub fetch_input {
    my( $self ) = @_;
    $self->throw("fetch_input not implemented");
}

sub fetch_output {
    my( $self ) = @_;
    $self->throw("fetch_output not implemented");
}

sub write_output {
    my( $self ) = @_;
    $self->throw("write_input not implemented");
}
