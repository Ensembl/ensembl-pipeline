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

Interface for running external programs

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::RunnableDBI;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::Root::Object;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI Bio::Root::Object);

sub dbobj {
    my ($self) = @_;

    $self->throw("dbobj not implemented in Bio::EnsEMBL::Pipeline::RunnableDBI");
}

sub fetch_input {
    my ($self) = @_;

    $self->throw("fetch_input not implemented in Bio::EnsEMBL::Pipeline::RunnableDBI");

}

