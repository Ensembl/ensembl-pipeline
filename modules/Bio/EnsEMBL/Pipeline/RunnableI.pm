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

Bio::EnsEMBL::Pipeline::Runnable::RunnableI

=head1 SYNOPSIS

=head1 DESCRIPTION

Interface for running external programs

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::Runnable::RunnableI;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object;

use Bio::Root::Object;

@ISA = qw(Bio::Root::Object);


sub run {
    my ($self) = @_;

    $self->throw("run not implemented in Bio::EnsEMBL::Pipeline::Runnable::RunnableI");
}

sub output {
    my ($self) = @_;

    $self->throw("output not implemented in Bio::EnsEMBL::Pipeline::Runnable::RunnableI");

}

