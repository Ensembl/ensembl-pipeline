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

Bio::EnsEMBL::Pipeline::RunnableI

=head1 SYNOPSIS

=head1 DESCRIPTION

Interface for running external programs

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the
object methods. Internal methods are usually
preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::RunnableI;

use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Analysis::Programs;

# Object preamble - inherits from Bio::Root::Object;

use Bio::Root::Object;

@ISA = qw(Bio::Root::Object);

=head1 ABSTRACT METHODS

These methods need to be implemented in any
module which implements
C<Bio::EnsEMBL::Pipeline::RunnableI>.

=head2 run

    $self->run();

Actually runs the analysis programs.  If the
analysis has fails, it should throw an
exception.  It should also remove any temporary
files created (before throwing the exception!).

=head2 output

    @output = $self->output();

Return a list of objects created by the analysis
run (eg: C<Bio::EnsEMBL::FeaturePair> objects).

=cut

sub run {
    my ($self) = @_;

    $self->throw("run not implemented");
}

sub output {
    my ($self) = @_;

    $self->throw("output not implemented");
}

#########################
# Added by MAQ 
# functions used by Runnable modules replacing hp.pl functions
#########################


sub locate_runnable {
    my ($self, $runnable) = @_;
    if ($runnable)
    {
        Bio::EnsEMBL::Programs->import($runnable);
        return Bio::EnsEMBL::Programs::Program_Paths{$runnable};
    }
}
