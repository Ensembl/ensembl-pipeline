#
# Mickeymouse implementation of RunnableI
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

Bio::EnsEMBL::Pipeline::Runnable::RevComp

=head1 SYNOPSIS

=head1 DESCRIPTION

Mickeymouse implementation of RunnableI that reverse complements a sequence

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::Runnable::RevComp;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::EnsEMBL::Root;

use Bio::EnsEMBL::Pipeline::RunnableI;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

sub new {
    my ($class,@args) = @_;

    my $self = $class->SUPER::new(@args);

    my ($seq) = $self->_rearrange([qw(SEQ)],@args);

    $seq->isa("Bio::Seq") || $self->throw("Input isn't a Bio::Seq");
    
    $self->seq($seq);
    
    return $self;
}

sub run {
    my ($self) = @_;

 
    my $str = $self->seq->seq;

    $str =~ tr/atgcATGC/tacgTACG/;
    $str = reverse($str);

    $self->{'_revseq'} = Bio::Seq->new(-seq => $str,
				     -id  => $self->seq->id . ".rc",
				     );
    
}

sub output {
    my ($self) = @_;

    return $self->{'_revseq'};
}


sub seq {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{'_seq'} = $arg;
    }
    return $self->{'_seq'};
}
