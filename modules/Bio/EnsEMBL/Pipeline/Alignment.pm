# Cared for by Dan Andrews <dta@sanger.ac.uk>
#
# Copyright EnsEMBL
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME
  
Bio::EnsEMBL::Pipeline::Alignment.pm
  
=head1 SYNOPSIS


=head1 DESCRIPTION

A very basic object that contains sequences of an alignment 
(as AlignmentSeq objects).  Stores (or at least will store) 
information about the alignment too.
  
=head1 CONTACT
  
Post general queries to B<ensembl-dev@ebi.ac.uk>

=cut

package Bio::EnsEMBL::Pipeline::Alignment;

use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

@ISA = qw();

sub new {
  my ($class, @args) = @_;

  my $self = bless {},$class;
  
  my ($name, 
      $align_seqs) = rearrange([qw(NAME
				   SEQS)],@args);

  if ($align_seqs) {
    foreach my $align_seq (@$align_seqs){
      $self->add_sequence($align_seq);
    }
  }

  $self->alignment_name($name) if $name;

  return $self;
}

sub alignment_name {
  my $self = shift;

  if (@_) {
    $self->{'_name'} = shift;
  }

  return $self->{'_name'};
}

sub add_sequence {
  my ($self, $align_seq) = @_;

  $self->throw("Must attach AlignmentSeq objects to Alignment, not [$align_seq]")
    unless $align_seq->isa("Bio::EnsEMBL::Pipeline::Alignment::AlignmentSeq");

  push (@{$self->{'_alignment_seqs'}}, $align_seq);

  return 1;
}

sub fetch_AlignmentSeqs {
  my $self = shift;

  return $self->{'_alignment_seqs'};
}

return 1;
