=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 AUTHORS

Dan Andrews <dta@sanger.ac.uk>
 
=head1 NAME

Bio::EnsEMBL::Pipeline::Alignment - 

=head1 SYNOPSIS


=head1 DESCRIPTION

A very basic object that contains sequences of an alignment 
(as AlignmentSeq objects).  Stores (or at least will store) 
information about the alignment too.

=cut

package Bio::EnsEMBL::Pipeline::Alignment;

use warnings ;
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
