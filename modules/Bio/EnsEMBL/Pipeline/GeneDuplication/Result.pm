=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2017] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Pipeline::GeneDuplication::Result -

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut

package Bio::EnsEMBL::Pipeline::GeneDuplication::Result;

use warnings ;
use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::Root;


@ISA = qw(Bio::EnsEMBL::Root);


sub new {
  my ($class, @args) = @_;

  my $self = bless {}, $class;

  my ($id,
      $distance_method) = $self->_rearrange([qw(ID
						DISTANCE_METHOD)],@args);

  $id              && $self->query_id($id);
  $distance_method && $self->distance_method($distance_method);

  return $self;
}

sub query_id {
  my $self = shift;

  if (@_) {
    $self->{_query_id} = shift;
  }

  return $self->{_query_id}
}

sub distance_method {
  my $self = shift;

  if (@_) {
    $self->{_distance_method} = shift;
  }

  return $self->{_distance_method}
}

sub matrix {
  my $self = shift;

  if (@_) {
    $self->{_matrix} = shift;
  }

  return $self->{_matrix}
}

sub otus {
  my $self = shift;

  if (@_) {
    $self->{_otus} = shift;
  }

  return $self->{_otus}
}


sub add_match {
  my ($self,
      $query_id,
      $match_id,
      $dn,
      $ds,
      $n,
      $s,
      $lnL) = @_;

  die "Dont have a match identifier to store with match."
    unless $match_id;

  die "Dont have nonsynonymous and synonymous values to store with match."
    unless defined $dn & defined $ds;

  my %match_hash;

  $match_hash{query_id} = $query_id;
  $match_hash{match_id} = $match_id;
  $match_hash{dN}       = $dn;
  $match_hash{dS}       = $ds;
  $match_hash{N}        = $n;
  $match_hash{S}        = $s;
  $match_hash{lnL}      = $lnL ? $lnL : 0; #Set to zero if no defined.

  push @{$self->{_matches}}, \%match_hash;

  return 1;
}

sub matches {
  my $self = shift;

  return $self->{_matches}
}



return 1;
