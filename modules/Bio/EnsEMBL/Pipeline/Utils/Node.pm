=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016] EMBL-European Bioinformatics Institute

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

Node -

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut

package Bio::EnsEMBL::Pipeline::Utils::Node;
use Bio::EnsEMBL::Root;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Root);


sub new {
  my ($class, @args) = @_;
  my $self = bless {}, $class;

  my ($anal) = $self->_rearrange([qw(ANAL)], @args);

  $self->throw("No anal") unless defined $anal;
  $self->anal($anal);

  $self->{_parents} = [];
  $self->{_children} = [];
  return $self;
}


sub anal {
  my ($self,$anal) = @_;
  if ($anal) {
    $self->{_anal} = $anal;
  }
  return $self->{_anal};
}


sub add_parent {
  my ($self, $parent) = @_;

  $self->throw("No parent argument to add_parent call\n") if (!defined($parent));

  foreach my $p (@{$self->{_parents}}) {
    #if ($p->anal->logic_name eq $parent->anal->logic_name) {
    if ($p == $parent) {
      print "Already added parent\n";
      return;
    }
  }

  push @{$self->{_parents}}, $parent;
}

sub add_child {
  my ($self, $child) = @_;

  $self->throw("No child argument to add_child call\n") if (!defined($child));

  foreach my $c (@{$self->{_children}}) {
    #if ($c->anal->logic_name eq $child->anal->logic_name) {
    if ($c == $child) {
      print "Already added child\n";
      return;
    }
  }

  push @{$self->{_children}}, $child;
}

sub parents {
  my ($self) = shift;

  return $self->{_parents};
}

sub children {
  my ($self) = shift;

  return $self->{_children};
}

1;
