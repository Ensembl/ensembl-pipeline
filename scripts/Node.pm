package Node;
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
