# holds a pmatch protein hit - simply the name(idenitifier) of the protein
# and a list of ContigHit objects

package ProteinHit;
use Bio::Root::Object;
#require ContigHit;
#BEGIN { print STDERR "\n\n*** new ProteinHit ***\n"; };

@ISA = qw(Bio::Root::Object);

sub new {
  my ($class, @args) = @_;
  my $self = bless {}, $class;

  my ($id) = $self->_rearrange(['ID'], @args);

  $self->throw("No id") unless defined $id;
  $self->id($id);
  
  # _contig_hits is a hash keyed by contig name for easy lookup
#  $self->{'_contig_hits'} = [];

  return $self;

}

sub id {
  my ($self,$id) = @_;
  if ($id) {
    $self->{'id'} = $id;
  }
  return $self->{'id'};
}

sub add_ContigHit {
  my ($self,$hit) = @_;
  $self->throw('No contig hit') unless defined $hit;
  $self->throw('$contig is not a ContigHit') unless $hit->isa("ContigHit");
  $self->{_contig_hits}{$hit->id()} = $hit;
}

sub each_ContigHit {
  my ($self) = @_;
  return values %{$self->{_contig_hits}};
}

sub get_ContigHit {
  my ($self,$contig) = @_;
#  print STDERR "wibble " . $self->{_contig_hits} . "\n";
  return ($self->{_contig_hits}{$contig}) if defined $contig;
  return undef;
}
1;

# holds a MergedHit - produced by munging together CoordinatePairs
# MergedHit knows which contig and protein it is pairing, strand, overall 
# coverage and details of component CoordinatePairs

package MergedHit;
use Bio::Root::Object;

#BEGIN { print STDERR "\n\n*** new CoordPair ***\n"; };
#require CoordPair;
@ISA = qw(Bio::Root::Object);

sub new {
  my ($class, @args) = @_;
  my $self = bless {}, $class;

  my ($query, $target, $strand,$coverage) = $self->_rearrange(['QUERY',
							       'TARGET',
							       'STRAND',
							       'COVERAGE'],@args);

  $self->throw("No query") unless defined $query;
  $self->query($query);

  $self->throw("No target") unless defined $target;
  $self->target($target);

  $self->throw("No strand") unless defined $strand;
  $self->strand($strand);

  $self->throw("No coverage") unless defined $coverage;
  $self->coverage($coverage);

  $self->{'_coord_pairs'} = [];

  return $self;

}

sub query {
  my ($self,$arg) = @_;
  if ($arg) {
    $self->{'query'} = $arg;
  }
  return $self->{'query'};
}


sub target {
  my ($self,$arg) = @_;
  if ($arg) {
    $self->{'target'} = $arg;
  }
  return $self->{'target'};
}

sub coverage {
  my ($self,$arg) = @_;
  if ($arg) {
    $self->{'coverage'} = $arg;
  }
  return $self->{'coverage'};
}

sub strand {
  my ($self,$arg) = @_;
  if ($arg) {
    $self->{'strand'} = $arg;
  }
  return $self->{'strand'};
}

sub add_CoordPair {
  my ($self,$pair) = @_;
  $self->throw('No coord pair') unless defined $pair;
  $self->throw('$pair is not a CoordPair') unless $pair->isa("CoordPair");
  push(@{$self->{_coord_pairs}},$pair);
}

sub each_CoordPair {
  my ($self) = @_;
  return @{$self->{_coord_pairs}};
}
1;

# holds a pmatch contig hit - simply the name(idenitifier) of the contig
# and a list of start-end positions

package ContigHit;
use Bio::Root::Object;
#require CoordPair;
# inherit from Bio::::Root::Object so can use throw, rearrange etc .

#BEGIN { print STDERR "\n\n*** new ContigHit ***\n"; };

@ISA = qw(Bio::Root::Object);

sub new {
  my ($class, @args) = @_;
  my $self = bless {}, $class;

  my ($id) = $self->_rearrange(['ID'], @args);

  $self->throw("No id") unless defined $id;
  $self->id($id);
 
  $self->{'_forward_pairs'} = [];
  $self->{'_reverse_pairs'} = [];

  return $self;

}

sub id {
  my ($self,$id) = @_;
  if ($id) {
    $self->{'id'} = $id;
  }
  return $self->{'id'};
}

sub add_CoordPair {
  my ($self,$pair) = @_;
  $self->throw('No coord pair') unless defined $pair;
  $self->throw('$pair is not a CoordPair') unless $pair->isa("CoordPair");
  if($pair->strand == 1) {
    push(@{$self->{_forward_pairs}},$pair);
  }
  else {
    push(@{$self->{_reverse_pairs}},$pair);
  }
}

sub each_ForwardPair {
  my ($self) = @_;
  return @{$self->{_forward_pairs}};
}

sub each_ReversePair {
  my ($self) = @_;
  return @{$self->{_reverse_pairs}};
}
1;

# holds a Coordinate Pair; name of the protein and contig that are hit (as a cross check) and 
# the start and end positions of the 

package CoordPair;
# inherit from Bio::::Root::Object so can use throw, rearrange etc .
use Bio::Root::Object;
#BEGIN { print STDERR "\n\n*** new CoordPair ***\n"; };

@ISA = qw(Bio::Root::Object);

sub new {
  my ($class, @args) = @_;
  my $self = bless {}, $class;

  my ($query, $target, $qstart, $qend, $tstart, $tend, $percent, $strand) = $self->_rearrange(['QUERY',
										     'TARGET',
										     'QSTART',
										     'QEND',
										     'TSTART',
										     'TEND', 
										     'PERCENT',
										     'STRAND'],@args);

  $self->throw("No query") unless defined $query;
  $self->query($query);

  $self->throw("No target") unless defined $target;
  $self->target($target);

  $self->throw("No qstart") unless defined $qstart;
  $self->qstart($qstart);

  $self->throw("No qend") unless defined $qend;
  $self->qend($qend);

  $self->throw("No tstart") unless defined $tstart;
  $self->tstart($tstart);

  $self->throw("No tend") unless defined $tend;
  $self->tend($tend);

  $self->throw("No percent") unless defined $percent;
  $self->percent($percent);

  $self->throw("No strand") unless defined $strand;
  $self->strand($strand);

  return $self;

}

sub query {
  my ($self,$arg) = @_;
  if ($arg) {
    $self->{'query'} = $arg;
  }
  return $self->{'query'};
}


sub target {
  my ($self,$arg) = @_;
  if ($arg) {
    $self->{'target'} = $arg;
  }
  return $self->{'target'};
}

sub qstart {
  my ($self,$arg) = @_;
  if ($arg) {
    $self->{'qstart'} = $arg;
  }
  return $self->{'qstart'};
}

sub qend {
  my ($self,$arg) = @_;
  if ($arg) {
    $self->{'qend'} = $arg;
  }
  return $self->{'qend'};
}

sub tstart {
  my ($self,$arg) = @_;
  if ($arg) {
    $self->{'tstart'} = $arg;
  }
  return $self->{'tstart'};
}

sub tend {
  my ($self,$arg) = @_;
  if ($arg) {
    $self->{'tend'} = $arg;
  }
  return $self->{'tend'};
}

sub percent {
  my ($self,$arg) = @_;
  if ($arg) {
    $self->{'percent'} = $arg;
  }
  return $self->{'percent'};
}

sub strand {
  my ($self,$arg) = @_;
  if ($arg) {
    $self->{'strand'} = $arg;
  }
  return $self->{'strand'};
}

1;
