# holds a Coordinate Pair; name of the protein and contig that are hit (as a cross check) and 
# the start and end positions of the 

package Bio::EnsEMBL::Pipeline::Tools::Pmatch::CoordPair;
use Bio::Root::Object;

@ISA = qw(Bio::Root::Object);

=head2 new

 Title   : new
 Usage   :
 Function: constructor
 Example :
 Returns : 
 Args    : 


=cut

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

  $self->throw("No query start") unless defined $qstart;
  $self->qstart($qstart);

  $self->throw("No query end") unless defined $qend;
  $self->qend($qend);

  $self->throw("No target start") unless defined $tstart;
  $self->tstart($tstart);

  $self->throw("No target end") unless defined $tend;
  $self->tend($tend);

  $self->throw("No percent") unless defined $percent;
  $self->percent($percent);

  $self->throw("No strand") unless defined $strand;
  $self->strand($strand);

  return $self;

}

=head2 query

 Title   : query
 Usage   :
 Function: get/set for query (contig name)
 Example :
 Returns : 
 Args    : 


=cut

sub query {
  my ($self,$arg) = @_;
  if ($arg) {
    $self->{'query'} = $arg;
  }
  return $self->{'query'};
}

=head2 target

 Title   : target
 Usage   :
 Function: get/set for target (protein name)
 Example :
 Returns : 
 Args    : 


=cut

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
