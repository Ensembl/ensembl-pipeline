# holds a MergedHit - produced by munging together CoordinatePairs
# MergedHit knows which contig and protein it is pairing, strand, overall 
# coverage and details of component CoordinatePairs

package Bio::EnsEMBL::Pipeline::Tools::Pmatch::MergedHit;
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

=head2 query

 Title   : query
 Usage   :
 Function: get/set for query (fpccontig name)
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

=head2 coverage

 Title   : coverage
 Usage   :
 Function: get/set for coverage (% of target covered by this hit
 Example :
 Returns : 
 Args    : 


=cut

sub coverage {
  my ($self,$arg) = @_;
  if ($arg) {
    $self->{'coverage'} = $arg;
  }
  return $self->{'coverage'};
}

=head2 strand

 Title   : strand
 Usage   :
 Function: get/set for (query) strand
 Example :
 Returns : 
 Args    : 


=cut

sub strand {
  my ($self,$arg) = @_;
  if ($arg) {
    $self->{'strand'} = $arg;
  }
  return $self->{'strand'};
}

=head2 add_CoordPair

 Title   : add_CoordPair
 Usage   :
 Function: adds a CoordPair to the supporting data for this MergedHit
 Example :
 Returns : 
 Args    : 


=cut


sub add_CoordPair {
  my ($self,$pair) = @_;
  $self->throw('No coord pair') unless defined $pair;
  $self->throw('$pair is not a Bio::EnsEMBL::Pipeline::Tools::Pmatch::CoordPair') unless $pair->isa("Bio::EnsEMBL::Pipeline::Tools::Pmatch::CoordPair");


  push(@{$self->{_coord_pairs}},$pair);
  
  # need to do some sorting?
  if ($self->strand == 1) {
    @{$self->{_coord_pairs}} = sort {$a->qstart <=> $b->qstart} @{$self->{_coord_pairs}};
  }
  else {
    @{$self->{_coord_pairs}} = sort {$b->qstart <=> $a->qstart} @{$self->{_coord_pairs}};
  }
  
  $self->tstart(@{$self->{_coord_pairs}}[0]->tstart);
  $self->qstart(@{$self->{_coord_pairs}}[0]->qstart);

  $self->tend(@{$self->{_coord_pairs}}[scalar(@{$self->{_coord_pairs}})-1]->tend);
  $self->qend(@{$self->{_coord_pairs}}[scalar(@{$self->{_coord_pairs}})-1]->qend);
  
}

=head2 subsume_MergedHit

 Title   : subsume_MergedHit
 Usage   :
 Function: Incorporate all the pairs from another MergedHit into this one
 Example :
 Returns : 
 Args    : 


=cut
sub subsume_MergedHit {
  my ($self,$hit) = @_;
  $self->throw('No hit') unless defined $hit;
  $self->throw('$hit is not a Bio::EnsEMBL::Pipeline::Tools::Pmatch::MergedHit') unless $hit->isa("Bio::EnsEMBL::Pipeline::Tools::Pmatch::MergedHit");

  push(@{$self->{_coord_pairs}},$hit->each_CoordPair);
  
  # need to do some sorting?
  if ($self->strand == 1) {
    @{$self->{_coord_pairs}} = sort {$a->qstart <=> $b->qstart} @{$self->{_coord_pairs}};
  }
  else {
    @{$self->{_coord_pairs}} = sort {$b->qstart <=> $a->qstart} @{$self->{_coord_pairs}};
  }
  
  $self->tstart(@{$self->{_coord_pairs}}[0]->tstart);
  $self->qstart(@{$self->{_coord_pairs}}[0]->qstart);

  $self->tend(@{$self->{_coord_pairs}}[scalar(@{$self->{_coord_pairs}})-1]->tend);
  $self->qend(@{$self->{_coord_pairs}}[scalar(@{$self->{_coord_pairs}})-1]->qend);
  
}

=head2 each_CoordPair

 Title   : each_CoordPair
 Usage   :
 Function: returns the CoordPairs contibuting to this MergedHit
 Example :
 Returns : 
 Args    : 


=cut

sub each_CoordPair {
  my ($self) = @_;
  # sorted by qstart, based on strand

  return @{$self->{_coord_pairs}};
}

=head2 tstart

 Title   : tstart
 Usage   :
 Function: returns overall start of hit in protein coords
 Example :
 Returns : 
 Args    : 


=cut

sub tstart {
  my ($self,$arg) = @_;
  if ($arg) {
    $self->{'tstart'} = $arg;
  }
  return $self->{'tstart'};
}

=head2 tend

 Title   : tend
 Usage   :
 Function: returns overall end of hit in protein coords
 Example :
 Returns : 
 Args    : 


=cut

sub tend {
  my ($self,$arg) = @_;
  if ($arg) {
    $self->{'tend'} = $arg;
  }
  return $self->{'tend'};
}

=head2 qstart

 Title   : qstart
 Usage   :
 Function: returns overall start of hit in fpc contig coords
 Example :
 Returns : 
 Args    : 


=cut

sub qstart {
  my ($self,$arg) = @_;
  if ($arg) {
    $self->{'qstart'} = $arg;
  }
  return $self->{'qstart'};
}

=head2 qend

 Title   : qend
 Usage   :
 Function: returns overall end of hit in fpc contig coords
 Example :
 Returns : 
 Args    : 


=cut

sub qend {
  my ($self,$arg) = @_;
  if ($arg) {
    $self->{'qend'} = $arg;
  }
  return $self->{'qend'};
}

1;
