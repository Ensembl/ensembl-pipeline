# holds a pmatch contig hit - simply the name(identifier) of the contig
# and a list of start-end positions

package Bio::EnsEMBL::Pipeline::Tools::Pmatch::ContigHit;
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

  my ($id) = $self->_rearrange(['ID'], @args);

  $self->throw("No id") unless defined $id;
  $self->id($id);
 
  $self->{'_forward_pairs'} = [];
  $self->{'_reverse_pairs'} = [];

  return $self;

}

=head2 id

 Title   : id
 Usage   :
 Function: get/set for contig id
 Example :
 Returns : 
 Args    : 


=cut

sub id {
  my ($self,$id) = @_;
  if ($id) {
    $self->{'id'} = $id;
  }
  return $self->{'id'};
}

=head2 add_CoordPair

 Title   : add_CoordPair
 Usage   :
 Function: adds a CoordPair to the list making up this hit
 Example :
 Returns : 
 Args    : 


=cut

sub add_CoordPair {
  my ($self,$pair) = @_;
  $self->throw('No coord pair') unless defined $pair;
  $self->throw('$pair is not a Bio::EnsEMBL::Pipeline::Tools::Pmatch::CoordPair') unless $pair->isa("Bio::EnsEMBL::Pipeline::Tools::Pmatch::CoordPair");
  if($pair->strand == 1) {
    push(@{$self->{_forward_pairs}},$pair);
  }
  else {
    push(@{$self->{_reverse_pairs}},$pair);
  }
}

=head2 each_ForwardPair

 Title   : each_ForwardPair
 Usage   :
 Function: returns CoordPairs represeting hits between a prtein and the forward strand of the contig
 Example :
 Returns : 
 Args    : 


=cut

sub each_ForwardPair {
  my ($self) = @_;
  return @{$self->{_forward_pairs}};
}

=head2 each_ReversePair

 Title   : each_Reverseair
 Usage   :
 Function: returns CoordPairs representing hits between a protein and the reverse strand of the contig
 Example :
 Returns : 
 Args    : 


=cut

sub each_ReversePair {
  my ($self) = @_;
  return @{$self->{_reverse_pairs}};
}
1;
