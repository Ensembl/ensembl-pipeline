# Cared for by Dan Andrews <dta@sanger.ac.uk>
#
# Copyright EnsEMBL
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME
  
Bio::EnsEMBL::Pipeline::Alignment::AlignmentSeq.pm
  
=head1 SYNOPSIS


=head1 DESCRIPTION


  
=head1 CONTACT
  
Post general queries to B<ensembl-dev@ebi.ac.uk>

=cut

package Bio::EnsEMBL::Pipeline::Alignment::AlignmentSeq;

use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Root;
use Bio::Seq;

# Object preamble - inherits from Bio::Root::Object;

@ISA = qw(Bio::EnsEMBL::Root);

sub new {
  my ($class, @args) = @_;

  my $self = bless {},$class;
  
  my ($name, 
      $seq, 
      $deletions,
      $start, 
      $exon,
      $type) = $self->_rearrange([qw(NAME
				     SEQ
				     DELETIONS
				     START
				     EXON
				     TYPE)],@args);


  # Optional arguments
  $self->seq($seq)     if $seq;
  $self->name($name)   if $name;
  $self->start($start) if $start;
  $self->exon($exon)   if $exon;
  $self->type($type)   if $type;

  return $self;
}

sub name {
  my $self = shift;

  if (@_) {
    $self->{'_name'} = shift;
  }

  return $self->{'_name'};
}


sub seq {
  my $self = shift;

  if (@_) {
    $self->{'_seq'} = shift;
  }

  $self->throw('No sequence attached to object')
    unless $self->{'_seq'};

  return $self->{'_seq'};
}

# Returns the sequence as an array, with one base per array element.

sub seq_array {
  my $self = shift;

  $self->throw("Cant retrieve seq_array when there is no sequence.") 
    unless $self->seq;

  my @array = split //, $self->seq; 

  return \@array;
}

# Slightly formatted output

sub seq_with_newlines {
  my ($self, $line_length) = @_;

  my $seq = $self->seq_array;

  my $concat = '';

  my $column_pos = 0;

  foreach my $element (@$seq){

   if (($column_pos%($line_length) == 0)&&($column_pos != 0)){
     $concat .= "\n";
     $column_pos = 0;
   }

   $concat .= $element;
   $column_pos++;
  }

  return $concat;
}


# Writes the sequence back from an array, where one base is
# resident in each array element.  Empty elements will be
# written as gaps '-'

sub store_seq_array {
  my $self = shift;
  my $input_array = shift;

  $self->throw("Array for storage contains nothing.")
    unless (scalar @$input_array > 0);

  my $concat_seq = '';

  foreach my $element (@$input_array) {
    $concat_seq .= $element;
  } 

  $self->seq($concat_seq);

  return 1;
}


sub fetch_base_at_position {
  my ($self, $base_position) = @_;

  $self->throw("Must specify a coordinate in order to retrieve a base from the sequence.")
    unless $base_position;

  $base_position--;

  return substr $self->seq, $base_position, 1;
}

# Specify the insertion position and length of
# gap for insertion into both the sequence and
# the deletion sequence.

sub insert_gap {
  my ($self, $insert_position, $gap_length) = @_;  

  unless ($insert_position && $gap_length){
    $self->throw("Need to specify gap insertion position [$insert_position] and length [$gap_length]");
  }

  $insert_position--;

  my $gap = '*' x $gap_length;
  my @gap = split //, $gap;

  # Stick the gap in the sequence

  my $sequence_as_array = $self->seq_array;
  splice (@$sequence_as_array, $insert_position, 0, @gap);
  $self->store_seq_array($sequence_as_array);

  # Take note of the gap in our set of deletion coordinates

  for (my $i = 0; $i < $gap_length; $i++) {
    $self->increment_deletions_above($insert_position+$i)
  }

print ">SEQUENCE\n" . $self->seq_with_newlines . "\n";

  return 1;
}

sub deletions {
  my $self = shift;

  my @sorted_array_of_deletion_coords = 
    sort {$a <=> $b} keys %{$self->{_deletion_hash}};

  return \@sorted_array_of_deletion_coords
}

sub add_deletion {
  my $self = shift;

  if (@_) {
    my $position = shift;

    $self->throw("Add deletion requires a numeric argument, not this [" .
		 ref($position) . "]") 
      if (ref($position) ne '');

    $self->_deletion_hash->{$position}++
  }

  return 1
}

sub is_a_deletion {
  my $self = shift;

  my $position = shift;

  if ($self->_deletion_hash->{$position}) {
    return 1
  }

  return 0
}

sub increment_deletions_above {
  my ($self, $coord) = @_;

print "Increment above this coord : " . $coord . "\n";

  my $coords = keys %{$self->_deletion_hash};


  for (my $i = 0; $i < scalar @$coords; $i++) {
    if ($self->_deletion_hash->{$coords->[$i]}){
      delete $self->_deletion_hash->{$coords->[$i]};
      $self->_deletion_hash->{$coords->[$i]+1} = 1;
    }
  }

  return 1
}

sub _deletion_hash {
  my $self = shift;

  unless ($self->{'_deletions'}){
    $self->{'_deletions'} = {};
  }

  return $self->{'_deletions'}
}


sub deletion_array {
  my $self = shift;

  $self->throw("deletion_array method removed due to" .
	       " the antisocial manner in which it hammers memory.");
}

sub store_deletion_array {
  my ($self, $input_array) = @_;

  $self->throw("store_deletion_array method removed due to" .
	       " the antisocial manner in which it hammers memory.");

}

sub fetch_deletion_at_position {
  my ($self, $base_position) = @_;

  $self->throw("fetch_deletion_at_position method removed due to" .
		 " the antisocial manner in which it hammers memory.");

}


# I'm not sure what this is used for...

sub start {
  my $self = shift;

  if (@_) {
    $self->{'_start'} = shift;
  }

  return $self->{'_start'}

}

# This stores the exon number.  Useful to
# Pipeline/Tools/AlignmentTool.pm

sub exon {
  my $self = shift;

  if (@_) {
    $self->{'_exon'} = shift;
  }

  return $self->{'_exon'};

}

# Sequence type store - can be 'nucleotide' or 'protein'

sub type {
  my $self = shift;

  if (@_) {

    my $type = shift;

    $self->throw("Unknown sequence type - needs to be either \'nucleotide\' or \'protein\'.")
      unless ($type eq 'nucleotide' || $type eq 'protein');

    $self->{'_type'} = $type;
    
  }

  return $self->{'_type'};
}

# Useful to know 

sub length {
  my $self = shift;

  return scalar @{$self->seq_array};
}


# Text dump sequence with identifier

sub fasta_string {
  my ($self, $line_length) = @_;

  $line_length = 60 unless $line_length;

  return '>' . $self->name . "\n" . $self->seq_with_newlines($line_length) . "\n";
}


return 1;


