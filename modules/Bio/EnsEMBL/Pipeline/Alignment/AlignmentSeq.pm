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

  # Mandatory argument
  $self->seq($seq);
  $self->deletions($deletions); # Even if there isn't a sequence
                                # an undefined value will cause
                                # a full sequence of '-' to be made.

  # Optional arguments
  $self->name($name)           if $name;
  $self->start($start)         if $start;
  $self->exon($exon)           if $exon;
  $self->type($type)           if $type;

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

  my $gap = '-' x $gap_length;
  my @gap = split //, $gap;
  
  my $sequence_as_array = $self->seq_array;
  my $deletions_as_array = $self->deletion_array;

  splice (@$sequence_as_array, $insert_position, 0, @gap);
  splice (@$deletions_as_array, $insert_position, 0, @gap);

  $self->store_seq_array($sequence_as_array);
  $self->store_deletion_array($deletions_as_array);

  return 1;
}

sub deletions {
  my $self = shift;

  if ((@_)||(!$self->{'_deletions'})) {
    my $value = shift;

    if (!$value) {
      $value = '-';
    }

    $self->throw("Deletion list needs to be a string, you have probably passed in an array.")
      if (ref $value);

    # If the deletion string is empty or shorter
    # than the sequence string, pad the end with
    # '-'.

    my @deletion_array = split //, $value;

    my $length_difference = (scalar @{$self->seq_array}) - (scalar @deletion_array);

    while ($length_difference) {
      $value .= '-';
      $length_difference--;
    }
    
    $self->{'_deletions'} = $value;
  }

  return $self->{'_deletions'};
}

sub deletion_array {
  my $self = shift;

  my @array = split //, $self->deletions;

  return \@array;
}

sub store_deletion_array {
  my ($self, $input_array) = @_;

  $self->throw("Array for storage contains nothing.")
    unless (scalar @$input_array > 0);

  my $concat_seq = '';

  foreach my $element (@$input_array) {
    $concat_seq .= $element;
  } 

  $self->deletions($concat_seq);

  return 1;

}

sub fetch_deletion_at_position {
  my ($self, $base_position) = @_;

  $self->throw("Must specify a base position to retrieve a specific position from the sequence.")
    unless $base_position;

  $base_position--;

  return substr $self->deletions, $base_position, 1;
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


