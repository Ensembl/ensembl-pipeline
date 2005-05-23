
# Cared for by Dan Andrews <dta@sanger.ac.uk>
#
# Copyright EnsEMBL
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod

=head1 NAME

Bio::EnsEMBL::Pipeline::Alignment::EvidenceAlignment

=head1 SYNOPSIS

This module allows a transcript to be displayed with its
attached supporting evidence.  It re-creates an alignment
from the database which is returned as an array of 
multiply-aligned sequences.  These can be printed as text 
and used to display the alignment with an alignment 
viewer/editor.

Quick start - use the following code if you want the 
alignment of an ensembl transcript with the evidence 
used to predict it:

my $evidence_alignment =
  Bio::EnsEMBL::Pipeline::Alignment::EvidenceAlignment->new(
      -transcript          => $transcript,
      -dbadaptor           => $db,
      -seqfetcher          => $seqfetcher,
      -padding             => 10);

my $alignment = 
  $evidence_alignment->retrieve_alignment(
     '-type'            => 'all',
   # Type can be 'all', 'nucleotide' or 'protein'
     '-remove_introns'  => 1);

my $align_seqs = $alignment->fetch_AlignmentSeqs;

foreach my $align_seq (@$align_seqs){
  print $align_seq->seq . "\n";
}

The '-padding' option specifies the amount of tailing
sequence upstream and downstream of the transcript.

Other alignment presentation options exist.  The intronic
regions of the transcript can be left in the alignment, but
in some cases can make it _very_ big:

my $alignment = 
  $align_tool->retrieve_alignment(
      '-type'            => 'all',
      '-remove_introns'  => 1);

If you are displaying protein sequences in your alignment, you
can use either single letter or three letter amino acid codes.
By default the alignment is generated with single letter codes,
which is fractionally faster than three letter codes.  To 
switch on the three letter codes specify the '-three_letter_aa'
argument when retrieving an alignment, like so:

my $alignment = 
  $evidence_alignment->retrieve_alignment(
     '-type'            => 'all',
     '-three_letter_aa' => 1,);

As it might be mind-bogglingly useful to some folks, it is possible
to display a transcript with a set of external supporting features
that can come from any source.  Where external features are used, 
the supporting features attached to the transcript are ignored and 
the external set used instead.  This feature is only going to work if 
the external features are mapped to the same assembly as the transcript.
It also usually helps if the external supporting features actually overlap
the transcript sequence :)  External features are attached to the evidence
alignment object at the time of object creation:

my $evidence_alignment = 
  Bio::EnsEMBL::Pipeline::Alignment::EvidenceAlignment->new(
     '-dbadaptor'           => $db,
     '-seqfetcher'          => $seqfetcher,
     '-transcript'          => $transcript,
     '-supporting_features' => \@supporting_features);


=head1 DESCRIPTION

Object for reconstructing alignments of predicted transcripts
and their associated supporting evidence.  Produces an alignment
that is an array of Bio::Seq objects that can be printed in
fasta format (for example) and displayed using an alignment 
viewer, such as SeaView or JalView.

=head1 CONTACT
  
Post general queries to B<ensembl-dev@ebi.ac.uk>

=cut

package Bio::EnsEMBL::Pipeline::Alignment::EvidenceAlignment;

use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Pipeline::Alignment;
use Bio::EnsEMBL::Pipeline::Alignment::AlignmentSeq;
use Bio::EnsEMBL::Utils::Exception qw(throw warning info);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

@ISA = qw();


### 'Public' methods ###


=head2 new

  Arg [db]         : a Bio::EnsEMBL::DBSQL::DBAdaptor object.
  Arg [transcript] : a Bio::EnsEMBL::Transcript object.
  Arg [seqfetcher] : a Bio::EnsEMBL::Pipeline::SeqFetcher object.
  Arg [padding]    : (optional) int, default 10, transcript slice padding.
  Arg [evidence_identity_cutoff] : (optional) percentage id cutoff to be
                                   applied to supporting features.
  Arg [supporting_features]      : (optional) a listref of 
                                   Bio::EnsEMBL::SeqFeature.  These will be
                                   displayed instead of the evidence attached
                                   to the transcript.
  Example    : my $evidence_alignment = 
                Bio::EnsEMBL::Pipeline::Alignment::EvidenceAlignment->new(
                 -transcript          => $transcript,
                 -dbadaptor           => $db,
                 -seqfetcher          => $seqfetcher,
                 -padding             => 10);
  Description: Constructs new EvidenceAlignment object.  If a listref
               of supporting features is provided, these will be displayed
               instead of the features attached to the transcript.
  Returntype : Bio::EnsEMBL::Pipeline::Alignment::EvidenceAlignment
  Exceptions : Will throw if isnt passed DBAdaptor/Transcript/SeqFetcher.
               Warns if transcript does not have evidence attached, or if all
               evidence is discarded or unfetchable.
  Caller     : General

=cut


sub new {
  my ($class, @args) = @_;

  my $self = bless {},$class;

  my ($db, 
      $transcript, 
      $seqfetcher, 
      $padding, 
      $evidence_identity_cutoff,
      $supporting_features) = rearrange([qw(DBADAPTOR
					    TRANSCRIPT
					    SEQFETCHER
					    PADDING
					    IDENTITY_CUTOFF
					    SUPPORTING_FEATURES
					   )],@args);

  # Throw an error if any of the below are undefined or
  # are the wrong thing.
  unless (defined $db && $db->isa("Bio::EnsEMBL::DBSQL::DBAdaptor")){
    throw("No database adaptor passed to AlignmentTool.  You passed a $db.");
  }
  unless (defined $transcript && $transcript->isa("Bio::EnsEMBL::Transcript")){
    throw("No transcript passed to AlignmentTool.  You passed a $transcript.");
  }
  unless (defined $seqfetcher && $seqfetcher->isa("Bio::DB::RandomAccessI")) {
    throw("No sequence fetcher passed to AlignmentTool.  You passed a $seqfetcher.");
  }

  $self->_db($db);

  # Set the amount of flanking sequence we want to include around
  # our slice.  This is quite important - many supporting features
  # extend past the beginnings and ends of predicted transcripts.  
  # Without padding, these sequences are truncated.  The default
  # of 50bp works OK, but you would want to set this manually
  # higher if you are interested in up- or down-stream goings-on.
  $self->_padding($padding) if (defined $padding);

  # Store our SeqFetcher

  $self->_seq_fetcher($seqfetcher);

  # Create the slice we will work on

  $self->_slice($transcript);

  # To take account of padding and/or possibly a reverse stranded transcript, it 
  # is necessary to transfer our transcript into the new slice coordinates.

  $transcript = $transcript->transfer($self->_slice);
  $self->_transcript($transcript);

  # Determine if our database contains translations.  If it doesn't
  # we'll have to skip adding a set of translated exons to our 
  # alignment.

  if ($self->_transcript->translation){
    $self->_translatable(1);
  } else {
    warning("Database doesn't contain translation.  Subsequently, ".
	    "wont be able to display translations of each exon.");
    $self->_type('nucleotide');
    $self->_translatable(0);
  }

  # When an external set of features are used, over-ride the
  # features from the transcript

  if (defined $supporting_features) {
    $self->_all_supporting_features($supporting_features)
  }

  # Optional evidence identity cut-off

  if ($evidence_identity_cutoff) {
    $self->{_evidence_identity_cutoff} = $evidence_identity_cutoff;
  } else {
    $self->{_evidence_identity_cutoff} = 0;
  }

  return $self;
}


=head2 retrieve_alignment

  Arg [type]            : String.  Alignment type (all, nucleotide, protein).
  Arg [remove_introns]  : 0 or 1.  Truncate intron sequences from alignment.
  Arg [three_letter_aa] : 0 or 1.  Display amino acid sequences with three-letter 
                          instead of single-letter codes.
  Example    : my $alignment = 
                $evidence_alignment->retrieve_alignment(
                 '-type'            => 'all',
                 '-remove_introns'  => 1,
                 '-three_letter_aa' => 1,);
  Description: Retrieves an alignment object of the transcript passed during object
               construction.  The alignment is not constructed prior to this call,
               hence this method drives the full alignment reconstruction process.
  Returntype : Bio::EnsEMBL::Pipeline::Alignment
  Exceptions : Warns if alignment generation fails.
  Caller     : General.

=cut

sub retrieve_alignment {
  my $self = shift @_;

  my ($type, 
      $remove_introns, 
      $three_letter_aa) = rearrange([qw(TYPE
					REMOVE_INTRONS
       					THREE_LETTER_AA
				       )],@_);

  $self->_type($type);
  $self->_three_letter_aa($three_letter_aa) if $three_letter_aa;

  unless ($self->_is_computed($type)){
    my $alignment_success = $self->_align($type, 
					  $remove_introns);
    unless ($alignment_success) {
      warning "Alignment generation failed.  There probably were" . 
	" not any displayable sequences.";
      return 0
    }
  }

  return $self->_create_Alignment_object($type);
}


### Main Internal Methods ###


=head2 _align

  Arg [1]    : String.  Alignment type (one of 'all', 'nucleotide' 
               or 'protein').
  Arg [2]    : int (1 or 0).  Set to 1 if introns are to be truncated.
  Example    :  my $status = $self->_align('all', 1);
  Description: Top level internal method that builds alignment from the
               database.  While other sudsidary methods perform tasks 
               such as building sequences from cigar strings, this 
               method does the difficult reconciliation of multiple
               pairwise alignments into a single multiple alignment.
  Returntype : int
  Exceptions : none.
  Caller     : retrieve_alignment

=cut

sub _align {
  my ($self, $type, $truncate_introns) = @_;

  my $evidence_sequences = $self->_corroborating_sequences;

  return 0 unless $evidence_sequences;

  # Insert deletions in the appropriate places in the genomic
  # sequence.  This will complete the alignment of the evidence
  # with the genomic sequence.  A bit of fiddling is needed to 
  # propagate gaps into the evidence sequences that need them.

    # Make some handy data structures

  my %deletion_sets;
  my %all_deletions;
  my %deletion_tracking;
  my %evidence_sequence_hash;

  foreach my $evidence_sequence (@$evidence_sequences) {
    # Note:  Purposely using the memory addresses as hash keys for evidence
    # sequence objects as they don't necessarily have unique names.

    $evidence_sequence_hash{$evidence_sequence} = $evidence_sequence;

    while (my ($deletion_coord, $length) = each %{$evidence_sequence->deletion_hash}){
      $deletion_sets{$evidence_sequence}->{$deletion_coord} = $length;

      $all_deletions{$deletion_coord} = $length
	unless ((defined $all_deletions{$deletion_coord} &&
		 $all_deletions{$deletion_coord} > $length));

      push @{$deletion_tracking{$deletion_coord}},
	[$length, $evidence_sequence->name];
    }
  }

  my @all_deletions = sort {$a <=> $b} keys %all_deletions;

  $self->{_total_inserted_deletions} = 0; # Initialise this method-less value.

  foreach my $deletion_coord (@all_deletions) {

    my $length = $all_deletions{$deletion_coord};

    # Add a gap to the genomic and exonic sequence

    $self->_genomic_sequence->insert_gap($deletion_coord, $length);
    $self->{_total_inserted_deletions}++; # Naughty, but fast.
    $self->_exon_nucleotide_sequence->insert_gap($deletion_coord, $length);
    $self->_exon_protein_translation->insert_gap($deletion_coord, $length);

    foreach my $evidence_name (keys %evidence_sequence_hash) {

      # Propagate any needed gaps into the aligned evidence
      # sequences

      if (! $deletion_sets{$evidence_name}->{$deletion_coord}) {
	$evidence_sequence_hash{$evidence_name}->insert_gap($deletion_coord, $length)
      } elsif (defined $deletion_sets{$evidence_name}->{$deletion_coord} && 
	       $deletion_sets{$evidence_name}->{$deletion_coord} < $length){

	my $insert_coord = 
	  $deletion_coord + $deletion_sets{$evidence_name}->{$deletion_coord};

	my $deletion_length_difference = 
	  $length - $deletion_sets{$evidence_name}->{$deletion_coord};

	$evidence_sequence_hash{$evidence_name}->insert_gap($insert_coord,
							    $deletion_length_difference)
      }
    }

    # Increment all deletion coords that are greater than this one.
    # There quite possibly is a better way to do this.

      # @all_deletions and %all_deletions

    my %new_all_deletions;

    for (my $i = 0; $i < scalar @all_deletions; $i++) {
      if ($all_deletions[$i] > $deletion_coord) {
        $all_deletions[$i] += $length;
	$new_all_deletions{$all_deletions[$i]} = 
	  $all_deletions{$all_deletions[$i] - $length};
      }
    }

    %all_deletions = %new_all_deletions;


      # %deletion_sets

    my %new_deletion_sets;
    foreach my $evidence_name (keys %evidence_sequence_hash) {
      my @coords = keys %{$deletion_sets{$evidence_name}};
      for (my $i = 0; $i < scalar @coords; $i++) {
	if ($deletion_sets{$evidence_name}->{$coords[$i]}){
	  $new_deletion_sets{$evidence_name}->{$coords[$i]+$length} = 
	    $deletion_sets{$evidence_name}->{$coords[$i]};
	} else {
	  $new_deletion_sets{$evidence_name}->{$coords[$i]} = 
	    $deletion_sets{$evidence_name}->{$coords[$i]};
	}
      }
    }
    %deletion_sets = %new_deletion_sets;


      # %deletion_tracking

    my %new_deletion_tracking;
    my @coords = keys %deletion_tracking;
    for (my $i = 0; $i < scalar @coords; $i++) {
      if ($coords[$i] > $deletion_coord){
	$new_deletion_tracking{$coords[$i]+$length} = $deletion_tracking{$coords[$i]}
      } else {
	$new_deletion_tracking{$coords[$i]} = $deletion_tracking{$coords[$i]}
      }
    }
    %deletion_tracking = %new_deletion_tracking;
  }
### Used by Dan for conserved indel finding project.  Please leave.
  # Print the locations of our deletions, handy for finding conserved gaps and
  # frameshifts
#
#  foreach my $tracked_deletion (sort {$a <=> $b} (keys %deletion_tracking)){
#    print STDOUT $tracked_deletion . "\t" . 
#      scalar @{$deletion_tracking{$tracked_deletion}} . "\t";
#    foreach my $deletion (@{$deletion_tracking{$tracked_deletion}}){
#      print STDOUT $deletion->[1] ." (" . $deletion->[0] . ")  ";
#    }
#    print STDOUT "\n";
#  }
###
  # Put our working alignments somewhere handy

  $self->_working_alignment('genomic_sequence', $self->_genomic_sequence);
  $self->_working_alignment('exon_nucleotide',  $self->_exon_nucleotide_sequence);

  if ($self->_translatable
      &&(($type eq 'protein')||($type eq 'all'))) {
    $self->_working_alignment('exon_protein', $self->_exon_protein_translation);
  }

  # Place our finalised evidence sequence alignments into the right place.

  foreach my $evidence_key (keys %evidence_sequence_hash) {
    $self->_working_alignment('evidence', 
			      $evidence_sequence_hash{$evidence_key});
  }

  # If introns are to be truncated, do this now.

  if ($truncate_introns) {
    $self->_truncate_introns;
  }

  # Set flag to indicate that the alignment has been computed.

  $self->_is_computed($type, 1);

  return 1;
}


=head2 _create_Alignment_object

  Arg [1]    : String.  Alignment type ('all', 'nucleotide' or 'protein').
  Example    : $self->_create_Alignment_object('all')
  Description: Builds the final Bio::EnsEMBL::Pipeline::Alignment object
               and and is responsible for adding each aligned sequence.
  Returntype : Bio::EnsEMBL::Pipeline::Alignment
  Exceptions : none
  Caller     : retrieve_alignment

=cut

sub _create_Alignment_object {
  my ($self, $type) = @_;

  my $alignment = Bio::EnsEMBL::Pipeline::Alignment->new(
			      '-name' => 'evidence alignment');

  $alignment->add_sequence($self->_working_alignment('genomic_sequence'));
  $alignment->add_sequence($self->_working_alignment('exon_nucleotide'));

  if ($self->_translatable
      &&(($type eq 'protein')||($type eq 'all'))) {
    $alignment->add_sequence($self->_working_alignment('exon_protein'));
  }

  foreach my $evidence_sequence (@{$self->_working_alignment('evidence')}){
    $alignment->add_sequence($evidence_sequence);
  }

  return $alignment;
}


=head2 _truncate_introns

  Args       : none
  Example    : $self->_truncate_introns
  Description: This method is responsible for removing intronic regions
               from the overall sequence alignment.
  Returntype : int
  Exceptions : Warns (via info) if the truncation of any intron causes
               part of the evidence sequence to be discarded.
  Caller     : _align

=cut

sub _truncate_introns {
  my ($self) = @_;

  # Get sequences from the working alignment.

  my @sequences;

  push (@sequences, $self->_working_alignment('genomic_sequence'));
  push (@sequences, $self->_working_alignment('exon_nucleotide'));

  if (($self->_type eq 'all')||($self->_type eq 'protein')){
    push (@sequences, $self->_working_alignment('exon_protein'));
  }

  foreach my $aligned_seq (@{$self->_working_alignment('evidence')}){
    push (@sequences, $aligned_seq);
  }

  if ($self->_working_alignment('unaligned')){
    foreach my $unaligned_seq (@{$self->_working_alignment('unaligned')}){
      push (@sequences, $unaligned_seq);
    }
  }

  # Determine where the introns are.

  my @coordinates;
  foreach my $exon (@{$self->_transcript->get_all_Exons}){

    my $intron_coord_1 = $exon->start - 1;
    my $intron_coord_2 = $exon->end + 1;

    push(@coordinates, $intron_coord_1, $intron_coord_2);

  }

  # Get locations of gaps in genomic sequence.

  my $genomic_gaps = $self->_genomic_sequence->all_gaps;

  # Work through the list of genomic gaps and increment all
  # coordinates that are greater than a gap position.

  foreach my $gap_coord (@$genomic_gaps) {
    for (my $i = 0; $i < scalar @coordinates; $i++) {
      $coordinates[$i]++ if $coordinates[$i] >= $gap_coord;
    }
  }

  # Sort in reverse, such that we splice out introns from the
  # 3-prime end first.  This means we dont have to adjust the
  # coordinates of the unspliced introns each time an upstream
  # intron is removed.

  @coordinates = sort {$b <=> $a} @coordinates;

  shift @coordinates;
  pop @coordinates;

  # Splice introns from each sequence.

  foreach my $align_seq (@sequences){

    my $seq = $align_seq->seq;

    my @working_coordinates = @coordinates;

  INTRON:
    while (@working_coordinates){
      my $intron_end = (shift @working_coordinates) - 1;
      my $intron_start = (shift @working_coordinates);

      next INTRON unless ($intron_start + $self->_padding + 22) < ($intron_end - $self->_padding - 22);

      my $offcut = substr($seq, $intron_start + $self->_padding - 1, 
			  ($intron_end - 2*$self->_padding));

      if ($self->_padding == 0) {
	$seq = 
	  substr($seq, 0, $intron_start - 1) . 
	    substr($seq, $intron_end + 1, length($seq));
      } else {
	$seq = 
	  substr($seq, 0, $intron_start + $self->_padding - 1) . 
	    '---intron-truncated---' . 
	      substr($seq, ($intron_end - 2*$self->_padding + 1), length($seq));
      }

      if (($offcut =~ /[atgc]/)&&($align_seq->name ne 'genomic_sequence')) {
	info("Truncating intron sequences has caused some aligned evidence\n" .
	     "to be discarded.  Try again with a higher amount of padding\n".
	     "around the exon sequences.\n");
      }
    }
    $align_seq->seq($seq);
  }

  return 1;
}


=head2 _is_computed

  Arg [1]    : String.  Alignment type ('all', 'nucleotide' or 'protein').
  Arg [2]    : int (1 or 0).  Set to 1 if alignment has been previously 
               computed.
  Example    : $self->_is_computed('nucleotide', 1)
  Description: Stores information about the type of alignment that may
               have been previously computed.
  Returntype : int
  Exceptions : Warns if alignment previous computed, but of another type.
               Throws if alignment type is unknown.
  Caller     : retrieve_alignment, _align

=cut

sub _is_computed {
  my ($self, $type, $value) = @_;

  # Quick check of whether alignment has been computed
  # - if it hasnt
  if ((!defined $type)&&(!$self->_type)) {
    return 0;
  }

  # - if it has
  if ((!defined $type)&&($self->_type)) {
    return 1;
  }

  # Paranoid initialisation

  if (!defined $self->{_is_computed}) {
    $self->{_is_computed} = 0;
  }

  # Check whether an alignment of a specific type
  # has been run.
  if ((!defined $value)&&($self->{_is_computed})&&($type ne $self->_type)) { 
    warning("Alignment has been previously computed, but was not\n" .
	    "of the same type.  The previously computed alignment\n" . 
	    "type was \'" . $self->_type . "\'.\n");

    return 0; 
  }

  if (defined $value && $value > 0) {

    if ((defined $type)
	&&(($type eq 'nucleotide')||
	   ($type eq 'protein')||
	   ($type eq 'all'))){
      $self->_type($type);
    } else {
      throw("Unknown alignment type.  Can be nucleotide, protein or all.\n");
      return 0;
    }

    $self->{_is_computed} = 1;
  }

  return $self->{_is_computed};
}


=head2 _type

  Arg [1]    : String.  Alignment type ('all', 'nucleotide' or 'protein').
  Example    : $self->_type('all')
  Description: Stores the type of alignment that has been computed.
  Returntype : String
  Exceptions : Throws if alignment type is unknown
  Caller     : _is_computed, _working_alignment, _corroborating_sequences,
               _truncate_introns, retrieve_alignment, new

=cut

sub _type {
  my ($self, $type) = @_;

  if (defined $type) {
    unless ($type eq 'all' ||
	    $type eq 'nucleotide' ||
	    $type eq 'protein') {
      throw("Type of alignment to retrieve not recognised.  Please use one " . 
	    "of \'all\', \'nucleotide\' or \'protein\'.");
    }

    $self->{_computed_type} = $type;
  }

  throw("Type of alignment to retrieve not specified.  Please use one " . 
	"of \'all\', \'nucleotide\' or \'protein\'.") 
    unless $self->{_computed_type};

  return $self->{_computed_type};
}


### Alignment information handling methods ###


=head2 _working_alignment

  Arg [1]    : String.  Alignment sequence identity (or 'slot').  Should 
               be one of 'genomic_sequence', 'exon_protein', 
               'exon_nucleotide', 'evidence', or 'unaligned'.

  Arg [2]    : Bio::EnsEMBL::Pipeline::Alignment::AlignmentSeq
               OR string 'empty'.  Pass an align seq to add it to the 
               slot or pass 'empty' to clear the slot.
  Example    : $self->('evidence', $align_seq) or perhaps 
               $self->('evidence', 'empty')
  Description: Stores the different types of alignment sequences
               while the alignment is being built.
  Returntype : Bio::EnsEMBL::Pipeline::Alignment::AlignmentSeq or
               a list of these objects ('evidence' is an array, while
               each of the other slots are single sequences).  Returns
               0 if slot is empty.
  Exceptions : Throws if slot type is unrecognised.
               Throws if sequence is not a 
                Bio::EnsEMBL::Pipeline::Alignment::AlignmentSeq.
  Caller     : _create_Alignment_object, _truncate_introns, _align

=cut

sub _working_alignment {
  my ($self, $slot, $align_seq) = @_;

  unless (defined $slot && 
      (($slot eq 'genomic_sequence')
       ||($slot eq 'exon_protein')
       ||($slot eq 'exon_nucleotide')
       ||($slot eq 'evidence')
       ||($slot eq 'unaligned'))){
    throw("Was trying to retrieve or write working alignments to "
	  . "a slot that isnt allowed ($slot)");
  }

  if (defined $slot && defined $align_seq && $align_seq eq 'empty'){
    undef $self->{_working_alignment_array}->{$slot};
    return 1
  }

  if (defined $slot && defined $align_seq){

    unless ($align_seq->isa("Bio::EnsEMBL::Pipeline::Alignment::AlignmentSeq")){
      throw("Sequence passed to _working alignment was not an " . 
	    "AlignmentSeq, it was a [$align_seq]")
    }

    push (@{$self->{_working_alignment_array}->{$slot}}, $align_seq);

  } elsif (defined $slot && !defined $align_seq) {

    my $slot_resident =  $self->{_working_alignment_array}->{$slot};

    if ((($slot eq 'genomic_sequence')||($slot eq 'exon_protein')||($slot eq 'exon_nucleotide')) 
	&& defined $slot_resident && scalar @$slot_resident == 1) {
      return $slot_resident->[0];
    }

    # Sort evidence so that nucleotide and protein sequences 
    # are returned in blocks of similar sequence type.
    if ($slot eq 'evidence' && $self->_type eq 'all') {

      throw "All the evidence for this gene has been thrown away or skipped.  Yuck."
	unless $slot_resident;

      @$slot_resident = sort {$a->type cmp $b->type} @$slot_resident;
    }

    return $slot_resident;
  }

  return 0;
}

##### Fiddlings with Slices and Transcripts #####


=head2 _transcript

  Arg [1]    :
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub _transcript {
  my $self = shift;

  if (@_){
    $self->{_transcript} = shift;
  }

  unless (defined $self->{_transcript} &&
	  $self->{_transcript}->isa("Bio::EnsEMBL::Transcript")){
    throw("Problem with transcript.  It is either unset or is not a " . 
	  "transcript.  It is a [".$self->{_transcript}."].");
  }

  return $self->{_transcript}
}


=head2 _slice

  Arg [1]    :
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub _slice {
  my ($self, $object) = @_;

  my $transcript;

  if (defined $object && 
      $object->isa("Bio::EnsEMBL::Slice")){
    $self->{_slice} = $object
  } elsif (defined $object && 
	   $object->isa("Bio::EnsEMBL::Transcript")) {
    $transcript = $object
  }

  if (! defined $self->{_slice} && defined $transcript) {

    my $slice_adaptor = $self->_db->get_SliceAdaptor;

    # Ideally, stable_ids should be used to fetch things - just in 
    # case the transcript comes from a different db to our current
    # db.  Otherwise, fall over to using dbID.

    if ($transcript->stable_id){

      $self->{_slice} = 
	$slice_adaptor->fetch_by_transcript_stable_id($transcript->stable_id, 
						      $self->_padding);
    } else {

      $self->{_slice} = 
	$slice_adaptor->fetch_by_transcript_id($transcript->dbID, 
					       $self->_padding);

    }

    # Check the orientation of our transcript.  If it is on the reverse
    # strand invert our slice.

    if ($transcript->strand == -1){
      $self->{_slice} = $self->{_slice}->invert;
    }
  }

  return $self->{_slice};
}


=head2 _strand

  Arg [1]    :
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub _strand {
  my $self = shift;

  return $self->{_transcript}->strand;
}


=head2 _seq_fetcher

  Arg [1]    :
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub _seq_fetcher {
  my ($self, $fetcher) = @_;

  if (defined $fetcher) {

    $self->{_seq_fetcher} = $fetcher;

    return 1;
  }

  return $self->{_seq_fetcher};
}


### Sequence handling methods ###


=head2 _corroborating_sequences

  Arg [1]    : 
  Example    : 
  Description: Uses all supporting features to generate a list of
               unique evidence sequences.
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub _corroborating_sequences {
  my $self = shift;

  # Work through supporting features attached to each exon and 
  # create a sequence for each.

  my %name_vs_type;

 FEATURE:
  foreach my $base_align_feature (@{$self->_all_supporting_features}){

    # Filter sequences by type, where necessary
    if ((($self->_type eq 'nucleotide')
	 &&($base_align_feature->isa("Bio::EnsEMBL::DnaPepAlignFeature")))
	||(($self->_type eq 'protein')
	   &&($base_align_feature->isa("Bio::EnsEMBL::DnaDnaAlignFeature")))){
      next FEATURE;
    }

    # Filter by evidence identity, if required.
    if ((defined $base_align_feature->percent_id)
	&&($base_align_feature->percent_id < $self->{_evidence_identity_cutoff})) {
      next FEATURE;
    }

    # Keep a hash of sequence names versus sequence type.
    unless (defined $name_vs_type{$base_align_feature->hseqname}){
      if ($base_align_feature->isa("Bio::EnsEMBL::DnaPepAlignFeature")){
	$name_vs_type{$base_align_feature->hseqname} = 'protein';
      } elsif ($base_align_feature->isa("Bio::EnsEMBL::DnaDnaAlignFeature")){
	$name_vs_type{$base_align_feature->hseqname} = 'nucleotide';
      }
    }

    # Hidden down here is the actual method call that constructs
    # a sequence fragment from the evidence.  This is just a single
    # ungapped feature sequence which will be merged with
    # other fragments from the same sequence (ie. ungapped sequences
    # to form a single gapped sequence).  The merging process is
    # managed from the _build_evidence_seq method.

    $self->_build_evidence_seq($base_align_feature);
  }

  # Retrieve all built evidence sequences and before returning
  # attach the correct sequence type to each.

  my $corroborating_sequences = $self->_feature_sequences;

  for (my $i = 0; $i < scalar @$corroborating_sequences; $i++) {
    $corroborating_sequences->[$i]->type($name_vs_type{$corroborating_sequences->[$i]->name});
  }

  unless (defined $corroborating_sequences &&
	  scalar @$corroborating_sequences > 0) {
    warning("There are no displayable supporting features for this transcript [" .
	    $self->_transcript->stable_id . "].  " .
	    "Try setting the -type attribute to 'all', instead of just " . 
	    "'protein' or 'nucleotide'.");
    return 0
  }

  return $corroborating_sequences;
}


=head2 _build_evidence_seq

  Arg [1]    :
  Example    : 
  Description: Takes a single align feature and constructs an aligned
               sequence of the correct portion of the hit sequence
               including all of the right gaps.  Creates a reference
               array that allows the correct insertion of gaps in the
               genomic sequence.
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub _build_evidence_seq {
  my ($self, $base_align_feature) = @_;

  # For convenience, localise a few coordinate variables.

  my $hstart  = $base_align_feature->hstart;
  my $hend    = $base_align_feature->hend;
  my $cigar   = $base_align_feature->cigar_string;
  my $hstrand = $base_align_feature->hstrand;
  my @cigar_instructions = $self->_cigar_reader($cigar);

  # Turn the cigar line around if the slice is on the
  # reverse strand.
#  if ($self->_slice->strand == -1){
#    @cigar_instructions = reverse @cigar_instructions;
#  }

  # Fetch our sequence from the cache.  If the sequence
  # is missing it means that it could not be fetched and
  # we will have to ignore it.

  my $fetched_seq = $self->_fetch_sequence($base_align_feature->hseqname);

  if ( ! $fetched_seq) {
    info("Error fetching sequence [" . 
	 $base_align_feature->hseqname . 
	 "].  Ignoring.");

    return 0;
  } elsif (($fetched_seq->seq eq '')||
	   ($fetched_seq->seq eq '-' x (length($fetched_seq->seq)))) {
    warning("Error fetching sequence [" . 
	    $base_align_feature->hseqname . 
	    "].  Sequence is an empty string.");

    return 0;
  }

  # If this is a protein align feature, splice out the relevant 
  # portion of the evidence sequence, then pad amino acids with gaps
  # to make them comparable to nucleotide coords.

  if ($base_align_feature->isa("Bio::EnsEMBL::DnaPepAlignFeature")){

    # Pop in a little sanity check here.  Protein features can not be
    # reverse complimented, hence they must have the same strand as 
    # the gene/slice (which is the same strand as the transcript).

    unless ($hstrand == $self->_strand) {
      print STDERR "Strand : " . $self->_strand . "\tHstrand : " . $hstrand . "\n";
      warning("Protein align feature [". $base_align_feature->hseqname . 
	      "] is reversed with respect to transcript - and I cant reverse " . 
	      "compliment an amino acid sequence.");
      return 0
    }

    # Continue with sequence fiddling.

    my $aa_seq = $fetched_seq->seq;

    if (($hstart > length($aa_seq))||
	($hend   > length($aa_seq))){
      warning("Evidence sequence coordinates lie outside " .
	      "the bounds of that sequence.  Data problem.");
    }

    $fetched_seq = substr($aa_seq, ($hstart - 1), ($hend - $hstart + 1));

    unless ($fetched_seq ne ''){
      throw("Error building amino acid sequence [" . $base_align_feature->hseqname 
	    . "].  Feature does not seem to lie within the bounds of the " . 
	    "evidence sequence.");
    }

    if ($self->_three_letter_aa) {
      $fetched_seq =~ s/(.)/$self->{_aa_names}{$1}/g;
    } else {
      $fetched_seq =~ s/(.)/$1\-\-/g;
    }
  }

  # If we have a dna align feature, extracting the correct portion
  # of the hit sequence is a bit easier than the method required
  # for a protein sequence.

  if ($base_align_feature->isa("Bio::EnsEMBL::DnaDnaAlignFeature")) {

    $fetched_seq = $fetched_seq->seq;

    # Splice out the matched region of our feature sequence
    $fetched_seq = substr($fetched_seq, ($hstart - 1), ($hend - $hstart + 1));

    # The coordinates for dna_align_features are stored with reference
    # to the way the sequence is was found in dbEST or whatever.  Hence,
    # if the EST is the reverse strand, it is necessary to chop out our 
    # fragment first, then reverse compliment it.

    if ($hstrand == -1) {
      my $backwards_seq = Bio::Seq->new('-seq' => $fetched_seq);

      my $forwards_seq = $backwards_seq->revcom;

      $fetched_seq = $forwards_seq->seq;
    }
  }

  # Add the needed insertion gaps to our supporting
  # feature sequence.  We only add gaps to the supporting
  # sequences at this stage.  Once we have all our supporting
  # sequences constructed and determined where all the 
  # 'deletion' gaps lie in these, we can transpose this 
  # information onto our genomic and exonic sequences - this is
  # done in $self->_align.

  my $hit_position = 1;
  my @deletions;
  my $feature_deletions  = 0;
  my $feature_insertions = 0;
  my $deletions_upstream_of_slice = 0;

    # Note to self - all genomic coordinates are actually the coordinates
    # within the slice.  This slice is not in chromosomal coordinates.  As
    # the slice was derived using a transcript, the slice will be oriented
    # in the same direction as the transcript.  As the transcript alignment 
    # will be displayed 5'-3' the slice will be in the correct orientation
    # with respect to the transcript and thus the genomic strand can be 
    # treated as being on the forward strand.

  my $slice_position = $base_align_feature->start;
  my $slice_rev_position = $self->_slice->length - $base_align_feature->end + 1;

  foreach my $instruction (@cigar_instructions) {

    if ($instruction->{type} eq 'I') {

      my $gap = '-' x $instruction->{length};

      if ($hit_position <= length $fetched_seq){
	# Splice in the insertion.
	$fetched_seq = substr($fetched_seq, 0, $hit_position - 1) . 
	  $gap . substr($fetched_seq, $hit_position - 1);

      } else {
	throw("Gap in sequence [" . $base_align_feature->hseqname . 
	      "] lies outside feature.  Code problem.")
      }

      $hit_position   += $instruction->{length};
      $slice_position += $instruction->{length};
      $slice_rev_position += $instruction->{length};
      $feature_insertions += $instruction->{length};

    } elsif ($instruction->{type} eq 'M') {

      $hit_position   += $instruction->{length};
      $slice_position += $instruction->{length};
      $slice_rev_position += $instruction->{length};

    } elsif ($instruction->{type} eq 'D') {

      # Deletions that fall upstream of the slice must
      # be counted and handled differently.
      if ($slice_position <= 0){
	$deletions_upstream_of_slice++;
	next
      }

      # Regular deletion tracking.
      if ($slice_position <= $self->_slice->length) {
	push @deletions, [$slice_position, $instruction->{length}];
	$feature_deletions += $instruction->{length};
      }

    }
  }

  if (scalar @deletions) {
    $self->_there_are_deletions(1);
  }

  # This little section of code handles any sequence that
  # overshoots the beginning or end of our slice.  Chop.

  if ($base_align_feature->start < 0 || 
      $base_align_feature->start > $self->_slice->length ||
      $base_align_feature->end > $self->_slice->length || 
      $base_align_feature->end < 0) {
    info("Feature [". $base_align_feature->hseqname . " start:" . 
	    $base_align_feature->start . " end:" . $base_align_feature->end 
	    ."] extends\npast the start or end of genomic slice.  Truncating\n" . 
	    "overhanging sequence.");

    if (($base_align_feature->start < 0 && 
	 $base_align_feature->end < 0)||
	($base_align_feature->start > $self->_slice->length &&
	 $base_align_feature->end > $self->_slice->length)) {

      info("Feature [" . $base_align_feature->hseqname . 
	   "] lies completely outside the bounds of Slice.  Chuck.");
      return 0

    } else {
      my $start_overshoot = 0;
      if ($base_align_feature->start < 0) {
	info("Feature [". $base_align_feature->hseqname .
	     "] extends past the beginning of the slice.  Truncating.\n");
	$start_overshoot = 
	  ($base_align_feature->start * -1) + $deletions_upstream_of_slice;

	$fetched_seq = substr($fetched_seq, $start_overshoot + 1);
      }

      if ($base_align_feature->end > $self->_slice->length) {
	info("Feature [". $base_align_feature->hseqname .
	     "] extends past the end of the slice.  Truncating.\n");
	my $end_overshoot = $base_align_feature->end - $self->_slice->length - 1;

	$fetched_seq = substr($fetched_seq, 0,
		(length $fetched_seq) - $start_overshoot - $end_overshoot - 1);
      }
    }
  }

  # Here we are actually building the sequence that will
  # align to our slice - each feature sequence is grafted onto
  # a sequence of the same length as the genomic slice.
  # Successive features of from the same source id are each
  # grafted into the same sequence, effectively merging these
  # features.  This does cause compulsory merging of features
  # with the same id, but by doing this here a very significant
  # performance boost is achieved.  Post-hoc merging was found
  # to be very slow.

  my $genomic_start;

  $genomic_start = 
    $base_align_feature->start > 0 ? $base_align_feature->start - 1 : 0;

  if ($genomic_start <= $self->_slice->length){

    # Use the _feature_sequence repository to retrieve this sequence.
    my $feat_align_seq = 
      $self->_feature_sequence($base_align_feature->hseqname);

    my $feature_sequence = $feat_align_seq->seq;

     # A messy and probably slow bit of code that figures out how
     # many deletions in the genomic sequence exist upstream
     # of the feature to be inserted.

    my $all_deletions   = $feat_align_seq->deletion_coords;

    my $upstream_deletion_count = 0;
    my $total_deletion_count    = $feature_deletions;

    foreach my $deletion (@$all_deletions) {

      if ($deletion <= $genomic_start) {
	$upstream_deletion_count += 
	  $feat_align_seq->deletion_hash->{$deletion};
      }

      $total_deletion_count += 
	$feat_align_seq->deletion_hash->{$deletion};
    }

    my $deletion_gap = '-' x $feature_deletions;

    my $insert_end = 
      $genomic_start + $upstream_deletion_count + length($fetched_seq);

    my $length_to_end = 
      $self->_slice->length + $total_deletion_count - $insert_end;

      # More mess to clean up the tail end of the sequence
      # where insertions in the feature cause the sequence to
      # align to genomic regions off the end of our slice.
      # These features are truncated to fit within the
      # alignment slice.

    unless (($insert_end + $feature_deletions) < 
	    ($self->_slice->length + $total_deletion_count)) {

      my $overrun = 
	$insert_end + $feature_deletions - 
	  $self->_slice->length - $total_deletion_count;

      $insert_end = 
	$self->_slice->length + $total_deletion_count - $feature_deletions;

      $length_to_end = 0;

      if (length($deletion_gap) > $overrun){ 
        $deletion_gap = '-' x (length($deletion_gap) - $overrun);
      } else {
        $overrun -= length($deletion_gap);
        $deletion_gap = '';

        $fetched_seq = 
	  substr($fetched_seq, 0, (length($fetched_seq) - $overrun));
      }
    }

      # Okay then, finally on with the seemingly simple task
      # of sticking all these bits of sequence together.

    $feature_sequence = 
      substr($feature_sequence, 0, ($genomic_start + $upstream_deletion_count)) 
      . $fetched_seq . $deletion_gap . 
	substr($feature_sequence, $insert_end, $length_to_end);

    if (length($feature_sequence) != 
	$self->_slice->length + $total_deletion_count){
      throw("Building a feature sequence [" . $base_align_feature->hseqname . 
            "] that is longer than our genomic slice.");
    }

    # Pass the merged sequence back to the _feature_sequence handler.
    $self->_feature_sequence($base_align_feature->hseqname,
			     $feature_sequence,
			     \@deletions);

  } else {
    info("Feature [". $base_align_feature->hseqname . " start:" . 
	 $base_align_feature->start . " end:" . $base_align_feature->end .
	 " strand:" . $base_align_feature->strand 
	 ."] starts beyond end of genomic slice.  This feature most probably " .
	 "aligns to an exon that is not part of this transcript.");
    return 0
  }

  return 1
}


=head2 _feature_sequence

  Arg [1]    :
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut


sub _feature_sequence {
  my ($self, $name, $sequence, $deletion_coords) = @_;

  unless (defined $name) {
    throw("Cant do anything without a sequence name.")
  }

  unless (defined $self->{_feature_sequences}){
    $self->{_feature_sequences} = {}
  }

  if (! defined $self->{_feature_sequences}->{$name}) {
    if (! defined $sequence) {
      $sequence = '-' x $self->_slice->length;
    }

    $self->{_feature_sequences}->{$name} = 
      Bio::EnsEMBL::Pipeline::Alignment::AlignmentSeq->new(
        '-name' => $name,
        '-seq'  => $sequence);

    $sequence = undef;

  } elsif (defined $sequence) {
    $self->{_feature_sequences}->{$name}->seq($sequence);

    if (defined $deletion_coords){
      foreach my $deletion (@$deletion_coords){
	$self->{_feature_sequences}->{$name}->add_deletion($deletion->[0], 
							   $deletion->[1])
      }
    }
  }

  return $self->{_feature_sequences}->{$name}
}


=head2 _feature_sequences

  Arg [1]    :
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub _feature_sequences {
  my $self = shift;

  my @feat_seqs = values %{$self->{_feature_sequences}};

  return \@feat_seqs
}


=head2 _print_tabulated_coordinates

  Arg [1]    :
  Example    :
  Description: A debugging subroutine for printing the coordinate locations 
               of our gene, exons, feature and evidence.
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub _print_tabulated_coordinates {
  my $self = shift;

  print STDERR 
    "Slice :      length    - " . $self->_slice->length           . "\n" .
    "             chr       - " . $self->_slice->seq_region_name  . "\n" .
    "             chr start - " . $self->_slice->start            . "\n" .
    "             chr end   - " . $self->_slice->end              . "\n" .
    "             strand    - " . $self->_slice->strand           . "\n";

  print STDERR "\n";

  my $exons = $self->_transcript->get_all_Exons;

  print STDERR 
    "Transcript : id         - " . $self->_transcript->stable_id . "\n" . 
    "             coding     - (start) " . $self->_transcript->coding_region_start . " (end) " . 
      $self->_transcript->coding_region_end . "\n" .
    "             exon count - " . scalar @$exons . "\n";
  foreach my $exon (@$exons) {
    print STDERR 
    "             exon       - (start) " . $exon->start . 
      " (end) " . $exon->end . " (strand) " . $exon->strand . "\n";
  }

  print STDERR "\n";

  my $supporting_features = $self->_all_supporting_features;

  print STDERR 
    "Evidence : total features - " . @$supporting_features . "\n";
  foreach my $feature (@$supporting_features){

    print STDERR
    "           Feature : " . $feature->hseqname . "\n" .
    "                     (start)  "  . $feature->start . "\t(end)  " . 
      $feature->end . "\t(strand)  "  . $feature->strand . "\n" . 
    "                    (hstart)  "  . $feature->hstart . "\t(hend) " . 
      $feature->hend . "\t(hstrand) " . $feature->hstrand . "\n" . 
    "                     (CIGAR)  "  . $feature->cigar_string . "\n";

  }

  return 1
}


=head2 _genomic_sequence

  Arg [1]    :
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub _genomic_sequence {
  my ($self) = @_;

  if (!defined $self->{_genomic_sequence}) {

    my $genomic_sequence;

    $genomic_sequence = $self->_slice->seq;

    $self->{_genomic_sequence} = 
      Bio::EnsEMBL::Pipeline::Alignment::AlignmentSeq->new(
	 '-seq'  => $genomic_sequence,
	 '-name' => 'genomic_sequence',
	 '-type' => 'nucleotide');

  }

  return $self->{_genomic_sequence};
}


=head2 _exon_nucleotide_sequence

  Arg [1]    :
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub _exon_nucleotide_sequence {
  my ($self) = @_;

  if (!defined $self->{_exon_nucleotide_sequence}) {

    # Make a sortable exon list.

    my @exon_list;

    foreach my $exon (@{$self->_transcript->get_all_Exons}) {
      my $exon_start;

      $exon_start = $exon->start - 1;

      push @exon_list, [$exon, $exon_start];
    }

    # Sort exon list.

    @exon_list = sort {$a->[1] <=> $b->[1]} @exon_list;


    # Build our string.

    my $genomic_exon_seq = '';
    my $gap_start = 1;

    foreach my $exon_item (@exon_list) {

      my $exon       = $exon_item->[0];
      my $exon_start = $exon_item->[1];
      my $exon_seq   = $exon->seq->seq;

      my $gap_length = $exon_start - $gap_start + 1;

      $genomic_exon_seq .= ('-' x $gap_length) . $exon_seq;

      $gap_start += $gap_length + length($exon_seq);

    }

      # Add final gap.

    my $final_gap_length = $self->_slice->length - $gap_start + 1;

    $genomic_exon_seq .= ('-' x $final_gap_length);


    # Build and store AlignmentSeq object for this sequence.

    $self->{_exon_nucleotide_sequence} 
      = Bio::EnsEMBL::Pipeline::Alignment::AlignmentSeq->new(
	   '-seq'  => $genomic_exon_seq,
	   '-name' => 'exon_sequence',
	   '-type' => 'nucleotide');
  }

  return $self->{_exon_nucleotide_sequence};
}


=head2 _exon_protein_translation

  Arg [1]    :
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub _exon_protein_translation {

  my ($self) = @_;

  if (! defined $self->{_exon_protein_translation}) {

    # Make a sorted list of exons, arranged in coding order.

    my @exon_list;

    foreach my $exon (@{$self->_transcript->get_all_Exons}){

      # Worry about parts of the exon that dont translate,
      # but only if they are in the final exon.
      my $end_coding_exon = 0;
      if (($exon->start < $self->_transcript->coding_region_end)&&
	  ($exon->end  >= $self->_transcript->coding_region_end)) {
	$end_coding_exon = 1;
      }

      my $start_exon = 0;
      if ($exon == $self->_transcript->start_Exon) {
	$start_exon = 1;
      }

      # Derive a translation of this exon peptide.

      my $seq     = $exon->peptide($self->_transcript);
      my $peptide = $seq->seq;

      next unless (length($peptide));

      $peptide =~ s/(.)/$1\-\-/g;


      # Worry about exon phases for a bit.  Make sure characters
      # from different exons are not included in this exon.

        # Remove preceeding characters belonging to the previous exon.

      if ($exon->phase == 1){
	$peptide = substr($peptide, 1);
      } elsif ($exon->phase == 2){
	$peptide = substr($peptide, 2);
      }

        # Remove trailing characters that dont belong in this exon.

      if ($exon->end_phase == 2){
	$peptide = substr($peptide, 0, length($peptide) - 1);

      } elsif ($exon->end_phase == 1){
	$peptide = substr($peptide, 0, length($peptide) - 2);
      }

      # Determining where to stick our peptide.

      my $peptide_genomic_start = $exon->end - length($peptide) + 1;

        # For all but the terminal exon it is possible to use the
        # end of the exon and work back.  A special treatment
        # it required for single exon genes.

      if ($start_exon && $end_coding_exon){ # Meaning, single exon gene.
	my $utr = $self->_transcript->five_prime_utr;
	my $utr_length = $utr ? length($utr->seq) : 0;
	$peptide_genomic_start = $exon->start + $utr_length;
      } elsif ($end_coding_exon) {
	$peptide_genomic_start = $exon->start;
      }

      push @exon_list, [$peptide_genomic_start, $peptide];

    }

    # Build the string representing the translated exons in a
    # genomic context.

    my $exon_translation_string = '';

    foreach my $exon_item (@exon_list){
      my $genomic_start = $exon_item->[0];
      my $peptide_seq   = $exon_item->[1];

      my $gap_length = $genomic_start - 1 - length($exon_translation_string);

      $exon_translation_string .= ('-' x $gap_length) . $peptide_seq;
    }

      # Add last gap.

    my $last_gap_length = $self->_slice->length - length($exon_translation_string);
    $exon_translation_string .= ('-' x $last_gap_length);

    # Convert to three letter AA names, if desired.

    if ($self->_three_letter_aa) {
      $exon_translation_string =~ s/([^\-])\-\-/$self->{_aa_names}{$1}/g;
    }

    # Build alignment seq from our sequences.

    $self->{_exon_protein_translation} = 
      Bio::EnsEMBL::Pipeline::Alignment::AlignmentSeq->new(
	 '-seq'  => $exon_translation_string,
	 '-name' => 'translated_exon_sequence',
	 '-type' => 'protein');
  }

  return $self->{_exon_protein_translation};
}


=head2 _all_supporting_features

  Arg [1]    :
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub _all_supporting_features {
  my $self = shift;

  # If there is a set of supporting features that are to be used
  # that are not already attached to the gene/transcript, they can
  # be added here.

  if (@_) {
    my $bafs = shift;

    foreach my $baf (@$bafs){
      throw("Evidence provided does not consist " . 
	    "of Bio::EnsEMBL::BaseAlignFeature")
	unless $baf->isa("Bio::EnsEMBL::BaseAlignFeature");

      $baf = $baf->transfer($self->_slice);

      push @{$self->{_all_supporting_features}}, $baf;
    }

  }

  # If there is not a set of external supporting features,
  # these can be yanked from the gene/transcript.

  unless ($self->{_all_supporting_features}){
    my $sfs = $self->_transcript->get_all_supporting_features;

    if (@$sfs) {
      push @{$self->{_all_supporting_features}}, @$sfs;
    } else {
      foreach my $exon (@{$self->_transcript->get_all_Exons}){
	push @{$self->{_all_supporting_features}}, 
	  @{$exon->get_all_supporting_features};
      }
    }

    $self->{_all_supporting_features} =
      $self->_remove_duplicate_features($self->{_all_supporting_features});
  }

  return $self->{_all_supporting_features}
}


=head2 _remove_duplicate_features

  Arg [1]    :
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub _remove_duplicate_features {
  my ($self, $features) = @_;

  unless (defined $features){
    throw("No features passed to filter.")
  }

  my %feature_lookup;
  my %feats_by_id;
  my @filtered_features;

 FEAT:
  foreach my $feat (@$features){

    my $feat_identifier 
      = $feat->hseqname . $feat->hstart . $feat->hend . $feat->hstrand;

    # Filter out identical matches.

    unless ($feature_lookup{$feat_identifier}){

      # Filter out overlapping matches - keep one, throw the
      # other away (although no consideration is given to
      # which might be the best to keep).

      foreach my $same_id_feat (@{$feats_by_id{$feat->hseqname}}) {
	if ($feat->end() >= $same_id_feat->start() &&
	    $feat->start() <= $same_id_feat->end()) {
	  next FEAT
	}
      }

      push @filtered_features, $feat;
      push(@{$feats_by_id{$feat->hseqname}}, $feat);
    }

    $feature_lookup{$feat_identifier}++
  }

  return \@filtered_features
}


### Methods that take care of sequence fetching and caching ###

=head2 _build_sequence_cache

  Arg [1]    :
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub _build_sequence_cache {
  my ($self) = @_;

  # Determine which sequences are likely to be needed.

  my %hash_of_accessions;

  foreach my $supporting_feature (@{$self->_all_supporting_features}){
    $hash_of_accessions{$supporting_feature->hseqname}++;
  }

  my @array_of_accessions = keys %hash_of_accessions;

  # Retrieve sequences.

  my $fetched_seqs;

  if ($self->_seq_fetcher->can("batch_fetch")){

    eval {
      $fetched_seqs = $self->_seq_fetcher->batch_fetch(@array_of_accessions);
    };

    if ($@){
      info("Not all evidence sequences could be pfetched.\n".
	      "Ignoring missing sequences.\n$@\n");
    }

  } else {

    foreach my $accession (@array_of_accessions){

      my $fetched_seq;

      eval {
	$fetched_seq = $self->_seq_fetcher->get_Seq_by_acc($accession);
      };

      if ($@) {
	warning("The seqfetcher is having trouble.\t$@\n");
      }

      push(@$fetched_seqs, $fetched_seq);

    }
  }

  # Build cache.

  foreach my $fetched_seq (@$fetched_seqs){

    next unless defined $fetched_seq;

    $self->{_fetched_seq_cache}->{$fetched_seq->accession_number} 
      = $fetched_seq;

  }

  $self->{_cache_is_built} = 1;
}


=head2 _fetch_sequence

  Arg [1]    :
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub _fetch_sequence {
  my ($self, $accession) = @_;

  $self->_build_sequence_cache 
    unless $self->{_cache_is_built};

  unless ($self->{_fetched_seq_cache}->{$accession}){
    info("Sequence $accession could not be retrieved from cache.");
  }

  return $self->{_fetched_seq_cache}->{$accession}; 
}


### Getter-Setters and Miscellaneous utilities ###


=head2 _db

  Arg [1]    :
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub _db {
  my $self = shift;


  if (@_) {
    $self->{_db_adaptor} = shift;
  }

  return $self->{_db_adaptor}
}


=head2 _padding

  Arg [1]    :
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut


sub _padding {
  my $self = shift;

  if (@_) {
    $self->{_padding} = shift;
  }

  return $self->{_padding} ? $self->{_padding} : 0;
}


=head2 _three_letter_aa

  Arg [1]    :
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut


sub _three_letter_aa {
  my $self = shift;

  if (@_) {
    $self->{_three_letter_aa} = shift;
  }
  if ($self->{_three_letter_aa}){
    $self->_set_aa_names;
  }


  return $self->{_three_letter_aa} ? $self->{_three_letter_aa} : 0;
}

=head2 _set_aa_names

  Arg [1]    :
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub _set_aa_names {
  my $self = shift;

  $self->{_aa_names} = {'A' => 'Ala',
                        'B' => 'Asx',
                        'C' => 'Cys',
                        'D' => 'Asp',
                        'E' => 'Glu',
                        'F' => 'Phe',
                        'G' => 'Gly',
                        'H' => 'His',
                        'I' => 'Ile',
                        'K' => 'Lys',
                        'L' => 'Leu',
                        'M' => 'Met',
                        'N' => 'Asn',
                        'P' => 'Pro',
                        'Q' => 'Gln',
                        'R' => 'Arg',
                        'S' => 'Ser',
                        'T' => 'Thr',
                        'V' => 'Val',
                        'W' => 'Trp',
                        'Y' => 'Tyr',
                        'Z' => 'Glx',
			'-' => '---',
			'X' => 'xxx'}
}


=head2 _translatable

  Arg [1]    :
  Example    : 
  Description: Toggle indicating whether translations are available.
  Returntype : 
  Exceptions : 
  Caller     : 

=cut


sub _translatable {
  my $self = shift;

  if (@_) {
    $self->{_translatable} = shift;
  }

  return $self->{_translatable};
}


=head2 _there_are_deletions

  Arg [1]    :
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut


sub _there_are_deletions {
  my $self = shift;

  if (@_) {
    $self->{_there_are_deletions} = shift;
  }

  return $self->{_there_are_deletions}
}




### CIGAR string handlers ###


=head2 _cigar_reader

  Arg [1]    :
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub _cigar_reader {
  my ($self, $cigar_string) = @_;

  my @cigar_array = split //, $cigar_string;

  my @cigar_elements;   # An array of hash references.

  my $current_digits;

  while (scalar @cigar_array) {

    my $next_char = shift @cigar_array;

    if ($next_char =~ /[MDI]/) {

      $current_digits = 1
	unless ($current_digits && $current_digits > 0);

      my %enduring_hash;
      $enduring_hash{type} = $next_char;
      $enduring_hash{length} = $current_digits;

      push (@cigar_elements, \%enduring_hash);

      $current_digits = '';

    } elsif ($next_char =~ /\d/) {

      $current_digits .= $next_char;

    } else {

      throw "There is an odd character in the CIGAR string retrieved from the database.\n" . 
	$cigar_string . "\n";
    }

  }

  return @cigar_elements;
}


1;
