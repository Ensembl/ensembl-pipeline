
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

# Quick start - use the following code if you want the 
# alignment of an ensembl transcript with the evidence 
# used to predict it:

my $alignment_tool = Bio::EnsEMBL::Pipeline::Alignment::EvidenceAlignment->new(
                           '-dbadaptor'         => $db,
			   '-seqfetcher'        => $seqfetcher,
			   '-transcript'        => $transcript);

my $alignment = $alignment_tool->retrieve_alignment('-type' => 'all');  
                             # Type can be 'all', 'nucleotide' or 'protein'

foreach my $align_seq (@$alignment){
  print $sequence->seq . "\n";
}


# Other alignment presentation options exist.  The intronic
# regions of the alignment can make it _very_ big.  Introns
# can be truncated thusly:

my $alignment = $align_tool->retrieve_alignment('-type'           => 'all',
						'-remove_introns' => 1,
						'-padding'        => 10);

# The '-padding' option specifies the amount of tailing
# sequence left attached to each 'exon'.  If you set this
# to zero be careful, as you could be truncating some
# sequence from the aligned evidence sequences.  If this 
# happens a warning is issued.

# Importantly, if you need to see non-aligned portions of 
# the evidence sequences use the following '-show_unaligned'
# option:

my $alignment = $align_tool->retrieve_alignment('-type'           => 'all',
						'-show_unaligned' => 1);


# This turns on an portion of code that generates sequences
# composed of all the fragments in the evidence sequence that
# havent been matched.  This is mainly useful to genebuilders 
# worrying about things that might be missed by various similarity 
# matching algorithms.  Note that for layout the fragments are 
# placed in the intron where they conceivably should go.  If 
# you are also trimming the intron sequences be careful that 
# you aren't unwittingly throwing these intronic fragments
# away.  A warning is raised if this happens.


# A few features of the actual alignment can be controled.
# The fasta line length and the amount of 5-prime and
# 3-prime padding can be stipulated:

my $alignment_tool = Bio::EnsEMBL::Pipeline::Alignment::EvidenceAlignment->new(
                           '-dbadaptor'         => $db,
			   '-seqfetcher'        => $pfetcher,
			   '-transcript'        => $transcript,
			   '-padding'           => 50,
			   '-fasta_line_length' => 60);


# Once the alignment is generated it is possible to calculate
# the identity of the best matching evidence for each exon.
# Coverage is also determined.

my $exon_identities = $alignment_tool->identity;

my $exon_counter = 1;
foreach my $exon_identity (@$exon_identities) {
  print "Exon $exon_counter\t" . 
    $exon_identity->[0] . "\t" .
      $exon_identity->[1] . "\t" .
	$exon_identity->[2] . "\t" .
	  $exon_identity->[3] ."\n";
  $exon_counter++;
}

# The identity scores are returned as a reference to an  array 
# of array references (I know, bleah - need a mini-object 
# perhaps).  Each referenced array represents the best identity 
# and coverage for an individual exon. The reference array has
# the format (nucleotide_identity, nucleotide_coverage,
# protein_identity, protein_coverage).  This could be easier...

# NOTE : The definition of identity used by this module ignores all
# gaps in the sequence.  Given than many of these alignments are
# gappy or fragmentary, including gaps in the identity score will
# dilute it somewhat according to coverage.


# Alternatively, you can check the coverage of each item of aligned 
# evidence.  This provides information about how well the aligned
# sequences are matched to the genomic sequence.

my $evidence_coverage = $alignment_tool->hit_coverage;
  
foreach my $hit (@$hit_coverage){

  print "Hit sequence             : " . $hit->{'name'} . "\n";
  print "Coverage                 : " . $hit->{'coverage'} . "\n";
  print "Number of presumed exons : " . $hit->{'exons'} . "\n";
  print "Unmatched 5prime bases   : " . $hit->{'unmatched_5prime'} . "\n";
  print "Unmatched 3prime bases   : " . $hit->{'unmatched_3prime'} . "\n";
  print "Unmatched internal bases : " . $hit->{'unmatched_internal'} . "\n\n";

}


# It is possible to determine the number of exons in an alignment 
# that have no evidence.

my $no_evidence_exons = $alignment_tool->rogue_exons;


=head1 DESCRIPTION
  
Object for dumping alignments of gene (transcript) predictions
and their associated supporting evidence.  Also produces key 
alignment statistics.  Pass in a transcript when the object 
is instantiated, use one of the methods to ask for the output 
you want, and there you are.

=head1 CONTACT
  
Post general queries to B<ensembl-dev@ebi.ac.uk>

=cut

package Bio::EnsEMBL::Pipeline::Alignment::EvidenceAlignment;

use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Pipeline::Alignment;
use Bio::EnsEMBL::Pipeline::Alignment::AlignmentSeq;

# Object preamble - inherits from Bio::Root::Object;

@ISA = qw(Bio::EnsEMBL::Root);


##### 'Public' methods #####


=head2 new

  Arg [1]    :
  Arg [2]    :
  Arg [3]    :
  Arg [4]    : (optional)
  Arg [5]    : (optional)
  Example    :
  Description: 
  Returntype : Bio::EnsEMBL::Pipeline::Alignment::EvidenceAlignment
  Exceptions : Will throw if isnt passed 
               DBAdaptor/Transcript/SeqFetcher
  Caller     : General

=cut


sub new {

  my ($class, @args) = @_;

  my $self = bless {},$class;
  
  my ($db, 
      $transcript, 
      $seqfetcher, 
      $padding, 
      $fasta_line_length, 
      $evidence_identity_cutoff) = $self->_rearrange([qw(DBADAPTOR
							 TRANSCRIPT
							 SEQFETCHER
							 PADDING
							 FASTA_LINE_LENGTH
							 IDENTITY_CUTOFF
							)],@args);
  

  # Throw an error if any of the below are undefined or
  # are the wrong thing.
  unless (defined $db && $db->isa("Bio::EnsEMBL::DBSQL::DBAdaptor")){
    $self->throw("No DB adaptor passed to AlignmentTool.  You passed a $db.");
  }
  unless (defined $transcript && $transcript->isa("Bio::EnsEMBL::Transcript")){
    $self->throw("No transcript passed to AlignmentTool.  You passed a $transcript.");
  }
  unless (defined $seqfetcher && $seqfetcher->isa("Bio::DB::RandomAccessI")) {
    $self->throw("No sequence fetcher passed to AlignmentTool.  You passed a $seqfetcher.");
  }

  $self->{'_db_adaptor'} = $db;

  # Set the amount of flanking sequence we want to include around
  # our slice.  This is quite important - many supporting features
  # extend past the beginnings and ends of predicted transcripts.  
  # Without padding, these sequences are truncated.  The default
  # of 50bp works OK, but you would want to set this manually
  # higher if you are interested in up- or down-stream goings-on.
  if ($padding) {
    $self->{'_transcript_padding'} = $padding;
  } else { 
    $self->{'_transcript_padding'} = 50;
  }

  # Store our SeqFetcher

  $self->_seq_fetcher($seqfetcher);

  # Create the slice we will work on
  
  $self->_slice($transcript);

  # Due to padding, it is necessary to re-construct our transcript
  # in proper slice coordinates.

  $self->_transcript($transcript, $self->_slice);


  # Determine if our database contains translations.  If it doesn't
  # we'll have to skip adding a set of translated exons to our 
  # alignment.

  if ($self->_transcript->translation){
    $self->_translatable(1);
  } else {
    $self->warn("Database doesn't contain translation.  Subsequently, ".
		"wont be able to display translations of each exon ".
		"or calculate protein identity scores.");
    $self->_type('nucleotide');
    $self->_translatable(0);
  }



  # The line length in the fasta alignment is set to a default
  # or a user specified value

  if ($fasta_line_length) {
    $self->_line_length($fasta_line_length);
  } else {
    $self->_line_length(60);
  }


  # Optional evidence identity cut-off

  if ($evidence_identity_cutoff) {
    $self->{'_evidence_identity_cutoff'} = $evidence_identity_cutoff;
  } else {
    $self->{'_evidence_identity_cutoff'} = 0;
  }

  return $self;
}


=head2 retrieve_alignment

  Arg [1]    :
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub retrieve_alignment {
  my $self = shift @_;

  my ($type, 
      $remove_introns, 
      $padding, 
      $show_missing_evidence,
      $merge_sequences) = $self->_rearrange([qw(TYPE
						REMOVE_INTRONS
						PADDING
						SHOW_UNALIGNED
						MERGE_SEQUENCES
					       )],@_);
  

  unless ($type) {
    $self->throw("Type of alignment to retrieve not specified.  Please use one " . 
		 "of \'all\', \'nucleotide\' or \'protein\'.");
  }

  $padding = 100 unless ($padding);

  unless ($self->_is_computed($type)){
    $self->_align($type, $show_missing_evidence, $remove_introns, 
		  $padding, $merge_sequences);
  }

  return $self->_create_Alignment_object($type, $show_missing_evidence);
}


=head2 identity

  Arg [1]    :
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub identity {

  my ($self) = @_;

  unless ($self->_is_computed){
    $self->_align('all');
  }

  return $self->_compute_identity;
}


=head2 hit_coverage

  Arg [1]    :
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub hit_coverage {
  my ($self) = @_;

  return $self->_compute_evidence_coverage;

}


=head2 rogue_exons

  Arg [1]    :
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub rogue_exons {
  my ($self) = @_;

  unless ($self->_is_computed){
    $self->_align('all');
  }

  unless ($self->_type eq 'all') {
    $self->warn("The alignment used to count rogue exons has\n".
		"not been created with both nucleotide and protein\n".
		"evidence.  Hence, it is quite likely that you\n".
		"will see rogue exons.");
  }

  my $evidence_alignments = $self->_working_alignment('evidence');

  my %seen_exons;

  foreach my $sequence (@$evidence_alignments){
    $seen_exons{$sequence->exon}++;
  }

  my $actual_exons = $self->_transcript->get_all_Exons;

  return ((scalar @$actual_exons) - (scalar keys %seen_exons))

}



##### Main Internal Methods #####


=head2 _align

  Arg [1]    :
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub _align {
  my ($self, $type, $show_missing_evidence, $truncate_introns, 
      $padding, $merge_sequences) = @_;

  # Collect all our sequences for aligning

  my $genomic_sequence         = $self->_genomic_sequence;
  my $exon_nucleotide_sequence = $self->_exon_nucleotide_sequence;
  
  my $exon_protein_sequence;
  if (($self->_translatable)
      &&(($type eq 'protein')||($type eq 'all'))){
    $exon_protein_sequence = $self->_exon_protein_translation;
  }
  
  my $evidence_sequences = $self->_corroborating_sequences($type);

  # Use information about 'deletions' (originally from the
  # cigar string and stored in the hashes returned by 
  # $self->_semialigned_sequences) to insert gaps into the genomic 
  # and exonic sequences.  This aligns the evidence with
  # the parent sequence (with a bit of fiddling to place
  # gaps into the evidence sequences that need them).

  for (my $i = 1; $i <= $self->_slice->length; $i++) {
    my $is_deletion = 0;

  DELETION_HUNT:
    foreach my $unaligned_sequence (@$evidence_sequences){

      if ($unaligned_sequence->fetch_deletion_at_position($i) eq 'D'){
	$is_deletion = 1;
	last DELETION_HUNT;
      }
    }
    
    if ($is_deletion) {

      $genomic_sequence->insert_gap($i, 1);
      $exon_nucleotide_sequence->insert_gap($i, 1);

      if ($self->_translatable
	  &&(($type eq 'protein')||($type eq 'all'))){
	$exon_protein_sequence->insert_gap($i, 1);  
      }


      for (my $j = 0; $j < scalar @$evidence_sequences; $j++) {
	unless ($evidence_sequences->[$j]->fetch_deletion_at_position($i) eq 'D') {

	  $evidence_sequences->[$j]->insert_gap($i, 1);

	}
      }
    }
    
  }

  # Put our working alignments somewhere handy

  $self->_working_alignment('genomic_sequence', $genomic_sequence);
  $self->_working_alignment('exon_nucleotide', $exon_nucleotide_sequence);

  if ($self->_translatable
      &&(($type eq 'protein')||($type eq 'all'))) {
    $self->_working_alignment('exon_protein', $exon_protein_sequence);
  }

  foreach my $evidence_sequence (@$evidence_sequences) {
    $self->_working_alignment('evidence', $evidence_sequence);
  }

  # If sequences are to be merged, do this now.
  if ($merge_sequences){
    $self->_merge_same_sequences;
  }

  # If unaligned fragments are to be shown, find these now.

  if ($show_missing_evidence) {
    $self->_compute_evidence_coverage;
    $self->_derive_unmatched_sequences;
  }

  # If introns are to be truncated, do this now.

  if ($truncate_introns) {
    $self->_truncate_introns($padding);
  }

  # Set flag to indicate that the alignment has been computed.

  $self->_is_computed($type, 1);
 
  return 1;
}


=head2 _create_Alignment_object

  Arg [1]    :
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub _create_Alignment_object {
  my ($self, $type, $show_missing_evidence) = @_;

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

  if ($show_missing_evidence && $self->_working_alignment('unaligned')) {
    foreach my $missing_sequence (@{$self->_working_alignment('unaligned')}){
      $alignment->add_sequence($missing_sequence);
    }
  }

  return $alignment;
}



=head2 _compute_identity

  Arg [1]    :
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub _compute_identity {

  my ($self) = @_;

  my $genomic_sequence = $self->_working_alignment('genomic_sequence');
  my $exon_protein_sequence = $self->_working_alignment('exon_protein');

  my $evidence = $self->_working_alignment('evidence');

  my @exon_identities;
  my %by_exon;

  foreach my $evidence_item (@$evidence) {
    push (@{$by_exon{$evidence_item->exon}}, $evidence_item);
  }

  my $exon_placemarker = 0;

  foreach my $exon (@{$self->_transcript->get_all_Exons}){

    my $highest_nucleotide_identity = 0;
    my $associated_nucleotide_coverage = 0;
    my $highest_protein_identity = 0;
    my $associated_protein_coverage = 0;

  EVIDENCE_ITEM:
    foreach my $evidence_item (@{$by_exon{$exon_placemarker}}){
      
      my $identity;
      my $coverage;

      # Here we are fetching the precent identity and coverage for
      # each evidence alignment.

      # We update the highest identity scores if the score just
      # calculated is higher AND has better than 80%
      # coverage OR better coverage than the present top identity 
      # match.

      # The top identities are grouped according to whether
      # they are protein or nucleotide sequences.

      if (($self->_translatable)
	  &&($evidence_item->type eq 'protein')
	  &&($self->_type ne 'nucleotide')){
	($identity, $coverage) = $self->_compare_to_reference($exon, 
							      $evidence_item, 
							      $exon_protein_sequence);

	if (($identity >= $highest_protein_identity)
	    &&(($coverage >= 80)
	       ||($coverage >= $associated_protein_coverage))) {
	  $highest_protein_identity = $identity;
	  $associated_protein_coverage = $coverage;
	}
      }

      elsif (($evidence_item->type eq 'nucleotide')
	     &&($self->_type ne 'protein')){
	($identity, $coverage) = $self->_compare_to_reference($exon, 
							      $evidence_item, 
							      $genomic_sequence);

      if (($identity >= $highest_nucleotide_identity)
	  &&(($coverage >= 80)
	     ||($coverage >= $associated_nucleotide_coverage))) {
	$highest_nucleotide_identity = $identity;
	$associated_nucleotide_coverage = $coverage;
      }

      } else {
	next EVIDENCE_ITEM;
      }

    }

    # Purely for neatness, some rounding
    $highest_nucleotide_identity    = sprintf("%.1f", $highest_nucleotide_identity);
    $associated_nucleotide_coverage = sprintf("%.1f", $associated_nucleotide_coverage);
    $highest_protein_identity       = sprintf("%.1f", $highest_protein_identity);
    $associated_protein_coverage    = sprintf("%.1f", $associated_protein_coverage);

    push (@exon_identities, [$highest_nucleotide_identity, 
			     $associated_nucleotide_coverage, 
			     $highest_protein_identity, 
			     $associated_protein_coverage]);

    $exon_placemarker++;
  }

  return \@exon_identities;

}

=head2 _compare_to_reference

  Arg [1]    :
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub _compare_to_reference { 
  my ($self, $exon, $evidence_align_seq, $reference_align_seq) = @_;

  # For nucleotide alignments each mismatch is counted
  # once.
  my $align_unit = 1;

  # If we are dealing this protein alignments we have to
  # multiply this by three.
  $align_unit *= 3 if ($evidence_align_seq->type eq 'protein');

  my $match_sequence = $evidence_align_seq->seq_array;
  my $reference_sequence = $reference_align_seq->seq_array;
  
  my $mismatches = 0;
  my $noncovered = 0;
  
  my $exon_start = $exon->start;
  my $exon_end = $exon->end;
  my $exon_length = $exon_end - $exon_start;

  if ($self->_strand == -1){
    $exon_start = $self->_slice->length - $exon_end + 1;
    $exon_end = $exon_start + $exon_length - 1;
  }

  for (my $i = $exon_start - 1; $i < $exon_end; $i++) {

    unless (defined ($match_sequence->[$i]) &&
	    defined ($reference_sequence->[$i]) &&
	    (($reference_sequence->[$i] eq $match_sequence->[$i])||
	     (($reference_sequence->[$i] eq '-')
	      ||($match_sequence->[$i] eq '-')))) {

      $mismatches += $align_unit;
    }
    
    if (($reference_sequence->[$i] ne '-')
	&&($match_sequence->[$i] eq '-')) {

      $noncovered += $align_unit;
    }
  }
  
  my $identity = (1 - ($mismatches/$exon_length))*100;
  
  # The next line gets around the problem of exon length not always
  # being a whole number of codons.  There can be cases where
  # there are more non-covered bases than there are bases in an exon.
  $noncovered = $exon_length if $noncovered > $exon_length;
  
  my $coverage = (1 - ($noncovered/$exon_length))*100;
  
  return ($identity, $coverage);
}


=head2 _compute_evidence_coverage

  Arg [1]    :
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub _compute_evidence_coverage {
  my ($self) = @_;

  my %coordinates_hash;
  my @evidence_coverage_stats;

  foreach my $exon (@{$self->_transcript->get_all_Exons}){

    foreach my $supporting_feature (@{$exon->get_all_supporting_features}){
      push (@{$coordinates_hash{$supporting_feature->hseqname}}, 
	    [$supporting_feature->hstart, $supporting_feature->hend, 
	     $supporting_feature->start, $supporting_feature->end]);
    }

  }  

 SEQUENCE:
  foreach my $sequence_identifier (keys %coordinates_hash) {

    my $length = $self->_fetch_sequence($sequence_identifier)->length;

    unless ($length){
      $self->warn("Evidence sequence with zero length found." 
		  . "  This sequence probably couldnt be retrieved.");
      next SEQUENCE;
    }

    my $presumed_exons = scalar @{$coordinates_hash{$sequence_identifier}};
	
    my $covered_length = 0;
    my $previous_end;
    my @sorted_matches = sort {$a->[0] <=> $b->[0]} @{$coordinates_hash{$sequence_identifier}};

    foreach my $coordinate_pair (@sorted_matches) {

                          # Hit end        minus hit start       plus one;
      $covered_length += $coordinate_pair->[1] - $coordinate_pair->[0] + 1;

      if ($previous_end && $coordinate_pair->[0] != $previous_end + 1) {
	$self->_add_unmatched_region($sequence_identifier, $coordinate_pair->[0], 
				     $coordinate_pair->[1], 'before',$coordinate_pair->[1]);
      }
    }

    if ($sorted_matches[0]->[0] != 1){ 
      $self->_add_unmatched_region($sequence_identifier, 1, ($sorted_matches[0]->[0] - 1), 
				   'before', $sorted_matches[0]->[2]);
    }

    if ($sorted_matches[-1]->[1] != $length){
      $self->_add_unmatched_region($sequence_identifier, ($sorted_matches[-1]->[1] + 1), 
				   $length, 'after', $sorted_matches[-1]->[3]);
    }

    my $uncovered_5prime_bases = $sorted_matches[0]->[0] - 1;
    my $uncovered_3prime_bases = $length - $sorted_matches[-1]->[1];
    my $uncovered_internal_bases = $length - $covered_length - 
      $uncovered_5prime_bases - $uncovered_3prime_bases;

    my %these_stats = ('name'               => $sequence_identifier,
		       'coverage'           => (($covered_length / $length) * 100),
		       'exons'              => $presumed_exons,
		       'unmatched_5prime'   => $uncovered_5prime_bases,
		       'unmatched_3prime'   => $uncovered_3prime_bases,
		       'unmatched_internal' => $uncovered_internal_bases
		      );

    push (@evidence_coverage_stats, \%these_stats);

  }

  return \@evidence_coverage_stats;

}

=head2 _add_unmatched_region

  Arg [1]    :
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub _add_unmatched_region {
  my ($self, $seqname, $start, $end, $before_or_after, $genomic_coord) = @_;

  if ($seqname && $start && $end && 
      ($before_or_after eq 'before' || $before_or_after eq 'after') && $genomic_coord) {
    
    push (@{$self->{'_unmatched_evidence_sequence'}->{$seqname}}, [$start, $end, 
								   $before_or_after, 
								   $genomic_coord]);

  } else {
    $self->throw("Incorrect arguments specified.");
  }

}

=head2 _derive_unmatched_sequences

  Arg [1]    :
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub _derive_unmatched_sequences {
  my ($self) = @_;

  my $spacing = 5;
  my @sequences;

  $self->warn("No unmatched sequences have been found")
    unless $self->{'_unmatched_evidence_sequence'};

  foreach my $seqname (keys %{$self->{'_unmatched_evidence_sequence'}}){

    my @sorted_missing_bits = sort {$a->[0] <=> $b->[0]} @{$self->{'_unmatched_evidence_sequence'}->{$seqname}};

    my $fetched_seq = $self->_fetch_sequence($seqname);
    my @fetched_seq_array = split //, $fetched_seq->seq;

    my $slice_seq = '-' x $self->_slice->length;
    my @slice_seq_array = split //, $slice_seq;

    foreach my $missed_fragment (@sorted_missing_bits){

      my $insert_start;

      if ($missed_fragment->[2] eq 'before'){

	$insert_start = $missed_fragment->[3] - 
	  ($missed_fragment->[1] - $missed_fragment->[0]) - $spacing;

      } elsif ($missed_fragment->[2] eq 'after') {

	$insert_start = $missed_fragment->[3] + $spacing;

      } else {
	$self->throw("Before or after not specified. Internal error.");
      }

      my $insert_point = $insert_start;

      for (my $i = $missed_fragment->[0]; $i <= $missed_fragment->[1]; $i++) {
	$slice_seq_array[$insert_point] = $fetched_seq_array[$i];
	$insert_point++;
      }      
    }

    my $sequence = '';

    foreach my $element (@slice_seq_array) {
      $sequence .= $element unless !$element;
    }

    my $display_seqname = "Unmatched " . $seqname;

    my $align_seq = Bio::EnsEMBL::Pipeline::Alignment::AlignmentSeq->new(
                              '-name'      => $display_seqname,
			      '-seq'       => $sequence,
			      '-deletions' => 0,
			      '-type'      => 'nucleotide');


    $self->_working_alignment('unaligned', $align_seq);

  }

  return 1;
}


=head2 _truncate_introns

  Arg [1]    :
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub _truncate_introns {
  my ($self, $padding) = @_;

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
    push(@coordinates, $exon->start);
    push(@coordinates, $exon->end);
  }

  @coordinates = sort {$b <=> $a} @coordinates;

  shift @coordinates;
  pop @coordinates;


  # Splice introns from each sequence.

  foreach my $align_seq (@sequences){

    my $seq_array = $align_seq->seq_array;

    my @working_coordinates = @coordinates;

  INTRON:
    while (@working_coordinates){
      my $intron_end = (shift @working_coordinates) - 1;
      my $intron_start = (shift @working_coordinates) + 1;

      next INTRON unless ($intron_start + $padding + 22) < ($intron_end - $padding - 22);

      my @offcut = splice (@$seq_array, $intron_start + $padding, 
			   ($intron_end - 2*$padding - $intron_start), 
			   '---intron-truncated---');

      my $offcut;
      foreach my $discarded_element (@offcut) {
	$offcut .= $discarded_element;
      }

      if (($offcut =~ /atgc/i)&&($align_seq->name ne 'genomic_sequence')) {
	$self->warn("Truncating intron sequences has caused some aligned evidence\n" .
		    "to be discarded.  Try again with a higher amount of padding\n".
		    "around the exon sequences.\n");
      }
      
    }
    
    
    $align_seq->store_seq_array($seq_array);
  }
  
  return 1;
}

=head2 _merge_same_sequences

  Arg [1]    :
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub _merge_same_sequences {
  my ($self) = @_;

  my @sequences = @{$self->_working_alignment('evidence')};

  # Blam!

  $self->_working_alignment('evidence', 'empty');

  my %same_sequences_hash;

  foreach my $sequence (@sequences) {
    
    push (@{$same_sequences_hash{$sequence->name}}, $sequence);
  }

  foreach my $sequence_name (keys %same_sequences_hash){

    my $chimera = '-' x $self->_slice->length;

    my @chimera = split (//, $chimera);

    foreach my $alike_sequence (@{$same_sequences_hash{$sequence_name}}){
      
      my $seq_array = $alike_sequence->seq_array;

      for (my $i = 0; $i <= scalar @$seq_array; $i++) {

	if ($seq_array->[$i] ne '-'){
	  if ($chimera[$i] && $chimera[$i] ne '-'){
	    $self->warn("Very bad - evidence from the same sequence overlaps between" .
			"exons.  You should not be merging sequences.  Set : " .
			"\$evidence_alignment->retrieve_alignment('-merge_sequences' => 0)");
	  } else {
	    $chimera[$i] = $seq_array->[$i] unless $seq_array->[$i] eq '-';
	  }
	}
      }
    }

    my $merged_align_seq = $same_sequences_hash{$sequence_name}->[0];

    $merged_align_seq->store_seq_array(\@chimera);

    $self->_working_alignment('evidence', $merged_align_seq);

  }

  return 1
}



=head2 _is_computed

  Arg [1]    :
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

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

  if (!defined $self->{'_is_computed'}) {
    $self->{'_is_computed'} = 0;
  }
  
  # Check whether an alignment of a specific type
  # has been run.
  if ((!defined $value)&&($self->{'_is_computed'})&&($type ne $self->_type)) { 
    $self->warn("Alignment has been previously computed, but was not\n" .
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
      warn "Unknown alignment type.  Can be nucleotide, protein or all.\n";
      return 0;
    }

    $self->{'_is_computed'} = 1;
  }

  return $self->{'_is_computed'};
}


=head2 _type

  Arg [1]    :
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub _type {
  my ($self, $type) = @_;

  if (defined $type) {
    $self->{'_computed_type'} = $type;
  }

  return $self->{'_computed_type'};
}



##### Alignment information handling methods #####


=head2 _working_alignment

  Arg [1]    :
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub _working_alignment {

  my ($self, $slot, $align_seq) = @_;

  unless (defined $slot && 
      (($slot eq 'genomic_sequence')
       ||($slot eq 'exon_protein')
       ||($slot eq 'exon_nucleotide')
       ||($slot eq 'evidence')
       ||($slot eq 'unaligned'))){
    $self->throw("Was trying to retrieve or write working alignments to "
		 . "a slot that isnt allowed ($slot)");
  }

  if (defined $slot && defined $align_seq && $align_seq eq 'empty'){
    undef $self->{'_working_alignment_array'}->{$slot};
    return 1
  }

  if (defined $slot && defined $align_seq){

    unless ($align_seq->isa("Bio::EnsEMBL::Pipeline::Alignment::AlignmentSeq")){
      $self->throw("Sequence passed to _working alignment was not an " . 
		   "AlignmentSeq, it was a [$align_seq]")
    }

    push (@{$self->{'_working_alignment_array'}->{$slot}}, $align_seq);

  } elsif (defined $slot && !defined $align_seq) {

    my $slot_resident =  $self->{'_working_alignment_array'}->{$slot};

    if ((($slot eq 'genomic_sequence')||($slot eq 'exon_protein')||($slot eq 'exon_nucleotide')) 
	&& defined $slot_resident && scalar @$slot_resident == 1) {
      return $slot_resident->[0];
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

# This method could be irradicated by adding some options
# to our object such that they can be instantiated with
# a slice and a bunch of transcripts from that slice.
# Something to think about.

sub _transcript {
  my ($self, $input_transcript, $slice) = @_;

  if (defined $input_transcript && defined $slice){

  # This is just a little horrible, but it
  # does make things much neater elsewhere.
  # We need to retrieve our transcript from
  # our slice.  We iterate through all genes on
  # the slice, then search the stable ids of the 
  # transcripts on this gene.  Yuck, but for 
  # most genes this is not much fiddling about.

  GENE:
    foreach my $candidate_gene (@{$slice->get_all_Genes}){

    TRANSCRIPT:
      foreach my $candidate_transcript (@{$candidate_gene->get_all_Transcripts}) {

	# Check whether the transcripts are the same, first using the stable_id
	# if it exists, otherwise with the dbID.
	unless (($candidate_transcript->stable_id 
		 && $input_transcript->stable_id 
		 && $candidate_transcript->stable_id eq $input_transcript->stable_id)
		||($candidate_transcript->dbID 
		   && $candidate_transcript->dbID == $input_transcript->dbID)){
	  next TRANSCRIPT;
	}
	
	$self->{'_transcript'} = $candidate_transcript;

	# Get the strand of our gene
	$self->_strand($candidate_gene->strand);

	last GENE;
      }
    }
    
    unless ($self->{'_transcript'}){
      $self->throw("Could not find transcript on Slice.  Very bad.");
    }
  }

  if (defined $self->{'_transcript'}){
    return $self->{'_transcript'}
  } else {
    $self->throw("Something has gone wrong.  A Transcript has not yet been".
		 " extracted from our new Slice.");
  }
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

  my ($self, $transcript) = @_;

  if (! defined $self->{'_slice'} && defined $transcript) {

    my $slice_adaptor = $self->{'_db_adaptor'}->get_SliceAdaptor;

    # Ideally, stable_ids should be used to fetch things - just in 
    # case the transcript comes from a different db to our current
    # db.  Otherwise, fall over to using dbID.

    if ($transcript->stable_id){

      $self->{'_slice'} = 
	$slice_adaptor->fetch_by_transcript_stable_id($transcript->stable_id, 
						      $self->{'_transcript_padding'});
    } else {
      
      $self->{'_slice'} = 
	$slice_adaptor->fetch_by_transcript_id($transcript->dbID, 
						      $self->{'_transcript_padding'});

    }
  }
  
  return $self->{'_slice'};
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
  my ($self, $value) = @_;

  # Potential gotcha.  The strand is most likely determined in
  # the _transcript method, as it is easiest to derive it then.
  # This method doesn't actually figure anything out.

  if (defined $value) {
    $self->{'_strand'} = $value;
  }

  if (! defined $self->{'_strand'}){
    $self->warn("No value for strand set.  Disaster awaits.");
  }

  return $self->{'_strand'};
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

    $self->{'_seq_fetcher'} = $fetcher;

    return 1;
  }

  return $self->{'_seq_fetcher'};
}


##### Sequence handling methods #####


=head2 _corroborating_sequences

  Arg [1]    :
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub _corroborating_sequences {
  my ($self, $type) = @_;

  # Possibly, it might be good to install some way of caching
  # these results, if ever there was a possibility that this
  # method could be called more than once per instantiation.

  # For each exon attached to our transcript, work through
  # each attached DnaXXXAlignFeature and create a sequence
  # for it.
  
  my $exons = $self->_transcript->get_all_Exons;

  # Put protein features separately, purely to keep the
  # nucleotide and protein parts of our alignment as
  # distinct blocks later.
  my @protein_features;

  my $exon_placemarker = 0;

  foreach my $exon (@{$exons}){
    # Work through each item of supporting evidence attached to our exon.

  FEATURE:
    foreach my $base_align_feature (@{$exon->get_all_supporting_features}){
 
      if ((($type eq 'nucleotide')
	   &&($base_align_feature->isa("Bio::EnsEMBL::DnaPepAlignFeature")))
	  ||(($type eq 'protein')
	     &&($base_align_feature->isa("Bio::EnsEMBL::DnaDnaAlignFeature")))){
	next FEATURE;
      }

      if ((defined $base_align_feature->percent_id)
	  &&($base_align_feature->percent_id < $self->{'_evidence_identity_cutoff'})) {
	  next FEATURE;
	}

      if ($base_align_feature->isa("Bio::EnsEMBL::DnaDnaAlignFeature")){

	my $align_seq = $self->_fiddly_bits($base_align_feature);
	next FEATURE unless $align_seq;

	$align_seq->exon($exon_placemarker);
	$align_seq->type('nucleotide');

	push (@{$self->{'_corroborating_sequences'}}, $align_seq);
	next FEATURE;
      }
      
      if ($base_align_feature->isa("Bio::EnsEMBL::DnaPepAlignFeature")){

	my $align_seq = $self->_fiddly_bits($base_align_feature);
	next FEATURE unless $align_seq;

	$align_seq->exon($exon_placemarker);
	$align_seq->type('protein');
	push (@protein_features, $align_seq);
      }
    }

    $exon_placemarker++;
  }
  # The order in which we add things to our array effects the order
  # in the final alignment.  Here, if they exist, protein features 
  # are added separately so that they are all together in the final 
  # alignment.  There could be a better way to do this, of course.
  push (@{$self->{'_corroborating_sequences'}}, @protein_features);

  return $self->{'_corroborating_sequences'};
}


=head2 _fiddly_bits

  Arg [1]    :
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub _fiddly_bits {
  my ($self, $base_align_feature) = @_;
  
  my $hstart = $base_align_feature->hstart;
  my $hend = $base_align_feature->hend;


  # Fetch our sequence from the cache.  If the sequence
  # is missing it means that it could not be fetched and
  # we will have to ignore it.

  my $fetched_seq = $self->_fetch_sequence($base_align_feature->hseqname);

  if ( ! $fetched_seq) {

    $self->warn("Error fetching sequence " 
		. $base_align_feature->hseqname . ".  Ignoring.");
    
    return 0;
  }

  my @fetched_seq;

  # If this is a protein align feature, pad amino acids with gaps
  # to make them comparable to nucleotide coords.  Then, make sure
  # splice out our sequence of interest using hit coordinates that
  # take account of the padding of our sequence.

  if ($base_align_feature->isa("Bio::EnsEMBL::DnaPepAlignFeature")){
    my $padded_aa_seq;
    ($padded_aa_seq = $fetched_seq->seq) =~ s/(.)/$1\-\-/g;
    
    my @full_seq = split //, $padded_aa_seq;
    
  
    # Splice out the matched region of our feature sequence
    my $first_aa = ($hstart - 1) * 3;
    my $last_aa = ($hend * 3) - 1;

    my $length = $last_aa - $first_aa + 1;

    @fetched_seq = splice(@full_seq, $first_aa, $length);
  } 

  # If we have a dna align feature, extracting the correct portion
  # of the hit sequence is a bit easier than the method required
  # for a protein sequence.

  if ($base_align_feature->isa("Bio::EnsEMBL::DnaDnaAlignFeature")) {

    @fetched_seq = split //, $fetched_seq->seq;
	
    # Splice out the matched region of our feature sequence
    @fetched_seq = splice(@fetched_seq, ($hstart - 1), ($hend -$hstart + 1));
    
  }

  # Add the needed insertion gaps to our supporting
  # feature sequence.  We only add gaps to the supporting
  # sequences at this stage.  Once we have all our supporting
  # sequences constructed and determined where all the 
  # 'deletion' gaps lie in these, we can transpose this 
  # information onto our genomic and exonic sequences - this is
  # done in $self->_align.

  my @cigar_instructions = $self->_cigar_reader($base_align_feature->cigar_string);
  
  my $added_gaps = 0;
  my $hit_sequence_position = $base_align_feature->start;
  my @deletion_sequence;

  foreach my $instruction (@cigar_instructions) {

    if ($instruction->{'type'} eq 'I') {
      my $gap = '-' x $instruction->{'length'};
      my @gap = split //, $gap;
      
      splice(@fetched_seq, $hit_sequence_position, 0, @gap);
      
      $hit_sequence_position += $instruction->{'length'};
      
    } elsif ($instruction->{'type'} eq 'M') {
      
      $hit_sequence_position += $instruction->{'length'};
      
    } elsif ($instruction->{'type'} eq 'D') {

      for (my $i = $hit_sequence_position; $i < ($hit_sequence_position + $instruction->{'length'});$i++){
	$deletion_sequence[$i] = 'D';
      }
      
    }
    
  }

  # Determine the point in the genomic sequence that the
  # features starts at.  Need to worry about strand and whether 
  # the feature starts or ends outside of our slice.

  my $genomic_start;  # Feature insertion point
  if ($self->_strand == 1) {
    $genomic_start = $base_align_feature->start - 1;
  } elsif ($self->_strand == -1) {
    $genomic_start = $self->_slice->length - $base_align_feature->end;
  }

  
  # This little section of code handles any sequence that
  # overshoots the beginning of our slice.  Chop.

  if ($genomic_start < 0) {
    $self->warn("Feature extends past the ends of genomic slice.  Truncating it to fit.");

    my $overshoot = $genomic_start * -1;      
        
    $genomic_start = 0;

    splice (@fetched_seq, 0, $overshoot);
  }
  
  # Here we are actually building the sequence that will
  # align to our slice

  my $feature_sequence = '-' x $self->_slice->length;
  my @feature_sequence = split //, $feature_sequence;

  splice (@feature_sequence, $genomic_start, (scalar @fetched_seq), @fetched_seq);

  $feature_sequence = '';
  
  foreach my $element (@feature_sequence) {
    $feature_sequence .= $element;
  }

  # Munch our array of deletion information into a string
  # This is not pretty, but will work for now.

  my $deletion_sequence = '';

  foreach my $element (@deletion_sequence){
 
    unless ($element){
      $deletion_sequence .= '-';
      next;
    }

    $deletion_sequence .= $element;
  }

  # Create a AlignmentSeq object with our valuable sequence:


  my $partially_aligned = Bio::EnsEMBL::Pipeline::Alignment::AlignmentSeq->new(
                                          '-name'      => $base_align_feature->hseqname,
					  '-seq'       => $feature_sequence,
					  '-deletions' => $deletion_sequence,
					  '-start'      => $base_align_feature->start);

  return $partially_aligned;
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

  if (!defined $self->{'_genomic_sequence'}) {

    my $genomic_sequence;

    if ($self->_strand == 1) {
      $genomic_sequence = $self->_slice->seq;
    } elsif ($self->_strand == -1) {
print STDERR "Reverse complimenting genomic sequence.\n";
      $genomic_sequence = $self->_slice->revcom->seq;
    }

    $self->{'_genomic_sequence'} = Bio::EnsEMBL::Pipeline::Alignment::AlignmentSeq->new(
					     '-seq'  => $genomic_sequence,
					     '-name' => 'genomic_sequence',
					     '-type' => 'nucleotide'
                                             );

  }

  return $self->{'_genomic_sequence'};
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

  if (!defined $self->{'_exon_nucleotide_sequence'}) {
    
    my @exon_only_sequence;

    foreach my $exon (@{$self->_transcript->get_all_Exons}) {

      # Add the exon sequence to our 'exons only' sequence.
      
      my @exon_seq = split //, $exon->seq->seq;

      my $exon_position;

      if ($self->_strand == 1) {
	$exon_position = $exon->start - 1;
      }elsif ($self->_strand == -1) {
	$exon_position = $self->_slice->length - $exon->end;
      }
      
      foreach my $exon_nucleotide (@exon_seq){
	if (!defined $exon_only_sequence[$exon_position]){
	  $exon_only_sequence[$exon_position] = $exon_nucleotide;
	  $exon_position++;
	} else {
	  die "Overlapping exons \!\?\!\n$@";
	}
      }
    }    
    
    # Fill in the blanks

    for (my $i = 0; $i < $self->_slice->length; $i++) {
      unless (defined $exon_only_sequence[$i]){
	$exon_only_sequence[$i] = '-';
      } 
    }

    # Convert back to a string

    my $exon_sequence = '';
    foreach my $element (@exon_only_sequence) {
      $exon_sequence .= $element;
    }


    $self->{'_exon_nucleotide_sequence'} = Bio::EnsEMBL::Pipeline::Alignment::AlignmentSeq->new(
					     '-seq'  => $exon_sequence,
					     '-name' => 'exon_sequence',
					     '-type' => 'nucleotide'
                                             );
  }

  return $self->{'_exon_nucleotide_sequence'};
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

  if (! defined $self->{'_exon_protein_translation'}) {

    my @exon_translation_sequence;

    my $exons = $self->_transcript->get_all_Exons;
    
    foreach my $exon (@{$exons}){
      # Add a translation of this exon peptide to our translated exon sequence.

      my $peptide_obj = $exon->peptide($self->_transcript);
      my $peptide = $peptide_obj->seq;

      $peptide =~ s/(.)/$1\-\-/g;

      my @peptide = split //, $peptide;

      # Whack off the first residue if it is only a partial 
      # codon (the internal rule is to:
      #   - include a whole residue for partial codons at ends
      #       of exons
      #   - ignore/remove residues from partial codons at 
      #       starts of exons)
      # By doing this, a complete undeleted/unrepeated sequence 
      # is displayed in the alignment.
      
      if ($exon->phase == 2){ 
	splice (@peptide, 0, 2);
      }

      if ($exon->phase == 1){
	shift @peptide;
      }

      # Darstardly, hidden in here is the coordinate shuffle
      # to turn reverse strand genes around (the protein sequence
      # is of course in the forward direction already, just the 
      # coordinates need to be reversed).
      
      my $exon_start = $exon->start;
      my $exon_end = $exon->end;
      my $exon_length = $exon_end - $exon_start;
      my $exon_phase = $exon->phase;
      my $exon_end_phase = $exon->end_phase;
      
      if ($self->_strand == -1) {
	
	$exon_start = $self->_slice->length - $exon_end;
	$exon_end = $exon_start + $exon_length - 1;

#	($exon_phase, $exon_end_phase) = ($exon_end_phase, $exon_phase);

      }

      # Jiggling the exons about to get frame right
      my $extra_length = 0;

#      $extra_length += 3 if (($exon->phase != 0) && ($exon->phase != -1));
      $extra_length += 3 if (($exon_end_phase != 0) && ($exon_end_phase != -1));

#      $extra_length -= 2 if $exon->phase == 2;
#      $extra_length -= 1 if $exon->phase == 1;

      $extra_length -= 2 if $exon_end_phase == 2;
      $extra_length -= 1 if $exon_end_phase == 1;


      my $peptide_genomic_start;

      if ($exon_end_phase != -1) {
	$peptide_genomic_start = $exon_end - (scalar @peptide) + $extra_length + 1 - 1;

      } else {
	$peptide_genomic_start = $exon_start - 1;
      }

      # This HAS to be removed, it is just necessary to make this work.
      # There is something about reversed sequences that puts them two
      # base positions upstream.
      if ($self->_strand == -1) {
	$peptide_genomic_start += 2;
      }
      
      my $insert_point = $peptide_genomic_start;
      
      foreach my $exon_aa (@peptide) {
	$exon_translation_sequence[$insert_point] = $exon_aa;

	$insert_point++;
      }
    }

    # Fill in the blanks

    for (my $i = 0; $i < $self->_slice->length; $i++) {
      unless (defined $exon_translation_sequence[$i]){
	$exon_translation_sequence[$i] = '-';
      } 
    }

    # Convert back to a string

    my $translated_exon_sequence = '';
    foreach my $element (@exon_translation_sequence) {
      $translated_exon_sequence .= $element;
    }

    $self->{'_exon_protein_translation'} = Bio::EnsEMBL::Pipeline::Alignment::AlignmentSeq->new(
					     '-seq'  => $translated_exon_sequence,
					     '-name' => 'translated_exon_sequence',
					     '-type' => 'protein'
                                             );
  }
  
  return $self->{'_exon_protein_translation'};
}

##### Methods that take care of sequence fetching and caching #####

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

  foreach my $exon (@{$self->_transcript->get_all_Exons}){

    foreach my $supporting_feature (@{$exon->get_all_supporting_features}){
      $hash_of_accessions{$supporting_feature->hseqname}++;
    }

  }  

  my @array_of_accessions = keys %hash_of_accessions;

  # Retrieve sequences.

  my $fetched_seqs;

  if ($self->_seq_fetcher->can("batch_fetch")){

    eval {
      $fetched_seqs = $self->_seq_fetcher->batch_fetch(@array_of_accessions);
    };
    
    if ($@){
      $self->warn("Not all evidence sequences could be pfetched.\n".
		  "Ignoring missing sequences.\n$@\n");
    }

  } else {

    foreach my $accession (@array_of_accessions){

      my $fetched_seq;

      eval {
	$fetched_seq = $self->_seq_fetcher->get_Seq_by_acc($accession);
      };

      if ($@) {
	$self->warn("The seqfetcher is having trouble.\t$@\n");
      }

      push(@$fetched_seqs, $fetched_seq);

    }
  }
  
  # Build cache.

  foreach my $fetched_seq (@$fetched_seqs){

    $self->{'_fetched_seq_cache'}->{$fetched_seq->accession_number} = $fetched_seq;

  }

  $self->{'_cache_is_built'} = 1;
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
    unless $self->{'_cache_is_built'};

  unless ($self->{'_fetched_seq_cache'}->{$accession}){
    $self->warn("Sequence $accession could not be retrieved from cache.");
  }

  return $self->{'_fetched_seq_cache'}->{$accession}; 
}


### Miscellaneous utilities ###


=head2 _line_length

  Arg [1]    :
  Example    : 
  Description: Getter/Setter for the line length in fasta output.
  Returntype : 
  Exceptions : 
  Caller     : 

=cut


sub _line_length {
  my $self = shift;

  if (@_) {
    $self->{'_fasta_line_length'} = shift;
  }

  return $self->{'_fasta_line_length'};
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
    $self->{'_translatable'} = shift;
  }

  return $self->{'_translatable'};
}




##### CIGAR string handlers #####


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

  while (my $next_char = shift @cigar_array) {

    if ($next_char =~ /[MDI]/) {

      my %enduring_hash;
      $enduring_hash{'type'} = $next_char;
      $enduring_hash{'length'} = $current_digits;

      push (@cigar_elements, \%enduring_hash);

      $current_digits = '';

    } elsif ($next_char =~ /\d/) {

      $current_digits .= $next_char;

    } else {

      die "There is an odd character in the CIGAR string retrieved from the database.\n" . 
	$cigar_string . "\n";

    }

  }
 
  return @cigar_elements;

}


return 1;
