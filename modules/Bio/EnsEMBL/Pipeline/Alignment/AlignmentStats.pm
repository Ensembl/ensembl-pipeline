
# Cared for by Dan Andrews <dta@sanger.ac.uk>
#
# Copyright EnsEMBL
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod

=head1 NAME

Bio::EnsEMBL::Pipeline::Alignment::AlignmentStats

=head1 SYNOPSIS

Give this object an EvidenceAlignment object, or pass it
a transcript and it will calculate all sorts of useful
information, such as:

* Identity score of best evidence for each exon in the
the transcript.
* The identity between the genomic sequence and each hit,
broken down by exon.
* Locations of parts of evidence sequences that do not align
to the genomic sequence.
* Exons without evidence.

And so on, blah, blah.

Importantly, if you need to see non-aligned portions of 
the evidence sequences use the following '-show_unaligned'
option:

my $alignment = 
  $align_tool->retrieve_alignment(
      '-type'           => 'all',
      '-show_unaligned' => 1);

This turns on an portion of code that generates sequences
composed of all the fragments in the evidence sequence that
havent been matched.  This is mainly useful to genebuilders 
worrying about things that might be missed by various similarity 
matching algorithms.  Note that for layout the fragments are 
placed in the intron where they conceivably should go.  If 
you are also trimming the intron sequences be careful that 
you arent unwittingly throwing these intronic fragments
away.  A warning is raised if this happens.

Once the alignment is generated it is possible to calculate
the identity of the best matching evidence for each exon.
Coverage is also determined.

my $exon_identities = $evidence_alignment->identity;

my $exon_counter = 1;
foreach my $exon_identity (@$exon_identities) {
  print "Exon $exon_counter\t" . 
    $exon_identity->[0] . "\t" .
      $exon_identity->[1] . "\t" .
	$exon_identity->[2] . "\t" .
	  $exon_identity->[3] ."\n";
  $exon_counter++;
}

The identity scores are returned as a reference to an  array 
of array references (I know, bleah - need a mini-object 
perhaps).  Each referenced array represents the best identity 
and coverage for an individual exon. The reference array has
the format (nucleotide_identity, nucleotide_coverage,
protein_identity, protein_coverage).  This could be easier...

NOTE : The definition of identity used by this module ignores all
gaps in the sequence.  Given than many of these alignments are
gappy or fragmentary, including gaps in the identity score will
dilute it somewhat according to coverage.


Alternatively, you can check the coverage of each item of aligned 
evidence.  This provides information about how well the aligned
sequences are matched to the genomic sequence.

my $evidence_coverage = $evidence_alignment->hit_coverage;

foreach my $hit (@$hit_coverage){

  print "Hit sequence             : " . $hit->{'name'} . "\n";
  print "Coverage                 : " . $hit->{'coverage'} . "\n";
  print "Number of presumed exons : " . $hit->{'exons'} . "\n";
  print "Unmatched 5prime bases   : " . $hit->{'unmatched_5prime'} . "\n";
  print "Unmatched 3prime bases   : " . $hit->{'unmatched_3prime'} . "\n";
  print "Unmatched internal bases : " . $hit->{'unmatched_internal'} . "\n\n";

}

It is possible to determine the number of exons in an alignment
that have no evidence.

my $no_evidence_exons = $evidence_alignment->rogue_exons;

print "Exons without evidence : " . $no_evidence_exons . "\n";


=head1 DESCRIPTION

Calculates a few useful numbers describing the quality of the evidence 
mapped to a given transcript.

=head1 CONTACT
  
Post general queries to B<ensembl-dev@ebi.ac.uk>

=cut

package Bio::EnsEMBL::Pipeline::Alignment::AlignmentStats;

use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Pipeline::Alignment;
use Bio::EnsEMBL::Pipeline::Alignment::AlignmentSeq;
use Bio::EnsEMBL::Utils::Exception qw(throw warning info);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

@ISA = qw();


##### 'Public' methods #####

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
    warning("The alignment used to count rogue exons has\n".
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

      # Here we are fetching the percent identity and coverage for
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

  # If we are dealing with protein alignments we have to
  # multiply this by three.
  $align_unit *= 3 if ($evidence_align_seq->type eq 'protein');

  my $match_sequence = $evidence_align_seq->seq_array;
  my $reference_sequence = $reference_align_seq->seq_array;

  my $mismatches = 0;
  my $noncovered = 0;

  my $exon_start = $exon->start;
  my $exon_end = $exon->end;
  my $exon_length = $exon_end - $exon_start;

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

  foreach my $supporting_feature (@{$self->_all_supporting_features}){
    push (@{$coordinates_hash{$supporting_feature->hseqname}}, 
	  [$supporting_feature->hstart, $supporting_feature->hend, 
	   $supporting_feature->start, $supporting_feature->end]);
  }

 SEQUENCE:
  foreach my $sequence_identifier (keys %coordinates_hash) {

    my $length = $self->_fetch_sequence($sequence_identifier)->length;

    unless ($length){
      info("Evidence sequence with zero length found." 
	      . "  This sequence probably couldnt be retrieved.");
      next SEQUENCE;
    }

    my $presumed_exons = scalar @{$coordinates_hash{$sequence_identifier}};
	
    my $covered_length = 0;
    my $previous_end;
    my @sorted_matches = sort {$a->[0] <=> $b->[0]} @{$coordinates_hash{$sequence_identifier}};

    my @genomic_starts = sort {$a->[2] <=> $b->[2]} @{$coordinates_hash{$sequence_identifier}};
    my @genomic_ends   = sort {$a->[3] <=> $b->[3]} @{$coordinates_hash{$sequence_identifier}};
    my $genomic_start  = $genomic_starts[0]->[2]; 
    my $genomic_end    = $genomic_ends[-1]->[3];

    foreach my $coordinate_pair (@sorted_matches) {
                          # Hit end                minus hit start    plus one;
      $covered_length += $coordinate_pair->[1] - $coordinate_pair->[0] + 1;

      if ($previous_end && $coordinate_pair->[0] != $previous_end + 1) {
	$self->_add_unmatched_region($sequence_identifier, 
				     $previous_end, 
				     $coordinate_pair->[0], 
				     'before', 
				     $coordinate_pair->[2]);
      }

      $previous_end = $coordinate_pair->[1];
    }

    if ($sorted_matches[0]->[0] != 1){ 
       $self->_add_unmatched_region($sequence_identifier, 1, ($sorted_matches[0]->[0] - 1), 
				    'before', $genomic_start);
    }

    if ($sorted_matches[-1]->[1] != $length){
	$self->_add_unmatched_region($sequence_identifier, ($sorted_matches[-1]->[1] + 1), 
				     $length, 'after', $genomic_end);
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
    throw("Incorrect arguments specified.");
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

  my $spacing = 0;
  my @sequences;

  warning("No unmatched sequences have been found")
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

	$insert_start = $missed_fragment->[3] - 1 - 
	  ($missed_fragment->[1] - $missed_fragment->[0] + 1) - $spacing;

      } elsif ($missed_fragment->[2] eq 'after') {

	$insert_start = $missed_fragment->[3] + $spacing;

      } else {
	throw("Before or after not specified. Internal error.");
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
