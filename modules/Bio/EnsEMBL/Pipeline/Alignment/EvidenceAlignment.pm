
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

Quick start - use the following code if you want the 
alignment of an ensembl transcript with the evidence 
used to predict it:

my $evidence_alignment = 
  Bio::EnsEMBL::Pipeline::Alignment::EvidenceAlignment->new(
      -transcript          => $transcript,
      -dbadaptor           => $db,
      -seqfetcher          => $seqfetcher,
      -padding             => 0);

my $alignment = 
  $evidence_alignment->retrieve_alignment(
     '-type'            => 'all',
   # Type can be 'all', 'nucleotide' or 'protein'
     '-remove_introns'  => 1,
     '-merge_sequences' => 1);

my $align_seqs = $alignment->fetch_AlignmentSeqs;

foreach my $align_seq (@$align_seqs){
  print $align_seq->seq . "\n";
}

The '-padding' option specifies the amount of tailing
sequence left attached to each 'exon'.  If you set this
to zero be careful, as you could be truncating some
sequence from the aligned evidence sequences.  If this 
happens a warning is issued.

Other alignment presentation options exist.  The intronic
regions of the alignment can make it _very_ big.  Introns
can be truncated thusly:

my $alignment = 
  $align_tool->retrieve_alignment(
      '-type'            => 'all',
      '-remove_introns'  => 1,
      '-merge_sequences' => 1);

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

It is possible to display a transcript with a set of external 
supporting features that can come from any source.  Where external
features are used, the supporting features attached to
the transcript are ignored and the external set used instead.  It
usually helps if the external supporting features actually overlap
the transcript sequence :)  This option is used by passing in 
supporting features at the time of object creation.

my $alignment_tool = 
  Bio::EnsEMBL::Pipeline::Alignment::EvidenceAlignment->new(
     '-dbadaptor'           => $db,
     '-seqfetcher'          => $seqfetcher,
     '-transcript'          => $transcript,
     '-supporting_features' => \@supporting_features);


A few features of the actual alignment can be controled.
The fasta line length and the amount of 5-prime and
3-prime padding can be stipulated:

my $alignment_tool = 
  Bio::EnsEMBL::Pipeline::Alignment::EvidenceAlignment->new(
     '-dbadaptor'         => $db,
     '-seqfetcher'        => $pfetcher,
     '-transcript'        => $transcript,
     '-padding'           => 50,
     '-fasta_line_length' => 60);


Once the alignment is generated it is possible to calculate
the identity of the best matching evidence for each exon.
Coverage is also determined.

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

my $evidence_coverage = $alignment_tool->hit_coverage;
  
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
use Bio::EnsEMBL::Utils::Exception qw(throw warning info);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

@ISA = qw();


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
      $evidence_identity_cutoff,
      $supporting_features) = rearrange([qw(DBADAPTOR
					    TRANSCRIPT
					    SEQFETCHER
					    PADDING
					    FASTA_LINE_LENGTH
					    IDENTITY_CUTOFF
					    SUPPORTING_FEATURES
					   )],@args);

  # Throw an error if any of the below are undefined or
  # are the wrong thing.
  unless (defined $db && $db->isa("Bio::EnsEMBL::DBSQL::DBAdaptor")){
    throw("No DB adaptor passed to AlignmentTool.  You passed a $db.");
  }
  unless (defined $transcript && $transcript->isa("Bio::EnsEMBL::Transcript")){
    throw("No transcript passed to AlignmentTool.  You passed a $transcript.");
  }
  unless (defined $seqfetcher && $seqfetcher->isa("Bio::DB::RandomAccessI")) {
    throw("No sequence fetcher passed to AlignmentTool.  You passed a $seqfetcher.");
  }

  $self->{'_db_adaptor'} = $db;

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

  # Due to padding, it is necessary to re-construct our transcript
  # in proper slice coordinates.

  $self->_transcript($transcript, $self->_slice);


  # Determine if our database contains translations.  If it doesn't
  # we'll have to skip adding a set of translated exons to our 
  # alignment.

  if ($self->_transcript->translation){
    $self->_translatable(1);
  } else {
    warning("Database doesn't contain translation.  Subsequently, ".
	    "wont be able to display translations of each exon ".
	    "or calculate protein identity scores.");
    $self->_type('nucleotide');
    $self->_translatable(0);
  }

  # When an external set of features are used, over-ride the
  # features from the transcript

  if (defined $supporting_features) {
    $self->_all_supporting_features($supporting_features)
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
      $show_missing_evidence,
      $merge_sequences,
      $three_letter_aa) = rearrange([qw(TYPE
					REMOVE_INTRONS
					SHOW_UNALIGNED
					MERGE_SEQUENCES
					THREE_LETTER_AA
				       )],@_);

  $self->_type($type);
  $merge_sequences = 1 unless defined $merge_sequences;
  $self->_three_letter_aa($three_letter_aa) if $three_letter_aa;

  unless ($self->_is_computed($type)){
    my $alignment_success = $self->_align($type, 
					  $show_missing_evidence, 
					  $remove_introns,
					  $merge_sequences);
    unless ($alignment_success) {
      warning "Alignment generation failed.  There probably were" . 
	" not any displayable sequences.";
      return 0
    }
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
      $merge_sequences) = @_;

  my $evidence_sequences = $self->_corroborating_sequences($type);

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
    # Note:  Purposely using the memory address for our evidence sequence objects
    # as they otherwise have non-unique names.

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
      }
    }

    # Increment all deletion coords that are greater than this one.
    # There quite possibly is a better way to do this.

      # @all_deletions and %all_deletions

    my %new_all_deletions;

    for (my $i = 0; $i < scalar @all_deletions; $i++) {
      if ($all_deletions[$i] > $deletion_coord) {
        $all_deletions[$i]++;
	$new_all_deletions{$all_deletions[$i]} = $all_deletions{$all_deletions[$i] - 1};
      }
    }

    %all_deletions = %new_all_deletions;


      # %deletion_sets

    my %new_deletion_sets;
    foreach my $evidence_name (keys %evidence_sequence_hash) {
      my @coords = keys %{$deletion_sets{$evidence_name}};
      for (my $i = 0; $i < scalar @coords; $i++) {
	if ($deletion_sets{$evidence_name}->{$coords[$i]}){
	  $new_deletion_sets{$evidence_name}->{$coords[$i]+1} = 
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
	$new_deletion_tracking{$coords[$i]+1} = $deletion_tracking{$coords[$i]}
      } else {
	$new_deletion_tracking{$coords[$i]} = $deletion_tracking{$coords[$i]}
      }
    }
    %deletion_tracking = %new_deletion_tracking;
  }

  # Print the locations of our deletions, handy for finding conserved gaps and
  # frameshifts

  foreach my $tracked_deletion (sort {$a <=> $b} (keys %deletion_tracking)){
#    print STDOUT $tracked_deletion . "\t" . 
      scalar @{$deletion_tracking{$tracked_deletion}} . "\t";
#    foreach my $deletion (@{$deletion_tracking{$tracked_deletion}}){
#      print STDOUT $deletion->[1] ." (" . $deletion->[0] . ")  ";
#    }
#    print STDOUT "\n";
  }


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
    $self->_truncate_introns;
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
	if ($self->_strand == 1){
	  $self->_add_unmatched_region($sequence_identifier, $previous_end, 
				       $coordinate_pair->[0], 'before', $coordinate_pair->[2]);
	  
	} elsif ($self->_strand == -1) {
	  $self->_add_unmatched_region($sequence_identifier, $previous_end, 
				       $coordinate_pair->[0], 'after', 
				       ($self->_slice->length - $coordinate_pair->[3]));
	}
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


=head2 _truncate_introns

  Arg [1]    :
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

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

    if ($self->_strand == -1) {
      $intron_coord_1 = $self->_slice->length - ($exon->start - 1) + 1;
      $intron_coord_2 = $self->_slice->length - ($exon->end + 1) + 1;
    }
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

  # Blam!  Obliterate evidence sequences in favour of the new
  # merged sequences that are about to be created.

  $self->_working_alignment('evidence', 'empty');

  # Group evidence by sequence id.

  my %same_sequences_hash;

  foreach my $sequence (@sequences) {

    push (@{$same_sequences_hash{$sequence->name}}, $sequence);
  }

  # Work through each distinct sequence merging all fragments into
  # a single string.

  foreach my $sequence_name (keys %same_sequences_hash){

    my @component_sequences;

    foreach my $unmerged_seq (@{$same_sequences_hash{$sequence_name}}){

      push @component_sequences, [$unmerged_seq, 
				  $unmerged_seq->first_base_coord, 
				  $unmerged_seq->last_base_coord];
    }

    @component_sequences = sort {$a->[1] <=> $b->[1]} @component_sequences;

    # Check for data problems.  Sequence from the same piece of evidence
    # should not overlap with itself.

    for (my $i = 0; $i < scalar @component_sequences - 1; $i++){
      if ($component_sequences[$i]->[2] >= $component_sequences[$i+1]->[1]) {
	warning("Parts of the same evidence sequence [" . 
		$component_sequences[$i]->[0]->name . "] have " . 
		"been matched to this gene twice.  Fiddling " . 
		"alignment coordinates to make things work.");

	  $component_sequences[$i+1]->[1] = $component_sequences[$i]->[2] + 1;
      }
    }

    # Merge our sequences

    my $gap_start = 1;
    my $merged_string = '';

    foreach my $component_sequence (@component_sequences){
      my $evidence_start  = $component_sequence->[1];
      my $evidence_length = $component_sequence->[2] - $evidence_start + 1;

      $merged_string .= 
	('-' x ($evidence_start - $gap_start)) . 
	substr($component_sequence->[0]->seq, 
	       $evidence_start - 1, 
	       $evidence_length);

      $gap_start = $component_sequence->[2] + 1;
    }

      # Add last gap

    $merged_string .= 
      '-' x ($component_sequences[-1]->[0]->length - 
	     $component_sequences[-1]->[2]);


    # Make new align seq and store

    my $merged_align_seq = 
      Bio::EnsEMBL::Pipeline::Alignment::AlignmentSeq->new(
       -name => $component_sequences[0]->[0]->name,
       -type => $component_sequences[0]->[0]->type,
       -seq  => $merged_string);

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
    unless ($type eq 'all' ||
	    $type eq 'nucleotide' ||
	    $type eq 'protein') {
      throw("Type of alignment to retrieve not recognised.  Please use one " . 
	    "of \'all\', \'nucleotide\' or \'protein\'.");
    }

    $self->{'_computed_type'} = $type;
  }

  throw("Type of alignment to retrieve not specified.  Please use one " . 
	"of \'all\', \'nucleotide\' or \'protein\'.") 
    unless $self->{'_computed_type'};

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
    throw("Was trying to retrieve or write working alignments to "
	  . "a slot that isnt allowed ($slot)");
  }

  if (defined $slot && defined $align_seq && $align_seq eq 'empty'){
    undef $self->{'_working_alignment_array'}->{$slot};
    return 1
  }

  if (defined $slot && defined $align_seq){

    unless ($align_seq->isa("Bio::EnsEMBL::Pipeline::Alignment::AlignmentSeq")){
      throw("Sequence passed to _working alignment was not an " . 
	    "AlignmentSeq, it was a [$align_seq]")
    }

    push (@{$self->{'_working_alignment_array'}->{$slot}}, $align_seq);

  } elsif (defined $slot && !defined $align_seq) {

    my $slot_resident =  $self->{'_working_alignment_array'}->{$slot};

    if ((($slot eq 'genomic_sequence')||($slot eq 'exon_protein')||($slot eq 'exon_nucleotide')) 
	&& defined $slot_resident && scalar @$slot_resident == 1) {
      return $slot_resident->[0];
    }

    # Sort evidence so that nucleotide and protein sequences 
    # are returned in blocks of similar sequence type.
    if ($slot eq 'evidence' && $self->_type eq 'all') {
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
      throw("Could not find transcript on Slice.  Very bad.");
    }
  }

  if (defined $self->{'_transcript'}){
    return $self->{'_transcript'}
  } else {
    throw("Something has gone wrong.  A Transcript has not yet been" .
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
						      $self->_padding);
    } else {

      $self->{'_slice'} = 
	$slice_adaptor->fetch_by_transcript_id($transcript->dbID, 
					       $self->_padding);

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
    warning("No value for strand set.  Disaster awaits.");
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

  # Work through supporting features and create a sequence for each.

  my @protein_features;

  my $exon_placemarker = 0;

  # Work through each item of supporting evidence attached to our exon.

 FEATURE:
  foreach my $base_align_feature (@{$self->_all_supporting_features}){

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

    my $align_seq = $self->_build_evidence_seq($base_align_feature);
    next FEATURE unless $align_seq;

    $align_seq->exon($exon_placemarker);

    if ($base_align_feature->isa("Bio::EnsEMBL::DnaDnaAlignFeature")){
      $align_seq->type('nucleotide')
    } elsif ($base_align_feature->isa("Bio::EnsEMBL::DnaPepAlignFeature")) {
      $align_seq->type('protein');

      if ($base_align_feature->hstrand == -1) {
	warning("It seems we have a reverse complimented " . 
		"protein sequence.  Surely not.");
      }
    }

    push (@{$self->{'_corroborating_sequences'}}, $align_seq);

    $exon_placemarker++;
  }

  unless (scalar @{$self->{'_corroborating_sequences'}} > 0) {
    warning("There are no displayable supporting features for this transcript [" .
	    $self->_transcript->stable_id . "].  " .
	    "Try setting the -type attribute to 'all', instead of just " . 
	    "'protein' or 'nucleotide'.") 
  }

  return $self->{'_corroborating_sequences'};
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

  # Create an AlignmentSeq object that will store our sequences
  # for this feature

  my $partially_aligned = 
    Bio::EnsEMBL::Pipeline::Alignment::AlignmentSeq->new(
          '-name'       => $base_align_feature->hseqname,
	  '-start'      => $base_align_feature->start);

  # For convenience, localise a few coordinate variables.

  my $hstart  = $base_align_feature->hstart;
  my $hend    = $base_align_feature->hend;
  my $cigar   = $base_align_feature->cigar_string;
  my $hstrand = $base_align_feature->hstrand;

    # Take care of unusual situation where genomic strand gets turned around.
  if (($self->_strand != $base_align_feature->strand)){
      # Force the hstrand around
    $hstrand = $hstrand * -1;
      # Reverse the order of things in the cigar
    my @cigar_items = $self->_cigar_reader($cigar);
    @cigar_items = reverse @cigar_items;
    $cigar = '';
    foreach my $item (@cigar_items){
      $cigar .= $item->{'length'} unless $item->{'length'} == 1;
      $cigar .= $item->{'type'};
    }
  }

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
	   ($fetched_seq->seq eq '-' x length($fetched_seq->seq))) {
    warning("Error fetching sequence [" . 
	    $base_align_feature->hseqname . 
	    "].  Sequence is an empty string.");

    return 0;
  }

  my @fetched_seq;

  # If this is a protein align feature, pad amino acids with gaps
  # to make them comparable to nucleotide coords.  Then, make sure
  # splice out our sequence of interest using hit coordinates that
  # take account of the padding of our sequence.

  if ($base_align_feature->isa("Bio::EnsEMBL::DnaPepAlignFeature")){
    my $padded_aa_seq;
    if ($self->_three_letter_aa) {
      ($padded_aa_seq = $fetched_seq->seq) =~ s/(.)/$self->{_aa_names}{$1}/g;
    } else {
      ($padded_aa_seq = $fetched_seq->seq) =~ s/(.)/$1\-\-/g;
    }

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

    # The coordinates for dna_align_features are stored with reference
    # to the way the sequence is was found in dbEST or whatever.  Hence,
    # if the EST is the reverse strand, it is necessary to chop out our 
    # fragment first, then reverse compliment it.

    if ($hstrand == -1) {
      my $wrong_way_around_exon = join '', @fetched_seq;

      my $backwards_seq = Bio::Seq->new('-seq' => $wrong_way_around_exon);

      my $forwards_seq = $backwards_seq->revcom;

      @fetched_seq = split //, $forwards_seq->seq;

    }
  }

  # Add the needed insertion gaps to our supporting
  # feature sequence.  We only add gaps to the supporting
  # sequences at this stage.  Once we have all our supporting
  # sequences constructed and determined where all the 
  # 'deletion' gaps lie in these, we can transpose this 
  # information onto our genomic and exonic sequences - this is
  # done in $self->_align.

  my @cigar_instructions = $self->_cigar_reader($cigar);

  my $hit_position = 1;

    # Note to self - all genomic coordinates are actually the coordinates
    # within the slice.  This slice is not in chromosomal coordinates.  As
    # the slice was derived using a transcript, the slice will be oriented
    # in the same direction as the transcript.  As the transcript alignment 
    # will be displayed 5'-3' the slice will be in the correct orientation
    # with respect to the transcript and thus the genomic strand can be 
    # treated as being on the forward strand.

  my $slice_position = $base_align_feature->start;
  my $slice_rev_position = $self->_slice->length - $base_align_feature->end + 1;

  my $deletions;
  my $deletions_upstream_of_slice = 0;

  foreach my $instruction (@cigar_instructions) {

    if ($instruction->{'type'} eq 'I') {

      my $gap = '-' x $instruction->{'length'};
      my @gap = split //, $gap;

      if ($hit_position <= scalar @fetched_seq){
	splice(@fetched_seq, $hit_position - 1, 0, @gap);
      } else {
	throw("Gap in sequence [" . $base_align_feature->hseqname . 
	      "] lies outside feature.  Code problem.")
      }

      $hit_position   += $instruction->{'length'};
      $slice_position += $instruction->{'length'};
      $slice_rev_position += $instruction->{'length'};

    } elsif ($instruction->{'type'} eq 'M') {

      $hit_position   += $instruction->{'length'};
      $slice_position += $instruction->{'length'};
      $slice_rev_position += $instruction->{'length'};

    } elsif ($instruction->{'type'} eq 'D') {
      $deletions++;

      if ($slice_position <= 0){
	$deletions_upstream_of_slice++;
	next
      }

      if ($self->_strand == 1) {
	$partially_aligned->add_deletion($slice_position, 
					 $instruction->{'length'});
      } elsif ($self->_strand == -1) {
	$partially_aligned->add_deletion($slice_rev_position, 
					 $instruction->{'length'});
      }

      $hit_position   += $instruction->{'length'};
      $slice_position += $instruction->{'length'};
      $slice_rev_position += $instruction->{'length'};
    }
  }

  if ($deletions) {
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
	    "overhanging sequence");

    if (($base_align_feature->start < 0 && 
	 $base_align_feature->end < 0)||
	($base_align_feature->start > $self->_slice->length &&
	 $base_align_feature->end > $self->_slice->length)) {
      info("Feature [" . $base_align_feature->hseqname . 
	   "] lies completely outside the bounds of Slice.  Chuck.");
      return 0
    } elsif ($self->_strand == 1) {
      my $start_overshoot = 0;
      if ($base_align_feature->start < 0) {
#print STDERR "Case 1.\n";
	$start_overshoot = ($base_align_feature->start * -1) + $deletions_upstream_of_slice;
	splice (@fetched_seq, 0, $start_overshoot + 1);
      }
      if ($base_align_feature->end > $self->_slice->length) {
#print STDERR "Case 2.\n";
	my $end_overshoot = $base_align_feature->end - $self->_slice->length - 1;
	splice (@fetched_seq, (scalar @fetched_seq) - $start_overshoot - 1, $end_overshoot);
      }
    } elsif ($self->_strand == -1) {
      my $start_overshoot = 0;
      if ($base_align_feature->end < 0) {
#print STDERR "Case 3.\n";
	$start_overshoot = ($base_align_feature->end * -1) + $deletions_upstream_of_slice;
	splice (@fetched_seq, 0, $start_overshoot + 1);
      }
      if ($base_align_feature->start > $self->_slice->length) {
#print STDERR "Case 4.\n";
	my $end_overshoot = $base_align_feature->start - $self->_slice->length - 1;
	splice (@fetched_seq, (scalar @fetched_seq) - $start_overshoot - 1, $end_overshoot);
      }
    }
  }

  # Here we are actually building the sequence that will
  # align to our slice

  my $feature_sequence = '-' x $self->_slice->length;
  my @feature_sequence = split //, $feature_sequence;

  my $genomic_start;

  if ($self->_strand == 1){
    $genomic_start = $base_align_feature->start > 0 ? $base_align_feature->start - 1 : 0;
  } elsif ($self->_strand == -1) {
    $genomic_start = $base_align_feature->end > 0 ? $base_align_feature->end - 1 : 0;
  }

  if ($genomic_start <= scalar @feature_sequence){
    if ($self->_strand == -1){
      $genomic_start = scalar @feature_sequence - $genomic_start - 1;
    }

    splice (@feature_sequence, $genomic_start, (scalar @fetched_seq), @fetched_seq)
  } else {
    info("Feature [". $base_align_feature->hseqname . " start:" . 
	 $base_align_feature->start . " end:" . $base_align_feature->end .
	 " strand:" . $base_align_feature->strand 
	 ."] starts beyond end of genomic slice.  This feature most probably " .
	 "aligns to an exon that is not part of this transcript.");
    return 0
  }

  $feature_sequence = join '', @feature_sequence;
  $partially_aligned->seq($feature_sequence);

  return $partially_aligned;
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

  if (!defined $self->{'_genomic_sequence'}) {

    my $genomic_sequence;

    if ($self->_strand == 1) {
      $genomic_sequence = $self->_slice->seq;
    } elsif ($self->_strand == -1) {
#print STDERR "Reverse complimenting genomic sequence.\n";
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
	  throw "Overlapping exons\n$@";
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

    my $exon_sequence = join('', @exon_only_sequence);

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
      }

      # Jiggling the exons about to get frame right
      my $extra_length = 0;

      $extra_length += 3 if (($exon_end_phase != 0) && ($exon_end_phase != -1));
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

    $self->{'_exon_protein_translation'} = 
      Bio::EnsEMBL::Pipeline::Alignment::AlignmentSeq->new(
	 '-seq'  => $translated_exon_sequence,
	 '-name' => 'translated_exon_sequence',
	 '-type' => 'protein');
  }

  return $self->{'_exon_protein_translation'};
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
    foreach my $exon (@{$self->_transcript->get_all_Exons}){
      push @{$self->{_all_supporting_features}}, 
	       @{$exon->get_all_supporting_features};
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
  my @filtered_features;

  foreach my $feat (@$features){
    my $feat_identifier 
      = $feat->hseqname . $feat->hstart . $feat->hend . $feat->hstrand;

    unless ($feature_lookup{$feat_identifier}){
      push @filtered_features, $feat;
    }

    $feature_lookup{$feat_identifier}++
  }

  return \@filtered_features
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
    info("Sequence $accession could not be retrieved from cache.");
  }

  return $self->{'_fetched_seq_cache'}->{$accession}; 
}


### Miscellaneous utilities ###

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
    $self->{'_padding'} = shift;
  }

  return $self->{'_padding'} ? $self->{'_padding'} : 0;
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
    $self->{'_three_letter_aa'} = shift;
  }
  if ($self->{'_three_letter_aa'}){
    $self->_set_aa_names;
  }


  return $self->{'_three_letter_aa'} ? $self->{'_three_letter_aa'} : 0;
}

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
                        'Z' => 'Glx'}
}

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
    $self->{'_there_are_deletions'} = shift;
  }

  return $self->{'_there_are_deletions'};
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

  while (scalar @cigar_array) {

    my $next_char = shift @cigar_array;

    if ($next_char =~ /[MDI]/) {

      $current_digits = 1
	unless ($current_digits && $current_digits > 0);

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
