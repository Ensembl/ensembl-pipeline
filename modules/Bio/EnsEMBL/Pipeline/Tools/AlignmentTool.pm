
# Cared for by Dan Andrews <dta@sanger.ac.uk>
#
# Copyright EnsEMBL
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME
  
Bio::EnsEMBL::Pipeline::Tools::AlignmentTool
  
=head1 SYNOPSIS

# Low stress (I just want an alignment, damnit):

my $alignment_tool = Bio::EnsEMBL::Pipeline::Tools::AlignmentTool->new(
                           '-dbadaptor'         => $db,
			   '-seqfetcher'        => $pfetcher,
			   '-transcript'        => $transcript);

my $alignment = $alignment_tool->retrieve_alignment('all');  
                             # Or just 'nucleotide' or 'protein'

foreach my $sequence (@$alignment){
  print $sequence . "\n";
}

# Fussy (set a few things to non-default, get worried about 
# alignment quality):

my $alignment_tool = Bio::EnsEMBL::Pipeline::Tools::AlignmentTool->new(
                           '-dbadaptor'         => $db,
			   '-seqfetcher'        => $pfetcher,
			   '-transcript'        => $transcript,
			   '-padding'           => 50,
			   '-fasta_line_length' => 60);


# Scrutinize each exon for the identity of best evidence:
my $exon_best_identities = $alignment_tool->identity('identity_by_exon');

### NOTE : The definition of identity used by this module ignores all
### gaps in the sequence.  Given than many of these alignments are
### gappy or fragmentary, including gaps in the identity score will
### dilute it somewhat according to coverage.

# ---> Everything below is presently unimplemented:
## Everything - an array of hashes containing match sequence id, coverage and 
## overall match identity.
#my $all_matches_stats = $alignment_tool->missing_method('all_evidence');
  
#print "Match ID         : " . $all_matches_stats->[0]->{'id'} . "\n" .
#      "Overall Identity : " . $all_matches_stats->[0]->{'identity'} . "\n" .
#      "Hit Coverage     : " . $all_matches_stats->[0]->{'coverage'} . "\n";

# Last, but certainly not least - find exons that have no evidence, or
# calculate the number of coding bases that have no evidence aligned
# to them.

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

package Bio::EnsEMBL::Pipeline::Tools::AlignmentTool;

use vars qw(@ISA);
use strict;

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
  Returntype : Bio::EnsEMBL::Pipeline::Tools::AlignmentTool
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
  
  $self->_slice($transcript->stable_id);

  # Due to padding, it is necessary to re-construct our transcript
  # in proper slice coordinates.

  $self->_transcript($transcript, $self->_slice);

  # The line length in the fasta alignment is set to a default
  # or a user specified value

  if ($fasta_line_length) {
    $self->{'_fasta_line_length'} = $fasta_line_length;
  } else {
    $self->{'_fasta_line_length'} = 60;
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
  my ($self, $type) = @_;

  unless ($type) {
    $self->throw("Type of alignment to retrieve not specified.  Please use one " . 
		 "of \'all\', \'nucleotide\' or \'protein\'.");
  }

  unless ($self->_is_computed($type)){
    $self->_align($type);
  }

  return $self->_alignment;
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
    $seen_exons{$sequence->{'exon'}}++;
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
  my ($self, $type) = @_;

  # Collect all our sequences for aligning

  my $genomic_sequence         = $self->_genomic_sequence;
  my $exon_nucleotide_sequence = $self->_exon_nucleotide_sequence;
  
  my $exon_protein_sequence;
  if (($type eq 'protein')||($type eq 'all')){
    $exon_protein_sequence = $self->_exon_protein_translation;
  }
  
  my $evidence_sequences = $self->_corroborating_sequences($type);

  # Insert use information about 'deletions' (originally from the
  # cigar string and stored in the hashes returned by 
  # $self->_semialigned_sequences) to insert gaps into the genomic 
  # and exonic sequences, therefore aligning the evidence with
  # the parent sequence (with a bit of fiddling to place
  # gaps into the evidence sequences that need them).

  for (my $i = 0; $i < $self->_slice->length; $i++) {
    my $is_deletion = 0;

  DELETION_HUNT:
    foreach my $unaligned_sequence (@$evidence_sequences){
      if ((defined $unaligned_sequence->{'deletions'}->[$i])
	  &&($unaligned_sequence->{'deletions'}->[$i] eq 'D')){
	$is_deletion = 1;
	last DELETION_HUNT;
      }
    }
    
    if ($is_deletion) {
      splice (@$genomic_sequence, $i, 0, '-');
      splice (@$exon_nucleotide_sequence, $i, 0, '-');

        # Rather than adding an if statement, just make the protein 
        # sequence and choose whether to use it later. 
      splice (@$exon_protein_sequence, $i, 0, '-');  
     


      for (my $j = 0; $j < scalar @$evidence_sequences; $j++) {
	unless ((defined $evidence_sequences->[$j]->{'deletions'}->[$i])
		&&($evidence_sequences->[$j]->{'deletions'}->[$i] eq 'D')) {
	  splice(@{$evidence_sequences->[$j]->{'seq'}}, $i, 0, '-');
	}
      }
    }
    
  }

  # Put our working alignments somewhere handy

  $self->_working_alignment('genomic_sequence', $genomic_sequence);
  $self->_working_alignment('exon_nucleotide', $exon_nucleotide_sequence);

  if (($type eq 'protein')||($type eq 'all')) {
    $self->_working_alignment('exon_protein', $exon_protein_sequence);
  }
  foreach my $evidence_sequence (@$evidence_sequences) {
    $self->_working_alignment('evidence', $evidence_sequence);
  }

  # Concatenate the sequences from their storage arrays and print

  my $genomic_string = '';
  my $exon_nucleotide_string = '';
  my $exon_protein_string = '';

  for (my $j = 0; $j < scalar @$evidence_sequences; $j++) {
    $evidence_sequences->[$j]->{'string_seq'} .= '';
  }

  for (my $i = 0; $i < $self->_slice->length; $i++) {
    my $multiple_of_line_length = 0;
    
    if (($i%($self->{'_fasta_line_length'}) == 0)&&($i != 0)){
      $genomic_string .= "\n";
      $exon_nucleotide_string .= "\n";
      $exon_protein_string .= "\n";
      $multiple_of_line_length = 1;
    }
    
    $genomic_string .= $genomic_sequence->[$i];
    
    if (defined $exon_nucleotide_sequence->[$i]){
      $exon_nucleotide_string .= $exon_nucleotide_sequence->[$i];
    } else {
      $exon_nucleotide_string .= '-';
    }
    
    if ((defined $exon_protein_sequence->[$i])
	&&(($type eq 'all')||($type eq 'protein'))){
      $exon_protein_string .= $exon_protein_sequence->[$i];
    } elsif (($type eq 'all')||($type eq 'protein')) {
      $exon_protein_string .= '-';
    }
    
    for (my $j = 0; $j < scalar @$evidence_sequences; $j++) {
      if ($multiple_of_line_length){
	$evidence_sequences->[$j]->{'string_seq'} .= "\n";
      }
      $evidence_sequences->[$j]->{'string_seq'} .= $evidence_sequences->[$j]->{'seq'}->[$i];
    }
    
  }

  # Place our aligned sequences in the alignment storage.  An
  # alignment object would be *so* useful here.

  my $genomic_alignment = '>' . $self->_transcript->stable_id. " genomic_sequence\n$genomic_string";
  $self->_alignment('add', $genomic_alignment);
  my $exon_alignment = '>' . $self->_transcript->stable_id. " exon_sequence\n$exon_nucleotide_string";
  $self->_alignment('add', $exon_alignment);

  if (($type eq 'protein')||($type eq 'all')){
    my $trans_alignment = '>' . $self->_transcript->stable_id . " translated_exons\n$exon_protein_string";
    $self->_alignment('add', $trans_alignment);
  }
  
  for (my $k = 0; $k < scalar @$evidence_sequences; $k++) {
    my $evidence_alignment =  '>' . $evidence_sequences->[$k]->{'name'} . "\n" . 
      $evidence_sequences->[$k]->{'string_seq'};
    $self->_alignment('add', $evidence_alignment);
  }
 
  # Set flag to indicate that the alignment has been computed.

  $self->_is_computed($type, 1);
 
  return $self->_alignment;

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
    push (@{$by_exon{$evidence_item->{'exon'}}}, $evidence_item);
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

      if (($evidence_item->{'type'} eq 'protein')
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

      elsif (($evidence_item->{'type'} eq 'nucleotide')
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
  my ($self, $exon, $evidence_item, $reference_sequence) = @_;

  # For nucleotide alignments each mismatch is counted
  # once.
  my $align_unit = 1;

  # If we are dealing this protein alignments we have to
  # multiply this by three.
  $align_unit *= 3 if ($evidence_item->{'type'} eq 'protein');

  my $match_sequence = $evidence_item->{'seq'};
  
  my $mismatches = 0;
  my $noncovered = 0;
  
  my $exon_start = $exon->start;
  my $exon_end = $exon->end;
  my $exon_length = $exon_end - $exon_start;

  if ($self->_strand == -1){
    $exon_start = $self->_slice->length - $exon_end + 1;
    $exon_end = $exon_start + $exon_length - 1;
  }

#print STDERR "Exon start : " . $exon_start . "\tExon end: " . $exon_end . "\n"; 
  for (my $i = $exon_start - 1; $i < $exon_end; $i++) {
#print STDERR $i . " " . $reference_sequence->[$i] . " " . $match_sequence->[$i] . "\n";
    unless (defined ($match_sequence->[$i]) &&
	    defined ($reference_sequence->[$i]) &&
	    (($reference_sequence->[$i] eq $match_sequence->[$i])||
	     (($reference_sequence->[$i] eq '-')
	      ||($match_sequence->[$i] eq '-')))) {
#print STDERR "MISMATCH\n";
      $mismatches += $align_unit;
    }
    
    if (($reference_sequence->[$i] ne '-')
	&&($match_sequence->[$i] eq '-')) {
#print STDERR "NONCOVERED\n";
      $noncovered += $align_unit;
    }
  }
  
  my $identity = (1 - ($mismatches/$exon_length))*100;
  
  # The next line gets around the problem of exon length not always
  # being a whole number of codons.  There can be cases where
  # there are more non-covered bases than there are bases in an exon.
  $noncovered = $exon_length if $noncovered > $exon_length;
  
  my $coverage = (1 - ($noncovered/$exon_length))*100;
  
#print STDERR "Identity : $identity\tCoverage : $coverage\tNoncovered : $noncovered\tExon length : $exon_length\n";
  
  return ($identity, $coverage);
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
    print "Alignment has been previously computed, but was not\n" .
      "of the same type.  The previously computed alignment\n" . 
	"type was \'" . $self->_type . "\'.\n";

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


=head2 _alignment

  Arg [1]    :
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub _alignment {

  my ($self, $action, $sequence) = @_;

  if (defined $action && $action eq 'add' && defined $sequence){

    # It would be nice to check the sequence for the 
    # correct format before adding it to our alignment

    push (@{$self->{'_alignment_array'}}, $sequence);
  }

  return $self->{'_alignment_array'};
}


=head2 _working_alignment

  Arg [1]    :
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub _working_alignment {

  my ($self, $slot, $seq_hash) = @_;

  unless (defined $slot && 
      (($slot eq 'genomic_sequence')
       ||($slot eq 'exon_protein')
       ||($slot eq 'exon_nucleotide')
       ||($slot eq 'evidence'))){
    $self->throw("Was trying to retrieve or write working alignments to "
		 . "a slot that isnt allowed ($slot)");
  }

  if (defined $slot && defined $seq_hash){

    push (@{$self->{'_working_alignment_array'}->{$slot}}, $seq_hash);

  } elsif (defined $slot && !defined $seq_hash) {

    my $slot_resident =  $self->{'_working_alignment_array'}->{$slot};

    if ((($slot eq 'genomic_sequence')||($slot eq 'exon_protein')||($slot eq 'exon_nucleotide')) 
	&& defined $slot_resident && scalar @$slot_resident == 1) {
      return $slot_resident->[0];
    }

    return $slot_resident;

  } 

  return 0;
}

##### Fiddlings with Slices and Transcripts. #####


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
	unless ($candidate_transcript->stable_id eq $input_transcript->stable_id){
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
    $self->throw("Transcript has not yet been extracted from new Slice.");
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

  my ($self, $transcript_stable_id) = @_;

  if (! defined $self->{'_slice'} && defined $transcript_stable_id) {

    my $slice_adaptor = $self->{'_db_adaptor'}->get_SliceAdaptor;

    $self->{'_slice'} = 
      $slice_adaptor->fetch_by_transcript_stable_id($transcript_stable_id, 
						    $self->{'_transcript_padding'});
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

    return
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

	my $half_aligned = $self->_fiddly_bits($base_align_feature);
	next FEATURE unless $half_aligned;

	$half_aligned->{'exon'} = $exon_placemarker;
	$half_aligned->{'type'} = 'nucleotide';

	push (@{$self->{'_corroborating_sequences'}}, $half_aligned);
	next FEATURE;
      }
      
      if ($base_align_feature->isa("Bio::EnsEMBL::DnaPepAlignFeature")){

	my $half_aligned = $self->_fiddly_bits($base_align_feature);
	next FEATURE unless $half_aligned;

	$half_aligned->{'exon'} = $exon_placemarker;
	$half_aligned->{'type'} = 'protein';

	push (@protein_features, $half_aligned);
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

  # Create a hash that has all the information we want to return.

  my %partially_aligned = ('name'      => $base_align_feature->hseqname,
			   'start'     => $base_align_feature->start,
			   'seq'       => \@feature_sequence,
			   'deletions' => \@deletion_sequence);
  
  return \%partially_aligned;
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

    my @genomic_sequence = split //, $genomic_sequence;

    # Fill in the blanks

    for (my $i = 0; $i < $self->_slice->length; $i++) {
      unless (defined $genomic_sequence[$i]){
	$genomic_sequence[$i] = '-';
      } 
    }

    $self->{'_genomic_sequence'} = \@genomic_sequence;

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
      
#      if ($exon_position != ($exon->end)){
#	$self->throw("Problem with exon length or sequence.\n$@");
#      }      
    }    
    
    # Fill in the blanks

    for (my $i = 0; $i < $self->_slice->length; $i++) {
      unless (defined $exon_only_sequence[$i]){
	$exon_only_sequence[$i] = '-';
      } 
    }

    $self->{'_exon_nucleotide_sequence'} = \@exon_only_sequence;
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
      
#      unless (($exon->phase == 0)||($exon->phase == -1)) {
#    }

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

    $self->{'_exon_protein_translation'} = \@exon_translation_sequence;
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

  eval {
    $fetched_seqs = $self->_seq_fetcher->batch_fetch(@array_of_accessions);
  };

  if ($@){
    $self->warn("Not all evidence sequences could be pfetched.\n".
		"Ignoring missing sequences.\n$@\n");
  }
  
  # Build cache.
print STDERR "Have " . scalar @$fetched_seqs . " fetched sequences.\n";
my $countthing = 0;
 BAILOUT:
  foreach my $fetched_seq (@$fetched_seqs){
print STDERR "Count : $countthing\t\tAccession : $array_of_accessions[$countthing]\n";
$countthing++;
unless ($fetched_seq) {
print STDERR "BAILING OUT\n";
  next BAILOUT;
}
#print STDERR "Count : $countthing\tFetched seq : " . $fetched_seq . "\n";
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
