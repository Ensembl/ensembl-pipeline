package Bio::EnsEMBL::Pipeline::Alignment::InformativeSites;

use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Pipeline::GeneDuplication::CodonBasedAlignment;
use Bio::EnsEMBL::Utils::Exception qw(throw warning info); 
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

# Object preamble - inherits from Bio::Root::Object;

@ISA = qw(Bio::EnsEMBL::Root);


sub new {
  my ($class, @args) = @_;

  my $self = bless {},$class;

  my ($translatable_seqs,
      $alignments
     ) = rearrange([qw(TRANSLATABLE_SEQS
		       ALIGNMENTS
                      )],@args);

  throw("No sequences to align.")
    unless $translatable_seqs;

  $self->_seqs($translatable_seqs);

  throw("No alignments to reconcile.")
    unless (defined $alignments);

  $self->_alignments($alignments);

  throw("Ids of translatable sequences do not fully match with " . 
	"alignment ids.  This sometimes happens because one of " . 
	"the clustered genes has zero EST matches and hence no " . 
	"alignment.")
    unless $self->_check_ids;

  return $self
}

sub informative_sites {
  my $self = shift;

  my $gene_informative_sites = $self->_informative_site_strings;

  my %inform_site_aligns;

  foreach my $seq (@{$self->_seqs}) {
    # Avoid making sequences where no informative sites exist
    # ie. identical sequences.
    next 
      unless $gene_informative_sites->{$seq->display_id} =~ /\*/; 

    $inform_site_aligns{$seq->display_id} = 
      $self->_derive_informative_site_alignment(
	 $seq, 
         $gene_informative_sites->{$seq->display_id});
  }

  return \%inform_site_aligns
}


sub _seqs {
  my $self = shift;

  if (@_) {
    $self->{_seqs} = shift;

    throw("Was expecting the translatable sequences to be Bio::Seq objects.")
      unless (scalar @{$self->{_seqs}} > 0 &&
	      defined $self->{_seqs}->[0] &&
	      $self->{_seqs}->[0]->isa("Bio::Seq")); 

    foreach my $seq (@{$self->{_seqs}}){
      throw("Gaps found in translatable sequence.  These should be gap-free.")
	if $seq->seq =~ /-/;
    }
  }

  return $self->{_seqs}
}

sub _alignments {
  my $self = shift;

  if (@_) {
    $self->{_aligns} = shift;

    throw("Was expecting alignments to be " .
	  "Bio::EnsEMBL::Pipeline::Alignment objects.")
      unless (scalar @{$self->{_aligns}} > 0 &&
	      defined $self->{_aligns}->[0] &&
	      $self->{_aligns}->[0]->isa("Bio::EnsEMBL::Pipeline::Alignment"));

    $self->{_align_vs_id} = undef;
  }

  return $self->{_aligns}
}

sub _retrieve_align_by_id {
  my ($self, $id) = @_;

  unless (defined $self->{_align_vs_id}) {
    foreach my $align (@{$self->_alignments}){
      $self->{_align_vs_id}->{$align->alignment_name} = $align
    }
  }

  warning ("Alignment with identifier [$id] does not exist.")
    unless $self->{_align_vs_id}->{$id};

  return $self->{_align_vs_id}->{$id}
}

sub _check_ids {
  my $self = shift;

  my %alignment_names;
  my $verdict = 1;

  foreach my $align (@{$self->_alignments}){
    $alignment_names{$align->alignment_name}++;
  }

  foreach my $seq (@{$self->_seqs}){
    $verdict = 0 
      unless $alignment_names{$seq->display_id}
  }

  return $verdict
}

sub _informative_site_strings {
  my $self = shift;

  # Do a codon based alignment on our translatable sequences

  my $cba = 
    Bio::EnsEMBL::Pipeline::GeneDuplication::CodonBasedAlignment->new(
      -genetic_code => 1);

  $cba->sequences($self->_seqs); 
  my $aligned_seqs = $cba->run_alignment;

  # Work through the length of this alignment marking sequence
  # variations/informative sites.  Informative site information
  # is stored as a long string, the same length as each UNALIGNED
  # sequence.  Informative sites are marked as '*' and uninformative
  # sites are marked as '-'

  my @seq_arrays;
  my @seq_name;
  foreach my $realigned_seq (@$aligned_seqs) {
    my @seq_array = split //, $realigned_seq->seq;
    push @seq_arrays, \@seq_array;
    push @seq_name, $realigned_seq->display_id;
  }

  my $alignment_length = length($aligned_seqs->[0]->seq);

  my %informative_site_string;
  for (my $i = 0; $i < $alignment_length; $i++) {
    my $informative_site = 0;
    my $reference_base = $seq_arrays[0]->[$i];
    for (my $j = 1; $j < scalar @seq_arrays; $j++) {
      if ($seq_arrays[$j]->[$i] ne $reference_base) {
	$informative_site = 1;
	last
      }
    }
    for (my $j = 0; $j < scalar @seq_arrays; $j++) {
      if ($informative_site && $seq_arrays[$j]->[$i] ne '-') {
	$informative_site_string{$seq_name[$j]} .= '*';
      } elsif ($seq_arrays[$j]->[$i] ne '-') {
	$informative_site_string{$seq_name[$j]} .= '-';
      }
    }
  }

  @seq_arrays = undef;

  return \%informative_site_string
}

sub _derive_informative_site_alignment {
  my ($self, $seq, $informative_sites_string) = @_;

  # Put gene sequence somewhere handy

  my $gene_string = $seq->seq;
  my $gene_length = length($gene_string);
#print STDERR ">gene seq\n" . $gene_string . "\n";
  # Take our existing external EvidenceAlignment and
  # put the evidence sequences somewhere accessible.

  my $ext_align = $self->_retrieve_align_by_id($seq->display_id);

  my %evidence_strings;
  foreach my $align_seq (@{$ext_align->fetch_AlignmentSeqs}){
    $evidence_strings{$align_seq->name} = $align_seq->seq;
  }

  my $align_length = length($evidence_strings{'genomic_sequence'});

  # Place the sequence for gene onto the identical sequence in the
  # external alignment.  Start by finding the 5 prime end of the 
  # gene and slide along the length of the alignment, matching identical
  # bases and hopping gaps in the sequence and intron.  Where
  # an informative site is apparently, splice out these bases from all
  # evidence sequences to create strings of informative sites.

  my %aligned_informative_site_strings;
  my $exon_string = $evidence_strings{'exon_sequence'};
  $exon_string =~ s/intron-truncated/----------------/g;
#print STDERR ">exon seq\n$exon_string\n";
  my $align_coord = 0;
  my $lost = 1;

  for (my $gene_coord = 0; $gene_coord < $gene_length; $gene_coord++) {
#print STDERR "Gene coord : " . $gene_coord . "\n";
    # This prone-to-breakage code loop shuffles along the alignment
    # and tries and place the unaligned gene.  As the gene sequence should be
    # base-for-base identical to the exonic sequence in the alignment (barring
    # missing UTR sequence), the code tries to place 10bp substrs of the gene
    # against the aligned exonic sequence.
    if ($lost) {
      my $gene_seed_string = substr($gene_string, $gene_coord, 10);
#print STDERR "  Gene seed : " . $gene_seed_string . "\n";
    LOST:
      while (1){
	my $exon_seed_string = substr($exon_string, $align_coord, 10);

	# Move along if we have starting gaps
	if ($exon_seed_string =~ /^(-+)/){
	  $align_coord += length($1);
	  next LOST;
	}

	# Make sure we havent run off the end of the sequence
	throw("Have run off the end of the alignment without matching the gene seed.  Bugger.")
	  if ($align_coord >= $align_length);

	# Get on with extending the exon seed, without gaps.

	$exon_seed_string =~ s/-//g;

#print STDERR "    Align coord : " . $align_coord . "  Seed : $exon_seed_string\n";
	my $increment = 0;

      EXTENTION:
	while (length ($exon_seed_string) < 8){
	  $increment++;
	  $exon_seed_string = substr($exon_string, $align_coord, 10 + $increment);

	  if ($exon_seed_string =~ /-$/){
	    $increment++;
	    $exon_seed_string =~ s/-//g;
	    next EXTENTION
	  }

	  $exon_seed_string =~ s/-//g;

	  throw("Have run off the end of the alignment without matching the gene seed.  Bugger.")
	    if ($align_coord + 10 + $increment >= $align_length);
#print STDERR "      Extending exon seed : " . $exon_seed_string . "\n";
	}
	
	if ($gene_seed_string =~ /^$exon_seed_string/) {
#print STDERR "Sequence FOUND\n";
	  $lost = 0;
	  last LOST
	}
	$align_coord++;
      }

    }

    # Sanity check
#print STDERR "Gene neighbourhood : " . substr($gene_string, $gene_coord, 10) . " Exon neighbourhood : " . substr($exon_string, $align_coord, 10) . "\n";
    throw("Gene and aligned exon sequence have unexpectedly become unaligned.")
      unless substr($gene_string, $gene_coord, 1) eq substr($exon_string, $align_coord, 1);

    # Do our work while we are still aligned
#print STDERR "Trying to get some work done.\n";
    if (substr($informative_sites_string, $gene_coord, 1) eq '*'){
#print STDERR "  Found an informative site.\n";
      foreach my $evidence_id (keys %evidence_strings) {
#	next if (($evidence_id eq 'genomic_sequence') ||
#	  ($evidence_id eq 'exon_sequence'));
	$aligned_informative_site_strings{$evidence_id} .= 
	  substr($evidence_strings{$evidence_id}, $align_coord, 1)
      }
    } 
else {
#print STDERR "Not an informative site, not interested and moving right along.\n";
}
    # Do a quick check that the next base of the gene will 
    # align to the alignment.  If there is gap, try hopping
    # over it, before assuming that we are lost again.

    my $hop = 1;
    if (($gene_coord + 1 <= $gene_length)&&
       (substr($exon_string, $align_coord + $hop, 1) ne '-')) {
#print STDERR "Next base is not a gap, hence not trying to hop it.\n";
      # Gene coord is about to be incremented in the next iteration 
      # of this loop, so had better do this to $align_coord as well.
      $align_coord++;
      next
    } else {
#print STDERR "Next base is a gap.  Trying to hop.\n";
      while (substr($exon_string, $align_coord + $hop, 1) eq '-'){
	$hop++;
#print STDERR "     Little hop.\n";
	while (substr($exon_string, $align_coord + $hop, 10) eq '-' x 10) {
	  $hop += 10;
#print STDERR "      Big hop.\n";
	  while (substr($exon_string, $align_coord + $hop, 100) eq '-' x 100) {
	    $hop += 100;
#print STDERR "      Really big hop.\n";
	  }
	}
      }
      $align_coord += $hop;
#print STDERR "Align coord after hop is $align_coord.\n";
#print STDERR "Gene neighbourhood : " . substr($gene_string, $gene_coord + 1, 10) . " Exon neighbourhood : " . substr($exon_string, $align_coord, 10) . "\n";
      my $exon_neighbourhood = substr($exon_string, $align_coord, 10);
      $exon_neighbourhood =~ s/-//g;
      if ((! substr($gene_string, $gene_coord + 1, 10) =~ /^$exon_neighbourhood/) 
	  && ($gene_coord < (length($gene_string) - 1))) {
	warning("Erk.  Have managed to become lost.");
	$lost = 1;
      }
    }
  }

#print STDERR "DONE.\n";

  # Finally, quickly make an alignment object out of the
  # informative site strings.

  my $informative_sites_alignment = 
    Bio::EnsEMBL::Pipeline::Alignment->new(
      -name => $seq->display_id);

  foreach my $evidence_id (keys %aligned_informative_site_strings) {
    my $new_is_seq = 
      Bio::EnsEMBL::Pipeline::Alignment::AlignmentSeq->new(
        -name => $evidence_id,
        -seq  => $aligned_informative_site_strings{$evidence_id});

    $informative_sites_alignment->add_sequence($new_is_seq);
  }

  return $informative_sites_alignment;
}
