package Bio::EnsEMBL::Pipeline::Runnable::ESTDescriminator;

use strict;
use Bio::EnsEMBL::Pipeline::Alignment::EvidenceAlignment;
use Bio::EnsEMBL::Pipeline::Alignment::InformativeSites;
use Bio::EnsEMBL::Pipeline::GeneDuplication::CodonBasedAlignment;
use Bio::EnsEMBL::Pipeline::Alignment;
use Bio::Seq;
use Bio::EnsEMBL::Utils::Exception qw(throw warning info); 
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

sub new {
  my ($class, @args) = @_;

  my $self = bless {},$class;

  my ($genes,
      $est_db,
      $est_coverage_cutoff,
      $seqfetcher,
     ) = rearrange([qw(GENES
                       EST_DB
                       EST_COVERAGE_CUTOFF
		       SEQUENCE_FETCHER
                      )],@args);

  throw("No genes to analyse")
    unless $genes;

  $self->_genes($genes);

  throw("Constructor needs adaptors to the database which " .
	"contains your ests.")
    unless (defined $est_db);

  $self->_est_db($est_db);

  throw("Cant proceed without a sequence fetcher.")
    unless (defined $seqfetcher);  

  $self->_seqfetcher($seqfetcher);

  if ($est_coverage_cutoff) {
    $self->_est_coverage_cutoff($est_coverage_cutoff)
  } else {
    $self->_est_coverage_cutoff(0.8)
  }

  return $self
}

sub run {
  my $self = shift;

  $self->_identify_shared_ESTs;

  my %hash_of_ests_and_gene_mapping = 
    %{$self->_cluster_ests_with_genes};

  return \%hash_of_ests_and_gene_mapping
}

sub print_shared_ESTs {
  my $self = shift;

  my %est_hash;

  foreach my $gene (@{$self->_genes}){
    my $overlapping_ests = $self->_find_overlapping_ests($gene);

    foreach my $est (@$overlapping_ests){
      $est_hash{$est->get_all_Exons->[0]->get_all_supporting_features->[0]->hseqname}++;
    }
  }

  my $shared_ests = 0;

  foreach my $est_name (keys %est_hash){
    if ($est_hash{$est_name} > 1) {
      print STDERR $est_name . "\tmapped to " . $est_hash{$est_name} . " genes.\n";
      $shared_ests++; 
    }
  }

  print STDERR "Had $shared_ests ESTs that mapped to more than one gene.\n";

  return 1
}

sub _cluster_ests_with_genes {
  my $self = shift;

  my %ests_clustered_to_genes;

  # Execute the subroutines that will generate the 
  # needed data.

#  $self->_calculate_gene_vs_est_distances;
#  $self->_calculate_intergene_distances;
  $self->_calculate_inf_site_gene_vs_est_distances;


  # Derive a list of unique EST ids

  my %all_shared_ests;
  foreach my $gene_id (keys %{$self->_gene_lookup}){
    my $est_ids = $self->_find_shared_ests_by_gene_id($gene_id);
    foreach my $est_id (@$est_ids){
      $all_shared_ests{$est_id}++;
    }
  }

  # Work through ESTs, grouping each with the most 
  # appropriate gene(s).

  foreach my $est_id (keys %all_shared_ests){

    # Make an array of simple hashes, containing
    # both EST vs gene and inter-gene distances
    # for this particular EST.

    my @est_gene_distances;

    my $genes_that_share_this_est = $self->_find_genes_by_est_id($est_id);

      # EST vs gene distances

    foreach my $gene_id (@$genes_that_share_this_est){

      my %pairwise_comparison = 
	('type'     => 'est_vs_gene',
	 'gene_ids'  => [$gene_id],
	 'distance' => $self->_pairwise_distance($gene_id, $est_id, undef, 1));
      push @est_gene_distances, \%pairwise_comparison;
    }

      # Inter-gene distances

    for (my $i = 0; $i < scalar @$genes_that_share_this_est; $i++) {
      for (my $k = $i + 1; $k < scalar @$genes_that_share_this_est; $k++) {

	my %pairwise_comparison = 
	  ('type'     => 'inter_gene',
	   'gene_ids' => [$genes_that_share_this_est->[$i], 
			  $genes_that_share_this_est->[$k]],
	   'distance' => $self->_pairwise_distance($genes_that_share_this_est->[$i],
						   $genes_that_share_this_est->[$k],
						   undef,
						   1));
	push @est_gene_distances, \%pairwise_comparison;
      }
    }


    # Sort the array (of hashes) by distance in ascending order.

    @est_gene_distances = 
      sort {$a->{distance} <=> $b->{distance}} @est_gene_distances;

    # Group the EST to specific genes according
    # to following rules:
    # * EST maps to one gene if EST vs gene distance
    #     is smaller than the smallest inter-gene 
    #     distance.
    # * Where genes are more closely related than the
    #     EST vs gene distances, this EST must be mapped
    #     to either both genes or neither gene.
    # * The inherant error attached to the distance
    #     must be considered when saying an EST is more
    #     closely related to one gene rather than
    #     another.  THIS IS PRESENTLY NOT IMPLEMENTED.


    my %closest_gene_ids;
    my $last_dist;
    foreach my $comparison (@est_gene_distances){
      if ($comparison->{type} eq 'est_vs_gene') {
	last if (defined $last_dist && 
		 $comparison->{distance} > $last_dist);

	$last_dist = $comparison->{distance};
	$closest_gene_ids{$comparison->{gene_ids}->[0]}++;

      } elsif ($comparison->{type} eq 'inter_gene') {
	last if (defined $last_dist && 
		 $comparison->{distance} > $last_dist);

	$closest_gene_ids{$comparison->{gene_ids}->[0]}++;
	$closest_gene_ids{$comparison->{gene_ids}->[1]}++;
      }
    }

    throw("Failed to find most closely related gene.")
      unless keys %closest_gene_ids;

    foreach my $close_gene_id (keys %closest_gene_ids){
      push @{$ests_clustered_to_genes{$est_id}}, $close_gene_id;
    }

    throw("Erroneously didn't map EST to any gene.")
      unless defined $ests_clustered_to_genes{$est_id};
  }

  return \%ests_clustered_to_genes
}

sub _calculate_gene_vs_est_distances {
  my $self = shift;

  my $genes_and_shared_ests = $self->_identify_shared_ESTs;
  my %gene_vs_est_distance_matrix;

  foreach my $gene_id (keys %$genes_and_shared_ests) {
    my $align = $self->_evidence_align($gene_id);

    my ($distances, $sequence_order) = $self->_compute_distances($align);

    # Store result of first column of result matrix (ie. gene vs ests)

    for (my $row = 1; $row < scalar @$distances; $row++) {
      next if $sequence_order->[$row-1] eq 'exon_sequence';
      next unless defined $distances->[$row];

      $self->_pairwise_distance($gene_id, 
				$sequence_order->[$row-1], 
				$distances->[$row]->[1])
    }
  }

  return 1
}

sub _calculate_intergene_distances {
  my $self = shift;

  # Derive the longest transcript from each gene.

  my %gene_transcript_lookup;
  my @transcript_seqs;

  foreach my $gene (@{$self->_genes}){
    my $longest_transcript = $self->_longest_transcript($gene);

    unless ($longest_transcript){
      throw("Gene has no translatable transcript [" . $gene->stable_id . "]")
    }

    my $seq = 
      Bio::Seq->new(-display_id => $gene->stable_id,
		    -seq        => $longest_transcript->translateable_seq);

    push @transcript_seqs, $seq;
    $gene_transcript_lookup{$longest_transcript->stable_id} = $gene->stable_id;
    $gene_transcript_lookup{$gene->stable_id} = $longest_transcript->stable_id;
  }

  # Obtain sequences for transcripts and perform an alignment

  my $cba = 
    Bio::EnsEMBL::Pipeline::GeneDuplication::CodonBasedAlignment->new(
      -genetic_code => 1); # Hard-code for universal genetic code.

  $cba->sequences(\@transcript_seqs);

  my $aligned_seqs = $cba->run_alignment;

  # Convert these aligned sequences to a proper Pipeline::Alignment object
  my $align = Bio::EnsEMBL::Pipeline::Alignment->new();

  foreach my $bioseq (@$aligned_seqs){
    my $align_seq = 
      Bio::EnsEMBL::Pipeline::Alignment::AlignmentSeq->new(
	-name => $bioseq->display_id,
	-seq  => $bioseq->seq);
    $align->add_sequence;
  }

  my ($distances, $sequence_order) = $self->_compute_distances($align);

  # Matrix of gene vs gene distances
  for (my $row = 1; $row < scalar @$distances; $row++) {
    next unless defined $distances->[$row];
    for (my $column = $row + 1; $column < scalar @{$distances->[$row]}; $column++){
      $self->_pairwise_distance($sequence_order->[$column - 1], 
				$sequence_order->[$row - 1], 
				$distances->[$row]->[$column])
    }
  }

  return 1
}

sub _calculate_inf_site_gene_vs_est_distances {
  my $self = shift;

  $self->_identify_shared_ESTs; # Just in case this has not been done.

  # Determine gene-pairs that share ESTs.  This is a bit
  # messy - genes are joined through sifting their EST
  # matches.

  my %gene_vs_gene;
  my @gene_pair;
  my %seen_est;

  foreach my $gene (@{$self->_genes}){
    my $gene_id = $gene->stable_id;

    foreach my $est_id (@{$self->_find_shared_ests_by_gene_id($gene_id)}){
      next if $seen_est{$est_id};

      foreach my $other_gene (@{$self->_find_genes_by_est_id($est_id)}){
	next if $other_gene eq $gene_id;
	$gene_vs_gene{$gene_id}->{$other_gene}++;
	$gene_vs_gene{$other_gene}->{$gene_id}++;
	if (($gene_vs_gene{$gene_id}->{$other_gene} + 
	     $gene_vs_gene{$other_gene}->{$gene_id}) == 2){
	  push @gene_pair, [$gene_id, $other_gene];
	}
      }
    }
  }

  # Foreach pair, find informative sites in gene and ESTs
  # and calculate the informative site distances.

  my %seen_genes;

  foreach my $gene_pair (@gene_pair){

    my @seqs = 
      (
      Bio::Seq->new(
       -display_id => $gene_pair->[0],
       -seq        => $self->_transcript_by_gene_lookup($gene_pair->[0])->translateable_seq)
      ,
      Bio::Seq->new(
       -display_id => $gene_pair->[1],
       -seq        => $self->_transcript_by_gene_lookup($gene_pair->[1])->translateable_seq)
      );

    my @aligns = ($self->_evidence_alignment($gene_pair->[0]),
		  $self->_evidence_alignment($gene_pair->[1]));

    my $inf_sites = 
      Bio::EnsEMBL::Pipeline::Alignment::InformativeSites->new(
        -translatable_seqs  => \@seqs,
        -alignments         => \@aligns);


    my $hr_inf_site_aligns = $inf_sites->informative_sites;

    # When calculating gene-gene distances from the informative
    # sites for a gene pair, the distance by definition must be 1,
    # meaning complete divergence.  Hence, we may as well just set
    # this value manually.

    $self->_pairwise_distance($gene_pair->[0], $gene_pair->[1], 1, 1);

    # Calculate and store gene-est distances from informative 
    # site alignments.

    foreach my $gene_id (@$gene_pair){

      next if $seen_genes{$gene_id};
      $seen_genes{$gene_id}++;

      my ($distances, $sequence_order) = 
	$self->_compute_distances($hr_inf_site_aligns->{$gene_id});

      for (my $row = 1; $row < scalar @$distances; $row++) {
	next unless defined $distances->[$row];
	for (my $column = 1; $column < scalar @{$distances->[$row]}; $column++){
	  my $est_id;

	  if ($sequence_order->[$row - 1] eq 'exon_sequence'){
	    $est_id = $sequence_order->[$column - 1];
	  } elsif ($sequence_order->[$column - 1] eq 'exon_sequence'){
	    $est_id = $sequence_order->[$row - 1];
	  } else {
	    next
	  }

	  $self->_pairwise_distance($gene_id, 
				    $est_id, 
				    $distances->[$row]->[$column],
				    1);
	}
      }

    }
  }

  return 1;
}

sub _compute_distances {
  my ($self, $align) = @_;

  my @sequence_order;
  my %seqs;

  foreach my $align_seq (@{$align->fetch_AlignmentSeqs}){

    # Skip the genomic sequence in our alignment - use the exon
    # sequence for calculating distance
    next if $align_seq->name eq 'genomic_sequence';

    $seqs{$align_seq->name} = $align_seq->seq;

    push @sequence_order, $align_seq->name;
  }

  my @distances;

  for (my $i = 0; $i < (scalar @sequence_order) - 1; $i++) {
    my @seq1 = split //, $seqs{$sequence_order[$i]};

    warning("Informative site alignments are length zero.  Try something different.")
      if scalar @seq1 < 1;

    for (my $j = $i + 1; $j < scalar @sequence_order; $j++) {
      my @seq2 = split //, $seqs{$sequence_order[$j]};

      my $matching_bases = 0;

      throw("Informative site sequences not the same length")
	unless scalar @seq1 == @seq2;

      for (my $k = 0; $k < scalar @seq1; $k++) {
	$matching_bases++
	  if $seq1[$k] eq $seq2[$k];
      }

      my $distance = 1 - ($matching_bases/(scalar @seq1));

      $distances[$j+1][$i+1] = sprintf("%4.2f", $distance);
    }
  }

  return (\@distances, \@sequence_order);
}

#use Bio::Align::DNAStatistics;
#use Bio::AlignIO;
#use Bio::LocatableSeq;
#
#sub _compute_distances {
#  my ($self, $align) = @_;
#
#  # Convert this alignment into a bioperl AlignI object.
#  my $bp_align = Bio::SimpleAlign->new();
#
#  my @sequence_order;
#  foreach my $align_seq (@{$align->fetch_AlignmentSeqs}){
#
#
#    # Skip the genomic sequence in our alignment - use the exon
#    # sequence for calculating distance
#    next if $align_seq->name eq 'genomic_sequence';
#
#    my $locatable_seq = Bio::LocatableSeq->new(
#			    -id       => $align_seq->name,
#                            -seq      => $align_seq->seq,
#                            -start    => 1,
#                            -end      => length($align_seq->seq),
#                            -alphabet => 'dna');
#
#    $bp_align->add_seq($locatable_seq);
#
#    push @sequence_order, $align_seq->name;
#  }
#
#  # Use bioperl object to calculate our pairwise sequence distances.
#  my $dna_stat = Bio::Align::DNAStatistics->new;
#
#  my $distances = $dna_stat->distance(-align  => $bp_align,
#                                      -method => 'Kimura');
#
#  return ($distances, \@sequence_order);
#}

sub _evidence_alignment {
  my ($self, $gene_id) = @_;

  # Return cached alignment if previously computed
  # for this gene id

  return $self->{_evidence_alignments}->{$gene_id}
    if defined $self->{_evidence_alignments}->{$gene_id};


  # If still here, proceed with retrieving evidence alignment.

  my $genes_and_shared_ests = $self->_identify_shared_ESTs;

  my $transcript = $self->_transcript_by_gene_lookup($gene_id);

  # Gather all EST features
  my @est_supporting_features;
  foreach my $est_id (keys %{$genes_and_shared_ests->{$gene_id}}){
    push @est_supporting_features, 
      @{$genes_and_shared_ests->{$gene_id}->{$est_id}};    
  }

  # Obtain an alignment of our transcript with the shared ESTs
  my $evidence_align = 
    Bio::EnsEMBL::Pipeline::Alignment::EvidenceAlignment->new(
        -transcript          => $transcript,
        -dbadaptor           => $self->_gene_lookup->{$gene_id}->adaptor->db,
        -seqfetcher          => $self->_seqfetcher,
        -padding             => 15,
        -supporting_features => \@est_supporting_features 
        );

  my $ens_align = $evidence_align->retrieve_alignment( 
                   -type            => 'nucleotide',
                   -remove_introns  => 0,
                   -merge_sequences => 1);

  # A horrible and hopefully temporary hack to 
  # ensure all sequences are of the same length.
  # This is due to a failing in Bio::EnsEMBL::Pipeline::Alignment::EvidenceAlignment.

  my @ens_align_seqs = @{$ens_align->fetch_AlignmentSeqs};
  my $longest_seq = 0;
  foreach my $ens_align_seq (@ens_align_seqs) {
    my $length = length($ens_align_seq->seq);
    if ($length > $longest_seq){
      $longest_seq = length($ens_align_seq->seq)
    }
  }
  for (my $i = 0; $i < scalar @ens_align_seqs; $i++) {
    my $length = length($ens_align_seqs[$i]->seq);
    if ($length < $longest_seq) {
      my $shortfall = $longest_seq - $length;
      my $tack_on = '-' x $shortfall;
      my $augmented_seq = $ens_align_seqs[$i]->seq . $tack_on;
      $ens_align_seqs[$i]->seq($augmented_seq); 
    }
  }

  my $new_align = 
    Bio::EnsEMBL::Pipeline::Alignment->new(
      -seqs => \@ens_align_seqs,
      -name => $gene_id);

  # Add this alignment to the pre-computed cache
  $self->{_evidence_alignments}->{$gene_id} = $new_align;

  return $new_align;
}

sub _identify_shared_ESTs {
  my $self = shift;

  # This is a fairly database heavy routine, hence
  # if it has been run previously return the cached results
  # rather than rerunning.

  return $self->{_shared_ests_cache}
    if defined $self->{_shared_ests_cache};

  # Retrieve all ESTs associated with our gene list and look 
  # for multiply-occurring EST ids.
  my %est_comap_counts;
  my %genes_with_their_ests;

  foreach my $gene (@{$self->_genes}){
    my $overlapping_ests = $self->_find_overlapping_ests($gene);

    my $gene_stable_id = $gene->stable_id;

    foreach my $est (@$overlapping_ests){
      my @est_supporting_features;

      foreach my $exon (@{$est->get_all_Exons}) {
        push @est_supporting_features, @{$exon->get_all_supporting_features} 
      }

      my $est_name = $est_supporting_features[0]->hseqname;
      $genes_with_their_ests{$gene_stable_id}->{$est_name} = \@est_supporting_features;
      $est_comap_counts{$est_name}++;
    }
  }

  # Purge our %genes_with_their_ests hash of singly-occuring ESTs

  foreach my $est_name (keys %est_comap_counts){
    if ($est_comap_counts{$est_name} == 1){
      foreach my $gene_name (keys %genes_with_their_ests){
        delete $genes_with_their_ests{$gene_name}->{$est_name};
      }
    }
  }

  $self->_build_id_lookup(\%genes_with_their_ests);
  $self->{_shared_ests_cache} = \%genes_with_their_ests;

  return \%genes_with_their_ests
}

sub _find_overlapping_ests {
  my ($self, $gene) = @_;
print STDERR "Looking at gene : " . $gene->stable_id . "\n";
  # First, move gene into chromosome coordinates.

  my $gene_slice_adaptor  = $gene->slice->adaptor;
  my $gene_chr_slice = $gene_slice_adaptor->fetch_by_region(
                           'chromosome',
                           $gene->slice->seq_region_name,
                           1,
                           $gene->slice->seq_region_length);
  $gene = $gene->transfer($gene_chr_slice);

  # Get on with finding EST genes that overlap this gene.

  my $est_sa = $self->_est_db->get_SliceAdaptor;
  my $est_ga = $self->_est_db->get_GeneAdaptor;

  my $est_slice =
    $est_sa->fetch_by_region('chromosome',
                             $gene->slice->seq_region_name,
                             $gene->start,
                             $gene->end);

  my $mapped_ests = $est_ga->fetch_all_by_Slice($est_slice);

  # Transform EST genes to chromosomal coordinates.

  my $est_chr_slice = $est_sa->fetch_by_region('chromosome',
                                               $gene->slice->seq_region_name, 
                                               1, 
                                               $gene->slice->seq_region_length);

  for (my $i = 0; $i < scalar @$mapped_ests; $i++) {
    $mapped_ests->[$i] = $mapped_ests->[$i]->transfer($est_chr_slice);
  }

  # Check mapped ESTs for good coverage against the exons of our gene.

    # Build a hash showing covered bases (an array with elements in
    # genomic coordinates is not a good thing).

  my %covered_bases;

  foreach my $exon (@{$gene->get_all_Exons}){
    my ($start, $end) = 
      sort {$a <=> $b} ($exon->start, $exon->end);

    for (my $base = $start; $base <= $end; $base++){
      $covered_bases{$base}++;
    }
  }

  my @overlapping_ests;
  my %seen_evidence;

  foreach my $mapped_est (@$mapped_ests){
    my $mapped_est_hseqname;

    my $covered_bases = 0;
    my $total_bases = 0;

    foreach my $exon (@{$mapped_est->get_all_Exons}){
      my ($start, $end) = 
          sort {$a <=> $b} ($exon->start, $exon->end);

      my $sf = $exon->get_all_supporting_features;
      $mapped_est_hseqname = $sf->[0]->hseqname;

      for (my $base = $start; $base <= $end; $base++){
        $total_bases++;
        $covered_bases++ if ($covered_bases{$base});
      }

    }

    unless ($total_bases) {
      throw("Somehow, found an gene that is zero bp long.")
    }

    my $coverage = $covered_bases/$total_bases;
    if (($coverage >= $self->_est_coverage_cutoff)&&
      (! defined($seen_evidence{$mapped_est_hseqname}))){
        push (@overlapping_ests, $mapped_est);
      }

    $seen_evidence{$mapped_est_hseqname}++;
  }

  print STDERR "Have " . scalar @overlapping_ests . " ESTs that map well to " . 
    $gene->stable_id . ".\n";

  return \@overlapping_ests
}

sub _pairwise_distance {
  my ($self, $id1, $id2, $distance, $inf) = @_;

  # If $inf is set, toggle to informative site distances,
  # otherwise just use the normal distances.

  my $dist_set = $inf ? 'inf' : 'norm';

  if (defined $distance){
    print "Storing [$dist_set][$id1][$id2][$distance]\n";
    warning("Trying to store an already stored value [" . $dist_set . 
	    "][" . $id1 . "][" . $id2 . "][" . $distance . "]")
      if defined $self->{_pairwise_distances}->{$dist_set}->{$id1}->{$id2};

    $self->{_pairwise_distances}->{$dist_set}->{$id1}->{$id2} = $distance;
    return $distance
  }

  $distance = undef;

  if (defined $self->{_pairwise_distances}->{$dist_set}->{$id1}->{$id2}){
    $distance = $self->{_pairwise_distances}->{$dist_set}->{$id1}->{$id2};
  } elsif (defined $self->{_pairwise_distances}->{$dist_set}->{$id2}->{$id1}) {
    $distance = $self->{_pairwise_distances}->{$dist_set}->{$id2}->{$id1};
  }

  return $distance
}

sub _build_id_lookup {
  my ($self, $genes_and_their_shared_ests) = @_;

  foreach my $gene_id (keys %$genes_and_their_shared_ests) {
    foreach my $est_id (keys %{$genes_and_their_shared_ests->{$gene_id}}){
      push @{$self->{_id_lookup}->{$gene_id}}, $est_id;
      push @{$self->{_id_lookup}->{$est_id}}, $gene_id;
    }
  }

  return 1
}

sub _find_shared_ests_by_gene_id {
  my ($self, $gene_id) = @_;

  my $shared_ests = [];

  if (defined ($self->{_id_lookup}->{$gene_id})) {
    $shared_ests = $self->{_id_lookup}->{$gene_id}
  }

  return $shared_ests
}

sub _find_genes_by_est_id {
  my ($self, $est_id) = @_;

  my $genes = [];

  if (defined ($self->{_id_lookup}->{$est_id})) {
    $genes = $self->{_id_lookup}->{$est_id}
  }

  return $genes
}

sub _longest_transcript {
  my ($self, $gene) = @_;

  my $longest = 0;
  my $longest_transcript;

  my $translatable = 0;

  foreach my $transcript (@{$gene->get_all_Transcripts}){

    eval {
      $transcript->translate;
    };

    next if $@;

    $translatable = 1;

     my $length = $transcript->translateable_seq =~ tr/[ATGCatgcNn]//;

     if ($length > $longest){
       $longest = $length;
       $longest_transcript = $transcript
     }
  }


  unless ($longest_transcript){
    warning ("Gene [" . $gene->stable_id . 
	 "] does not have a transcript with a translation.");
    $longest_transcript = 0;
  }

  return $longest_transcript
}

sub _genes {
  my $self = shift;

  if (@_) {
    $self->{_genes} = shift;

    throw("Gene is not a Bio::EnsEMBL::Gene, it is a [" . $self->{_genes}->[0] . "]")
      unless defined $self->{_genes}->[0] &&
        $self->{_genes}->[0]->isa("Bio::EnsEMBL::Gene");

    foreach my $gene (@{$self->{_genes}}){
      throw("Gene is not a Bio::EnsEMBL::Gene, it is a [$gene]") 
       unless $gene->isa("Bio::EnsEMBL::Gene")
    }
  }

  return $self->{_genes}
}

sub _gene_lookup {
  my $self = shift;

  unless ($self->{_gene_lookup}) {
    foreach my $gene (@{$self->_genes}) {
      $self->{_gene_lookup}->{$gene->stable_id} = $gene;
    }
  }

  return $self->{_gene_lookup}
}

sub _transcript_by_gene_lookup {
  my ($self, $gene_id) = @_;

  unless (defined $self->{_transcript_lookup}){
    foreach my $gene (@{$self->_genes}) {
      $self->{_transcript_lookup}->{$gene->stable_id}
	= $self->_longest_transcript($gene);
    }
  }

  return $self->{_transcript_lookup}->{$gene_id};
}


sub _est_db {
  my $self = shift;

  if (@_) {
    $self->{_est_db} = shift;
    throw("Gene db is not a Bio::EnsEMBL::DBSQL::DBAdaptor")
      unless defined $self->{_est_db} && 
        $self->{_est_db}->isa("Bio::EnsEMBL::DBSQL::DBAdaptor");
  }

  return $self->{_est_db}
}

sub _est_coverage_cutoff {
  my $self = shift;

  if(@_){
    $self->{_est_coverage_cutoff} = shift;
    if ($self->{_est_coverage_cutoff} > 1){
      $self->{_est_coverage_cutoff} = $self->{_est_coverage_cutoff}/100;
    }
  }

  return $self->{_est_coverage_cutoff}
}

sub _seqfetcher {
  my $self = shift;

  if (@_) {
    $self->{_seqfetcher} = shift;
    throw("Seqfetcher is not a Bio::EnsEMBL::Pipeline::SeqFetcher")
      unless (defined $self->{_seqfetcher} &&
        ($self->{_seqfetcher}->isa("Bio::EnsEMBL::Pipeline::SeqFetcher")||
         $self->{_seqfetcher}->can("get_Seq_by_acc"))); 
  }

  return $self->{_seqfetcher}
}

1;
