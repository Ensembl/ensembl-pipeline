package Bio::EnsEMBL::Pipeline::Runnable::ESTDescriminator;

use strict;
use Bio::EnsEMBL::Pipeline::Alignment::EvidenceAlignment;
use Bio::EnsEMBL::Pipeline::GeneDuplication::CodonBasedAlignment;
use Bio::Align::DNAStatistics;
use Bio::AlignIO;
use Bio::LocatableSeq;
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

  my %hash_of_ests_and_gene_mapping = 
    %{$self->_cluster_ests_with_genes};

  ### This is were I fell asleep


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

  $self->_calculate_gene_vs_est_distances;
  $self->_calculate_intergene_distances;

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
	 'distance' => $self->_pairwise_distance($gene_id, $est_id));
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
						   $genes_that_share_this_est->[$k]));
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
    #     another. ----------------------------------------THIS IS PRESENTLY NOT IMPLEMENTED.


    my $closest_gene_id;
    my @closely_related_genes;
    foreach my $comparison (@est_gene_distances){
      if ($comparison->{type} eq 'est_vs_gene') {
	$closest_gene_id = $comparison->{gene_ids}->[0];
	last
      } elsif ($comparison->{type} eq 'inter_gene') {
	push @closely_related_genes, $comparison;
      }
    }

    throw("Failed to find most closely related gene.")
      unless defined $closest_gene_id;

    if (scalar @closely_related_genes == 0){
      push @{$ests_clustered_to_genes{$est_id}}, $closest_gene_id;
      next
    }

    foreach my $gene_pair (@closely_related_genes){
      if ($gene_pair->{gene_ids}->[0] eq $closest_gene_id) {
	push @{$ests_clustered_to_genes{$est_id}}, $gene_pair->{gene_ids}->[1]
      } elsif ($gene_pair->{gene_ids}->[1] eq $closest_gene_id) {
	push @{$ests_clustered_to_genes{$est_id}}, $gene_pair->{gene_ids}->[0]
      }
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
    # THERE IS A SMARTER WAY TO DO THIS-------------------
    my $transcript = $self->_gene_lookup->{$gene_id}->get_all_Transcripts->[0];

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
    # This is due to a failing in Bio::EnsEMBL::Pipeline::Alignmenti::EvidenceAlignment.

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

    # Convert this alignment into a bioperl AlignI object.
    my $bp_align = Bio::SimpleAlign->new();

    my @sequence_order;
#    foreach my $align_seq (@{$ens_align->fetch_AlignmentSeqs}){
    foreach my $align_seq (@ens_align_seqs){

      # Skip the genomic sequence in our alignment - use the exon
      # sequence for calculating distance
      next if $align_seq->name eq 'genomic_sequence';

      my $locatable_seq = Bio::LocatableSeq->new(
                            -id       => $align_seq->name,
                            -seq      => $align_seq->seq,
                            -start    => 1,
                            -end      => length($align_seq->seq),
                            -alphabet => 'dna');

      $bp_align->add_seq($locatable_seq);

      push @sequence_order, $align_seq->name;
    }

    # Use bioperl object to calculate our pairwise sequence distances.
    my $dna_stat = Bio::Align::DNAStatistics->new;

    my $distances = $dna_stat->distance(-align  => $bp_align,
					-method => 'Kimura');

    # Store result of first column of result matrix (ie. gene vs ests)
print STDERR "Pairwise distances between gene sequence [$gene_id] and ESTs : \n";
    for (my $row = 1; $row < scalar @$distances; $row++) {
      next if $sequence_order[$row-1] eq 'exon_sequence';
      next unless defined $distances->[$row];
print STDERR "Pairwise distance : [$row] : " . $distances->[$row]->[1] . "\n";
      $self->_pairwise_distance($gene_id, 
				$sequence_order[$row-1], 
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

  # Make a bioperl align object and calculate pairwise distances

  my $bp_align = Bio::SimpleAlign->new();

  my @sequence_order;
  foreach my $align_seq (@$aligned_seqs){
    my $locatable_seq = Bio::LocatableSeq->new(
                          -id    => $align_seq->display_id,
                          -seq   => $align_seq->seq,
                          -start => 1,
                          -end   => $align_seq->length);
    $bp_align->add_seq($locatable_seq);
    push @sequence_order, $align_seq->display_id;
  }

  # Use bioperl object to calculate our pairwise sequence distances.

  my $dna_stat = Bio::Align::DNAStatistics->new;

  my $distances = $dna_stat->distance(-align  => $bp_align,
				      -method => 'Kimura');
  # Matrix of gene vs gene distances
print STDERR "Gene-gene pairwise distances : \n";
  for (my $row = 1; $row < scalar @$distances; $row++) {
    next unless defined $distances->[$row];
    for (my $column = $row + 1; $column < scalar @{$distances->[$row]}; $column++){
print STDERR "Pairwise distance : [$row] [$column] : " . $distances->[$row]->[$column] . "\n";
      $self->_pairwise_distance($sequence_order[$column - 1], 
				$sequence_order[$row - 1], 
				$distances->[$row]->[$column])
    }
  }

  return 1
}

sub _identify_shared_ESTs {
  my $self = shift;

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

  return \%genes_with_their_ests
}

sub _find_overlapping_ests {
  my ($self, $gene) = @_;

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
      throw("Somehow, found an exon that is zero bp long.")
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
  my ($self, $id1, $id2, $distance) = @_;

  if (defined $distance){
    $self->{_pairwise_distances}->{$id1}->{$id2} = $distance;
    return $distance
  }

  $distance = undef;

  if (defined $self->{_pairwise_distances}->{$id1}->{$id2}){
    $distance = $self->{_pairwise_distances}->{$id1}->{$id2};
  } elsif (defined $self->{_pairwise_distances}->{$id2}->{$id1}) {
    $distance = $self->{_pairwise_distances}->{$id2}->{$id1};
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

  my $shared_ests = undef;

  if (defined ($self->{_id_lookup}->{$gene_id})) {
    $shared_ests = $self->{_id_lookup}->{$gene_id}
  }

  return $shared_ests
}

sub _find_genes_by_est_id {
  my ($self, $est_id) = @_;

  my $genes = undef;

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

#     my $length = $transcript->spliced_seq =~ tr/[ATGCatgcNn]//;
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
