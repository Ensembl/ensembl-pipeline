package Bio::EnsEMBL::Pipeline::Runnable::ESTDescriminator;

use strict;
use Bio::EnsEMBL::Pipeline::Alignment::EvidenceAlignment;
use Bio::Align::DNAStatistics;
use Bio::AlignIO;
use Bio::EnsEMBL::Utils::Exception qw(throw warning info); 
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

sub new {
  my ($class, @args) = @_;

  my $self = bless {},$class;

  my ($genes,
      $est_db,
      $est_coverage_cutoff,
     ) = rearrange([qw(GENES
                       EST_DB
                       EST_COVERAGE_CUTOFF
                      )],@args);

  throw("No genes to analyse")
    unless $genes;

  $self->_genes($genes);

  throw("Constructor needs adaptors to the database which " .
	"contains your ests.")
    unless $est_db;

  $self->_est_db($est_db);

  if ($est_coverage_cutoff) {
    $self->_est_coverage_cutoff($est_coverage_cutoff)
  } else {
    $self->_est_coverage_cutoff(0.8)
  }

  return $self
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

sub _calculate_gene_vs_est_distances {
  my $self = shift;

  my $genes_and_shared_ests = $self->_identify_shared_ests;
  my %gene_vs_est_distance_matrix;

  foreach my $gene_id (keys %$genes_and_shared_ests) {
    # THERE IS A SMARTER WAY TO DO THIS-------------------
    my $transcript = $self->_gene_lookup->{$gene_id}->get_all_Transcripts->[0];

    my $evidence_align = 
      Bio::EnsEMBL::Pipeline::Alignment::EvidenceAlignment->new(
        -transcript          => $transcript,
        -dbadaptor           => $self->_gene_lookup->{$gene_id}->db,
        -seqfetcher          => $self->_seqfetcher,
        -padding             => 15,
        -supporting_features => $genes_and_shared_ests->{$gene_id});
        );

    my $ens_align = $evidence_align->retrieve_alignment( 
                   -type            => 'nucleotide',
                   -remove_introns  => 0,
                   -merge_sequences => 1);

    $align_filename = $self->_write_align($ens_align);  ################## Must implement

    my $dna_stat = Bio::Align::DNAStatistics->new;
    my $alignio  = Bio::AlignIO('-format' => 'fasta',
                                -file     => $align_filename);

    my $bp_align = $alignio->next_aln;

    my $kimura_distances = $stats->distance(-align  => $bp_align,
                                            -method => 'Kimura');


    foreach my $d ( @$kimura_distances )  {
      print "\t";
      foreach my $r ( @$d ) {
        print "$r\t";
      }
      print "\n";
    }
################# Fell asleep here.
  }


}

sub _identify_shared_ESTs {
  my $self = shift;

  # Retrieve all ESTs associated with our gene list and look 
  # for multiply-occurring EST ids.
  my %est_comap_counts;
  my %genes_with_their_ests;

  foreach my $gene (@{$self->_genes}){
    my $overlapping_ests = $self->_find_overlapping_ests($gene);

    my $genes_stable_id = $gene->stable_id;

    foreach my $est (@$overlapping_ests){
      my @est_supporting_features;

      foreach my $exon (@{$est->get_all_Exons}) {
        push @est_supporting_features, @{$exon->get_all_supporting_features} 
      }

      my $est_name = $est_supporting_features->[0]->hseqname;
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
    foreach my $gene (@$self->_genes) {
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
    $self->{_seqfetche} = shift;
    throw("Seqfetcher is not a Bio::EnsEMBL::Pipeline::SeqFetcher")
      unless defined $self->{_seqfetcher} &&
        $self->{_est_db}->isa("Bio::EnsEMBL::Pipeline::SeqFetcher");
  }

  return $self->{_seqfetcher}
}




1;
