#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::BlastMiniBuilder;

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 CONTACT

Dan Andrews <dta@sanger.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::Runnable::BlastMiniBuilder;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Pipeline::Runnable::MultiMiniGenewise;
use Bio::EnsEMBL::Pipeline::Runnable::Blast;
use Bio::EnsEMBL::Pipeline::Runnable::BlastDB;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::PrimarySeqI;
use Bio::SeqIO;
use Bio::DB::RandomAccessI;
use Bio::EnsEMBL::Pipeline::GeneConf qw (
					 GB_INPUTID_REGEX
					);
use Bio::EnsEMBL::Pipeline::Config::Blast;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);


=head2 build_runnables

  Args [1]   : @features - list of FeaturePair objects 
               constructed from blast align features.  Within
               each FeaturePair, feature1 and feature2 should
               correspond to the genomic sequence and the
               protein/est/cdna sequence.
  Example    : $self->build_runnables(@features);
  Description: Builds the runnable objects for conducting
               BlastMini(Genewise/Est2Genome) runs.  In 
               the simplest case this
               method creates a single runnable object.  The
               method itself was written to cope with
               situations where the presence of repeated genes
               means that more than one whole gene of high 
               homology is present in a minigenomic sequence.
               This module separates possible repeated genes
               into clusters of exons such that multiple 
               minigenomic sequence fragments can be
               generated for passing to genewise/est2genome
               - when this happens this method returns
               multiple runnables that examine smaller
               portions of genomic DNA.
  Returntype : A list of Bio::EnsEMBL::Pipeline::Runnable
               - the exact identity of the runnable sub-class
               is determined by the calling class (as the 
               actual final runnable construction is done in 
               the inheriting class.)
  Exceptions : None.
  Caller     : $self->run

=cut

sub build_runnables {
  my ($self, @features) = @_;

  my @runnables;

  my %unfiltered_partitioned_features;

  foreach my $raw_feat (@features){
      push (@{$unfiltered_partitioned_features{$raw_feat->hseqname}}, $raw_feat);
  }

  # We partition the features by sequence id such that we can
  # iterate for each id.  Most ids  will fall through
  # to the simplest case where we make a single runnable with one 
  # unique gene per minigenomic sequence fragment.  In cases
  # where multiple similar genes are likely to be present in
  #  the minigenomic sequence we switch to the
  # special case algorithm.  When this happens we cluster the 
  # likely exons into groups of similar exons, which we then
  # use make a set of distinct genes.  For each likely gene
  # a minigenomic sequence is made that tightly spans the 
  # likely location of the gene - and most probably excludes
  # any similar genes nearby.  With these smaller minigenomic
  # sequence fragments a runnable is made using ALL of the
  # original blast-derived feature pairs.

  my %partitioned_features;

  foreach my $raw_feat (@features){

    # Imporantly, the features are filtered for low identity matches.
    # This is just for the clustering process.  The full set of features
    # are still passed to the runnable.
    if($raw_feat->percent_id >= 80){
      push (@{$partitioned_features{$raw_feat->hseqname}}, $raw_feat);
    }
  }

  my %runnable_featids;

  foreach my $seqname (keys %partitioned_features){

    # Preliminary check to see whether our filtered features
    # give enough coverage to make it worth continuing.

    my $coverage = $self->check_coverage($seqname, $partitioned_features{$seqname});

    # If not enough coverage, go home early.
    unless ($coverage > 0.75) {
      my $runnable = $self->make_object($self->genomic_sequence, $unfiltered_partitioned_features{$seqname});
      push (@runnables, $runnable);

      $runnable_featids{$seqname}++;
      next
    }

    # Now, we attempt to cluster our features into groups of
    # similar exons.

     my $clustered_features = $self->cluster_features(\@{$partitioned_features{$seqname}});

    # We have to do a little test to see if the clusters
    # represent the bulk of the features present in the 
    # genomic region.  This test helps avoid the situation
    # where a non-repeated gene has exons that are very 
    # similar to one another.  If the flag $clusters_seem_real
    # is not set to 1, the the analysis will default to the  
    # non-repeated gene way of doing things.

    my $clusters_seem_real = 0;

    if (@$clustered_features) {
      my $number_of_clustered_features = 0;
      foreach my $cluster (@$clustered_features) {
	foreach my $clust_feat (@$cluster){
	  $number_of_clustered_features++;
	}
      }

      if ($number_of_clustered_features  >= (0.75 * (scalar @{$partitioned_features{$seqname}}))) {
	$clusters_seem_real = 1;
      }
    }

    # If we have clusters in such numbers that it looks like
    # there might be repeated genes, we break the mini genomic
    # sequence into fragments and construct multiple runnables.
    # Otherwise, we default to the normal way of building our
    # runnable.
    
    if ((@$clustered_features)&&($clusters_seem_real > 0)) {

      print STDERR "Minigenomic sequence could contain a number "
	. "of highly similar genes.  Fragmenting the minigenomic "
	  . "sequence to try to resolve these genes individually.\n";

      my $gene_clusters = $self->form_gene_clusters($clustered_features);

      # Determine the extreme start and ends of each gene
      # cluster.  Decide on positions to split the genomic
      # sequence for genomewise - mid-way between 
      # ends and starts of neighboring clusters.  The
      # fragments from clusters that don't have the maximum
      # number of members are a problem - these presumably 
      # will fall at the ends of the genomic fragment and 
      # should just be discarded.

    GENE:      
      foreach my $gene_cluster (@$gene_clusters){

	my @sorted_gene_cluster = sort {$a->{_gsf_start} <=> $b->{_gsf_start};} @$gene_cluster;
       	my $cluster_start = $sorted_gene_cluster[0]->{_gsf_start} - 1000;
	my $cluster_end = $sorted_gene_cluster[-1]->{_gsf_end} + 1000;
        unless ($cluster_start > 0) {
	  if ($sorted_gene_cluster[0]->{_gsf_end} > 0) {
	    $cluster_start = 1;
	  } else {
	    next GENE;
	  }
	}
        unless ($cluster_end <= $self->genomic_sequence->length) {
	  if ($sorted_gene_cluster[0]->{_gsf_start} <= $self->genomic_sequence->length) {
	    $cluster_end = $self->genomic_sequence->length;
	  } else {
	    next GENE;
	  }
	}        

	# Here we create our genomic fragment for passing to 
	# our runnable.  This is the fragment of genomic sequence
	# that we have calculated should only contain the 
	# exons of one single gene.  It is padded with Ns to
        # maintain the coordinate system of the features that
	# may be returned.  Another way of looking at this is
	# that we make a genomic sequence where the regions 
	# flanking our gene are masked.
	my $string_seq = ('N' x ($cluster_start - 1)) . 
	         $self->genomic_sequence->subseq($cluster_start, $cluster_end)
		 . ('N' x ($self->genomic_sequence->length - ($cluster_end + 1)));

	my $genomic_subseq = Bio::Seq->new(-seq => $string_seq,
					   -id  => $self->genomic_sequence->id );

	my $runnable = $self->make_object($genomic_subseq, $unfiltered_partitioned_features{$seqname});
	push (@runnables, $runnable);

	$runnable_featids{$seqname}++;
      }
    
    }else{
      # This is what we do when we dont have multiple genes
      # in our genomic fragment.  This is "Normal mode".

      my $runnable = $self->make_object($self->genomic_sequence, $unfiltered_partitioned_features{$seqname});
      push (@runnables, $runnable);

      $runnable_featids{$seqname}++;
    }

  }

  # Gah, what if we havent constructed a runnable for each
  # feature id by this late stage?
  # It means that the homology of the blast features for
  # an id was so low that we threw all the features 
  # away.  Had better just make a default runnable...

  my %feature_ids;
  foreach my $feature (@features){
    $feature_ids{$feature->hseqname}++;
  }

  foreach my $feature_id (keys %feature_ids){
    unless ($runnable_featids{$feature_id} > 0){
      my $runnable = $self->make_object($self->genomic_sequence, $unfiltered_partitioned_features{$feature_id});
      push (@runnables, $runnable);
    }
  }

  return \@runnables;
}

=head2 cluster_features

  Args [1]   : $features - reference to a list of 
               FeaturePairs derived from the align feature
               objects resulting from a blast run.
               Ideally, these features should have been
               filtered to remove garbage - say, filter by
               percent_id > 80.
  Example    : $self->cluster_features($features);
  Description: Takes a list of features from a blast run
               and attempts to cluster these into groups 
               of high homology.  If the blast features
               result from matches to a region with several
               highly similar copies of the same gene this
               clustering should produce a series of exon
               clusters.  If this is not the case only a
               couple of grubby clusters will be produced.
  Returntype : A reference to an array of clusters.  Each
               cluster in the array is an array of similar
               features.
  Exceptions : none
  Caller     : self->build_runnables

=cut

sub cluster_features {

  my ($self, $features) = @_;

  # Sort the list of features according to their percent id.
  my @sorted_features = sort { $b->percent_id <=> $a->percent_id;} @$features;

  # With the sorted features, we iteratively check each     
  # highest scoring hit against each lesser scoring
  # sequence.  Features that correspond to the same part of 
  # the hit sequence (ie, are overlapping) are clustered
  # in an array of 'exons' (and removed from the 
  # sorted_features array). 
  
  my @all_feature_clusters;
  
  while (@sorted_features) {

    my $top_feature = shift @sorted_features;
    
    my @similar_features;
    push (@similar_features, $top_feature);

    my @subtracted_sorted_features;
    foreach my $lesser_feature (@sorted_features) {

      if ($self->check_overlap($top_feature,$lesser_feature)){
	push (@similar_features, $lesser_feature);

      } else {
	push (@subtracted_sorted_features, $lesser_feature);
      }
    }
    @sorted_features = @subtracted_sorted_features;
    unless ((scalar @similar_features) == 1){
      push (@all_feature_clusters, \@similar_features);
    }
  }

  return \@all_feature_clusters;
}


=head2 form_gene_clusters

  Args [1]   : Array of exon clusters generated by
               $self->cluster_features (an array of arrays,
               each sub-array containing similar blast 
               features)
  Example    : $self->for_gene_clusters($exon_clusters);
  Description: Using an array of exons clusters (similar
               blast features) to build genes with a feature
               from each cluster.  Genes are formed from
               exons that lie closest together.
  Returntype : A reference to an array of arrays, each sub-
               array containing features that should 
               correspond to exons. 
  Exceptions : none
  Caller     : $self->build_runnables

=cut



sub form_gene_clusters {

  my ($self, $exon_clusters) = @_;

  # Choose a cluster with the maximum number of members.
  # Bump the other clusters onto a new array.


  my @sorted_by_start = sort {$a->[0]->{_hstart} <=> $b->[0]->{_hstart}} @$exon_clusters;  

  my $big_cluster = shift @sorted_by_start;
  my @other_clusters = @sorted_by_start;

  my @final_gene_clusters;      
  while (@$big_cluster){

    my $seed_exon = shift @$big_cluster;
	
    my @gene_cluster;
    push (@gene_cluster, $seed_exon);
    
    foreach my $other_cluster (@other_clusters) {
      if (@$other_cluster){
	
	my @subtracted_other_clusters;
	my $closest = 1000000;
	my $closest_feature;
	foreach my $candidate_exon (@$other_cluster){
	  
	  my $distance = abs($candidate_exon->{_gsf_start} 
			     - $seed_exon->{_gsf_start});
	  
	  if(($distance < $closest)&&($candidate_exon->{_hstart} > $seed_exon->{_hend})){
	    if (defined $closest_feature){
	      push (@subtracted_other_clusters, $closest_feature);
	    }
	    $closest = $distance;
	    $closest_feature = $candidate_exon;
	  } else {
	    push (@subtracted_other_clusters, $candidate_exon);
	  }
	}
	@$other_cluster = @subtracted_other_clusters;

	push (@gene_cluster, $closest_feature) if (defined $closest_feature);
	undef $closest_feature;

      }
    }
    push (@final_gene_clusters, \@gene_cluster);
  }	

  return \@final_gene_clusters;
}

sub make_object {
  my ($self) = @_;

  $self->throw("make_object method must be implemented in the inheriting class.");

  # Sample implementation for Runnable::MultiMiniGenewise:
  #
  #  my ($self, $miniseq, $features) = @_;
  #
  #  my $mg      = new Bio::EnsEMBL::Pipeline::Runnable::MultiMiniGenewise(
  #	       		       '-genomic'    => $miniseq,
  #	       		       '-features'   => \@newf,
  #	       		       '-seqfetcher' => $self->seqfetcher,
  #	       		       '-endbias'    => $self->endbias
  #	       		      );
  #
  # return $mg;

}


=head2 check_overlap

  Args [1]   : query_feature - any kind of align feature.
  Args [2]   : other feature - another align feature that is
               to be checked for substantial overlap with the
               query feature.
  Example    : $self->check_overlap($query_feature, $other_feature);
  Description: Checks two features for greater than 90% hit overlap.
               Features must be blast-generated align features 
               resultant from a blast run against a common hit
               sequence.
  Returntype : 0 or 1
  Exceptions : none
  Caller     : $self->cluster_features

=cut


sub check_overlap {
  my ($self, $query_feature, $other_feature) = @_;

  my $query_start = $query_feature->{_hstart};
  my $query_end   = $query_feature->{_hend};

  if ($query_start > $query_end){($query_start, $query_end) = ($query_end, $query_start)}
  
  my $other_start = $other_feature->{_hstart};
  my $other_end   = $other_feature->{_hend};
  
  if ($other_start > $other_end) {($other_start, $other_end) = ($other_end, $other_start)}
  
  unless ($query_end < $other_start || $query_start > $other_end){
    # We have overlapping features
    
    my $query_length = $query_end - $query_start;
    my $other_length = $other_end - $other_start;
    
    my ($start1, $start2, $end1, $end2) = sort {$a <=> $b} ($query_start, $other_start, $query_end, $other_end);
    
    my $overlap = $end1 - $start2;
    
    if (($overlap >= (0.9 * $query_length))&&($overlap >= (0.9 * $other_length))){
      return 1
    }
  }
    
  return 0;    
}
  
=head2 check_coverage

  Args [1]   : protein id - id used by seqfetcher to retrieve
               protein object.
  Args [2]   : feature listref - reference to list of align features
               that align to parts of the protein specified in Arg [1]
  Example    : $self->check_coverage($seqname, $features);
  Description: Using a protein and a list of corresponding align
               features, calculates the percentage of the protein
               that the features cover.
  Returntype : A real (e.g. a percentage like 0.80)
  Exceptions : none
  Caller     : $self->make_runnables

=cut


sub check_coverage {
  my ($self, $seqname, $features) = @_;

  my $protein = $self->get_Sequence($seqname);  
  my $protein_length = $protein->length;

  my @voider;
  my @score_keeper;

  foreach my $feature (@$features){

    my $start = $feature->{_hstart};
    my $end = $feature->{_hend};

    ($start, $end) = sort {$a <=> $b} ($start, $end);

    for (my $j = $start; $j <= $end; $j++){
      $score_keeper[$j] = 1;
    }
  }

  my $tally = 0;

  foreach my $score (@score_keeper){
    $tally++ if $score;
  }

  return ($tally/$protein_length);
}


=head2 check_repeated

    Title   :   check_repeated
    Usage   :   $self->check_repeated(1)
    Function:   Get/Set method for check_repeated
    Returns :   0 (False) or 1 (True)
    Args    :   0 (False) or 1 (True)

=cut

sub check_repeated {
  my( $self, $value ) = @_;    

  if ($value) {
    $self->{'_check_repeated'} = $value;
  }

  return $self->{'_check_repeated'};
}



return 1;
