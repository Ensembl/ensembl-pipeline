#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise->new
    ('-genomic'        => $genseq,
     '-features'       => $features,
     '-protein'        => $protein,
     '-seqfetcher'     => $seqfetcher,
     '-check_repeated' => 1);

    
    $obj->run

    my @newfeatures = $obj->output;

(where $protein and $genseq are Bio::Seq objects, 
 $features are X objects and $seqfetcher is a 
 SeqFetcher object.)


=head1 DESCRIPTION

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise;

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

require "Bio/EnsEMBL/Pipeline/pipeConf.pl";

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

sub new {
  my ($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  
  my( $genomic, $ids, $seqfetcher, $endbias, $check_repeated) = $self->_rearrange([qw(GENOMIC
										      IDS
										      SEQFETCHER
										      ENDBIAS
										      CHECK_REPEATED)],
										  @args);
  
  $self->throw("No genomic sequence input")            unless defined($genomic);
  $self->throw("No seqfetcher provided")               unless defined($seqfetcher);
  $self->throw("No ids arrary ref provided")           unless defined($ids);

  $self->throw("[$genomic] is not a Bio::PrimarySeqI") 
    unless $genomic->isa("Bio::PrimarySeqI");
	
  $self->ids($ids)                                     if defined($ids);
  $self->genomic_sequence($genomic)                    if defined($genomic);
  $self->endbias($endbias)                             if defined($endbias);
  $self->seqfetcher($seqfetcher)                       if defined($seqfetcher);

  if (defined $check_repeated){
    $self->check_repeated($check_repeated);
  }else {
    $self->check_repeated(0);
  }

  return $self;
}

sub ids {
  my ($self,$ids) = @_;

	if (!defined($self->{_idlist})) {
		$self->{_idlist} = [];
	}
	if (defined($ids)) {
    if (ref($ids) eq "ARRAY") {
      push(@{$self->{'_idlist'}},@$ids);
    } else {
      $self->throw("[$ids] is not an array ref.");
    }
  }
	return @{$self->{_idlist}};
}

=head2 genomic_sequence

    Title   :   genomic_sequence
    Usage   :   $self->genomic_sequence($seq)
    Function:   Get/set method for genomic sequence
    Returns :   Bio::Seq object
    Args    :   Bio::Seq object

=cut

sub genomic_sequence {
    my( $self, $value ) = @_;    
    if ($value) {
        #need to check if passed sequence is Bio::Seq object
        $value->isa("Bio::PrimarySeqI") || $self->throw("Input isn't a Bio::PrimarySeqI");
        $self->{'_genomic_sequence'} = $value;
    }
    return $self->{'_genomic_sequence'};
}

=head2 endbias

    Title   :   endbias
    Usage   :   $self->endbias($endbias)
    Function:   Get/set method for genewise endbias
    Returns :   
    Args    :   

=cut

sub endbias {
    my ($self,$arg) = @_;

    if (defined($arg)) {
        $self->{'_endbias'} = $arg;
    }

    if (!defined($self->{'_endbias'})) {
      $self->{'_endbias'} = 0;
    }    

    return $self->{'_endbias'};
}

=head2 seqfetcher

    Title   :   seqfetcher
    Usage   :   $self->seqfetcher($seqfetcher)
    Function:   Get/set method for SeqFetcher
    Returns :   Bio::EnsEMBL::Pipeline::SeqFetcher object
    Args    :   Bio::EnsEMBL::Pipeline::SeqFetcher object

=cut

sub seqfetcher {
  my( $self, $value ) = @_;    
  if ($value) {
    $self->{'_seqfetcher'} = $value;
  }
  return $self->{'_seqfetcher'};
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


=head2 run

  Title   : run
  Usage   : $self->run()
  Function: 
  Returns : none
  Args    : 

=cut

sub run {
  my ($self) = @_;
  
  my @features = $self->run_blast;

  #print STDERR "There are ".@features." remaining features after reblasting.\n";
  unless (@features) {
    print STDERR "Contig has no associated features.  Finishing run.\n";
    return;
  }

  my $mg_runnables;

  if ($self->check_repeated){ 
    $mg_runnables = $self->build_runnables(@features);
  } else {
    my $runnable = $self->make_mmgw($self->genomic_sequence, \@features);
    push (@$mg_runnables, $runnable); 
  }

  #print STDERR "Running the MultiMiniGenewise runnables.\n";
  foreach my $mg (@$mg_runnables){
    $mg->run;
    my @f = $mg->output;
    #print STDERR "There were " . scalar @f . " $f[0]  " 
    #  . " features after the MiniGenewise run.\n";

    push(@{$self->{'_output'}},@f);
  }
  
  return 1;

}

sub run_blast {
  my ($self) = @_;
  
  my @seq         = $self->get_Sequences;
  my @valid_seq   = $self->validate_sequence(@seq);
  #print STDERR "there are ".@valid_seq." valid sequences\n";

  my $blastdb     = new Bio::EnsEMBL::Pipeline::Runnable::BlastDB(
					 -sequences => [$self->genomic_sequence],
					 -type      => 'DNA');
  #print STDERR "\n";
  $blastdb->run;
  #print STDERR "\n";
  my @features;
  my $dbname = $blastdb->dbname;
  my @sorted_seqs = sort {$a->id cmp $b->id} @valid_seq;
  foreach my $seq (@sorted_seqs) {
    # First sort out the header parsing. Blergh! cb25.NA_057.31208-61441 Slice, no descrtipion 
     
    #print STDERR "ID ".$self->genomic_sequence->id."\n";
    if($GB_INPUTID_REGEX && $self->genomic_sequence->id =~ /$GB_INPUTID_REGEX/){
      $::fasta_header_re{$dbname} = $GB_INPUTID_REGEX;
    }elsif ($self->genomic_sequence->id =~ /^(.*)\|(.*)\|(.*)/) {
      $::fasta_header_re{$dbname} = '^.*\|(.*)\|.*';
    } elsif ($self->genomic_sequence->id =~ /^..\:(.*)/) {
      $::fasta_header_re{$dbname} = '^..\:(.*)';
    }else {
      $::fasta_header_re{$dbname} = '^(\w+)\s+';
    }

    my $run = new Bio::EnsEMBL::Pipeline::Runnable::Blast(-query    => $seq,
							  -program  => 'wutblastn',
							  -database => $blastdb->dbfile,
							  -filter => 1,
							 );
    $run->run;

    push(@features,$run->output);
  }

  $blastdb->remove_index_files;
  unlink $blastdb->dbfile;

  return @features;
}


=head2 build_runnables

  Args [1]   : @features - list of features generated from a 
               blast run.
  Example    : $self->build_runnables(@features);
  Description: Builds the runnable objects for conducting
               minigenewise runs.  In the simplest case this
               method creates a single runnable object.  The
               method itself was written to cope with
               situations where the presence of repeated genes
               means that more than one whole gene of high 
               homology is present in a minigenomic sequence.
               This module separates possible repeated genes
               into clusters of exons such that minigenomic
               sequence fragments can be generated for passing
               to genewise - when this happens this method
               returns multiple minigenewise runnables that
               examine smaller portions of genomic DNA.
  Returntype : A list of Bio::EnsEMBL::Pipeline::Runnable::MultiMiniGenewise
  Exceptions : SOME
  Caller     : $self->run

=cut

sub build_runnables {
  my ($self, @features) = @_;

  my @mg_runnables;

  # We partition the features by protein id such that we can
  # iterate for each protein id.  Most ids  will fall through
  # to the simplest case where we are running genewise on one 
  # unique gene per minigenomic sequence fragment.  In cases
  # where multiple similar genes are present we switch to the
  # special case algorithm.  When this happens we cluster the 
  # likely exons of each individual gene and make small 
  # minigenomic sequences than span just one cluster.

  # Hmm, Im not entirely sure that partitioning the ids is necessary.
  
  my %partitioned_features;

  foreach my $raw_feat (@features){
    # Imporantly, the features are filtered for low identity matched.
    # This is just for the clustering process.  The full set of features
    # are still passed to genewise.
    if($raw_feat->percent_id >= 80){
      push (@{$partitioned_features{$raw_feat->hseqname}}, $raw_feat);
    }
  }

  foreach my $hseqname (keys %partitioned_features){

     # Now, we attempt to cluster our features into groups of
     # similar exons.

     my $clustered_features = $self->cluster_features(\@{$partitioned_features{$hseqname}});

    # We have to do a little test to see if the clusters
    # represent the bulk of the features present in the 
    # genomic region.  This test helps avoid the situation
    # where a non-repeated gene has exons that are very 
    # similar to one another.  If the flag $clusters_seem_real
    # is not set to one then the analysis will default to the  
    # non-repeated gene way of doing things.

    my $clusters_seem_real = 0;

    if (@$clustered_features) {
      my $number_of_clustered_features = 0;
      foreach my $cluster (@$clustered_features) {
	foreach my $clust_feat (@$cluster){
	  $number_of_clustered_features++;
	}
      }

      if ($number_of_clustered_features  > (0.9 * (scalar @{$partitioned_features{$hseqname}}))) {
	$clusters_seem_real = 1;
      }
    }

    # With the sets of clustered blast hits, find those
    # clusters with the highest number of members.  Choose
    # one of these clusters and choose the centre-most 
    # cluster member (in genomic DNA fragment terms).  With
    # this blast feature, select from each of the other
    # clusters the nearest member blast feature.  Create a
    # new array with these features.  Iterate this process
    # until the original cluster has no more members.  With
    # the arrays of features created check for clusters that
    # overlap.  If an overlap is found, discard a feature
    # from one of the arrays, choosing the feature furthest
    # from the core of its cluster.  If doing this still
    # doesn't fix the overlap, throw a nasty-looking warning 
    # and keep going.
    
    if ((@$clustered_features)&&($clusters_seem_real > 0)) {

      print STDERR "Minigenomic sequence could contain as many as " 
	. scalar @$clustered_features . " highly similar genes." 
	  . "  Fragmenting the minigenomic sequence to try to " 
	    . " resolve these genes individually.\n";

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

	my @sorted_gene_cluster = sort {$a->{_hstart} <=> $b->{_hstart};} @$gene_cluster;
       	my $cluster_start = $sorted_gene_cluster[0]->{_hstart} - 1000;
	my $cluster_end = $sorted_gene_cluster[-1]->{_hend} + 1000;
        unless ($cluster_start > 0) {
	  if ($sorted_gene_cluster[0]->{_hend} > 0) {
	    $cluster_start = 1;
	  } else {
	    next GENE;
	  }
	}
        unless ($cluster_end <= $self->genomic_sequence->length) {
	  if ($sorted_gene_cluster[0]->{_hstart} <= $self->genomic_sequence->length) {
	    $cluster_end = $self->genomic_sequence->length;
	  } else {
	    next GENE;
	  }
	}        

	# Here we create our genomic fragment for passing to 
	# genewise.  This is the fragment of genomic sequence
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

	my $mmgw = $self->make_mmgw($genomic_subseq, \@features);

	push (@mg_runnables, $mmgw);
      }
    
    }else{
      # This is what we do when we dont have multiple genes
      # in our genomic fragment.  This is "Normal mode".

      my $mmgw = $self->make_mmgw($self->genomic_sequence, \@features);

      push (@mg_runnables, $mmgw);
    }
      
  }
    return \@mg_runnables;
}

=head2 cluster_features

  Args [1]   : $features - reference to a list of 
               PepDnaAlignFeatures generated by a blast run.
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
  Caller     : $self->build_runnables

=cut

sub cluster_features {

  my ($self, $features) = @_;

  # Sort the list of features according to their percent id.
  my @sorted_features = sort { $b->percent_id <=> $a->percent_id;} @$features;

  # With sorted features, sequentially blast the 
  # highest scoring hit against the lesser scoring
  # sequences.  Hits that have high identity are clustered
  # in an array of features and these hits removed from 
  # the blast db.  Next, the  remaining highest hit is 
  # blasted against the remaining lesser scoring hits, and
  # so on.  If a feature doesn't hit anything it is
  # discarded.  
  
  # NOTE: At the moment each lesser scoring feature is 
  # separately blasted with the top scoring feature.  If all
  # features are blasted together at the same time there is
  # no way to keep track of which feature is which when then
  # are returned.  This multiple blasting may be slow way of
  # finding similar features.  An alternative might be to
  # pad the non-matching nucleotides of a genomic fragment
  # with Ns - this will still give the right  matches, in a
  # recognisable coordinate system.  Yuck though.
  
  my @all_feature_clusters;
  
  while (@sorted_features) {
    
    my $top_feature = shift @sorted_features;
    
    my @similar_features;
    push (@similar_features, $top_feature);

    my @subtracted_sorted_features;
    for(my $i = 0; $i < scalar @sorted_features; $i++) {
      my $lesser_feature = @sorted_features[$i];
      
      my @matched_feature = $self->blast_features($top_feature,
						  [$lesser_feature]);
      
      if ((@matched_feature)&&($matched_feature[0]->percent_id > 95)){

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

  # If more than one cluster has the biggest number of 
  # members, should try to choose the cluster that houses
  # a central exon (rather than a terminal exon).
  # Perhaps could add all the start coordinates per cluster and
  # take the central value...  UNIMPLEMENTED AS YET.

  my $biggest = 0;
  my $big_cluster;
  my @other_clusters;
  foreach my $cluster (@$exon_clusters){

    if ((scalar @$cluster > $biggest)) {
      if (defined $big_cluster){
	push (@other_clusters, $big_cluster);
      }
      $big_cluster = $cluster;
      $biggest = scalar @$cluster;
    } else {
      push (@other_clusters, $cluster);
    }
  }
  
  my @final_gene_clusters;
      
  while (@$big_cluster){

    # Choose the centre-most feature (genomic-location speaking).
    my @sorted_cluster = sort { $a->{_hstart} <=> $b->{_hstart};} @$big_cluster;
    my $centre_most_element = int(((scalar @$big_cluster)-1)/2);
    my $centre_most_feature = $big_cluster->[$centre_most_element];
    
    splice (@$big_cluster, $centre_most_element, 1);
	
    my @gene_cluster;
    push (@gene_cluster, $centre_most_feature);
    
    foreach my $other_cluster (@other_clusters) {
      if (@$other_cluster){
	
	my @subtracted_other_clusters;
	my $closest = 1000000;
	my $closest_feature;
	foreach my $other_feature (@$other_cluster){
	  
	  my $distance = abs($other_feature->{_hstart} 
			     - $centre_most_feature->{_hstart});
	  
	  if($distance < $closest){
	    if (defined $closest_feature){
	      push (@subtracted_other_clusters, $closest_feature);
	    }
	    $closest = $distance;
	    $closest_feature = $other_feature;
	  } else {
	    push (@subtracted_other_clusters, $other_feature);
	  }
	}
	@$other_cluster = @subtracted_other_clusters;

	push (@gene_cluster, $closest_feature) if (defined $closest_feature);
	$closest_feature = <UNDEF>;
      }
    }
    push (@final_gene_clusters, \@gene_cluster);
  }	

  return \@final_gene_clusters;
}



=head2 make_mmgw

  Args [1]   : $miniseq - a Bio::Seq object representing the
               target sequence for the genewise run.
  Args [2]   : $features - reference to a list of 
               PepDnaAlignFeatures generated by a blast run.
  Example    : $self->makemmgw($miniseq, $features);
  Description: Takes a genomic sequence and builds a
               MultiMiniGenewise runnable object using the 
               list of PepDnaAlignFeatures.
  Returntype : A list of 
               Bio::EnsEMBL::Pipeline::Runnable::MultiMiniGenewise
  Exceptions : none
  Caller     : $self->build_runnables

=cut


sub make_mmgw {
  my ($self, $miniseq, $features) = @_;

  # Before we pass our blast generated features to 
  # MultiMiniGenewise we must first convert them from 
  # PepDnaAlignFeatures to FeaturePairs.

  my @newf;
  foreach my $f (@$features){
    my $newf = new Bio::EnsEMBL::FeaturePair(-feature1 => $f->feature2,
					     -feature2 => $f->feature1);
    push(@newf,$newf);
  }

  # Create a MultiMiniGenewise object with the features weve
  # just converted.

  my $mg      = new Bio::EnsEMBL::Pipeline::Runnable::MultiMiniGenewise(
				       '-genomic'    => $miniseq,
				       '-features'   => \@newf,
				       '-seqfetcher' => $self->seqfetcher,
				       '-endbias'    => $self->endbias
				      );

  return $mg;
}




=head2 blast_features

  Args [1]   : query_feature - any kind of align feature.
  Args [2]   : other features - a reference to a list of 
               features that will be used as the db sequences.
  Example    : $self->blast_features($query_feature, $db_features);
  Description: Runs a blast job of a query blast feature
               against an array of other blast features.
               Currently is used to determine they level of
               identity between blast features.  This is by
               no means the best way to this, but for now it
               works well.  This method is sufficiently
               different to $self->run_blast that it can't be
               merged in a neat way.
  Returntype : A list of matched features.              
  Exceptions : none
  Caller     : *ONLY* $self->make_runnables

=cut


sub blast_features {
  my ($self, $query_feature, $blast_features) = @_;

  # Obtain a sequence from each blast feature in order to
  # create a blast-able database.  We are only concerned
  # with strand so as to hand seq->subseq start and end
  # coords in the right order.  The database may be made with
  # sequence from the wrong strand, but this doesn't matter
  # for our purposes.

  my @db_seqs;

  foreach my $blast_feature (@$blast_features){

    my ($feat_start, $feat_end) = sort {$a <=> $b;} $blast_feature->{_hstart}, $blast_feature->{_hend};

    my $feat_seq = Bio::Seq->new(
                     -seq => $self->genomic_sequence->subseq($feat_start, $feat_end),
		     -id  => $blast_feature->seqname);

    push (@db_seqs, $feat_seq);
  }


  my $blast_db     = new Bio::EnsEMBL::Pipeline::Runnable::BlastDB(
					 -sequences => \@db_seqs,
					 -type      => 'DNA');

  $blast_db->run;

  my $db_name = $blast_db->dbname;

  # With a database now at hand, bam, make a query sequence object and time for the run.
  my ($query_start, $query_end) = sort {$a <=> $b;} $query_feature->{_hstart}, $query_feature->{_hend};

  my $query_seq = Bio::Seq->new(
	           -seq => $self->genomic_sequence->subseq($query_start, $query_end),
		   -id  => $query_feature->seqname);

  # BlastConf.pl magic that sets the correct regular expression
  # for parsing the header of the query sequence.

  if ($query_seq->id =~ /\W*\w+/){
    $::fasta_header_re{$db_name} = '\W*(\w+)';
  } else {
    $self->warn('Problem with regex - unlikely to be parsing correct sequence ids.');
  }

  my $blast = new Bio::EnsEMBL::Pipeline::Runnable::Blast(
			      -query    => $query_seq,
			      -program  => 'wublastn',
			      -database => $blast_db->dbfile,
			      -filter => 0,
			      );
  $blast->run;
  
  $blast_db->remove_index_files;
  unlink $blast_db->dbfile;
  
  return $blast->output;
}



sub get_Sequences {
    my ($self) = @_;

    my @seq;

    foreach my $id ($self->ids) {
        my $seq = $self->get_Sequence($id);

        if (defined($seq) && $seq->length > 0) {
            push(@seq,$seq);
        } else {
            print STDERR "Invalid sequence for $id - skipping\n";
        }
    }

    return @seq;

}

sub validate_sequence {
    my ($self,@seq) = @_;
    my @validated;

    foreach my $seq (@seq) {

        my $sequence = $seq->seq;

        if ($sequence !~ /[^acgtn]/i) {
            push (@validated, $seq);
        } else {
            $_ = $sequence;
            my $len = length ($_);
            my $invalidCharCount = tr/bB/xX/;

            if ($invalidCharCount / $len > 0.05) {
                $self->warn("Ignoring ".$seq->display_id()
                    ." contains more than 5% ($invalidCharCount) "
                    ."odd nucleotide codes ($sequence)\n Type returns "
                    .$seq->moltype().")\n");
            } else {
                $seq->seq($_);
                push (@validated, $seq);
            }
        }
    } 
    return @validated;  
}

=head2 get_Sequence

  Title   : get_Sequence
  Usage   : my $seq = get_Sequence($id)
  Function: Fetches sequences with id $id
  Returns : Bio::PrimarySeq
  Args    : none

=cut
    
sub get_Sequence {
    my ($self,$id) = @_;
    my $seqfetcher = $self->seqfetcher;
    my $seq;

    if (!defined($id)) {
      $self->warn("No id input to get_Sequence");
    }  
    
    eval {
      $seq = $seqfetcher->get_Seq_by_acc($id);
    };

    if($@) {
      $self->warn("Problem fetching sequence for id [$id] $@\n");
      return undef;
    }
    
    if(!defined($seq)){
      $self->warn("Could not find sequence for [$id]");
    }

    return $seq;
	}

1;
