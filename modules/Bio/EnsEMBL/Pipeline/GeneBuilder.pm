#
# BioPerl module for GeneBuilder
#
# Cared for by EnsEMBL <ensembl-dev@ebi.ac.uk>
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::GeneBuilder

=head1 SYNOPSIS

# This is the main analysis database

my $db = new Bio::EnsEMBL::DBSQL::Obj(-host   => 'obi-wan',
				      -user   => 'ensro',
				      -dbname => 'ens500',
				      );

# Fetch a clone and its contigs from the database
my $clone       = $db   ->get_Clone($clone);
my @contigs     = $clone->get_all_Contigs;

# The genebuilder object will fetch all the features from the contigs
# and use them to first construct exons, then join those exons into
# exon pairs.  These exon apris are then made into transcripts and
# finally all overlapping transcripts are put together into one gene.


my $genebuilder = new Bio::EnsEMBL::Pipeline::GeneBuilder
    (-contigs => \@contigs);

my @genes       = $genebuilder->build_Genes;

# After the genes are built they can be used to order the contigs they
# are on.

my @contigs     = $genebuilder->order_Contigs;


=head1 DESCRIPTION

This module reads your favourite annotations (genewise, combined_genes, est2genome, genomewise, ... )
on the one hand, and ab initio predictions plus features on the other hand. Ab initio predictions and features
are passed to Bio::EnsEMBL::Pipeline::Runnable::PredictionGeneBuilder which generates putative transcripts
from supported prediction exons (see documentation in that module for details).
The product of Bio::EnsEMBL::Pipeline::Runnable::PredictionGeneBuilder is combined with
all the other annotations and redundant transcripts are eliminated in the method prune_Transcripts().
The resulting transcripts are combined into genes. For more details, follow the list of methods called
by build_Genes() method and the description in each one.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::GeneBuilder;

use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::Runnable::PredictionGeneBuilder;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster;

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Blessed   qw (
							       GB_BLESSED_GENETYPES
							      );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Targetted   qw (
							       GB_TARGETTED_GW_GENETYPE
							      );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Combined    qw (
							       GB_COMBINED_GENETYPE
							      );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Similarity  qw (
							       GB_SIMILARITY_GENETYPE
							      );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::General     qw (
							       GB_INPUTID_REGEX
							      );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::GeneBuilder qw (
							       GB_MIN_GENSCAN_EXONS
							       GB_GENSCAN_MAX_INTRON
							       GB_MIN_FEATURE_SCORE
							       GB_MIN_FEATURE_LENGTH
							       GB_ABINITIO_TYPE
							       GB_ABINITIO_SUPPORTED_TYPE
							       GB_MAXSHORTINTRONLEN
							       GB_MINSHORTINTRONLEN
							       GB_MAX_TRANSCRIPTS_PER_GENE
							       GB_USE_ABINITIO
							       GB_CONFIRM_PFAM
							      );

use vars qw(@ISA);
use strict;

@ISA = qw(Bio::EnsEMBL::Root);


############################################################

sub new {
    my ($class,@args) = @_;

    my $self = $class->SUPER::new(@args);

    my ($slice,$input_id) = $self->_rearrange([qw(SLICE INPUT_ID)],
					      @args);

    $self->throw("Must input a slice to GeneBuilder") unless defined($slice);
    $self->{_final_genes} = [];
    $self->{_gene_types}  = [];

    $self->query($slice);
    $self->gene_types($GB_COMBINED_GENETYPE);
    $self->gene_types($GB_TARGETTED_GW_GENETYPE);
    $self->gene_types($GB_SIMILARITY_GENETYPE);
    $self->gene_types("KnownUTR");
    foreach my $bgt(@{$GB_BLESSED_GENETYPES}){
      $self->gene_types($bgt->{'type'});
    }

    unless ( $input_id =~ /$GB_INPUTID_REGEX/ ){
      $self->throw("format of the input is not defined in Config::GeneBuild::General::GB_INPUTID_REGEX = $GB_INPUTID_REGEX");
    }
    $self->input_id($input_id);

    return $self;
}

############################################################

=head2 input_id

 Function: get/set for input id
 Returns : string
 Args    : string (it expects a string of the format chr_name.start_coord-end_coord

=cut
  
sub input_id {
  my ($self,$id) = @_;
  
  if (defined($id)) {
    $self->{_input_id} = $id;
  }
  return $self->{_input_id};
}

############################################################

=head2 build_Genes

 Example    : my @genes = $self->build_Genes
 Description: builds genes. It is like the run method in Runnables. It calls everything that needs to be done.
 Returns    : none
 Args       : none
 Caller     : Bio::EnsEMBL::Pipeline::RunnableDB::Gene_Builder

=cut

sub build_Genes{
  my ($self) = @_;
  
  print STDERR "Building genes...\n";
  
  # get all genes of type defined in gene_types() on this slice
  $self->get_Genes;
  print STDERR "After checks: Number of genewise and combined transcripts " . scalar($self->combined_Transcripts) . "\n";
  
  #test
#  foreach my $t ( $self->combined_Transcripts ){
  #  Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($t);

#  }

  
  if ( $GB_USE_ABINITIO ){
    # get all Genscan predictions on this slice
    $self->get_Predictions;
    print STDERR "Number of ab initio predictions ". scalar($self->predictions)  . "\n";
    
    # get all the dna/protein align features from the pre-computes pipeline on this slice
    $self->get_Similarities;
    print STDERR "Number of similarity features ". scalar($self->features) . "\n";
  }
  
  my @supported_predictions;
  my @annotations = $self->combined_Transcripts;
  
  if ($GB_USE_ABINITIO ){
    # process PredictionTranscripts using the features and the annotations:
      
      my @predictions = $self->predictions;
      my @features    = $self->features;
      if(@predictions && @features){
	  my $genecooker  = Bio::EnsEMBL::Pipeline::Runnable::PredictionGeneBuilder->new(
											 -predictions => \@predictions,
											 -features    => \@features,
											 -annotations => \@annotations,
											 );
	  
	  $genecooker->query($self->query);

	  if ($GB_CONFIRM_PFAM) {
	      @supported_predictions = $genecooker->run_pfam;
	  }
	  else {
	       @supported_predictions = $genecooker->run;
	   }
      }
  }
  
  # cluster all the transcripts according to mere genomic overlap
  my @all_transcripts;
  push ( @all_transcripts, @annotations );
  
  if (@supported_predictions){
    push( @all_transcripts, @supported_predictions );
  }
  
  unless( @all_transcripts ){
      print STDERR "no transcripts left to cook. Exiting...\n";
      return;
  }

  print STDERR "clustering transcripts...\n";
  my @transcript_clusters = $self->cluster_Transcripts(@all_transcripts);
  print STDERR scalar(@transcript_clusters)." clusters formed\n";
  
  # prune the redundant transcripts for each cluster
  print STDERR "pruning transcripts...\n";
  my @pruned_transcripts = $self->prune_Transcripts(@transcript_clusters);
  print STDERR scalar(@pruned_transcripts)." transcripts obtained\n";
  
  # cluster transcripts into genes
  print STDERR "clustering into genes...\n";
  
  # do a preliminary clustering
  my @preliminary_genes = $self->cluster_into_Genes(@pruned_transcripts);

  # Remove CDS from gene where transcript has no UTR and has same CDS as one with UTR
  $self->_prune_redundant_CDS(\@preliminary_genes);

  # select the best ones per gene (use the GB_MAX_TRANSCRIPTS_PER_GENE )
  my @best_transcripts = $self->_select_best_transcripts( @preliminary_genes );
  
  # recluster the chosen transcripts into genes
  my @tmp_genes = $self->cluster_into_Genes(@best_transcripts);
  
  # make shared exons unique objects
  my @genes =  $self->_make_shared_exons_unique( @tmp_genes );
  
  print STDERR scalar(@genes)." genes built\n";
  
  print STDERR "Final_results:\n";
  my $count = 0;
  foreach my $gene ( @genes ){
    $count++;
    print STDERR "Gene $count:\n";
    foreach my $tran ( @{$gene->get_all_Transcripts} ){
      Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($tran);
    }
  }    
  print STDERR scalar( @genes )." final genes\n";
  $self->final_genes( @genes );
}

sub _prune_redundant_CDS {
  my ( $self, $genes ) = @_;

  my $nremoved = 0;

  my %blessed_genetypes;
  foreach my $bgt(@{$GB_BLESSED_GENETYPES}){
    $blessed_genetypes{$bgt->{'type'}} = 1;
  }

# For each gene
  foreach my $gene (@$genes) {
    my @trans_with_utrs;
    my @trans_without_utrs;

# Separate into ones with UTRs and ones without
    foreach my $trans (@{$gene->get_all_Transcripts}) {
      $trans->sort;
      my @exons = @{$trans->get_all_Exons};
      if ($trans->translation) {
        if ($trans->translation->start_Exon == $exons[0] &&
            $trans->translation->start == 1 &&
            $trans->translation->end_Exon == $exons[$#exons] &&
            $trans->translation->end == $exons[$#exons]->length) {
          push @trans_without_utrs, $trans;
        } else {
          push @trans_with_utrs, $trans;
        }
      } else {
        $self->warn("No translation for transcript");
      }
    }

# Generate CDS exons just once (saves CPU time)
    my %cds_exon_hash;
    foreach my $trans (@{$gene->get_all_Transcripts}) {
      # Shouldn't need sort but paranoid about this
      $trans->sort;
      if ($trans->translation) {
        my @cds_exons = @{$trans->get_all_translateable_Exons};
        $cds_exon_hash{$trans} = \@cds_exons;
      }
    }
      
# Compare CDSs of ones with UTRs to CDSs of ones without
    foreach my $utr_trans (@trans_with_utrs) {
      my $utr_trans_cds_exons = $cds_exon_hash{$utr_trans};

      CDS_TRANS: foreach my $cds_trans (@trans_without_utrs) {
        my $cds_trans_cds_exons = $cds_exon_hash{$cds_trans};

        if (scalar(@$cds_trans_cds_exons) != scalar(@$utr_trans_cds_exons)) {
          # print "Different numbers of exons\n";
          next CDS_TRANS;
        }

        my @cds_exons = @$cds_trans_cds_exons;
        foreach my $utr_trans_exon (@$utr_trans_cds_exons) {
          my $cds_trans_exon = shift @cds_exons;
          if ($cds_trans_exon->start     != $utr_trans_exon->start  ||
              $cds_trans_exon->end       != $utr_trans_exon->end    ||
              $cds_trans_exon->strand    != $utr_trans_exon->strand) {
            # print "Exon difference for " . $cds_trans_exon->start . "-" . $cds_trans_exon->end . " and " . $utr_trans_exon->start . "-" . $utr_trans_exon->end . "\n";
            next CDS_TRANS;
          }

        }
# Remove non UTR one if CDS is same  (only get to here if all exons match 
        # print "Removing transcript as redundant CDS\n";
        # Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($cds_trans);

        if (  exists $blessed_genetypes{$cds_trans->type} && 
            ! exists $blessed_genetypes{$utr_trans->type}) {
# Hack to make sure transcript gets through if its like a blessed one
          # print "Hacking transcript type from " . $utr_trans->type . " to " . $cds_trans->type . "\n";
          $utr_trans->type($cds_trans->type);
        }

        $nremoved++;
        $self->_remove_transcript_from_gene($gene,$cds_trans);
      }
    }
  }
  print "Removed $nremoved transcripts because of redundant CDS\n";
}

sub _remove_transcript_from_gene {
  my ($self, $gene, $trans_to_del)  = @_;

  my @newtrans;
  foreach my $trans (@{$gene->get_all_Transcripts}) {
    if ($trans != $trans_to_del) {
      push @newtrans,$trans;
    }
  }

# The naughty bit!
  $gene->{_transcript_array} = [];

  foreach my $trans (@newtrans) {
    $gene->add_Transcript($trans);
  }

  return scalar(@newtrans);
}






############################################################

sub _select_best_transcripts{
  my ( $self, @genes ) = @_;
  my @selected_transcripts;


  # make sure we don't discard blessed transcripts
  my %blessed_genetypes;
  foreach my $bgt(@{$GB_BLESSED_GENETYPES}){
    $blessed_genetypes{$bgt->{'type'}} = 1;
  }

 GENE:
  foreach my $gene ( @genes ){
    # sort the transcripts, get the longest CDS + UTR first (like in prune_Transcripts(); 
    # blessed transcripts come at the top )
    my @sorted_transcripts = $self->_bin_sort_transcripts( @{$gene->get_all_Transcripts} );
    my $count = 0;
  TRAN:
    foreach my $transcript( @sorted_transcripts ){
      $count++;
      unless (exists $blessed_genetypes{$transcript->type}){
	next GENE if ($count > $GB_MAX_TRANSCRIPTS_PER_GENE);
      }
      push ( @selected_transcripts, $transcript );
    }
  }
  return @selected_transcripts;
}

############################################################

sub _make_shared_exons_unique{
  my ( $self, @genes ) = @_;
  my @pruned_genes;
  foreach my $gene ( @genes ){
    
    # make different exon objects that are shared between transcripts 
    # ( regarding attributes: start, end, etc )
    # into unique exon objects 
    my $new_gene = $self->prune_Exons($gene);
    push( @pruned_genes, $new_gene );
  }
  return @pruned_genes;
}

############################################################


=head2 get_Genes

 Description: retrieves genewise and combined gene annotations with supporting evidence. 
              Splits transcripts with very long introns, discards transcripts with strand problems etc.
 ReturnType : none, but $self->combined_Transcripts is filled
 Args       : none

=cut

sub get_Genes {
  my ($self) = @_;
  my @transcripts;
  my $db = $self->genes_db;
  my $sa = $db->get_SliceAdaptor;
  
  my $input_id = $self->input_id;
  $input_id =~/$GB_INPUTID_REGEX/;
  my $chr   = $1;
  my $start = $2;
  my $end   = $3;
  my $slice = $sa->fetch_by_chr_start_end($chr,$start,$end);
  
  my @unchecked_genes;

  foreach my $type ($self->gene_types) {
    my @genes = @{$slice->get_all_Genes_by_type($type)};
    print STDERR "Retrieved ".scalar(@genes)." genes of type ".$type."\n";
    foreach my $gene ( @genes ){
    
    TRANSCRIPT:
      foreach my $tran (@{$gene->get_all_Transcripts}) {
	
	# do NOT check intron sizes - they're already tightly controlled by Targetted & Similarity stages
	unless ( Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_check_Transcript( $tran,$self->query ) && 
		 Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_check_Translation($tran)
	       ){
	  next TRANSCRIPT;
	}
	$tran->type($type);
	push(@transcripts, $tran);
      }
    }
  }
  $self->combined_Transcripts(@transcripts);
}


###########################################################c

=head2 cluster_Transcripts

 Description : It separates transcripts according to strand and then clusters 
               each set of transcripts by calling _cluster_Transcripts_by_genomic_range()
  Args       : Array of Bio::EnsEMBL::Transcript
  Return     : Array of Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster

=cut

sub cluster_Transcripts {
  my ($self,@transcripts) = @_;
 
  my @forward_transcripts;
  my @reverse_transcripts;
 
  foreach my $transcript (@transcripts){
    my @exons = @{ $transcript->get_all_Exons };
    if ( $exons[0]->strand == 1 ){
      push( @forward_transcripts, $transcript );
    }
    else{
      push( @reverse_transcripts, $transcript );
    }
  }
  
  my @forward_clusters;
  my @reverse_clusters;
  
  if ( @forward_transcripts ){
    @forward_clusters = $self->_cluster_Transcripts_by_genomic_range( @forward_transcripts );
  }
  if ( @reverse_transcripts ){
    @reverse_clusters = $self->_cluster_Transcripts_by_genomic_range( @reverse_transcripts );
  }
  my @clusters;
  if ( @forward_clusters ){
    push( @clusters, @forward_clusters);
  }
  if ( @reverse_clusters ){
    push( @clusters, @reverse_clusters);
  }
  return @clusters;
}

############################################################

=head2 _cluster_Transcripts_by_genomic_range

 Description : It clusters transcripts according to genomic overlap
  Args       : Array of Bio::EnsEMBL::Transcript
  Return     : Array of Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster

=cut

sub _cluster_Transcripts_by_genomic_range{
  my ($self,@mytranscripts) = @_;
  # first sort the transcripts
  my @transcripts = sort { my $result = ( $self->transcript_low($a) <=> $self->transcript_low($b) );
				 if ($result){
				     return $result;
				 }
				 else{
				     return ( $self->transcript_high($b) <=> $self->transcript_high($a) );
				 }
			     } @mytranscripts;

  # create a new cluster 
  my $cluster=Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster->new();
  my $count = 0;
  my @cluster_starts;
  my @cluster_ends;
  my @clusters;
  
  # put the first transcript into these cluster
  $cluster->put_Transcripts( $transcripts[0] );

  $cluster_starts[$count] = $self->transcript_low($transcripts[0]);
  $cluster_ends[$count]   = $self->transcript_high($transcripts[0]);
  
  # store the list of clusters
  push( @clusters, $cluster );
  
  # loop over the rest of the transcripts
 LOOP1:
  for (my $c=1; $c<=$#transcripts; $c++){
    #print STDERR "\nIn cluster ".($count+1)."\n";
    #print STDERR "start: $cluster_starts[$count] end: $cluster_ends[$count]\n";
    #print STDERR "comparing:\n";
    #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript( $transcripts[$c] );
    
    if ( !( $transcripts[$c]->end < $cluster_starts[$count] ||
	    $transcripts[$c]->start > $cluster_ends[$count] ) ){
      $cluster->put_Transcripts( $transcripts[$c] );
      
      # re-adjust size of cluster
      if ($self->transcript_low($transcripts[$c]) < $cluster_starts[$count]) {
	$cluster_starts[$count] = $self->transcript_low($transcripts[$c]);
      }
      if ( $self->transcript_high($transcripts[$c]) > $cluster_ends[$count]) {
	$cluster_ends[$count] =   $self->transcript_high($transcripts[$c]);
      }
    }
    else{
      # else, create a new cluster with this feature
      $count++;
      $cluster = Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster->new();
      $cluster->put_Transcripts( $transcripts[$c] );
      $cluster_starts[$count] = $self->transcript_low( $transcripts[$c]);
      $cluster_ends[$count]   = $self->transcript_high($transcripts[$c]);
      
      # store it in the list of clusters
      push(@clusters,$cluster);
    }
  }
  return @clusters;
}

############################################################


=head2 prune_Transcripts

 Example    : my @pruned_transcripts = $self->prune_Transcripts(@transcript_clusters)
 Description: rejects duplicate transcripts, transfers supporting feature data from the rejected transcripts
               to the accepted ones
 Returns    : array of Bio::EnsEMBL::Transcript
  Args      : array of Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster

=cut

sub prune_Transcripts {
  my ($self, @transcript_clusters) = @_;
  my @newtran;
  my %blessed_genetypes;
  foreach my $bgt(@{$GB_BLESSED_GENETYPES}){
    $blessed_genetypes{$bgt->{'type'}} = 1;
  }

  my $cluster_count = 0;
 CLUSTER:
  foreach my $transcript_cluster ( @transcript_clusters ){
    $cluster_count++;
    print STDERR "Cluster $cluster_count\n";
    my @mytranscripts = @{$transcript_cluster->get_Transcripts};

    ########################
    #
    # sort the transcripts
    #
    ########################

    print STDERR "sorting transcripts in cluster $cluster_count...\n";
    my @transcripts = $self->_bin_sort_transcripts( @mytranscripts );

    ##############################
    #
    # deal with single exon genes
    #
    ##############################
#    my @maxexon = @{$transcripts[0]->get_all_Exons};

    # do we really just want to take the first transcript only? What about supporting evidence from other transcripts?
    # also, if there's a very long single exon gene we will lose any underlying multi-exon transcripts
    # this may increase problems with the loss of valid single exon genes as mentioned below. 
    # it's a balance between keeping multi exon transcripts and losing single exon ones
    #if ($#maxexon == 0 && $max_num_exons == 1) {
    #  push(@newtran, $transcripts[0] );
    ## we are done with this cluster
    #  next CLUSTER;
    #}

    my $maxexon_number = 0;
    foreach my $t (@transcripts){
      if ( scalar(@{$t->get_all_Exons}) > $maxexon_number ){
	$maxexon_number = scalar(@{$t->get_all_Exons});
      }
    }
    if ($maxexon_number == 1){
      # take the longest:
      @transcripts = map { $_->[1] } sort { $b->[0]->length <=> $a->[0]->length } map{ [ $_->start_Exon, $_ ] } @transcripts;
      my $tran = shift( @transcripts );
      push (@newtran, $tran);
#      print STDERR "found single_exon_transcript\n";
#      Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($tran);
      my @es = @{$tran->get_all_Exons};
      my $e  = $es[0];
      foreach my $transcript (@transcripts){
	# make sure we keep it if it's blessed
	if(exists $blessed_genetypes{$transcript->type}){

	  push(@newtran, $transcript);
	}
	else{
	  foreach my $exon ( @{$transcript->get_all_Exons} ){
	    $self->transfer_supporting_evidence($exon, $e);
	  }
	}
      }
      next CLUSTER;
    }
    # otherwise we need to deal with multi exon transcripts and reject duplicates.

    # links each exon in the transcripts of this cluster with a hash of other exons it is paired with
    my %pairhash;

    # allows retrieval of exon objects by exon->id - convenience
    my %exonhash;

    # keep track of single exon transcripts (automatically rejected if the "top" transcript is multi exon), 
    # so we can check their supporting evidence at the very end.
    my %single_exon_rejects;    

    ##############################
    #
    # prune redundant transcripts
    #
    ##############################
  TRANSCRIPT:
    foreach my $tran (@transcripts) {

      $tran->sort;
      my @exons = @{$tran->get_all_Exons};

      #print STDERR "\ntranscript: ".$tran->dbID."\n";
      #foreach my $exon ( @exons ){
      #  print STDERR $exon->start."-".$exon->end." ";
      #}
      #print STDERR "\n";

      my $i     = 0;
      my $found = 1;

      # if this transcript has already been seen, this
      # will be used to transfer supporting evidence
      my @evidence_pairs;

      # 10.1.2002 VAC we know there's a potential problem here - single exon transcripts which are in a 
      # cluster where the longest transcriopt has > 1 exon are not going to be considered in 
      # this loop, so they'll always be marked "transcript already seen"
      # How to sort them out? If the single exon overlaps an exon in a multi exon transcript then 
      # by our rules it probably ought to be rejected the same way transcripts with shared exon-pairs are.
      # 27/5/2003 VAC But the supporting evidence should be transferred if we can! 
      # Tough one.
      # more of a problem is that if the transcript with the largest number of exons is really a
      # single exon with frameshifts, it will get rejected here based on intron size but in addition
      # any valid non-frameshifted single exon transcripts will get rejected - which is definitely not right.
      # We need code to represent frameshifted exons more sensibly so the frameshifted one doesn't 
      # get through the check for single exon genes above.

    EXONS:
      for ($i = 0; $i < $#exons; $i++) {
	my $foundpair = 0;
	my $exon1 = $exons[$i];
	my $exon2 = $exons[$i+1];
		
	# Only count introns > 50 bp as real introns
	my $intron;
	if ($exon1->strand == 1) {
	  $intron = abs($exon2->start - $exon1->end - 1);
	}
	else {
	  $intron = abs($exon1->start - $exon2->end - 1);
	}
	
	if ($intron < $GB_MAXSHORTINTRONLEN && $intron > $GB_MINSHORTINTRONLEN ) {
	  print STDERR "Intron too short: $intron bp. Transcript will be rejected\n";
	  $foundpair = 1;	# this pair will not be compared with other transcripts
	}
	else {
	
	  # go through the exon pairs already stored in %pairhash. 
	  # If there is a pair whose exon1 overlaps this exon1, and 
	  # whose exon2 overlaps this exon2, then these two transcripts are paired
	
	  foreach my $first_exon_id (keys %pairhash) {
	    my $first_exon = $exonhash{$first_exon_id};

	    foreach my $second_exon_id (keys %{$pairhash{$first_exon}}) {
	      my $second_exon = $exonhash{$second_exon_id};

	      if ( $exon1->overlaps($first_exon) && $exon2->overlaps($second_exon) ) {
		$foundpair = 1;
		
		# eae: this method allows a transcript to be covered by exon pairs
		# from different transcripts, rejecting possible
		# splicing variants. Needs rethinking
		
		# we put first the exon from the transcript being tested:
		push( @evidence_pairs, [ $exon1 , $first_exon  ] );
		push( @evidence_pairs, [ $exon2 , $second_exon ] );
		
		# transfer evidence between exons, assuming the suppfeat coordinates are OK.
		# currently not working as the supporting evidence is not there - 
		# can get it for genewsies, but why not there for genscans?
		$self->transfer_supporting_evidence($exon1, $first_exon);
		$self->transfer_supporting_evidence($first_exon, $exon1);
		$self->transfer_supporting_evidence($exon2, $second_exon);
		$self->transfer_supporting_evidence($second_exon, $exon2);
	      }
	    }
	  }
	}
	
	if ($foundpair == 0) {	# ie this exon pair does not overlap with a pair yet found in another transcript
	  $found = 0;		# ie currently this transcript is not paired with another

	  # store the exons so they can be retrieved by id
	  $exonhash{$exon1} = $exon1;
	  $exonhash{$exon2} = $exon2;

	  # store the pairing between these 2 exons
	  $pairhash{$exon1}{$exon2} = 1;
	}
      }				# end of EXONS

      # decide whether this is a new transcript or whether it has already been seen

      # if it's blessed, we keep it and there's nothing more to do
      if(exists $blessed_genetypes{$tran->type}){
	push (@newtran, $tran);
      }

      elsif ($found == 0) {
	#print STDERR "found new transcript:\n";
	#Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript( $tran );
	
	push(@newtran,$tran);
	@evidence_pairs = ();
      }
      elsif ($found == 1 && $#exons == 0){
	# save the transcript and check though at the end to see if we can transfer supporting
	# evidence; if we try it now we may not (yet) have any exons that overlap in %exonhash
	$single_exon_rejects{$tran} = $tran;
      }
      else {
	# print STDERR "transcript already seen:\n";
	if ( $tran == $transcripts[0] ){
	  print STDERR "Strange, this is the first transcript in the cluster!\n";
	}
	# Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript( $tran );
	
	## transfer supporting feature data. We transfer it to exons
	foreach my $pair ( @evidence_pairs ){
	  my @pair = @$pair;
	
	  # first in the pair is the 'already seen' exon
	  my $source_exon = $pair[0];
	  my $target_exon = $pair[1];
	
	  #print STDERR "\n";
	  $self->transfer_supporting_evidence($source_exon, $target_exon)
	}
      }
    } # end of this transcript

    # check to see if we can transfer evidence from rejected single exon transcripts
    foreach my $reject(values %single_exon_rejects){
      my @exons = @{$reject->get_all_Exons};

      # is there any supporting evidence?
      my @sf = @{$exons[0]->get_all_supporting_features};

      foreach my $stored_exon(values %exonhash){
	
	if($exons[0]->overlaps($stored_exon)){
	  # note that we could end up with bizarre situations of a single exon transcript overlapping 
	  # two exons in a multi exon transcript, so the supporting evidence would be transferred in 
	  # entirety to both exons.
	  $self->transfer_supporting_evidence($exons[0], $stored_exon);
	
	}
      }
    }

  } #end CLUSTER

  return @newtran;
}

############################################################

sub _bin_sort_transcripts{

  my ($self,@transcripts) = @_;

  my %lengths;
  my $max_num_exons = 0;

  # sizehash holds transcript length - based on sum of exon lengths
  my %sizehash;

  # orfhash holds orf length - based on sum of translateable exon lengths
  my %orfhash;

  my %tran2orf;
  my %tran2length;

  # keeps track of to which transcript(s) this exon belong
  my %exon2transcript;

  foreach my $tran (@transcripts) {

    # keep track of number of exons in multiexon transcripts
    my @exons = @{ $tran->get_all_Exons };
    if(scalar(@exons) > $max_num_exons){
      $max_num_exons = scalar(@exons);
    }

    # total exon length
    my $length = 0;
    foreach my $e ( @exons ){
      $length += $e->end - $e->start + 1;
      push ( @{ $exon2transcript{ $e } }, $tran );
    }
    $sizehash{$tran} = $length;
    $tran2length{ $tran } = $length;

    # now for ORF length
    $length = 0;
    foreach my $e(@{$tran->get_all_translateable_Exons}){
      $length += $e->end - $e->start + 1;
    }
    $tran2orf{ $tran } = $length;
    push(@{$orfhash{$length}}, $tran);
  }

  # VAC 15/02/2002 sort transcripts based on total exon length - this
  # introduces a problem - we can (and have) masked good transcripts
  # with long translations in favour of transcripts with shorter
  # translations and long UTRs that are overall slightly longer. This
  # is not good.
  # better way? hold both total exon length and length of translateable exons. Then sort:
  # long translation + UTR > long translation no UTR > short translation + UTR > short translation no UTR
  ##########
  @transcripts = ();
  ##########

  # sort first by orfhash{'length'}
  my @orflengths = sort {$b <=> $a} (keys %orfhash);

  # strict sort by translation length is just as wrong as strict sort by UTR length
  # bin translation lengths - 4 bins (based on 25% length diff)? 10 bins (based on 10%)?
  my %orflength_bin;
  my $numbins = 4;
  my $currbin = 1;

  foreach my $orflength(@orflengths){
    last if $currbin > $numbins;
    my $percid = ($orflength*100)/$orflengths[0];
    if ($percid > 100) { 
      $percid = 100; 
    }
    my $currthreshold = $currbin * (100/$numbins);
    $currthreshold    = 100 - $currthreshold;

    if($percid <$currthreshold) { 
      $currbin++; 
    }
    my @tmp = @{$orfhash{$orflength}};
    push(@{$orflength_bin{$currbin}}, @{$orfhash{$orflength}});
  }

  # now, foreach bin in %orflengthbin, sort by exonlength
  $currbin = 1;
 EXONLENGTH_SORT:
  while( $currbin <= $numbins){
    if(!defined $orflength_bin{$currbin} ){
      $currbin++;
      next EXONLENGTH_SORT;
    }	
    my @sorted_transcripts = sort { $sizehash{$b} <=> $sizehash{$a} } @{$orflength_bin{$currbin}};
    push(@transcripts, @sorted_transcripts);

    $currbin++;
  }

  # post filter - we want blessed transcripts to be at the top of the list so that 
  # identical build transcripts get properly pruned.
#  my %blessed_genetypes;
#  foreach my $bgt(@{$GB_BLESSED_GENETYPES}){
#    $blessed_genetypes{$bgt->{'type'}} = 1;
#  }
#
#  my @blessed_transcripts;
#  my @heathen_transcripts;
#  foreach my $transcript(@transcripts){
#
#    if(exists $blessed_genetypes{$transcript->type}){
#      push(@blessed_transcripts, $transcript);
#    }
#    else{
#      push(@heathen_transcripts, $transcript);
#    }
#  }
#    @transcripts = ();
#    push(@transcripts, @blessed_transcripts);
#    push(@transcripts, @heathen_transcripts);


  #test
  my $debug = 0;
  if ($debug == 1 ){

    print STDERR "2.- sorted transcripts:\n";

    foreach my $tran (@transcripts){
      if ( $tran->dbID ){
	print STDERR $tran->dbID." ";
      }
      print STDERR "orf_length: $tran2orf{$tran}, total_length: $tran2length{$tran}\n";
      my @exons = @{$tran->get_all_Exons};
      if ( $exons[0]->strand == 1 ){
	@exons = sort { $a->start <=> $b->start } @exons;
      }
      else{
	@exons = sort { $b->start <=> $a->start } @exons;
      }
      foreach my $exon ( @{$tran->get_all_Exons} ){
	print "  ".$exon->start."-".$exon->end." ".( $exon->end - $exon->start + 1)." phase: ".$exon->phase." end_phase ".$exon->end_phase." strand: ".$exon->strand."\n";
      }
    }
  }
  return @transcripts;
}

############################################################

=head2 cluster_into_Genes

    Example :   my @genes = $self->cluster_into_Genes(@transcripts);
Description :   it clusters transcripts into genes according to exon overlap.
                It will take care of difficult cases like transcripts within introns.
                It also unify exons that are shared among transcripts.
    Returns :   a beautiful list of geen objects
    Args    :   a list of transcript objects

=cut

sub cluster_into_Genes{
  my ($self, @transcripts_unsorted) = @_;
  
  my $num_trans = scalar(@transcripts_unsorted);
  print STDERR "clustering $num_trans transcripts into genes\n";

  
  # flusold genes
  #$self->flush_Genes;

  my @transcripts = sort { $a->start <=> $b->start ? $a->start <=> $b->start  : $b->end <=> $a->end } @transcripts_unsorted;
  my @clusters;

  # clusters transcripts by whether or not any exon overlaps with an exon in 
  # another transcript (came from original prune in GeneBuilder)
  foreach my $tran (@transcripts) {

    my @matching_clusters;
  CLUSTER: 
    foreach my $cluster (@clusters) {
      foreach my $cluster_transcript (@$cluster) {
        if ($tran->end  >= $cluster_transcript->start &&
            $tran->start <= $cluster_transcript->end) {

          foreach my $exon1 (@{$tran->get_all_Exons}) {
  	  foreach my $cluster_exon (@{$cluster_transcript->get_all_Exons}) {
              if ($exon1->overlaps($cluster_exon) && $exon1->strand == $cluster_exon->strand) {
                push (@matching_clusters, $cluster);
                next CLUSTER;
              }
            }
          }

        }
      }
    }
    
    if (scalar(@matching_clusters) == 0) {
      my @newcluster;
      push(@newcluster,$tran);
      push(@clusters,\@newcluster);
    } 
    elsif (scalar(@matching_clusters) == 1) {
      push @{$matching_clusters[0]}, $tran;
      
    } 
    else {
      # Merge the matching clusters into a single cluster
      my @new_clusters;
      my @merged_cluster;
      foreach my $clust (@matching_clusters) {
        push @merged_cluster, @$clust;
      }
      push @merged_cluster, $tran;
      push @new_clusters,\@merged_cluster;
      # Add back non matching clusters
      foreach my $clust (@clusters) {
        my $found = 0;
      MATCHING: 
	foreach my $m_clust (@matching_clusters) {
          if ($clust == $m_clust) {
            $found = 1;
            last MATCHING;
          }
        }
        if (!$found) {
          push @new_clusters,$clust;
        }
      }
      @clusters =  @new_clusters;
    }
  }
  
  # safety and sanity checks
  $self->check_Clusters(scalar(@transcripts), \@clusters);
  
  # make and store genes
  print STDERR scalar(@clusters)." created, turning them into genes...\n";
  my @genes;
  foreach my $cluster(@clusters){
    my $count = 0;
    my $gene = new Bio::EnsEMBL::Gene;
    foreach my $transcript (@$cluster){
      $gene->add_Transcript($transcript);
    }
    push( @genes, $gene );
  }
  return @genes;
}

############################################################

sub check_Clusters{
  my ($self, $num_transcripts, $clusters) = @_;
  #Safety checks
  my $ntrans = 0;

my $cluster_num = 0;

  my %trans_check_hash;
  foreach my $cluster (@$clusters) {
    $ntrans += scalar(@$cluster);

    foreach my $trans (@$cluster) {

      if (defined($trans_check_hash{$trans})) {
        $self->throw("Transcript " . $trans->dbID . " added twice to clusters\n");
#        $self->warn("Transcript " . $trans->dbID . " added twice to clusters\n");
      }
      $trans_check_hash{$trans} = 1;
    }
    if (!scalar(@$cluster)) {
      $self->throw("Empty cluster");
    }
  }
  if ($ntrans != $num_transcripts) {
    $self->throw("Not all transcripts have been added into clusters $ntrans and " . $num_transcripts. " \n");
  } 
  #end safety checks
  return;
}


############################################################

sub transcript_high{
  my ($self,$tran) = @_;
  my $high;
  #$tran->sort;
  if ( $tran->start_Exon->strand == 1){
    $high = $tran->end_Exon->end;
  }
  else{
    $high = $tran->start_Exon->end;
  }
  return $high;
}

############################################################

sub transcript_low{
  my ($self,$tran) = @_;
  my $low;
  #$tran->sort;
  if ( $tran->start_Exon->strand == 1){
    $low = $tran->start_Exon->start;
  }
  else{
    $low = $tran->end_Exon->start;
  }
  return $low;
}

############################################################

sub by_transcript_high {
  my $alow;
  my $blow;

  my $ahigh;
  my $bhigh;
  
  # alow and ahigh are the left most and right most coordinates for transcript $a 
  if ($a->start_Exon->strand == 1) {
    $alow  = $a->start_Exon->start;
    $ahigh = $a->end_Exon->end;
  } 
  else {
    $alow  = $a->end_Exon->start;
    $ahigh = $a->start_Exon->end;
  }

  # blow and bhigh are the left most and right most coordinates for transcript $b 
  if ($b->start_Exon->strand == 1) {
    $blow  = $b->start_Exon->start;
    $bhigh = $b->end_Exon->end;
  } 
  else {
    $blow  = $b->end_Exon->start;
    $bhigh = $b->start_Exon->end;
  }

  # return the ascending comparison of the right-most coordinates if they're different
  if ($ahigh != $bhigh) {
    return $ahigh <=> $bhigh;
  } 
  # if they'r equal, return the ascending comparison of the left most coordinate
  else {
    return $alow <=> $blow;
  }
}



############################################################

sub prune_Exons {
  my ($self,$gene) = @_;
  
  my @unique_Exons; 
  
  # keep track of all unique exons found so far to avoid making duplicates
  # need to be very careful about translation->start_Exon and translation->end_Exon
  
  foreach my $tran (@{$gene->get_all_Transcripts}) {
    my @newexons;
    foreach my $exon (@{$tran->get_all_Exons}) {
      my $found;
      #always empty
    UNI:foreach my $uni (@unique_Exons) {
	if ($uni->start  == $exon->start  &&
	    $uni->end    == $exon->end    &&
	    $uni->strand == $exon->strand &&
	    $uni->phase  == $exon->phase  &&
	    $uni->end_phase == $exon->end_phase
	   ) {
	  $found = $uni;
	  last UNI;
	}
      }
      if (defined($found)) {
	push(@newexons,$found);
	if ($exon == $tran->translation->start_Exon){
	  $tran->translation->start_Exon($found);
	}
	if ($exon == $tran->translation->end_Exon){
	  $tran->translation->end_Exon($found);
	}
      } else {
	push(@newexons,$exon);
	push(@unique_Exons, $exon);
      }
    }          
    $tran->flush_Exons;
    foreach my $exon (@newexons) {
      $tran->add_Exon($exon);
    }
  }
  return $gene;
}

############################################################

=head2 get_Predictions

Description:  gets your favourite ab initio predictions (genscan,genefinder,fgenesh,genecooker...8-)
Returns : none, but $self->predictions is filled
Args    : none

=cut
  
sub get_Predictions {
  my ($self) = @_;
  my @checked_predictions;
  foreach my $prediction ( @{ $self->query->get_all_PredictionTranscripts } ){
    $prediction->type("ab-initio");
    Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Peptide( $prediction );
    unless ( Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_check_Transcript( $prediction, $self->query ) ){
      $self->warn("We let in a prediction with wrong phases!");
    }
    
    # now need to explicitly check intron sizes if required
    unless ( Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_check_introns( $prediction, $self->query ) ){
      $self->warn("Rejecting prediction with long introns!");
      next;
    }
    push ( @checked_predictions, $prediction );
  }
  $self->predictions(@checked_predictions);
}

############################################################

=head2 get_Similarities

 Title   : get_Similarities
 Usage   : $self->get_Similarities
 Function: gets similarity features for this region
 Returns : none, but $self->feature is filled
 Args    : none

=cut

sub get_Similarities {
  my ($self) = @_;

  my @features = @{ $self->query->get_all_SimilarityFeatures('',$GB_MIN_FEATURE_SCORE) };
  
  print STDERR "retrieved ".scalar(@features)." features\n";
  my %idhash;
  my @other_features;
  
  foreach my $feature (@features) {
      if ($feature->length > $GB_MIN_FEATURE_LENGTH) {
	  unless ( $idhash{ $feature->hseqname } ){
	      $idhash{ $feature->hseqname } = [];
	  }
	  if ($feature->isa("Bio::EnsEMBL::BaseAlignFeature")) {
	      push( @{ $idhash{ $feature->hseqname } }, $feature );
	  }
      }
      # any other feature is rejected
  }
  
  # test:
  #foreach my $id ( keys %idhash ){
  #    print STDERR "Id: $id ".scalar(@{$idhash{$id}})." features\n";
  #}
 
  my @pruned_features = $self->prune_features(\%idhash);

  $self->features(@pruned_features);
}
   

############################################################


=head2 merge
  Note(eae): very suspicious, I do not why it is doing that, and if it
             should be doing what I think, it is not doung correctly, so
             do not use for the time being, subjet to future revision.

 Description: wicked method that merges two or more homol features into one 
              if they are close enough together
  Returns   : nothing
  Args      : none

=cut

sub merge {
    my ($self,$feature_hash,$overlap,$query_gap,$homol_gap) = @_;
    
    $overlap   = 20  unless $overlap;
    $query_gap = 15  unless $query_gap;
    $homol_gap = 15  unless $homol_gap;
    
    my @mergedfeatures;
    
  ID:
    foreach my $id (keys %{ $feature_hash }) {
	
	my $count = 0;
	my @newfeatures;
	my @features = @{$feature_hash->{$id}};
	
	@features = sort { $a->start <=> $b->start} @features;
	unless ( @features ){
	    print STDERR "No features here for id: $id\n";
	    next ID;
	}
	while ( @features && !defined $features[0] ){
	    print STDERR "jumping an undefined feature\n";
	    shift @features;
	}
	
	# put the first feature in the new array;
	print STDERR "pushing feature[0]: ".$features[0]->gffstring."\n";
	push(@newfeatures,$features[0]);

	
	for (my $i=0; $i < $#features; $i++) {
	    my $id  = $features[$i]  ->id;
	    my $id2 = $features[$i+1]->id;
	    
	    print STDERR "Comparing:\n";
	    print STDERR "feature[$i]: ".$features[$i]->gffstring."\n";
	    print STDERR "feature[".($i+1)."]: ".$features[$i+1]->gffstring."\n";
	    
	    # First case is if start of next hit is < end of previous
	    if ( $features[$i]->end > $features[$i+1]->start && 
		 ($features[$i]->end - $features[$i+1]->start + 1) < $overlap) {
		print STDERR "overlap\n";
		
		if ($features[$i]->strand == 1) {
		    $newfeatures[$count]-> end($features[$i+1]->end);
		    $newfeatures[$count]->hend($features[$i+1]->hend);
		} 
		else {
		    $newfeatures[$count]-> end($features[$i+1]->end);
		    $newfeatures[$count]->hend($features[$i+1]->hstart);
		}
		
		# Take the max score
		if ($features[$i+1]->score > $newfeatures[$count]->score) {
		    $newfeatures[$count]->score($features[$i+1]->score);
		}
		
		if ($features[$i+1]->hstart == $features[$i+1]->hend) {
		    $features[$i+1]->strand($features[$i]->strand);
		}
		
		print STDERR "newfeature[$count]: ".$newfeatures[$count]->gffstring."\n";
		
	    }
	    # Allow a small gap if < $query_gap, $homol_gap
	    elsif (($features[$i]->end < $features[$i+1]->start) &&
		   abs($features[$i+1]->start - $features[$i]->end) <= $query_gap) {
		
		if ($features[$i]->strand eq "1") {
		    $newfeatures[$count]->end($features[$i+1]->end);
		    $newfeatures[$count]->hend($features[$i+1]->hend);
		} 
		else {
		    $newfeatures[$count]->end($features[$i+1]->end);
		    $newfeatures[$count]->hstart($features[$i+1]->hstart);
		}
		
		if ($features[$i+1]->score > $newfeatures[$count]->score) {
		    $newfeatures[$count]->score($features[$i+1]->score);
		}
		
		if ($features[$i+1]->hstart == $features[$i+1]->hend) {
		    $features[$i+1]->strand($features[$i]->strand);
		}
		print STDERR "newfeature[$count]: ".$newfeatures[$count]->gffstring."\n";
		
	    } 
	    else {
		# we can't extend the merged homologies so start a
		# new homology feature
		
		# first do the coords on the old feature
		if ($newfeatures[$count]->hstart > $newfeatures[$count]->hend) {
		    my $tmp = $newfeatures[$count]->hstart;
		    $newfeatures[$count]->hstart($newfeatures[$count]->hend);
		    $newfeatures[$count]->hend($tmp);
		}
		
		$count++;
		$i++;
		
		print STDERR "pushing feature[$i]: ".$features[$i]->gffstring."\n";
		push(@newfeatures,$features[$i]);
		$i--;
	    }
	}
	
	if ( @newfeatures ){
	    # Adjust the last new feature coords
	    if ($newfeatures[$#newfeatures]->hstart > $newfeatures[$#newfeatures]->hend) {
		my $tmp = $newfeatures[$#newfeatures]->hstart;
		$newfeatures[$#newfeatures]->hstart($newfeatures[$#newfeatures]->hend);
		$newfeatures[$#newfeatures]->hend($tmp);
	    }
	    
	    my @pruned = $self->prune_features(@newfeatures);
	    
	    push(@mergedfeatures,@pruned);
	}
	else{
	    print STDERR "No features merged\n";
	}
    }
    
    return @mergedfeatures;
   
}

############################################################

=head2 prune_features

 Description: prunes out duplicated features
 Returntype : array of Bio::EnsEMBL::SeqFeature
 Args       : array of Bio::EnsEMBL::SeqFeature

=cut
    
sub prune_features {
    my ($self,$feature_hash)  = @_;
    my @pruned;
 
  ID:
    foreach my $id (keys %{ $feature_hash }) {
	my @features = @{$feature_hash->{$id}};
	@features = sort {$a->start <=> $b->start} @features;
	
	unless ( @features ){
	    print STDERR "No features here for id: $id\n";
	    next ID;
	}
	while ( @features && !defined $features[0] ){
	    print STDERR "jumping an undefined feature\n";
	    shift @features;
	}
	
	my $prev = -1;
	
      FEATURE: 
	foreach  my $f (@features) {
	    if ($prev != -1 && $f->hseqname eq $prev->hseqname &&
		$f->start   == $prev->start &&
		$f->end     == $prev->end   &&
		$f->hstart  == $prev->hstart &&
		$f->hend    == $prev->hend   &&
		$f->strand  == $prev->strand &&
		$f->hstrand == $prev->hstrand) 
	    {
		#keep the one with highest score
		if ( $f->score > $prev->score ){
		    $prev->score( $f->score );
		}
		#print STDERR "pruning duplicated feature\n";
		#print STDERR "previous: ".$prev->gffstring."\n";
		#print STDERR "thisone : ".$f->gffstring."\n";
		next FEATURE;
	    } 
	    else {
		push(@pruned,$f);
		$prev = $f;
	    }
	}
    }
    return @pruned;
}

############################################################


############################################################
#
# GETSET METHODS
#
############################################################

# get/set method holding a reference to the db with genewise and combined genes
# this reference is set in Bio::EnsEMBL::Pipeline::RunnableDB::Gene_Builder

sub genes_db{
 my ($self,$genes_db) = @_;
 if ( $genes_db ){
   $self->{_genes_db} = $genes_db;
 }
 return $self->{_genes_db};
}

############################################################

sub combined_Transcripts {
    my ($self,@transcripts) = @_;

    if (!defined($self->{_genewise_andthelike_transcripts})) {
        $self->{_genewise_andthelike_transcripts} = [];
    }

    if (scalar @transcripts > 0) {
	push(@{$self->{_genewise_andthelike_transcripts}},@transcripts);
    }

    return @{$self->{_genewise_andthelike_transcripts}};
}

############################################################

=head2 my_genes

 Description: this holds and returns the genes that are produced after putting together genewise, combined and
              processed_supporte_ab_initio predictions and removing the redundant set, giving priority
              to long CDSs + UTR

=cut


sub my_genes {
  my ($self,@genes) = @_;
  
  unless($self->{_my_genes}){
    $self->{_my_genes} = [];
  }

  if (@genes){
    push(@{$self->{_my_genes}},@genes);
  }
  return @{$self->{_my_genes}};
}

############################################################

=head2 final_genes

 Descripton: this holds/returns the final genes produced after clustering transcripts and sharing common exons

=cut

sub final_genes{
  my ($self, @genes) = @_;
  
  if ( @genes ){
    push( @{$self->{_final_genes}}, @genes );
  }
  return @{$self->{_final_genes}};
}

############################################################

=head2 gene_types

 Description: get/set for the type(s) of genes (usually TGE_gw, similarity_genewise and combined_e2g genes) 
              to be used in the genebuilder they get set in new()
              Does not include the ab inition predictions
=cut

sub gene_types {
  my ($self,$type) = @_;

  if (defined($type)) {
     push(@{$self->{_gene_types}},$type);
  }

  return @{$self->{_gene_types}};
}
############################################################

=head2 predictions

 Description: get/set for the PredictionTranscripts. It is  set in new()

=cut

sub predictions {
  my ($self,@predictions) = @_;

  if(!$self->{_predictions}){
    $self->{_predictions} = [];
  }
  if ( @predictions ) {
     push(@{$self->{_predictions}},@predictions);
  }
  return @{$self->{_predictions}};
}

############################################################

sub features {
  my ($self,@features) = @_;
  
  if (!defined($self->{_feature})) {
    $self->{_feature} = [];
  }
  if ( scalar @features ) {
    push(@{$self->{_feature}},@features);
  }
  return @{$self->{_feature}};
}

############################################################

sub query {
  my ($self,$slice) = @_;
  
  if (defined($slice)) {
    $self->{_query} = $slice;
  }
  return $self->{_query};
}

#############################################################################
# 
# Printing routines
#
#############################################################################

sub print_ExonPairs {
  my ($self) = @_;
  
  foreach my $pair ($self->get_all_ExonPairs) {
    $self->print_ExonPair($pair);
  }
}

############################################################

sub print_ExonPair {
  my ($self,$pair) = @_;
  
  print STDERR $pair->exon1->gffstring."\n";
  print STDERR $pair->exon2->gffstring."\n";
  
  print(STDERR "\nExon Pair (splice - " . $pair->splice_seq->seq . ")\n");
  
  foreach my $ev ($pair->get_all_Evidence) {
    print(STDERR "   -  " . $ev->hseqname . "\t" . $ev->hstart . "\t" . $ev->hend . "\t" . $ev->strand . "\n");
  }
}

############################################################


=head2 split_transcript

 Title   : split_transcript 
 Usage   : my @splits = $self->split_transcript($transcript)
 Function: splits a transcript into multiple transcripts at long introns. Rejects single exon 
           transcripts that result. 
 Returns : @Bio::EnsEMBL::Transcript
 Args    : Bio::EnsEMBL::Transcript

=cut


sub split_transcript{
  my ($self, $transcript) = @_;
  $transcript->sort;
  my @split_transcripts   = ();

  if(!($transcript->isa("Bio::EnsEMBL::Transcript"))){
    $self->warn("[$transcript] is not a Bio::EnsEMBL::Transcript - cannot split");
    return (); # empty array
  }
  
  my $prev_exon;
  my $exon_added = 0;
  my $curr_transcript = new Bio::EnsEMBL::Transcript;
  my $translation     = new Bio::EnsEMBL::Translation;
  $curr_transcript->type($transcript->type);
  $curr_transcript->translation($translation);


EXON:   
  foreach my $exon( @{$transcript->get_all_Exons} ){
    
    $exon_added = 0;
      # is this the very first exon?
    if($exon == $transcript->start_Exon){
      $prev_exon = $exon;      
      # set $curr_transcript->translation start and start_exon
      $curr_transcript->add_Exon($exon);
      $exon_added = 1;
      $curr_transcript->translation->start_Exon($exon);
      $curr_transcript->translation->start($transcript->translation->start);
      push(@split_transcripts, $curr_transcript);
      next EXON;
    }
    
    if ($exon->strand != $prev_exon->strand){
      return (); # empty array
    }

    # We need to start a new transcript if the intron size between $exon and $prev_exon is too large
    my $intron = 0;
    if ($exon->strand == 1) {
      $intron = abs($exon->start - $prev_exon->end + 1);
    } else {
      $intron = abs($exon->end   - $prev_exon->start + 1);
    }
    
    if ($intron > $GB_GENSCAN_MAX_INTRON) {
      $curr_transcript->translation->end_Exon($prev_exon);
      # need to account for end_phase of $prev_exon when setting translation->end
      $curr_transcript->translation->end($prev_exon->end - $prev_exon->start + 1 - $prev_exon->end_phase);
      
      # start a new transcript 
      my $t  = new Bio::EnsEMBL::Transcript;
      my $tr = new Bio::EnsEMBL::Translation;
      $t->type($transcript->type);
      $t->translation($tr);

      # add exon unless already added, and set translation start and start_exon
      $t->add_Exon($exon) unless $exon_added;
      $exon_added = 1;

      $t->translation->start_Exon($exon);

      if ($exon->phase == 0) {
	$t->translation->start(1);
      } elsif ($exon->phase == 1) {
	$t->translation->start(3);
      } elsif ($exon->phase == 2) {
	$t->translation->start(2);
      }

      # start exon always has phase 0
      $exon->phase(0);      

      # this new transcript becomes the current transcript
      $curr_transcript = $t;

      push(@split_transcripts, $curr_transcript);
    }

    if($exon == $transcript->end_Exon){
      # add it unless already added
      $curr_transcript->add_Exon($exon) unless $exon_added;
      $exon_added = 1;

      # set $curr_transcript end_exon and end
      $curr_transcript->translation->end_Exon($exon);
      $curr_transcript->translation->end($transcript->translation->end);
    }

    else{
      # just add the exon
      $curr_transcript->add_Exon($exon) unless $exon_added;
    }
    
    # this exon becomes $prev_exon for the next one
    $prev_exon = $exon;

  }

  # discard any single exon transcripts
  my @t = ();
  my $count = 1;
  
  foreach my $st(@split_transcripts){
    $st->sort;
    my @ex = @{$st->get_all_Exons};
    if(scalar(@ex) > 1){
      $count++;
      push(@t, $st);
    }
  }

  return @t;

}

############################################################

=head2 validate_transcript

 Title   : validate_transcript 
 Usage   : my @valid = $self->validate_transcript($transcript)
 Function: Validates a transcript - rejects if mixed strands, splits if long introns
 Returns : @Bio::EnsEMBL::Transcript
 Args    : Bio::EnsEMBL::Transcript

=cut

sub validate_transcript{
  my ($self, $transcript) = @_;
  my @valid_transcripts;

  my $valid = 1;
  my $split = 0;

  # check exon phases:
  my @exons = @{$transcript->get_all_Exons};
  $transcript->sort;
  for (my $i=0; $i<(scalar(@exons-1)); $i++){
    my $end_phase = $exons[$i]->end_phase;
    my $phase    = $exons[$i+1]->phase;
    if ( $phase != $end_phase ){
      $self->warn("rejecting transcript with inconsistent phases( $phase <-> $end_phase) ");
      return undef;
    }
  }
  

  my $previous_exon;
  foreach my $exon (@exons){
    if (defined($previous_exon)) {
      my $intron;
      
      if ($exon->strand == 1) {
	$intron = abs($exon->start - $previous_exon->end + 1);
      } else {
	$intron = abs($exon->end   - $previous_exon->start + 1);
      }
      
      if ($intron > $GB_GENSCAN_MAX_INTRON) {
	print STDERR "Intron too long $intron  for transcript " . $self->transcript_id($transcript) . "\n";
	$split = 1;
	$valid = 0;
      }
      
      if ($exon->strand != $previous_exon->strand) {
	print STDERR "Mixed strands for gene " . $self->transcript_id($transcript) . "\n";
	$valid = 0;
	return;
      }
    }
    $previous_exon = $exon;
  }
  
  if ($valid) {
    push(@valid_transcripts,$transcript);
  }
  elsif ($split){
    # split the transcript up.
    my @split_transcripts = $self->split_transcript($transcript);
    push(@valid_transcripts, @split_transcripts);
  }
  return @valid_transcripts;
}

############################################################

=head2 transfer_supporting_evidence

 Title   : transfer_supporting_evidence
 Usage   : $self->transfer_supporting_evidence($source_exon, $target_exon)
 Function: Transfers supporting evidence from source_exon to target_exon, 
           after checking the coordinates are sane and that the evidence is not already in place.
 Returns : nothing, but $target_exon has additional supporting evidence

=cut

sub transfer_supporting_evidence{
  my ($self, $source_exon, $target_exon) = @_;
  
  my @target_sf = @{$target_exon->get_all_supporting_features};
  #  print "target exon sf: \n";
  #  foreach my $tsf(@target_sf){ print STDERR $tsf; $self->print_FeaturePair($tsf); }
  
  #  print "source exon: \n";
 
  # keep track of features already transferred, so that we do not duplicate
  my %unique_evidence;
  my %hold_evidence;

 SOURCE_FEAT:
  foreach my $feat ( @{$source_exon->get_all_supporting_features}){
    next SOURCE_FEAT unless $feat->isa("Bio::EnsEMBL::FeaturePair");
    
    # skip duplicated evidence objects
    next SOURCE_FEAT if ( $unique_evidence{ $feat } );
    
    # skip duplicated evidence 
    if ( $hold_evidence{ $feat->hseqname }{ $feat->start }{ $feat->end }{ $feat->hstart }{ $feat->hend } ){
      #print STDERR "Skipping duplicated evidence\n";
      next SOURCE_FEAT;
    }

    #$self->print_FeaturePair($feat);
    
  TARGET_FEAT:
    foreach my $tsf (@target_sf){
      next TARGET_FEAT unless $tsf->isa("Bio::EnsEMBL::FeaturePair");
      
      if($feat->start    == $tsf->start &&
	 $feat->end      == $tsf->end &&
	 $feat->strand   == $tsf->strand &&
	 $feat->hseqname eq $tsf->hseqname &&
	 $feat->hstart   == $tsf->hstart &&
	 $feat->hend     == $tsf->hend){
	
	#print STDERR "feature already in target exon\n";
	next SOURCE_FEAT;
      }
    }
    #print STDERR "from ".$source_exon->dbID." to ".$target_exon->dbID."\n";
    #$self->print_FeaturePair($feat);
    $target_exon->add_supporting_features($feat);
    $unique_evidence{ $feat } = 1;
    $hold_evidence{ $feat->hseqname }{ $feat->start }{ $feat->end }{ $feat->hstart }{ $feat->hend } = 1;
  }
}

############################################################

sub print_FeaturePair{
  my ($self, $fp) = @_;
  return unless $fp->isa("Bio::EnsEMBL::FeaturePair");
  print STDERR $fp;
  print STDERR $fp->seqname . " " .
    $fp->start . " " .
      $fp->end . " " .
	$fp->strand . " " .
	  $fp->hseqname . " " .
	    $fp->hstart . " " .
	      $fp->hend . "\n";
}

############################################################

1;
