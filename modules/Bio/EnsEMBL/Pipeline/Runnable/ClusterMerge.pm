#
# Written by Eduardo Eyras
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::ClusterMerge

=head1 SYNOPSIS

    my $cluster_merge = Bio::EnsEMBL::Pipeline::Runnable::CusterMerge->new(
									   -transcripts => \@transcripts,
									   -comparison_level => 3,
									   -intron_mismatch  => 10,
									   -splice_mismatch  => 40,
									   -minimum_order    => 2,
									   );
   $cluster_merge->run;

    one can retrieve the lists of transcripts that can merge:
    my @lists = $cluster_merge->sub_clusters;
   
    where @lists is a list of listrefs, each one cointaining a list of transcript objects

    one can retrieve the transcripts already merged:
    my @merged_transcripts = $cluster_merge->output;



=head1 DESCRIPTION

ClusterMerge takes a set of transcripts, it clusters them first and them
create sets of transcripts that can all merge with each other. This merging can be
exact, so that exon boundaries must match exactily, or fuzzy where 
exon boundaries do not necessarily coincide. The latter case could be useful for
ests, whereas the former case is advisable with full length cdnas.
The output is given in the form of the resulting transcripts.


Options:

exon_match = BOOLEAN ------------> TRUE if we want both transcripts to match 1-to-1 all their exons 
                                   Appropriate for comparison levels 1, 2 and 3

splice_mismatch = INT -----------> Maximum number of bases mismatch that we allow in the internal splice sites
                                   Alignment programs are sometimes not able to resolve some ambiguities in
                                   the splice sites, this might help to identify truly equivalent splice sites.

intron_mismatch = INT -----------> Maximum number of bases that we consider for an intron to be non-real.
                                   Any intron of this size or smaller will be allowed to be 
                                   merged with exons covering it.
                                   The reason for this is that we do not expect very very small intron to be
                                   real

internal_splice_overlap  --------> (default = 0 ) 
                                   number of base pairs (N) we allow an external exon overlap
                                   an intron in another transcript:

                                                     |--N--|
                                    ######-------##########
                                   #######-------####-------------#######

minimum_order -------------------> positive INT, it is the minimum number of ESTs that must be in a list
                                   to use that list for merging. The default is 1.

comparison_level = INT ----------> There are currently 5 comparison levels:

 
 1 --> strict: exact exon matching (unrealistic). This one does not use any other parameteres passed in

                 #####-----#####-----#####
                 #####-----#####-----#####

   2 --> allow edge exon mismatches. This one will use the parameter 'exon_match' if is defined
   
                 #####-----#####-----#######
       #####-----#####-----#####-----#####


   3 ---> allow internal mismatches. This one can use the parameters 'exon_match' and 'splice_mismatch' if they are defined

                 #####-----######----#######
       #####-----######----####------#####

   4 ---> allow intron mismatches. This one can use all the parameters if they have been defined.

                 ################----#######
       #####-----######----####------#####

  
   5 ---> loose match. It allows intron mismatches if so desired. There is no limitation on
          the number of mismatches at the splice-sites.

                #################----#######   OR                 #######------####----#######  
       #####-----######----####------#####                #####-----######----####------#####



=head1 AUTHOR

eae@sanger.ac.uk

=head1 CONTACT

ensembl-dev@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::Runnable::ClusterMerge;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster;
use Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptComparator;
use Bio::EnsEMBL::Pipeline::GeneComparison::ScoreModel;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Bio::EnsEMBL::Pipeline::Node;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  
  my( $transcripts, $comparison_level, $splice_mismatch, $intron_mismatch, $exon_match, $minimum_order, $internal_splice_overlap, $use_score, $speed_up) 
      = $self->_rearrange([qw(
			      TRANSCRIPTS
			      COMPARISON_LEVEL
			      SPLICE_MISMATCH
			      INTRON_MISMATCH
			      EXON_MATCH
			      MINIMUM_ORDER
			      INTERNAL_SPLICE_OVERLAP
			      USE_SCORE
			      SPEED_UP
			      )], 
			  @args);
  
  unless( $transcripts ){
    $self->warn("No transcripts passed it. We cannot go further");
    exit(0);
  }
  $self->input_transcripts(@{$transcripts});
  
  # default values are
  # $comparison_level = 4,
  # $splice_mismatch  = 10,
  # $intron_mismatch  = 7,
  # $exon_match       = 0,
  # $internal_splice_overlap = 0,
  
  # to get some blah,blah
  $self->verbose(0);

  $use_score = 1;
  if ( $use_score ){
      $self->use_score($use_score);
  }
  

  if ( $comparison_level ){
    $self->_comparison_level($comparison_level);
    print STDERR "comparison_level = ".$self->_comparison_level."\n";
  }
  else{
    print STDERR "defaulting comparison_level to 2\n";
    $self->_comparison_level(2 );
  }
  
  if ($minimum_order){
    $self->_minimum_order($minimum_order);
  }
  else{
    $self->_minimum_order(1);
  }

  if ( defined $speed_up ){
    $self->_speed_up( $speed_up );
  }
  else{
    $self->_speed_up( 0 );
  }
  
  if ( defined $splice_mismatch ){
    $self->_splice_mismatch( $splice_mismatch );
    #print STDERR "splice mismatch ".$self->_splice_mismatch."\n";
  }
  else{
    $self->_splice_mismatch( 10 );
    #print STDERR "defaulting splice mismatch to ".$self->_splice_mismatch."\n";
  }

  if ( defined $intron_mismatch ){
    $self->_intron_mismatch( $intron_mismatch );
    #print STDERR "intron mismatch ".$self->_intron_mismatch."\n";
  }
  else{
    $self->_intron_mismatch( 7 );
    #print STDERR "defaulting intron mismatch to ".$self->_intron_mismatch."\n";
  }

  if ( defined $internal_splice_overlap ){
    $self->_internal_splice_overlap( $internal_splice_overlap );
    #print STDERR "internal_splice_overlap ".$self->_internal_splice_overlap."\n";
  }
  else{
    $self->_internal_splice_overlap( 0 );
    #print STDERR "defaulting internal_splice_overlap = 0\n";
  }

  if ( defined $exon_match ){
    $self->_exon_match( $exon_match );
    #print STDERR "exon match ".$self->_exon_match."\n";
  }
  else{
    $self->_exon_match( 0 );
    #print STDERR "defaulting exon match to ".$self->_exon_match."\n";
  }
  
  return $self;
}


############################################################
#
# METHODS FOR COLLECTING ALL THE NON-REDUNDANT LISTS
#
############################################################

sub solutions{
  my ($self,$node) = @_;
  my @all_solutions;

  my $verbose = $self->verbose;

  print STDERR "solution at node: ".$node->transcript->dbID."\n" if $verbose;
  if ( $node->is_candidate ){
    print STDERR "node: ".$node->transcript->dbID." is a candidate\n" if $verbose;
  }
  
  ############################################################
  # if node has extension parents
  my @extension_parents;
  if ($node->is_candidate){
    @extension_parents = @{$node->candidate_extension_parents};
  }
  else{
    @extension_parents = @{$node->extension_parents};
  }
  if ( @extension_parents ){
      
      if ( scalar(@extension_parents) > 1 ){
	  my @paps = map { $_->transcript->dbID } @extension_parents;
	  print STDERR "node ".$node->transcript->dbID." has ".scalar( @extension_parents)." parents: @paps\n" if $verbose;
      }
      
    PARENT:
      foreach my $extension_parent ( @extension_parents ){
	  
	  #if  ( $extension_parent->collected ){
	  #	next PARENT;
	  #}
	  my $solutions = $self->solutions( $extension_parent );
	  
	  foreach my $solution ( @$solutions ){
	    # add itself to the list
	    push ( @$solution, $node );
	  }
	  
	  # if the node $node has inclusion children
	  if ( @{$node->inclusion_children} ){
	    my %added;
	    foreach my $inclusion_child ( @{$node->inclusion_children} ){
	      
	      ############################################################
	      # add the inclusion children recursively
	      # from this child unless this is included as well in the extension parent
	      # We only need checking in the first generation
	      ############################################################
	      unless ( $added{$inclusion_child} 
		       ||
		       $self->compare( $inclusion_child, $extension_parent) eq 'inclusion' 
		     ){
		my $list = [];
		$self->collect_inclusion_children( $inclusion_child , $list);
		$added{$inclusion_child}++;
		foreach my $solution ( @$solutions ){
		  push ( @$solution, @$list );
		}
	      }
	    }
	  }
	  
	  push (@all_solutions, @$solutions);
	}
    }
  else{
    print STDERR "no ext-parent - collecting inclusion-children and itself\n" if $verbose;
    my $list = [];
    $self->collect_inclusion_children($node,$list);
    #$node->collected(1);
    return [$list];
  }

  
  if ( @all_solutions ){
      return \@all_solutions;
  }
  else{
      my $list = [ $node ];
      #$node->collected(1);
      return [ $list ];
  }
}

############################################################

sub collect_inclusion_children{
  my ($self, $node, $list) = @_;

  my $verbose = $self->verbose();
  my %seen;
  print STDERR "collecting inclusion tree from node ".$node->transcript->dbID."\n" if $verbose;
  my $generation = 1;
  if ( @{$node->inclusion_children} ){
    my @next_generation = @{$node->inclusion_children};
    
    while( @next_generation ){
      push( @$list, @next_generation);
      print STDERR "generation $generation: ".scalar(@next_generation)."\n" if $verbose;
      $generation++;
      my @new_generation = ();
      while( @next_generation ){
	my $child = shift @next_generation;
	next if $seen{$child};
	$seen{$child} = 1;
	if ( @{$child->inclusion_children} ){
	  push( @new_generation, @{$child->inclusion_children} );
	}
      }
      @next_generation = @new_generation;
      @new_generation = ();
    }
  }
  push ( @$list, $node );
  return;
}			       
############################################################

sub get_transcript_string{
  my ($self, $t) = @_;
  my $string;
  my @exons = sort { $a->start <=> $b->start } @{$t->get_all_Exons};
  while (@exons){
    my $exon = shift @exons;
    $string .= $exon->start."_".$exon->end."_";
  }
  return $string;
}

############################################################
#
# LINKING ALGORITHM
#
############################################################

sub link_Transcripts{
  my ($self,$transcript_clusters) = @_;
  
  my @all_solutions;  

  my $verbose  = $self->verbose;
  my %seen_transcripts;

  # look in each cluster
 CLUSTER:
  foreach my $cluster ( @{ $transcript_clusters} ){
    
    print STDERR "~~~ CLUSTER ~~~\n" if $verbose;

    # matrix for cacheing the comparison results
    my %overlap_matrix;
    $self->matrix(\%overlap_matrix);
    
    # get the transcripts in this cluster
    my @transcripts = @{$cluster->get_Transcripts};
    
    ############################################################
    # sort the transcripts by the left most coordinate in ascending order
    # for equal value, order by the right most coordinate in descending order
    ############################################################
    @transcripts = sort { my $result = ( $self->transcript_low($a) <=> $self->transcript_low($b) );
			  if ($result){
			      return $result;
			  }
			  else{
			      return ( $self->transcript_high($b) <=> $self->transcript_high($a) );
			  }
		      } @transcripts;
    
    # just to track the sets
    # take 10 in each cluster
    #my @reduced_cluster;
    #my $count = 0;
    #while ( @transcripts && $count<10 ){
    #push(@reduced_cluster,shift @transcripts);
    #$count++;
    #}
    #@transcripts = @reduced_cluster;
    

    ############################################################
    # search all the trees
    ############################################################
    my $all_the_leaves;  
    my @candidates;
  
  TRANSCRIPT:
    foreach my $transcript ( @transcripts ){
	
      if ( $self->_speed_up ){
	my $string = $self->get_transcript_string($transcript);
	if ( $seen_transcripts{$string} ){
	  print STDERR "seen string $string\n";
	  next TRANSCRIPT;
	}
	$seen_transcripts{ $string } = 1;
      }
      
	# put this in a node
	my $newnode = Bio::EnsEMBL::Pipeline::Node->new();
	$newnode->transcript($transcript);
	print STDERR "newnode: ".$newnode->transcript->dbID."\n" if $verbose;
	
	# initialize the tags
	$self->flush_visited_tags;
	
	############################################################
	# the first element is always a new leaf
	unless ( $all_the_leaves){
	    push ( @$all_the_leaves, $newnode );
	    next TRANSCRIPT;
	}
	
	############################################################
	# check the node against all the leaves
	my $result = $self->check_tree($newnode,$all_the_leaves);
	unless ( $result eq 'placed'){
	  if ( $newnode->is_included ){
	    print STDERR "####### WARNING: adding leaf that is included in other node\n" if $verbose;
	  }
	  print STDERR "adding new leaf: ".$newnode->transcript->dbID."\n" if $verbose;
	  push( @$all_the_leaves, $newnode);
	  
	}
	############################################################
	# check whether this node is a candidate for a missing link:
	#
	#                 zN <- ...<- z1 <- x
	#                  |                |
	#                 \|/              \|/
	#                 wN <= ...<= w1 <= y   
	# 
	# wN is a potential missing link because it is not leaf but
	# it could generate a new solution
	if ( $result eq 'placed' ){
	  if ( $newnode->is_included && @{$newnode->extension_parents} ){
	    
	    ############################################################
	    # every extension parent must account for a triangle, otherwise we have a candidate
	    unless( $self->check_triangle( $newnode ) ){
	      
	      ############################################################
	      # check that $newnode is not already in the leaves
	      my $is_leaf = 0;
	      foreach my $leaf (@{$all_the_leaves}){
		if ( $leaf == $newnode ){
		  $is_leaf++;
		  last;
		}
	      }
	      unless( $is_leaf ){
		push( @candidates, $newnode );
	      }
	    }
	  }	    
	}
	############################################################
	# check the leaves - this is a bit redundant but it is a sanity-check
	$self->_check_leaves( $all_the_leaves );
      
      } # end of TRANSCRIPT
    # clean up some memory
    @transcripts = ();
    
    # check #   
    #if ($verbose){
    #      foreach my $leaf ( @$all_the_leaves ){
    #	print STDERR "leaf:\n";
    #	Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript($leaf->transcript);
    #	foreach my $parent ( @{$leaf->extension_parents} ){
    #	  print STDERR "extension parent:\n";
    #	  Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript($parent->transcript);
    #	}
    #      }
    #    }
    
    ############################################################
    # get rid of the candidates which have been extended:
    my @final_candidates;
    foreach my $candidate ( @candidates ){
      unless ( $candidate->is_extended ){
	push( @final_candidates, $candidate );
      }
    }

    ############################################################
    # recover the lists from the trees in this cluster:
    ############################################################
    my @solutions;
    my @t = map {$_->transcript->dbID} @$all_the_leaves;
    print STDERR "**** LEAVES: @t\n" if $verbose;
    #print STDERR "******************************************************\n";
    
    my @final_leaves = ( @$all_the_leaves, @final_candidates );
    
    # clean up some memory;
    @$all_the_leaves = ();
    @final_candidates = ();

    foreach my $leaf ( @final_leaves ){
	print STDERR "finding solutions from leaf ".$leaf->transcript->dbID."\n" if $verbose;
	my @sol = @{$self->solutions($leaf)};
	print STDERR scalar(@sol)." solutions found\n" if $verbose;
	push (@solutions, @{$self->solutions($leaf)} );
    }
    # @solutions is a list of listrefs holding the different non-redundant lists
    push ( @all_solutions, @solutions );
    
    # clean up some memory
    @solutions = ();
    
  } # end of CLUSTER  

  # clean up some memory
  @{ $transcript_clusters} = ();
  
  ############################################################
  # convert the Bio::EnsEMBL::Pipeline::Node objects
  # back into Bio::EnsEMBL::Transcript objects
  ############################################################
  my @lists;
  foreach my $solution ( @all_solutions ){
    my @tmp = map { $_->transcript } @$solution;
    push (@lists, \@tmp );
  }  
  # clean up some memory
  @all_solutions = ();
  return @lists;
}

############################################################

# this method is only for the main leaves
# it rejects those that are either included or
# extended

sub _check_leaves{
    my ($self,$leaves) = @_;
    my $copy_leaves;
    my $verbose = $self->verbose;

    foreach my $leaf (@$leaves){
	if ( $leaf->is_included || $leaf->is_extended ){
	    print STDERR "rejecting leaf ".$leaf->transcript->dbID." for being included/extended\n" if $verbose;
	    next;
	}
	push( @$copy_leaves, $leaf );
    }
    @$leaves = @$copy_leaves;
    $copy_leaves = [];
}


############################################################
# Every extension parent must account for a triangle, otherwise we have a candidate
#
# This method basically checks for situations like this:
#
#                     z       s              s     z <=== w
#                     |\      |              |    /     //
#                     |  \    |      or      |  /     //
#                     |    \  |              |/     //
#                     x ====> y              y <==//
#
# if every extension parent of y can be put into a triangle of either of these
# two forms, the node y cannot be the leaf of a 'hidden' branch.
# In the cartoons above, node s represents a hidden branch, and y has no triangle
# with this extension parent. As a consequence, y is labelled as candidate
# and later one, when retrieving the solutions, the branch y --- s ---... 
# will be included.
#

sub check_triangle{
    my ($self, $node ) = @_;

    my $verbose = $self->verbose();
    print STDERR "checking for triangle for node ".$node->transcript->dbID."\n" if $verbose;
    unless ( @{$node->extension_parents} ){
      print STDERR "No extension parents!!?\n";
      return 1;
    }
    my $is_candidate = 0;
  EXT_PARENT:
    foreach my $extension_parent ( @{$node->extension_parents} ){
      print STDERR "extension parent: ".$extension_parent->transcript->dbID."\n" if $verbose;
      my $found_triangle = 0;
   
    INCL_PARENT:
      foreach my $inclusion_parent ( @{$node->inclusion_parents} ){	
	my $incl_id = $inclusion_parent->transcript->dbID;
	print STDERR "inclusion parent: ".$incl_id."\n" if $verbose;
	
      OTHER_EXT_PARENT:
	foreach my $other_parent ( @{$inclusion_parent->extension_parents} ){
	  my $other_ext_id = $other_parent->transcript->dbID;
	  print STDERR "extension parent of $incl_id: ".$other_ext_id."\n" if $verbose;
	  if ( $other_parent == $extension_parent ){
	    print STDERR "found triangle\n" if $verbose;
	    $found_triangle = 1;
	    last INCL_PARENT;
	  }
	}
	
	unless ( $found_triangle ){
	  
	OTHER_EXT_PARENT2:
	  foreach my $other_parent ( @{$extension_parent->inclusion_parents} ){
	    my $other_ext_incl_id = $other_parent->transcript->dbID;
	    print STDERR "inclusion parent of ".$extension_parent->transcript->dbID." = ".$other_ext_incl_id."\n" if $verbose;
	    if ( $other_parent == $inclusion_parent ){
	      print STDERR "found triangle\n" if $verbose;
	      $found_triangle = 1;
	      last INCL_PARENT;
	    }
	  }
	}
      }
      if ( $found_triangle == 0 ){
	print STDERR "id candidate\n" if $verbose;
	$is_candidate = 1;
	$node->is_candidate(1);
	$node->add_candidate_extension_parent($extension_parent);
      }
    }
    if ( $is_candidate ){
      return 0;
    }
    else{
      return 1;
    }
}

############################################################

sub check_tree{
    my ($self, $newnode, $tree , $root) = @_;
    
    my $verbose = $self->verbose;
    
    # tree is a ref to a list of leaves
    my @ids = map { $_->transcript->dbID } @$tree;
    print STDERR "checking tree @ids:\n" if $verbose;
    my @nodes_to_check = @$tree;
    
    # as we check we regenerate the leaves as needed
    my $newleaves;
    
    # keep track of which ones are leaves 
    my %is_leaf;
    foreach my $node (@nodes_to_check){
	$is_leaf{$node} = 1;
    }
    
    # and which ones are added to the re-generated list of leaves
    my %added;
    
    # keep track of parents being added to the @nodes_to_check
    my %added_parent;
    
    # keep track whether the new element has been placed at all or not
    my $placed=0;
    
  NODE:
    while( @nodes_to_check ){
	my $node = shift @nodes_to_check;
	
	############################################################
	# skip if the node has been tagged as visited
	############################################################
	if ( $self->is_visited($node) ){
	    if ( $is_leaf{$node} && !$added{$node} ){
		push( @$newleaves, $node );
	    }
	    print STDERR "skipping visited node\n" if $verbose;
	    next NODE;
	}
	print STDERR "check_tree(): node = ".$node->transcript->dbID."\n" if $verbose;
	
	my $result = $self->check_node($newnode,$node);
	
	############################################################
	# if no-overlap, we jump to the next leaf
	############################################################
	if ( $result eq 'no-overlap' ){
	    print STDERR "no-overlap\n" if $verbose;
	    # recover the leaf if it hasn't been added before
	    if ( $is_leaf{$node} &&  !$added{$node} ){
		push(@$newleaves,$node);
		$added{$node} = 1;
	    }	  
	    next NODE;
	}
	############################################################
	# if 'continue' that means that we did not manage to place it
	# we might be in a worst case scenario - we jump to the next leaf
	############################################################
	elsif ( $result eq 'continue' ){
	    # recover the leaf if it hasn't been added before
	    print STDERR "check_tree(): CONTINUE\n" if $verbose;
	    if ( $is_leaf{$node} &&  !$added{$node} ){
		push(@$newleaves,$node);
		$added{$node} = 1;
	    }
	    next NODE;
	}
	############################################################
	# if placed, 
	############################################################
	elsif ( $result eq 'placed' ){
	    $placed++;
	    
	    print STDERR "placed\n" if $verbose;
	    # check whether newnode is a new leaf
	    if ( !$root && !$newnode->is_extended && !$newnode->is_included ){
	      if ( $added{$newnode} ){
		print STDERR "newnode is already a leaf\n" if $verbose;
	      }
	      unless ( $added{$newnode} ){
		print STDERR "adding newnode as leaf\n" if $verbose;
		push(@$newleaves,$newnode);
		$added{$newnode} = 1;
	      }
	    }
	    
	    # recover the leaf if it hasn't been added before, and if it is still a leaf
	    if ( $is_leaf{ $node } && !$added{$node} && !$node->is_extended ){
	      push(@$newleaves,$node);
	      $added{$node} = 1;
	    }
	    next NODE;
	    #return 'placed';
	  }
	############################################################
	# place here if it extends a node. 
	# 'Delete' leaf and add new leaf.
	# Tag all the extension parents
	############################################################
	elsif ($result eq 'extension'){
	    
	    ############################################################
	    # before placing here, check that it is a valid extension
	    if ( $self->_check_for_closed_loops( $newnode, $node ) ){
		print STDERR "extension: adding new extension\n" if $verbose;
		$newnode->add_extension_parent($node);
		$node->is_extended(1);
		
		############################################################
		# add to the root if it is included in it 
		if ( $root && $self->compare( $newnode, $root) eq 'inclusion' ){
		    print STDERR "adding inclusion of ".$newnode->transcript->dbID." into ".$root->transcript->dbID."\n" if $verbose;
		    $root->add_inclusion_child( $newnode );
		    $newnode->add_inclusion_parent( $root);
		    $newnode->is_included(1);
		    $self->tag_extension_ancestors($root);
		    
		}

		############################################################
		# this is a new leaf of this inclusion sub-tree 
		# ( and of the main tree if it has not been included or extended before
		unless( $added{$newnode} ){
		  my $label;
		  if ($root){
		    $label = "subtree";
		    print STDERR "adding newnode as leaf of this $label\n" if $verbose;
		    push(@$newleaves,$newnode);
		    $added{$newnode} = 1;
		  }
		  elsif( !$newnode->is_included && !$newnode->is_extended ){
		    $label = "tree";
		    print STDERR "adding newnode as leaf of this $label\n" if $verbose;
		    push(@$newleaves,$newnode);
		    $added{$newnode} = 1;
		  }
		  # tag this node and all the extension ancestors recursively
		  $self->tag_extension_ancestors($node);
		  $placed++;
		  next NODE;
		}
	      }
	    ############################################################
	    # else, recover the leaf if it has not been added before and continue
	    else{
		if ( $is_leaf{$node} && !$added{$node} && !$node->is_extended ){
		    push(@$newleaves,$node);
		    $added{$node} = 1;
		}
		next NODE;
	    }
	}
	############################################################
	# if there is a clash, 
	# check first the inclusion trees
	############################################################
	elsif ( $result eq 'clash' ){
	    
	    print STDERR "recursing inclusion tree\n" if $verbose;
	    my $result2 = $self->_recurse_inclusion_tree($newnode, $node);
	    
	    ############################################################
	    # if we manage to place it,
	    if ( $result2 eq 'placed' ){
		$placed++;
		
		############################################################
		# if we are in the tree without root ( the main list of leaves)
		# and our placed node is not extended nor included, it is a new leaf
		if ( !$root && !$newnode->is_extended && !$newnode->is_included ){
		    if ( $added{$newnode} ){
			print STDERR "newnode is already a leaf\n" if $verbose;
		    }
		    unless ( $added{$newnode} ){
			print STDERR "adding newnode as leaf\n" if $verbose;
			push(@$newleaves,$newnode);
			$added{$newnode} = 1;
		    }
		}
		
		# recover the leaf if it hasn't been added before, and if it is still a leaf
		if ( $is_leaf{ $node } && !$added{$node} && !$node->is_extended ){
		    push(@$newleaves,$node);
		    $added{$node} = 1;
		}
		next NODE;
	    }
	    ############################################################
	    # if clash, do breadth-first on the parents
	    elsif ($result2 eq 'clash' ){
		if ( @{$node->extension_parents} ){
		    my @parents = map { $_->transcript->dbID }  @{$node->extension_parents};
		    print STDERR "adding extension parents ( @parents) to the list to check\n" if $verbose;
		    foreach my $parent ( @{$node->extension_parents} ){
			unless ( $added_parent{ $parent } ){
			    push( @nodes_to_check, $parent );
			    $added_parent{$parent}++;
			}
		    }
		}
		# recover the leaf if it hasn't been added before
		if ( $is_leaf{$node} &&  !$added{$node} && !$node->is_extended){
		    push(@$newleaves,$node);
		    $added{$node} = 1;
		}
		next NODE;
	    }
	    ############################################################
	    # else simply skip to the next node to check
	    elsif( $result2 eq 'continue' ){
		# recover the leaf if it hasn't been added before
		if ( $is_leaf{$node} &&  !$added{$node} && !$node->is_extended){
		    push(@$newleaves,$node);
		    $added{$node} = 1;
		}
		next NODE;
	    }
	    
	    
	}#end of elsif(clash)
	    
	    
	} # end of NODE
    
    # this is our new tree
    @$tree = @$newleaves;
    $newleaves = [];
    
    if ($placed){
	return 'placed';
    }
    else{
	return 'continue';
    }
}  


############################################################
# this method checks for cases like this:
#   x ==> z ==> w
#   |          /
#   |        /
#   |      /
#   |    /
#   |  /
#    y
#
# which give rise to redundant solutions

sub _check_for_closed_loops{
    my ($self,$newnode,$node) = @_;
    
    my $verbose = $self->verbose;

    ############################################################
    # $newnode extends $node, but we need to check that $node
    # is not included in any of the extension parents of $newnode.

    ############################################################
    # the difficult question is how far we should go
    # In principle let's check the inmediate parents only
    if ( @{$newnode->extension_parents} ){
	foreach my $parent ( @{$newnode->extension_parents} ){
	    next if $parent == $node;
	    if ( $self->compare( $node, $parent ) eq 'inclusion' ){
		print STDERR "potential extension parent ".$node->transcript->dbID.
		    " is included in an extension parent ".$parent->transcript->dbID."\n" if $verbose;
		return 0;
	    }
	}
	return 1;
    }
    else{
	return 1;
    }
}

############################################################

sub _recurse_inclusion_tree{
  my ($self, $newnode, $node ) = @_;
  
  my $verbose = $self->verbose;

  ############################################################
  # start checking the inclusion tree, if there is any
  if ( defined $node->inclusion_tree ){
      
      print STDERR "clash -> checking inclusion-tree\n" if $verbose;
      return $self->check_tree( $newnode, $node->inclusion_tree, $node);
  }
  ############################################################
  # if the inclusion tree is not valid, skip to the next node to check
  elsif( !defined $node->inclusion_tree && @{$node->inclusion_children} ){
      print STDERR "not-valid inclusion tree -> continue\n" if $verbose;
      return 'continue';
  }
  ############################################################
  # else do breadth-first on the parents
  else{ 
      print STDERR "no inclusion tree - returning clash\n" if $verbose;
      return 'clash';
  }
}


############################################################


sub check_node{
  my ($self,$newnode,$node) = @_;

  #test#
  if ( $newnode == $node || $newnode->transcript == $node->transcript ){
    $self->throw( "<<<<<< Comparing node with its self - something has gone wrong - exiting >>>>>");
  }
  
  my $verbose = $self->verbose;
  
  #print STDERR "check_node(): newnode: $newnode, node: $node\n" if $verbose;
  my $result = $self->compare( $newnode, $node );
  print STDERR "check_node(); result = $result\n" if $verbose;

  ############################################################
  # inclusion is recursive in the extension tree
  ############################################################
  if( $result eq 'inclusion' ){
      print STDERR "inclusion:\n" if $verbose;
      my $result2 = $self->_recurse_extension_branch( $newnode, $node, 'inclusion');
      if ( $result eq 'placed' ){
	  return 'placed';
      }
      else{
	  print STDERR "returning ".$result2."\n" if $verbose;
	  return $result2;
      }
  }
  else{
      return $result;
  }
}

############################################################


sub _recurse_extension_branch{
  my ($self, $newnode,$node) = @_;
  
  my $verbose = $self->verbose;
  
  ############################################################
  # recurse along the extension branch if there is any
  if ( @{$node->extension_parents} ){
    print STDERR "inclusion -> go to extension parent\n" if $verbose;
      my $placed = 0;
    
  PARENT:
    foreach my $parent_node ( @{$node->extension_parents} ){
	  
      # skip if already visited (extended)
      if ( $self->is_visited($parent_node) ){
	print STDERR "skipping this node - already visited\n" if $verbose;
	next PARENT;
      }
      
      my $result2 = $self->check_node( $newnode, $parent_node );
      if ( $result2 eq 'extension' || $result2 eq 'no-overlap' ){
	
	#print STDERR "is there an inclusion tree on this node: ".$node->transcript->dbID."\n" if $verbose;
	
	############################################################
	# if there is an inclusion tree and is valid, then check
	if ( defined $node->inclusion_tree ){
	  my @tree = map { $_->transcript->dbID } @{$node->inclusion_tree};
	  print STDERR "checking the inclusion tree : @tree\n" if $verbose;
	  
	  my $result3 = $self->check_tree( $newnode, $node->inclusion_tree , $node);
	  
	  print STDERR "result = $result3\n" if $verbose;
	  if ( $result3 eq 'placed' ){
	    $placed++;
	  }
	  ############################################################
	  # place here if no overlap, we know it is included
	  elsif( $result3 eq 'continue' || $result3 eq 'no-overlap' ){
	    print STDERR "adding inclusion of ".$newnode->transcript->dbID." into ".$node->transcript->dbID."\n" if $verbose;
	    $node->add_inclusion_child($newnode);
	    $newnode->add_inclusion_parent($node);
	    $newnode->is_included(1);
	    $self->tag_extension_ancestors($node);
	    if ( $self->compare( $newnode, $parent_node ) eq 'extension' ){
	      print STDERR "adding extension of ".
		$newnode->transcript->dbID." to ".$parent_node->transcript->dbID."\n" if $verbose;
	      $newnode->add_extension_parent($parent_node);
	      $parent_node->is_extended(1);
	      $self->tag_extension_ancestors($parent_node);
	    }
	    $placed++;}
	  elsif( $result3 eq 'placed' ){
	    $placed++;
	    next PARENT;
	  }
	}
	############################################################
	# else place here
	else{
	  print STDERR "adding inclusion of ".$newnode->transcript->dbID." into ".$node->transcript->dbID."\n" if $verbose;
	  $node->add_inclusion_child( $newnode );
	  $newnode->add_inclusion_parent( $node );
	  $newnode->is_included(1);
	  $self->tag_extension_ancestors($node);
	  if ( $self->compare( $newnode, $parent_node ) eq 'extension' ){
	    print STDERR "adding extension of ".
	      $newnode->transcript->dbID." to ".$parent_node->transcript->dbID."\n" if $verbose;
	    $newnode->add_extension_parent($parent_node);
	    $parent_node->is_extended(1);
	    $self->tag_extension_ancestors($parent_node);
	  }
	  $placed++;
	}
      }
      elsif( $result2 eq 'placed'){
	$placed++;
      }
    } # end of PARENT
    
    ############################################################
    # continue if we didn't place it
    if ($placed){
      return 'placed';
    }
    else{
      print STDERR "could not place it - continue\n" if $verbose;
      return 'continue';
    }
  }
  ############################################################
  # else, check the inclusion branch if any
  elsif ( defined $node->inclusion_tree ){
      print STDERR "inclusion -> go to inclusion tree\n" if $verbose;
      return $self->check_tree( $newnode, $node->inclusion_tree , $node);
  }
  #elsif( @{$node->inclusion_children} && !defined $node->inclusion_tree ){
  #    #my @list = @{$node->inclusion_children};
  #    #print STDERR "inclusion children = @list\n";
  #    print STDERR "non-valid inclusion tree -> continue\n";
  #    return 'continue';
  #}
  ############################################################
  # else, we know we have inclusion, we place it here
  # and continue with the next global leaf
  else{
      
      # place here
      print STDERR "adding inclusion of ".$newnode->transcript->dbID." into ".$node->transcript->dbID."\n" if $verbose;
      $node->add_inclusion_child( $newnode );
      $newnode->add_inclusion_parent( $node );
      $newnode->is_included(1);
      $self->tag_extension_ancestors($node);
      return 'placed';
  }
}


############################################################


=head2 compare

  Function: encapsulates the comparison method between two transcripts.

=cut

sub compare{
  my ($self,$newnode,$node) = @_;
  my $newtrans = $newnode->transcript;
  my $trans    = $node->transcript;


  #print STDERR "comparirson_level: ".$self->_comparison_level."\n";
  my $comparator = 
    Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptComparator
      ->new(
	    -comparison_level                 => $self->_comparison_level, 
	    -exon_match                       => $self->_exon_match,
	    -splice_mismatch                  => $self->_splice_mismatch,
	    -intron_mismatch                  => $self->_intron_mismatch,
	    -internal_splice_overlap          => $self->_internal_splice_overlap,
	   );
  
  my %overlap_matrix = %{$self->matrix};
  
  if ( defined( $overlap_matrix{$newtrans}{$trans} ) ){
    #print STDERR "Using cached matrix\n";
    return $overlap_matrix{$newtrans}{$trans};
  }
  else{
    #my $t1 = time;
    $overlap_matrix{$newtrans}{$trans} = $comparator->discrete_compare($newtrans, $trans);
    #my $t2 = time;
    #print STDERR "COMPARATOR TIME (".$newtrans->dbID.",".$trans->dbID.") = ".($t2-$t1)."\n";
    $self->matrix( \%overlap_matrix );
    return $overlap_matrix{$newtrans}{$trans};
  }
}


############################################################

sub tag_extension_ancestors{
    my ($self,$node) = @_;
    unless ( $self->is_visited($node ) ){
	$self->is_visited( $node, 1 );
    }
    foreach my $parent (@{$node->extension_parents} ){
	unless ( $self->is_visited($parent) ){
	    $self->tag_extension_ancestors($parent);
	} 
    }
}

############################################################

sub is_visited{
  my ($self,$node,$boolean) = @_;
  #unless( $self->{_visited} ){
  #  $self->{_visited} = {};
  #}
  my $verbose = $self->verbose;

  if ( $node && defined $boolean){
    $self->{_visited}{$node} = $boolean;
    print STDERR "tagging node ".$node->transcript->dbID." as visited (".$self->{_visited}{$node}.")\n" if $verbose;
    return $self->{_visited}{$node};
  }
  return $self->{_visited}{$node};
}

############################################################

sub flush_visited_tags{
  my ($self) = @_;
  delete $self->{_visited};
}

############################################################
#
# RUN METHOD
#
############################################################

sub run {
  my $self = shift;
  my @transcripts = $self->input_transcripts;
	
  ############################################################
  # cluster the transcripts
  my @transcript_clusters = $self->_cluster_Transcripts(\@transcripts);
  print STDERR scalar(@transcript_clusters)." clusters returned from _cluster_Transcripts\n";
	 
  ############################################################
  # restrict the number of ESTs per cluster, for the time being, to avoid deep recursion 
  #my @transcript_clusters = $self->_filter_clusters(\@transcript_clusters1,100);

  ############################################################
  # find the lists of clusters that can be merged together according to consecutive exon overlap
  my @lists = $self->link_Transcripts( \@transcript_clusters );
  print STDERR scalar(@lists)." lists returned from link_Transcripts\n";

  ############################################################
  # merge each list into a putative transcript
  my @merged_transcripts  = $self->_merge_Transcripts(\@lists);
  print STDERR scalar(@merged_transcripts)." transcripts returned from _merge_Transcripts\n";
  
  ############################################################
  # we can score the predictions:
  my @final_transcripts;
  if ( $self->use_score ){
    print STDERR "calculating scores\n";
      my $score_model = 
	Bio::EnsEMBL::Pipeline::GeneComparison::ScoreModel
	    ->new(
		  -hold_list   => $self->hold_list,
		  -transcripts => \@merged_transcripts,
		  );
    @final_transcripts = $score_model->score_Transcripts;
  }
  else{
    @final_transcripts = @merged_transcripts;
  }
  
  @merged_transcripts = ();
  #foreach my $tran ( @more_merged_transcripts ){
  #  Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript($tran);
  #}
  
  $self->output(@final_transcripts);
}

############################################################

sub _filter_clusters{
  my ( $self, $clusters, $threshold ) = @_;
  my @new_clusters;
 
 CLUSTER:
  foreach my $cluster ( @$clusters ){
    my $new_cluster = Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster->new();
    my @transcripts = 
      sort { scalar( @{$b->get_all_Exons}) <=> scalar(@{$a->get_all_Exons}) } @{$cluster->get_Transcripts};
    my $count = 0;
  
    my @to_be_added =();
  TRAN:
    foreach my $t ( @transcripts ){
      $count++;
      if ( $count > $threshold){
	last TRAN;
      }
      push( @to_be_added, $t );
    }
    $new_cluster->put_Transcripts( @to_be_added );
    push( @new_clusters, $new_cluster );
  }
  return @new_clusters;
}
      
############################################################
#
# CLUSTERING OF TRANSCRIPTS BY GENOMIC EXTENT
#
############################################################

=head2 _cluster_Transcripts

 Description : It separates transcripts according to strand and then clusters 
               each set of transcripts by calling _cluster_Transcripts_by_genomic_range()
  Args       : Array of Bio::EnsEMBL::Transcript
  Return     : Array of Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster

=cut

sub _cluster_Transcripts {
  my ($self,$transcripts) = @_;
 
  my $forward_transcripts;
  my $reverse_transcripts;
 
  while ( @$transcripts){
    my $transcript = shift @$transcripts;
    my @exons = @{ $transcript->get_all_Exons };
    if ( $exons[0]->strand == 1 ){
      push( @$forward_transcripts, $transcript );
    }
    else{
      push( @$reverse_transcripts, $transcript );
    }
  }
  my @clusters;
  if ( $forward_transcripts && @$forward_transcripts ){
    my @forward_clusters = $self->_cluster_Transcripts_by_genomic_range( $forward_transcripts );
    if ( @forward_clusters){
      print STDERR scalar( @forward_clusters )." clusters in forward strand\n";
      push( @clusters, @forward_clusters);
    }
  }
  if ( $reverse_transcripts && @$reverse_transcripts ){
    my @reverse_clusters = $self->_cluster_Transcripts_by_genomic_range( $reverse_transcripts );
    if ( @reverse_clusters ){
      print STDERR scalar( @reverse_clusters )." clusters in reverse strand\n";
      push( @clusters, @reverse_clusters);
    }
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
  my ($self,$mytranscripts) = @_;

  unless (@$mytranscripts){
    return;
  }

  # filter some small transcripts:
  #my @transcripts = sort { scalar( @{$b->get_all_Exons} ) <=> scalar( @{$a->get_all_Exons} ) } @mytranscripts;
  
  # take only the longest 50;
  #@mytranscripts = ();
  #my $c = 0;
  #foreach my $t ( @transcripts ){
  #  $c++;
  #  push( @mytranscripts, $t);
  #  last if $c == 50;
  #}
  
  # first sort the transcripts
  my @transcripts = sort { my $result = ( $self->transcript_low($a) <=> $self->transcript_low($b) );
			if ($result){
			  return $result;
			}
			else{
			  return ( $self->transcript_high($b) <=> $self->transcript_high($a) );
			}
		      } @$mytranscripts;
  
  @$mytranscripts = ();
  print STDERR "clustering ".scalar(@transcripts)." transcripts\n";
  
  # create a new cluster 
  my $cluster=Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster->new();
  my $count = 0;
  my @cluster_starts;
  my @cluster_ends;
  my @clusters;
  
  # put the first transcript into these cluster
  $cluster->put_Transcripts( $transcripts[0] );
  #print STDERR "first cluster:\n";
  #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript( $transcripts[0] );
  
  $cluster_starts[$count] = $self->transcript_low( $transcripts[0]);
  $cluster_ends[$count]   = $self->transcript_high($transcripts[0]);
  
  # store the list of clusters
  push( @clusters, $cluster );
  
  # loop over the rest of the transcripts
 LOOP1:
  for (my $c=1; $c<=$#transcripts; $c++){
    
    #print STDERR "\nIn cluster ".($count+1)."\n";
    #print STDERR "start: $cluster_starts[$count] end: $cluster_ends[$count]\n";
    #print STDERR "comparing:\n";
    #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript( $transcripts[$c] );
    
    if ( !( $self->transcript_high($transcripts[$c]) < $cluster_starts[$count] ||
	    $self->transcript_low($transcripts[$c])  > $cluster_ends[$count] ) ){
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
  
  my @new_clusters;
  my $cutoff = 0;
  if ( $cutoff ){
    my $ccount = 0;
    foreach my $cluster ( @clusters ){
      $ccount++;
      my $new_cluster = Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster->new();
      push( @new_clusters, $new_cluster );
      my @transcripts = @{$cluster->get_Transcripts};
      print STDERR "cluster $ccount: ".scalar(@transcripts)." transcripts\n";
      @transcripts = sort { scalar( @{$b->get_all_Exons} ) <=> scalar( @{$a->get_all_Exons} ) } @transcripts;
      
      #take only the longest 50;
      my $c = 0;
      foreach my $trans ( @transcripts ){
	$c++;
        $new_cluster->put_Transcripts( $trans );
	last if $c == 50;
      }
      print STDERR "new_cluster $ccount: ".scalar( @{$new_cluster->get_Transcripts} )." transcripts\n";
    }
    return @new_clusters;
  }
  else{
    return @clusters;
  }
}

############################################################



=head2 by_transcript_high

Description:  it returns the comparison ( $a <=> $b ) of the right-most coordinates if they
              are different, else it returns the comparison ( $a <=> $b ) 
              of the left most coordinate.
=cut

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

=head2 transcript_high

Description: it returns the highest coordinate of a transcript

=cut

sub transcript_high{
  my ($self,$tran) = @_;
  my @exons = sort { $a->start <=> $b->start } @{$tran->get_all_Exons};
  return $exons[-1]->end;
}

############################################################

=head2 transcript_low

Description: it returns the lowest coordinate of a transcript

=cut

sub transcript_low{
  my ($self,$tran) = @_;
  my @exons = sort { $a->start <=> $b->start } @{$tran->get_all_Exons};
  return $exons[0]->start;
}

############################################################

sub matrix{
  my($self,$matrix) = @_;
  if (defined $matrix){
    $self->{_matrix} = $matrix;
  }
  return $self->{_matrix};
}

############################################################


############################################################
#
# METHODS TO MERGE TRANSCRIPTS
#
############################################################

=head2 _merge_Transcripts

description: make a transcript for every list built above in link_Transcripts().
             the procedure goes by clustering the exons, finding the splice ends
             and then producing the transcripts.
        Arg: a listref with al the listrefs produced in link_Transcripts()

=cut

sub _merge_Transcripts{
    my ($self,$lists) = @_;
    print STDERR "<<<<<<<<<< merging transcripts >>>>>>>>>>\n";
	
    my $verbose = 1;#$self->verbose;
    
    
    # $list is an arrayref of the ests/cdnas that we can merge
    my @merged_transcripts;
    
    my $count = 0;
  
  LIST:
    foreach my $list ( @$lists ){
      $count++;
      
      if ($verbose){
	print STDERR "list $count:\n";
	foreach my $t ( @$list ){
	  Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript( $t );
	}
      }
      my $order = scalar( @$list );
      if ( $order < $self->_minimum_order ){
	print STDERR "list has only $order elements and the minimum allowed is ".$self->_minimum_order." - rejecting\n";
	next LIST;
      }
      
      my @allexons;
      my %exon2transcript;			
      my %is_first;
      my %is_last;			
      
      # collect all exons
      my %seen_already;
      my $newlist;
    EST:
      foreach my $tran (@{ $list }){
	next EST if $seen_already{$tran};
	$seen_already{$tran} = 1;
	push ( @$newlist, $tran );
	my @exons = @{$tran->get_all_Exons};
	
	# sort them in genomic order regardless of the strand
	@exons    = sort { $a->start <=> $b->start } @exons;
	
	push ( @allexons, @exons );
	
	# keep track of whether the exons is first or last and the transcript it belongs to
	for (my $i = 0; $i< scalar( @exons ); $i++){
	  if ( 0 == $i ){
	    $is_first{$exons[$i]} = 1;
	  }
	  else{
	    $is_first{$exons[$i]} = 0;
	  }
	  if ( $#exons == $i ){
	    $is_last{$exons[$i]} = 1;
	  }
	  else{
	    $is_last{$exons[$i]} = 0;
	  }
	  $exon2transcript{$exons[$i]} = $tran;
	}
      }
      $list = $newlist;
      
      # cluster the exons
      my $first_cluster_list = $self->_cluster_Exons( @allexons );
      
      # set start and end of the clusters (using the info collected above)
      my $cluster_list = $self->_set_splice_Ends($first_cluster_list,\%exon2transcript,\%is_first,\%is_last);
      
      # we turn each cluster into an exon and create a new transcript with these exons
      my $transcript    = Bio::EnsEMBL::Transcript->new();
      $transcript->type( "merged_".$count );
      
      my @exon_clusters = $cluster_list->sub_SeqFeature;
      
      my $exon_count = 0;
      foreach my $exon_cluster (@exon_clusters){
	$exon_count++;

	my $new_exon = Bio::EnsEMBL::Exon->new();
	
	$new_exon->start($exon_cluster->start );
	$new_exon->end  ($exon_cluster->end   );
	$new_exon->seqname( $transcript->type."_".$exon_count);
	$new_exon->score(100);

	############################################################
	# they should not have a translation
	############################################################
	$new_exon->phase(0);
	$new_exon->end_phase(0);
	
	# set the strand to be the same as the exon_cluster
	$new_exon->strand($exon_cluster->strand);
	
	my %evidence_hash;
	my %evidence_obj;
	foreach my $exon ( $exon_cluster->sub_SeqFeature ){
	  $self->_transfer_Supporting_Evidence($exon,$new_exon);
	}
	
	$transcript->add_Exon($new_exon);
      }
      
      if ($verbose){
	  print STDERR "Produced transcript:\n";
	Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript( $transcript );
      }
      
      if ( $self->use_score ){
	my %duplicated;
	foreach my $est ( @$list ){
	  $duplicated{$est}++;
	  if ( $duplicated{$est}==2 ){
	    print STDERR "CHECK1 est $est duplicated\n";
	  }
	}


	$self->hold_list( $transcript, $list );
      }
      push ( @merged_transcripts, $transcript );
      
  } # end of LIST
    return @merged_transcripts;
    
}

############################################################

=head2 _cluster_Exons
 
 Function: it cluster exons according to overlap,
           it returns a Bio::EnsEMBL::SeqFeature, where the sub_SeqFeatures
           are exon_clusters, which are at the same time Bio::EnsEMBL::SeqFeatures,
           whose sub_SeqFeatures are exons
=cut

sub _cluster_Exons{
  my ($self, @exons) = @_;

  #print STDERR "EST_GeneBuilder: clustering exons...\n";
  
  # no point if there are no exons!
  return unless ( scalar( @exons) > 0 );   

  # keep track about in which cluster is each exon
  my %exon2cluster;
  
  # main cluster feature - holds all clusters
  my $cluster_list = new Bio::EnsEMBL::SeqFeature; 
  
  # sort exons by start coordinate
  @exons = sort { $a->start <=> $b->start } @exons;

  # Create the first exon_cluster
  my $exon_cluster = new Bio::EnsEMBL::SeqFeature;
  
  # Start off the cluster with the first exon
  $exon_cluster->add_sub_SeqFeature($exons[0],'EXPAND');

  $exon_cluster->strand($exons[0]->strand);    
  $cluster_list->add_sub_SeqFeature($exon_cluster,'EXPAND');
  
  # Loop over the rest of the exons
  my $count = 0;
  
 EXON:
  foreach my $exon (@exons) {
    if ($count > 0) {
      
      # Add to cluster if overlap AND if strand matches
      if ( $exon_cluster->overlaps($exon) && ( $exon->strand == $exon_cluster->strand) ) { 
	$exon_cluster->add_sub_SeqFeature($exon,'EXPAND');
      }  
      else {
	# Start a new cluster
	$exon_cluster = new Bio::EnsEMBL::SeqFeature;
	$exon_cluster->add_sub_SeqFeature($exon,'EXPAND');
	$exon_cluster->strand($exon->strand);
		
	# and add it to the main_cluster feature
	$cluster_list->add_sub_SeqFeature($exon_cluster,'EXPAND');	
      }
    }
    $count++;
  }
  return $cluster_list;
}

############################################################

=head2 _set_splice_Ends

  Function: Resets the ends of the clusters to the most frequent coordinate
            not allowing to contribute exons which are at one end of its
            original transcript. When the cluster represents the first or the lasst exon
            of the merged transcript, we always take the external exon coordinate that
            extends the transcript the most, to allow for UTRs to build up.
    Arg[1]: listref with the exon_clusters
    Arg[2]: hashref assigning to each exon, the transcripts it comes from
    Arg[3]: hashref holding a BOOLEAN for each exon to test whether it is the first in the transcript
    Arg[4]: hashref holding a BOOLEAN for each exon to test whether it is the last in the transcript
            Note - order in transcripts is here considered low-coord to high-coord, regardles of
            the strand
=cut

sub _set_splice_Ends {
  my ($self, $cluster_list, $ref_exon2transcript_hash, $ref_is_first, $ref_is_last) = @_;
  
  # $cluster_list is a arrayref of exon-clusters, each cluster 
  # representing all the exons that can be merged

  # hash having exons as keys and mother-transcript as value
  my %exon2transcript = %$ref_exon2transcript_hash;

  # keep track of whether the exon is first or last in the transcript
  my %is_first = %$ref_is_first;
  my %is_last  = %$ref_is_last;

  #print STDERR "EST_GeneBuilder: setting common ends...\n";

  # get the exon clusters
  my @exon_clusters = $cluster_list->sub_SeqFeature;

  # sort clusters according to their start coord.
  @exon_clusters = sort { $a->start <=> $b->start  } @exon_clusters;
  my $count =  0;

  # check whether a cluster has fused two separate exons from the same transcript
  my $position_is = 0;
  my @exon_list;
  my $need_to_recluster = 0;		      

 CLUSTERS:		      
 foreach my $cluster ( @exon_clusters ){
   $position_is++;
   
   # keep track of the transcripts used in this cluster
   my %seen_transcript;
   my %exons_from_transcript;
   
   # each cluster is a SeqFeature with sub_SeqFeatures
   foreach my $exon ( $cluster->sub_SeqFeature ){
     
     push ( @{ $exons_from_transcript{ $exon2transcript{ $exon } } }, $exon );  
     
     #if ( $seen_transcript{ $exon2transcript{$exon} } && 1 == $seen_transcript{ $exon2transcript{$exon} } ){
     #  #print STDERR "There is more than one exon from transcript $exon2transcript{ $exon }\n"; 
     #  $fused = 1;
     #  $need_to_recluster = 1;
     #}
     $seen_transcript{ $exon2transcript{ $exon } } = 1;
   } # end of exon loop
   
   # isimply collect the exons
   push( @exon_list, $cluster->sub_SeqFeature);
      
   ## we don't check any more for fusion. Whether one wants to allow fusion or not should happen during
   ## the comparison between transcripts.
   
   #   elsif ( $fused == 1 ){
      
   #     my @exons = sort{ $b->length <=> $a->length } $cluster->sub_SeqFeature;
   
   #   EXONS:
   #     foreach my $exon ( @exons ){
   
   #     TRANS_IN_CLUSTER:
   #       foreach my $tran ( keys( %exons_from_transcript ) ){
   #	 if ( $exon2transcript{$exon} eq $tran ){
   #	   next;
   #	 }	  
   #	 my $overlap_count = 0;
   #	  foreach my $exon2 ( @{ $exons_from_transcript{ $tran } } ){
   #	    if ( $exon->overlaps($exon2) ){
   #	      $overlap_count++;
   #	    }
   #	  }
   #	  # if  $exon overlaps 2 or more exons from the same transcript, in this cluster ...
   #	  if ( $overlap_count >= 2 ){
   
   #	    # ... and $exon is the only one from its transcript in this cluster
   #	    if ( scalar( @{ $exons_from_transcript{ $exon2transcript{$exon} } } ) == 1 ){
   
   #	      # ... and $exon is at the edge of its transcript
   #	      if ( $is_first{ $exon } == 1 || $is_last{ $exon } == 1 ){
   
   #		#print STDERR "ditching one exon\n";
   #		# then we ditch it and continue ...
   #		next EXONS;
   #	      }
   #	    }
   #	  }
   #	} # end of TRANS_IN_CLUSTER
   
   #	# if it is good, it gets here
   #        push ( @exon_list, $exon );
   
   #      }   # end of EXONS
   
   #    }     # end of 'if fused == 1'
   
   
    
  }       # end of CLUSTERS
	
  ## we don't recluster, let's keep things simple
  ## if needed, re-cluster
  #if ( $need_to_recluster == 1){
  #  @exon_clusters   = ();
  #  #print STDERR " *** Reclustering ***\n";
  #  $cluster_list  = $self->_cluster_Exons( @exon_list );
  #  @exon_clusters = $cluster_list->sub_SeqFeature;
  #}

  # keep track of the cluster position
  my $position = 0;

CLUSTER:		     
  foreach my $cluster (@exon_clusters) {
    $position++;
    my %start;
    my $new_start;
    my $max_start;
    
  #### set first the start-coordinates ####
  
  EXON:
    foreach my $exon ($cluster->sub_SeqFeature){

	# for a start-coord in the middle 
	if ( $position > 1 && $position <= scalar( @exon_clusters ) ){
	
	# don't use the exon if it is the first of a transcript
	if ( $is_first{ $exon } == 1 ){
	  next EXON;
	}
      }
      $start{$exon->start}++;
    }
    # we can as well sort them:
    my @starts = sort{ $start{ $b } <=> $start{ $a } } keys( %start );

    # take the most common start (note that we do not resolve ambiguities here)
    # at some point we could resolve them by checking for the splice site consensus sequence

    ## test:
    #print STDERR "starts: ";
    #foreach my $s (@starts){
    #  print STDERR $s."\t";
    #}
    #print STDERR "\n";

    $new_start = shift( @starts );
    if ( defined $new_start ){
	$max_start = $start{ $new_start };
    }
    # if we have too little exons to obtain the start, take the original value
    if ( !defined $max_start ){
	$new_start = $cluster->start;
    }
    
    # the first cluster is a special case - potential UTRs, take the longest one.
    if( $position == 1) {
      $new_start = $cluster->start;  
    }

    #### now set the end coordinate ####

    my %end;
    my $new_end;
    my $max_end;
    
  EXON:
    foreach my $exon ($cluster->sub_SeqFeature){
      
      # for an end-coord in the middle 
      if ( $position >= 1 && $position < scalar( @exon_clusters ) ){
	
	# don't use the exon if it is the last of a transcript
	if ( $is_last{ $exon } == 1 ){
	  next EXON;
	}
      }
      $end{$exon->end}++;
    }
    my @ends = sort{ $end{ $b } <=> $end{ $a } } keys( %end );
    
    # take the most common end (note that we do not resolve ambiguities here)
    
    $new_end = shift( @ends );
    if ( defined $new_end ){
      $max_end = $end{ $new_end };
    }
    # if we have too little exons to obtain the end, take the original value
    if ( !defined $max_end ){
      #print STDERR "In last position, cluster end wins!\n";
      $new_end = $cluster->end;
    }
    
    # the last cluster is a special case - potential UTRs, take the longest one.
    if( $position == scalar( @exon_clusters ) ) {
      $new_end = $cluster->end; # recall that the cluster is expanded  
    }
    
    ######  if new_start> new_end change them in turns until we get start < end ######
    # this can happen when there are many exons in a cluster but they vary greatly in size
    my $other_start = $new_start;
    my $other_end   = $new_end;
    my $trouble     = 0;
    my $stop_start  = 0;
    my $stop_end    = 0;

    while ( $other_start >= $other_end ){
      #print STDERR "*** Trouble: $new_start >= $new_end ***\n";
      $trouble = 1;
      
      # take the next one
      if ( $stop_start == 0 ){
	my $re_start = shift( @starts );
	if ( $re_start ){
	  $other_start = $re_start;
	}
	else{
	  $stop_start = 1;
	}
      }
      if ( $other_start >= $other_end ){
	if ( $stop_end == 0 ){
	  my $re_end = shift( @ends );
	  if ( $re_end ){
	    $other_end = $re_end;
	  }
	  else{
	    $stop_end = 1;
	  }
	}
      }
      if ( $stop_start == 1 && $stop_end ==1 ){
	last;
      }
    }

    if ($trouble == 1 ){
      if ($other_end && $other_start && $other_start < $other_end ){
	$new_start = $other_start;
	$new_end   = $other_end;
	#print STDERR "Reset: new_start: $new_start\t new_end: $new_end\n";
      }
      else{
	## ok, you tried to change, but you got nothing, what can we do about it?
	#print STDERR "Could not reset the start and end coordinates\n";
	if ( $other_end <= $other_start ){
	  #print STDERR "Sorry will have to put the end = end of cluster\n";
	  $new_end  = $cluster->end;
	  if ( $new_start >= $new_end ){
	    #print STDERR "Last resort: I'm afraid we'll also have to put start = start of cluster, good luck\n";
	  }
	}
      }
    }
    
    # reset the cluster start (that's the coordinate used in _produce_Transcript() )
    $cluster->start($new_start);

    # reset the cluster end (that's the coordinate used in _produce_Transcript() )
    $cluster->end($new_end);
    
  } # end of CLUSTER

  return $cluster_list;
}

############################################################

=head2 _transfer_Supporting_Evidence

 Usage   : $self->transfer_supporting_evidence($source_exon, $target_exon)
 Function: Transfers supporting evidence from source_exon to target_exon, 
           after checking the coordinates are sane and that the evidence is not already in place.
 Returns : nothing, but $target_exon has additional supporting evidence

=cut

sub _transfer_Supporting_Evidence{
  my ($self, $source_exon, $target_exon) = @_;
  
  # this method fails when first called in a new exon without any evidence
  my $target_sf;
  eval{
    $target_sf = $target_exon->get_all_supporting_features;
  };   

  # keep track of features already transferred, so that we do not duplicate
  my %unique_evidence;
  my %hold_evidence;

 SOURCE_FEAT:
  foreach my $feat (@{$source_exon->get_all_supporting_features}){
    next SOURCE_FEAT unless $feat->isa("Bio::EnsEMBL::FeaturePair");
    
    # skip duplicated evidence objects
    next SOURCE_FEAT if ( $unique_evidence{ $feat } );
    
    # skip duplicated evidence 
    if ( $hold_evidence{ $feat->hseqname }{ $feat->start }{ $feat->end }{ $feat->hstart }{ $feat->hend } ){
      #print STDERR "Skipping duplicated evidence\n";
      #$self->print_FeaturePair($feat);
      next SOURCE_FEAT;
    }

    #$self->print_FeaturePair($feat);
    
    if ( (scalar @$target_sf) > 0){
   TARGET_FEAT:
     foreach my $tsf (@$target_sf){
       next TARGET_FEAT unless $tsf->isa("Bio::EnsEMBL::FeaturePair");
      
       if($feat->start    == $tsf->start &&
	  $feat->end      == $tsf->end &&
	  $feat->strand   == $tsf->strand &&
	  $feat->hseqname eq $tsf->hseqname &&
	  $feat->hstart   == $tsf->hstart &&
	  $feat->hend     == $tsf->hend){
	
 	#print STDERR "feature already in target exon\n";
        #$self->print_FeaturePair($feat);
	next SOURCE_FEAT;
      }
     }
    }
    #print STDERR "transfering evidence\n";
    #$self->print_FeaturePair($feat);
    #$feat->analysis($self->analysis);
    $target_exon->add_supporting_features($feat);
    $unique_evidence{ $feat } = 1;
    $hold_evidence{ $feat->hseqname }{ $feat->start }{ $feat->end }{ $feat->hstart }{ $feat->hend } = 1;
  }
}

############################################################

sub _print_List{
  my ($self,$list) = @_;
  foreach my $t ( @$list ){
    Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript($t);
  }
}


############################################################    
#
# GET/SET METHODS
#
############################################################

############################################################
# $list is an arrayref with the transcripts that make up the transcript $tran
sub hold_list{
    my ($self,$tran,$list) = @_;

    if ($tran){
	unless ( $self->{_hold_list}{$tran} ){
	    $self->{_hold_list}{$tran} = [];
	}
	if ( $list ){
	    $self->{_hold_list}{$tran} = $list;
	}
	return $self->{_hold_list}{$tran};
    }
    else{
	return $self->{_hold_list};
    }
}

############################################################

sub _minimum_order{
  my ($self,$minimum_order) = @_;
  if ($minimum_order) {
    $self->{_minimum_order} = $minimum_order;
  }
  return $self->{_minimum_order};
}

############################################################

sub _speed_up{
  my ($self,$boolean) = @_;
  if ( defined $boolean ){
    $self->{_speed_up} = $boolean;
  }
  return $self->{_speed_up};
}

############################################################

sub _comparison_level{
  my ($self, $level ) = @_;
  if ($level){
    $self->{_comparison_level} = $level;
  }
  return $self->{_comparison_level};
}
############################################################

sub _splice_mismatch{
  my ($self, $mismatch ) = @_;
  if (defined $mismatch){
    $self->{_splice_mismatch} = $mismatch;
  }
  return $self->{_splice_mismatch};
}

############################################################

sub _exon_match{
  my ($self,$exon_match) =@_;
       if ( defined $exon_match){
       $self->{_exon_match} = $exon_match;
       }
     return $self->{_exon_match};
     }

############################################################

sub _intron_mismatch{
  my ($self, $mismatch ) = @_;
  if (defined $mismatch){
    $self->{_intron_mismatch} = $mismatch;
  }
  return $self->{_intron_mismatch};
}

############################################################

sub _internal_splice_overlap{
  my ($self, $value) = @_;
  if (defined $value){
    $self->{_internal_splice_overlap} = $value;
  }
  return $self->{_internal_splice_overlap};
}

############################################################

sub _mismatch_allowed{
  my ($self, $mismatch) = @_;
  if (defined $mismatch){
    $self->{_mismatch_allowed} = $mismatch;
  }
  return $self->{_mismatch_allowed};
}

############################################################

sub input_transcripts{
    my ($self, @transcripts) = @_;
    if (@transcripts){
	push( @{$self->{_input_transcripts} }, @transcripts );
    }
    return @{$self->{_input_transcripts} };
}

############################################################

sub output {
    my ($self, @transcripts ) = @_;
    if (@transcripts){
	push( @{$self->{_output_transcripts} }, @transcripts );
    }
    return @{$self->{_output_transcripts} };
}

############################################################

# this holds a lists of lists, where each list
# contains the transcripts that can be merged with each other

sub sub_clusters {
  my ($self,@lists) = @_;
  if ( @lists ){
    push ( @{ $self->{_sub_clusters} }, @lists );
  }
  return @{$self->{_sub_clusters}};
}

############################################################

sub verbose{
  my ($self,$verbose) = @_;
  if ( $verbose ){
    $self->{_verbose} = $verbose;
  }
  return $self->{_verbose};
}
############################################################

sub use_score{
    my ($self,$boolean) = @_;
    if ( $boolean ){
	$self->{_use_score} = $boolean;
    }
    return $self->{_use_score};
}

############################################################

#sub check_tree{
#  my ($self, $newnode, $tree , $root) = @_;

#  my $verbose = $self->verbose;
  
#  # tree is a ref to a list of leaves
#  my @ids = map { $_->transcript->dbID } @$tree;
#  print STDERR "checking tree @ids:\n";
#  my @nodes_to_check = @$tree;

#  # as we check we regenerate the leaves as needed
#  my $newleaves;

#  # keep track of which ones are leaves 
#  my %is_leaf;
#  foreach my $node (@nodes_to_check){
#    $is_leaf{$node} = 1;
#  }
  
#  # and which ones are added to the re-generated list of leaves
#  my %added;

#  # keep track of parents being added to the @nodes_to_check
#  my %added_parent;

#  # keep track whether the new element has been placed at all or not
#  my $placed=0;
  
# NODE:
#  while( @nodes_to_check ){
#    my $node = shift @nodes_to_check;
    
#    ############################################################
#    # skip if the node has been tagged as visited
#    ############################################################
#    if ( $self->is_visited($node) ){
#      if ( $is_leaf{$node} && !$added{$node} ){
#	push( @$newleaves, $node );
#      }
#      next NODE;
#    }
#    print STDERR "check_tree(): node = ".$node->transcript->dbID."\n" if $verbose;
    
#    my $result = $self->check_node($newnode,$node);
    
#    ############################################################
#    # if no-overlap, we jump to the next leaf
#    ############################################################
#    if ( $result eq 'no-overlap' ){
#      print STDERR "no-overlap\n" if $verbose;
#      # recover the leaf if it hasn't been added before
#      if ( $is_leaf{$node} &&  !$added{$node} ){
#	push(@$newleaves,$node);
#	$added{$node} = 1;
#      }	  
#      next NODE;
#    }
#    ############################################################
#    # if no-overlap, we jump to the next leaf
#    ############################################################
#    elsif ( $result eq 'continue' ){
#      # recover the leaf if it hasn't been added before
#      print STDERR "check_tree(): CONTINUE - not sure what to do here\n" if $verbose;
#      if ( $is_leaf{$node} &&  !$added{$node} ){
#	push(@$newleaves,$node);
#	$added{$node} = 1;
#      }	  
#      next NODE;
#    }
    
#    ############################################################
#    # if placed, 
#    ############################################################
#    elsif ( $result eq 'placed' ){
#      $placed++;
      
#      print STDERR "placed\n";
#      # check whether newnode is a new leaf
#      if ( !$root && !$newnode->is_extended && !$newnode->is_included ){
#	if ( $added{$newnode} ){
#	  print STDERR "newnode is already a leaf\n";
#	}
#	unless ( $added{$newnode} ){
#	  print STDERR "adding newnode as leaf\n" if $verbose;
#	  push(@$newleaves,$newnode);
#	  $added{$newnode} = 1;
#	}
#      }
      
#      # recover the leaf if it hasn't been added before, and if it is still a leaf
#      if ( $is_leaf{ $node } && !$added{$node} && !$node->is_extended ){
#	push(@$newleaves,$node);
#	$added{$node} = 1;
#      }
#      next NODE;
#      #return 'placed';
#    }

#    ############################################################
#    # place here if it extends a node. 
#    # 'Delete' leaf and add new leaf.
#    # Tag all the extension parents
#    ############################################################
#    elsif ($result eq 'extension'){
      
#      ############################################################
#      # place here
#      print STDERR "extension: adding new extension\n" if $verbose;
#      $newnode->add_extension_parent($node);
#      $node->is_extended(1);
      
#      ############################################################
#      # add to the root if it is included in it 
#      if ( $root && $self->compare( $newnode, $root) eq 'inclusion' ){
#	print STDERR "adding inclusion of ".$newnode->transcript->dbID." into ".$root->transcript->dbID."\n";
#	$root->add_inclusion_child( $newnode );
#	$newnode->is_included(1);
#	$self->tag_extension_ancestors($root);
#      }
#      # this is a new leaf of this inclusion sub-tree (or main tree, where there is no root)
#      unless ( $added{$newnode} ){
#	print STDERR "adding newnode as leaf\n" if $verbose;
#	push(@$newleaves,$newnode);
#	$added{$newnode} = 1;
#      }
#      # tag this node and all the extension ancestors recursively
#      $self->tag_extension_ancestors($node);
#      $placed++;
#      next NODE;
#    }
    
#    ############################################################
#    # if there is a clash, 
#    # check first the inclusion trees
#    ############################################################
#    elsif ( $result eq 'clash' ){
      
#      if ( defined $node->inclusion_tree ){
	
#	print STDERR "clash -> checking inclusion-tree\n" if $verbose;
#	my $result = $self->check_tree( $newnode, $node->inclusion_tree );
	
#	############################################################
#	# if we are in the tree without root ( the main list of leaves)
#	# and our placed node is not extended nor included, it is a new leaf
#	if ( $result eq 'placed' ){
#	  $placed++;
	  
#	  # check whether newnode is a new leaf
#	  if ( !$root && !$newnode->is_extended && !$newnode->is_included ){
#	    if ( $added{$newnode} ){
#	      print STDERR "newnode is already a leaf\n";
#	    }
#	    unless ( $added{$newnode} ){
#	      print STDERR "adding newnode as leaf\n" if $verbose;
#	      push(@$newleaves,$newnode);
#	      $added{$newnode} = 1;
#	    }
#	  }
	  
#	  # recover the leaf if it hasn't been added before, and if it is still a leaf
#	  if ( $is_leaf{ $node } && !$added{$node} && !$node->is_extended ){
#	    push(@$newleaves,$node);
#	    $added{$node} = 1;
#	  }
#	  next NODE;
#	}
#	############################################################
#	# do breadth-first on the parents
#	elsif( $result eq 'clash' ){
	  
#	  if ( @{$node->extension_parents} ){
#	    my @parents = map { $_->transcript->dbID }  @{$node->extension_parents};
#	    print STDERR "adding extension parents ( @parents) to the list to check\n";
#	    foreach my $parent ( @{$node->extension_parents} ){
#	      unless ( $added_parent{ $parent } ){
#		push( @nodes_to_check, $parent );
#		$added_parent{$parent}++;
#	      }
#	    }
#	  }
#	  # recover the leaf if it hasn't been added before
#	  if ( $is_leaf{$node} &&  !$added{$node} && !$node->is_extended){
#	    push(@$newleaves,$node);
#	    $added{$node} = 1;
#	  }	  
#	}
#	############################################################
#	# nothing else should happen here - so this is a check
#	else{
#	  print STDERR "result: $result\n";
#	  print STDERR "<<< STRANGE - didn't expect to get an extension here \n" if $verbose;
#	}
	
#      }
#      ############################################################
#      # if the inclusion tree is not valid, skip to the next node to check
#      elsif( !defined $node->inclusion_tree && @{$node->inclusion_children} ){
#	next NODE;
#      }
#      ############################################################
#      # else do breadth-first on the parents
#      else{ 
	
#	if ( @{$node->extension_parents} ){
#	  my @parents = map { $_->transcript->dbID }  @{$node->extension_parents};
#	  print STDERR "adding extension parents ( @parents ) to the list to check\n";
#	  foreach my $parent ( @{$node->extension_parents} ){
#	    unless ( $added_parent{ $parent } ){
#	      push( @nodes_to_check, $parent );
#	      $added_parent{$parent}++;
#	    }
#	  }
#	}
#	# recover the leaf if it hasn't been added before
#	if ( $is_leaf{$node} &&  !$added{$node} ){
#	  push(@$newleaves,$node);
#	  $added{$node} = 1;
#	}	  
#	next NODE;
#      }
#    }
    
#  } # end of NODE
  
#  # this is our new tree
#  @$tree = @$newleaves;
#  $newleaves = [];
  
#  if ($placed){
#    return 'placed';
#  }
#  else{
#    return 'continue';
#  }
#}


1;
