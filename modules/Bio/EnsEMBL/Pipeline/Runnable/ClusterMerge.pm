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

    my $gene_machine = Bio::EnsEMBL::Pipeline::Runnable::CusterMerge->new(
									  -transcripts => \@transcripts,
									  -exact_merge => $value,
									  );
    $gene_machine->run;

    my @new_transcripts = $gene_machine->output;


=head1 DESCRIPTION

ClusterMerge takes a set of transcripts, it clusters them first and them
create sets of transcripts that can all merge with each other. This merging can be
exact, so that exon boundaries must match exactily, or fuzzy where 
exon boundaries do not necessarily coincide. The latter case could be useful for
ests, whereas the former case is advisable with full length cdnas.
The output is given in the form of the resulting transcripts.

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

# config file; parameters searched for here if not passed in as @args
use Bio::EnsEMBL::Pipeline::ESTConf;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  
  my( $transcripts, $exact_merge ) = $self->_rearrange([qw(
							   TRANSCRIPTS
							   EXACT_MERGE
							  )], 
						       @args);
  
  
  unless( $transcripts ){
    $self->warn("No transcripts passed it. We cannot go further");
    exit(0);
  }
  $self->input_transcripts(@{$transcripts});
  
  if (defined($exact_merge)){
    $self->exact_merge($exact_merge);
  }
  
  return $self;
}

############################################################

=head2 run

Usage   :   $self->run
 Function:   main magic and witchcraft on the transcripts. 
  it fills up the holder $self->output with transcript objects
  
=cut

sub run {
  my $self = shift;
  my @transcripts = $self->input_transcripts;
	 
  # cluster the transcripts
 my @transcript_clusters = $self->_cluster_Transcripts(@transcripts);
	 print STDERR scalar(@transcript_clusters)." clusters returned from _cluster_Transcripts\n";
	 
	 # merge the transcripts in each cluster according to consecutive exon overlap
	 my @merged_transcripts  = $self->_merge_Transcripts(\@transcript_clusters);
	 print STDERR scalar(@merged_transcripts)." transcripts returned from _merge_Transcripts\n";
	 
	 $self->output(@merged_transcripts);
	
}

############################################################
#
# METHODS CALLED BY THE RUN METHOD
#
############################################################

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
  
  my @forward_clusters = $self->_cluster_Transcripts_by_genomic_range( @forward_transcripts );
  my @reverse_clusters = $self->_cluster_Transcripts_by_genomic_range( @reverse_transcripts );
  my @clusters;
  push( @clusters, @forward_clusters);
  push( @clusters, @reverse_clusters);
  
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
  my @transcripts = sort by_transcript_high @mytranscripts;

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
  my $high;
  $tran->sort;
  if ( $tran->start_Exon->strand == 1){
    $high = $tran->end_Exon->end;
  }
  else{
    $high = $tran->start_Exon->end;
  }
  return $high;
}

############################################################

=head2 transcript_low

Description: it returns the lowest coordinate of a transcript

=cut

sub transcript_low{
  my ($self,$tran) = @_;
  my $low;
  $tran->sort;
  if ( $tran->start_Exon->strand == 1){
    $low = $tran->start_Exon->start;
  }
  else{
    $low = $tran->end_Exon->start;
  }
  return $low;
}

############################################################





#########################################################################

=head2 _merge_Transcripts

 Function: reads all the transcripts per cluster and merges those that
           are redundant with each other, producing brand new transcripts from that,
           see below the description of the algorithm

Algorithm:
    
    for each cluster of transcripts T={t}
      for each transcript t1
        loop over the rest t2
           loop over those in the t1-set and only add t2 to this set if: 
             at least merges to one in the set and with those that does not 
             merge it should not overlap at all (to avoid merging things that shouldn't merge)
     
           allow for transcripts to be used twice
      proceed as above with the next transcript in the list, 


    once we've gone through the whole lot, merge the chosen transcripts and put the resulting 
    transcript in an array.
    
    check that the resulting transcripts do not merge among themselves, 
    we merge and create new transcripts out of those that still can merge
    
=cut

sub _merge_Transcripts{
  my ($self,$ref_transcript_clusters,$strand) = @_;
    
  my @total_merged_transcripts;
  
  # look in each cluster
  my $count =0;

 CLUSTER:
  foreach my $cluster ( @{ $ref_transcript_clusters } ){
    
    $count++;
    
    # keep track of the transcripts originating each newly created one
    # each $origin_list{$new_transcript} is an array of transcripts
    my %origin_list;

    # get the transcripts in this cluster
    my @transcripts = $cluster->get_Transcripts;
    
    # sort the transcripts by the number of exons in descending order
    # for equal number of exons, order according to start coordinate
    # this is crucial!!
    @transcripts = sort { my $result = ( scalar( @{$b->get_all_Exons} ) <=> scalar( @{$a->get_all_Exons} ) );
			  if ($result){
			      return $result;
			  }
			  else{
			      return ( $a->start_Exon->start <=> $b->start_Exon->start );
			  }
		      } @transcripts;
    
    ##############################
    #print STDERR "\nNew Cluster:\n";
    #foreach my $tran (@transcripts){
    #  print STDERR "transcript: $tran\n";
    #  foreach my $exon ( $tran->get_all_Exons ){
    #	print STDERR $exon->start.":".$exon->end."  ";
    #  }
    #  print STDERR "\n";
    #}
    ##############################

    # now we loop over the transcripts, 
    my %is_merged;
    
    # store the transcripts we create
    my %created_transcripts;
    my %chosen_list;
    
    # for each transcript
  TRAN1:
    for (my $u=0; $u<scalar(@transcripts); $u++){

      # store the transcripts that can merge
      my @current_list = ();
      push (@current_list, $transcripts[$u]);
      
      # go over the rest
    TRAN2:
      for (my $v=$u+1; $v<scalar(@transcripts); $v++){
	
	my $overlap_ifnot_merged   = 0;
	my $merge_to_current       = 0;
	
	# loop over the accepted transcripts
	for (my $w=0; $w<scalar(@current_list); $w++){
	  
	  # in order to merge a transcript...
	  #print STDERR "comparing $current_list[$w] ($w) and $transcripts[$v] ($v)\n";
	  
	  my ($merge,$overlaps) = $self->_test_for_Merge( $current_list[$w], $transcripts[$v] );
	  
	  # ...there must be at least one in @current_list to which $transcripts[$v] merges
	  #print STDERR "merge = $merge, overlaps = $overlaps\n";
	  if ( 1 == $merge ){
	    $merge_to_current = 1;
	  }

	  # ...and it should not overlap with those to which it does not merge
	  unless (1 == $merge){
	    $overlap_ifnot_merged += $overlaps;
	  }
        }
	
	# ... and then it can merge to the list in @current_list
	if ( 1 == $merge_to_current && 0 == $overlap_ifnot_merged ){
	  push (@current_list, $transcripts[$v]);
	}
	
      } # end of TRAN2
      
      # keep track of the transcripts that are going to be used
      $chosen_list{ $u } = \@current_list;

    }   # end of TRAN1

    # create a new transcript with those merged above
    for (my $u=0; $u<scalar(@transcripts); $u++){
      $created_transcripts{$u} = $self->_produce_Transcript( $chosen_list{ $u }, $strand );
    }
 
    # at this point we have $created_transcripts{u=0,1,...,$#transcripts} transcripts
    my @created_transcripts = values ( %created_transcripts );
    
#    # print created transcripts
#    @created_transcripts  = sort { my $result = ( scalar( $b->get_all_Exons ) <=> scalar( $a->get_all_Exons ) );
#				     if ($result){
#				       return $result;
#				     }
#				     else{
#				       return ( $a->start_Exon->start <=> $b->start_Exon->start );
#				     }
#				   } @created_transcripts; 

#    print STDERR "\ntranscripts created:\n";
#    for (my $u=0; $u<scalar(@transcripts); $u++){
#      print STDERR "tran ".$created_transcripts{$u}.":";
#      foreach my $e1 ( $created_transcripts{$u}->get_all_Exons ){
#	print STDERR $e1->start.":".$e1->end."\t";
#      }
#      print STDERR "\nFrom the following ones:\n";
      
#      foreach my $tran ( @{ $chosen_list{$u} } ){
#	print STDERR "   tran ".$tran.":";
#	foreach my $e1 ( $tran->get_all_Exons ){
#	  print STDERR $e1->start.":".$e1->end."\t";
#	}
#	print STDERR "\n";
#      }
#    }

    # we process again the transcripts for possible residual merges
    my @merged_transcripts = $self->_check_for_residual_Merge(\@created_transcripts,$strand);

    ############################################################
    # why do we still get residual merge?
    # maybe we need to refine the algorithm so that it comes
    # clean first time around
    ############################################################

    # We put here yet another check, to make sure that we really get rid of redundant things
    
    # if we have produce more than one transcript, check again
    
    if ( scalar( @merged_transcripts ) > 1 ){
	my @merged_transcripts2 = $self->_check_for_residual_Merge(\@merged_transcripts,$strand);
	
	# replace then the previous list by this new one:
	@merged_transcripts = ();
	push ( @merged_transcripts, @merged_transcripts2);
      } 
  
    # check for the single exon transcripts ( dandruff ), they are no good
    print STDERR "Eliminating single exon transcripts this shouldn't be in ClusterMerge, it should be in EST_GeneBuilder\n";
    my @final_merged_transcripts;
    foreach my $tran ( @merged_transcripts ){
	if ( scalar ( @{$tran->get_all_Exons} ) == 1 ){
	    next;
	}
	else{
	    push( @final_merged_transcripts, $tran );
	}
    }
    
    # keep track of everything it is created for every cluster
    push (@total_merged_transcripts, @final_merged_transcripts);
    
    # empty arrays to save memory
    @merged_transcripts       = ();
    @final_merged_transcripts = ();
    
  }     # end of CLUSTER		       
		       
  return @total_merged_transcripts;
		      
}

#########################################################################

=head2 _check_for_residual_Merge
 Function: this function is called from _merge_Transcripts and checks whether after merging the transcript
           there is any residual merge
 Returns : It returns two numbers ($merge,$overlaps), where
           $merge = 1 (0) when they do (do not) merge,
           and $overlaps is the number of exon-overlaps.

=cut

sub _check_for_residual_Merge{

   my ($self, $created_tran_ref, $strand) = @_;
   my @created_transcripts = @$created_tran_ref;

   # first of all, sort them again, this is CRUCIAL!!!
   @created_transcripts = sort { my $result = ( scalar( @{$b->get_all_Exons} ) <=> scalar( @{$a->get_all_Exons} ) );
				     if ($result){
				       return $result;
				     }
				     else{
				       return ( $a->start_Exon->start <=> $b->start_Exon->start );
				     }
				   } @created_transcripts;
    

    # we process the transcripts, the results go in this hash
    my %new_transcripts;
    my $max_label = 0;
    my @labels;

    # for each created transcript
  TRAN: 
    foreach my $tran (@created_transcripts){
      my $found = 0;
      my $label = 0;

      # loop over the existing new_transcripts
    NEW_TRAN: 
      foreach my $k ( @labels ){

	# take each new_transcript
	my $new_tran = $new_transcripts{ $k };
	
#	# test for merge
#	print STDERR "comparing\n";
#	print STDERR $new_tran.":";
#	foreach my $e1 ( $new_tran->get_all_Exons ){
#	  print STDERR $e1->start.":".$e1->end."\t";
#	}
#	print STDERR "\n".$tran.":";
#	foreach my $e1 ( $tran->get_all_Exons ){
#	  print STDERR $e1->start.":".$e1->end."\t";
#	}
#	print STDERR "\n";
      
	my ($merge,$overlaps) = $self->_test_for_Merge( $new_tran, $tran );

	# if they can merge, merge them and put the new transcript in this place 
	if ( $merge == 1 ){
	  $found = 1;
	  my @list = ( $new_tran, $tran );
	  #print STDERR "They MERGE,\n";
	  #print STDERR "adding it to new_transcripts[ $label ] = $new_transcripts{ $label }\n";
	  my $new_transcript = $self->_produce_Transcript( \@list, $strand );
	  $new_transcripts{ $label } = $new_transcript;
	  #print STDERR "producing a new  new_transcripts[ $label ] = $new_transcripts{ $label }\n\n";
	  next TRAN;
	}
	else{
	  # go to the next cluster
	  $label++;
	}
      } # end of NEW_TRAN
      
      # if it does not merge to any transcript, create a new new_tran
      if ($found == 0) {
	$max_label = $label;
	push ( @labels, $label);
	#print STDERR "*** LABEL: $label ***\n";
	$new_transcripts{ $label } = $tran;
      }    
    
    }   # end of TRAN

    my @merged_transcripts = values( %new_transcripts );

   # ### print out the results
#    print STDERR "Transcripts created:\n";
#    foreach my $tran ( @merged_transcripts ){
#      print STDERR "tran ".$tran.":";
#      foreach my $e1 ( $tran->get_all_Exons ){
#	print STDERR $e1->start.":".$e1->end."\t";
#      }
#      print STDERR "\n";
#    }

  return @merged_transcripts;
}



#########################################################################

=head2 _test_for_Merge
 Function: this function is called from _merge_Transcripts and actually checks whether two transcripts
           inputs merge.
 Returns : It returns two numbers ($merge,$overlaps), where
           $merge = 1 (0) when they do (do not) merge,
           and $overlaps is the number of exon-overlaps.

=cut

sub _test_for_Merge{
  my ($self,$tran1,$tran2) = @_;
  my @exons1 = @{$tran1->get_all_Exons};
  my @exons2 = @{$tran2->get_all_Exons};	
 
  my $foundlink = 0; # flag that gets set when starting to link exons
  my $start     = 0; # start looking at the first one
  my $overlaps  = 0; # independently if they merge or not, we compute the number of exon overlaps
  my $merge     = 0; # =1 if they merge

  #print STDERR "comparing ".$tran1->dbID." ($tran1)  and ".$tran2->dbID." ($tran2)\n";


EXON1:
  for (my $j=0; $j<=$#exons1; $j++) {
  
  EXON2:
    for (my $k=$start; $k<=$#exons2; $k++){
    #print STDERR "comparing ".($j+1)." and ".($k+1)."\n";
	    
      # if exon 1 is not the first, check first whether it matches the previous exon2 as well, i.e.
      #                        ____     ____        
      #              exons1 --|____|---|____|------ etc... $j
      #                        ____________  
      #              exons2 --|____________|------ etc...  $k
      #
      if ($foundlink == 1 && $j != 0){
	if ( $k != 0 && $exons1[$j]->overlaps($exons2[$k-1]) ){
	  #print STDERR ($j+1)." <--> ".($k)."\n";
	  $overlaps++;
          next EXON1;
	}
      }
      
      # if texons1[$j] and exons2[$k] overlap go to the next exon1 and  next $exon2
      if ( $exons1[$j]->overlaps($exons2[$k]) ){
	#print STDERR ($j+1)." <--> ".($k+1)."\n";
        $overlaps++;
	
        # in order to merge the link always start at the first exon of one of the transcripts
        if ( $j == 0 || $k == 0 ){
          $foundlink = 1;
        }
      }          
      else {  
	# if you haven't found an overlap yet, look at the next exon 
	if ( $foundlink == 0 ){
	  next EXON2;
	}
	# leave if we stop finding links between exons before the end of transcripts
	if ( $foundlink == 1 ){
	  $merge = 0;
	  last EXON1;
	}
      }
      
      # if foundlink = 1 and we get to the end of either transcript, we merge them!
      if ( $foundlink == 1 && ( $j == $#exons1 || $k == $#exons2 ) ){
	
	# and we can leave
        $merge = 1;
	last EXON1;
      }
      # if foundlink = 1 but we're not yet at the end, go to the next exon 
      if ( $foundlink == 1 ){
	
	# but first check whether in exons2 there are further exons overlapping exon1, i.e.
        #                       ____________        
	#             exons1 --|____________|------ etc...
	#                       ____     ___  
	#             exons2 --|____|---|___|------ etc...
	# 
	my $addition = 0;
	while ( $k+1+$addition < scalar(@exons2) && $exons1[$j]->overlaps($exons2[$k+1+$addition]) ){
	  #print STDERR ($j+1)." <--> ".($k+2+$addition)."\n";
	  $overlaps++;
          $addition++;
	}      
	$start = $k+1+$addition;
	next EXON1;
      }    
    
    } # end of EXON2 
  
    if ($foundlink == 0){
      $start = 0;
    }
 
  }   # end of EXON1      

  # if we haven't returned at this point, they don't merge, thus
  return ($merge,$overlaps);
}
  


############################################################

=head2 _produce_Transcript
 Function: reads all the est2genome transcripts that can be merged and make a single transcript
           out of them
=cut

sub _produce_Transcript{
  my ($self,$merged,$strand) = @_;

  my @allexons;
  my %exon2transcript;			
  my %is_first;
  my %is_last;			
	       
  # collect all exons
  foreach my $tran (@{ $merged }){

    my @exons = @{$tran->get_all_Exons};

    # sort them in genomic order
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

  # cluster the exons
  my $first_cluster_list = $self->_cluster_Exons( @allexons );

  # set start and end of the clusters (using the info collected above)
  my $cluster_list = $self->_set_splice_Ends($first_cluster_list,\%exon2transcript,\%is_first,\%is_last,$strand);

  # we turn each cluster into an exon and create a new transcript with these exons
  my $transcript    = Bio::EnsEMBL::Transcript->new();
  my @exon_clusters = $cluster_list->sub_SeqFeature;

  foreach my $exon_cluster (@exon_clusters){

    my $new_exon = Bio::EnsEMBL::Exon->new();
    $new_exon->start ($exon_cluster->start );
    $new_exon->end   ($exon_cluster->end   );
    
    ###  dont't set strand yet, genomewise cannot handle that ###
    # we put it to 1 anyway, so that minigenomewise does not complain?
    # the real strand will be dealt with in the run method
    $new_exon->strand(1);

    my %evidence_hash;
    my %evidence_obj;
    foreach my $exon ( $exon_cluster->sub_SeqFeature ){
      $self->_transfer_Supporting_Evidence($exon,$new_exon);
    }

    $transcript->add_Exon($new_exon);
  }
 			
  return $transcript;
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
      my @overlap = $self->match($exon, $exon_cluster);    
      
      # Add to cluster if overlap AND if strand matches
      if ( $overlap[0] && ( $exon->strand == $exon_cluster->strand) ) { 
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

    Usage   :   $cluster_list = $self->_set_splice_Ends($cluster_list)
    Function:   Resets the ends of the clusters to the most frequent coordinate
   
=cut

sub _set_splice_Ends {
  my ($self, $cluster_list, $ref_exon2transcript_hash, $ref_is_first, $ref_is_last, $strand) = @_;

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
    my $fused  = 0;

    # keep track of the transcripts used in this cluster
    my %seen_transcript;
    my %exons_from_transcript;
    
    foreach my $exon ( $cluster->sub_SeqFeature ){

      push ( @{ $exons_from_transcript{ $exon2transcript{ $exon } } }, $exon );  

      if ( $seen_transcript{ $exon2transcript{$exon} } && 1 == $seen_transcript{ $exon2transcript{$exon} } ){
	#print STDERR "There is more than one exon from transcript $exon2transcript{ $exon }\n"; 
	$fused = 1;
	$need_to_recluster = 1;
      }
      $seen_transcript{ $exon2transcript{ $exon } } = 1;
    } # end of exon loop

    # if it is not fused, simply collect the exons
    if ( $fused != 1 ){
      push( @exon_list, $cluster->sub_SeqFeature);
    }    
    # if there is fussion going on (be-bop?), get rid of the bad guys and re-cluster
    elsif ( $fused == 1 ){

      my @exons = sort{ $b->length <=> $a->length } $cluster->sub_SeqFeature;
      
    EXONS:
      foreach my $exon ( @exons ){
	
      TRANS_IN_CLUSTER:
	foreach my $tran ( keys( %exons_from_transcript ) ){
	  if ( $exon2transcript{$exon} eq $tran ){
	    next;
	  }	  
	  my $overlap_count = 0;
	  foreach my $exon2 ( @{ $exons_from_transcript{ $tran } } ){
	    if ( $exon->overlaps($exon2) ){
	      $overlap_count++;
	    }
	  }
	  # if  $exon overlaps 2 or more exons from the same transcript, in this cluster ...
	  if ( $overlap_count >= 2 ){
	    
	    # ... and $exon is the only one from its transcript in this cluster
	    if ( scalar( @{ $exons_from_transcript{ $exon2transcript{$exon} } } ) == 1 ){
	      
	      # ... and $exon is at the edge of its transcript
	      if ( $is_first{ $exon } == 1 || $is_last{ $exon } == 1 ){
	         
		#print STDERR "ditching one exon\n";
		# then we ditch it and continue ...
		next EXONS;
	      }
	    }
	  }
	} # end of TRANS_IN_CLUSTER
			   
	# if it is good, it gets here
        push ( @exon_list, $exon );
      
      }   # end of EXONS
      
    }     # end of 'if fused == 1'
    
  }       # end of CLUSTERS

  # if needed, re-cluster
  if ( $need_to_recluster == 1){
    @exon_clusters   = ();
    #print STDERR " *** Reclustering ***\n";
    $cluster_list  = $self->_cluster_Exons( @exon_list );
    @exon_clusters = $cluster_list->sub_SeqFeature;
  }

  # at this point we (hopefully) have got rid of the bad fussion, we can set the splice ends on the new clusters
      
  # there is a GOTCHA, however. Once we have reclustered the exons, exons that were together before may be 
  # separated in different clusters, which means they will be part of different potential exons, for the same
  # transcript! However, it is possible that we have lost the actual connection between those two exons,
  # creating then a fake exon.
  # SOLUTION: to introduce here a check for links between the exons clusters, if every two neighbouring
  # clusters share a est2genome transcript, then it is o.k., otherwise, we should then split the transcript.
  # This means that we need to return two (or maybe more!) differentiated sets of clusters which will
  # be then used to create two ( or more ) new transcripts.

  # keep track of the cluster position
  my $position = 0;

CLUSTER:		     
  foreach my $cluster (@exon_clusters) {
    $position++;
    my %start;
    my $new_start;
    my $max_start = 0;
    
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
    $max_start = $start{ $new_start };

    # if we have too little exons to obtain the start, take the original value
    if ( $max_start == 0 ){
      $new_start = $cluster->start;
    }

    # the first cluster is a special case - potential UTRs, take the longest one.
    if( $position == 1) {
      $new_start = $cluster->start;  
    }

    #### now set the end coordinate ####

    my %end;
    my $new_end;
    my $max_end = 0;
    
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
    
    ## test:
    #print STDERR "ends: ";
    #foreach my $e (@ends){
    #  print STDERR $e."\t";
    #}
    #print STDERR "\n";
    
    $new_end = shift( @ends );
    $max_end = $end{ $new_end };
    
    # if we have too little exons to obtain the end, take the original value
    if ( $max_end == 0 ){
      print STDERR "In last position, cluster end wins!\n";
      $new_end = $cluster->end;
    }
    
    # the last cluster is a special case - potential UTRs, take the longest one.
    if( $position == scalar( @exon_clusters ) ) {
      $new_end = $cluster->end;  
    }
    
    #    print STDERR "new_start: $new_start\n";
    #    print STDERR "start array:\n";
    #    foreach my $s ( sort{ $start{ $b } <=> $start{ $a } } keys( %start ) ){
    #      print STDERR "start: $s\tnum_times: ".$start{ $s }."\t";
    #    }
    #    print STDERR "\n";
    #    print STDERR "new_end: $new_end\n";
    #    print STDERR "end array:\n";
    #    foreach my $e ( sort{ $end{ $b } <=> $end{ $a } } keys( %end ) ){
    #      print STDERR "end: $e\tnum_times: ".$end{ $e }."\t";
    #    }
    #    print STDERR "\n";
    

    ######  if new_start> new_end change them in turns until we get start < end ######
    # this can happen when there are many exons in a cluster but they vary greatly in size
    my $other_start = $new_start;
    my $other_end   = $new_end;
    my $trouble     = 0;
    my $stop_start  = 0;
    my $stop_end    = 0;

    while ( $other_start >= $other_end ){
      print STDERR "*** Trouble: $new_start >= $new_end ***\n";
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
	print STDERR "Reset: new_start: $new_start\t new_end: $new_end\n";
      }
      else{
	## ok, you tried to change, but you got nothing, what can we do about it?
	print STDERR "Could not reset the start and end coordinates\n";
	if ( $other_end <= $other_start ){
	  print STDERR "Sorry will have to put the end= ens of cluster\n";
	  $new_end   = $cluster->end;
	  if ( $new_start >= $new_end ){
	    print STDERR "Last resort: I'm afraid we'll also have to put start = start of cluster, good luck\n";
	  }
	}
      }
    }

    # reset the cluster start (that's the coordinate used in _produce_Transcript() )
    #print STDERR "reseting cluster start to: $new_start\n";
    $cluster->start($new_start);

    # reset the cluster end (that's the coordinate used in _produce_Transcript() )
    #print STDERR "reseting cluster end to: $new_end\n";
    $cluster->end($new_end);
    
  }
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
#
# GET/SET METHODS
#
############################################################

sub exact_merge{
 my ($self, $value);
 if (defined($value)){
     $self->{_exact_merge} = $value;
 }
 return $self->{_exact_merge};
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


1;
