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
									   -exact_merge => $value,
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
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;

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

  print STDERR "we've passed ".scalar(@$transcripts)." trans to ClusterMerge\n";
  $self->input_transcripts(@{$transcripts});
  
  if (defined($exact_merge)){
    $self->exact_merge($exact_merge);
  }

  if ( $EST_GENEBUILDER_MERGE eq 'fuzzy_semiexact_merge' ){ 
    $self->_merge_type('fuzzy_semiexact_merge');
    $self->_mismatch_allowed(7);
  }
  elsif( $EST_GENEBUILDER_MERGE eq 'semiexact_merge' ){
    $self->_merge_type('semiexact_merge');
    $self->_mismatch_allowed(2);
  }
  elsif( $EST_GENEBUILDER_MERGE eq 'merge_allow_gaps' ){
    $self->_merge_type('merge_allow_gaps');
    $self->_mismatch_allowed(7);
  }
  else{
    $self->throw("this won't work, you must define ESTConf::EST_GENEBUILDER_MERGE for your ests/cdnas to be compared!")
  }

  return $self;
}

############################################################
#
# RUN METHOD
#
############################################################

=head2 run

 Function:   main magic and witchcraft on the transcripts. 
  it fills up the holder $self->output with transcript objects
  
=cut

sub run {
  my $self = shift;
  my @transcripts = $self->input_transcripts;
	 
  # cluster the transcripts
  my @transcript_clusters = $self->_cluster_Transcripts(@transcripts);
  print STDERR scalar(@transcript_clusters)." $transcript_clusters[0] clusters returned from _cluster_Transcripts\n";
	 
  # find the lists of clusters that can be merged together according to consecutive exon overlap
  my @lists = $self->link_Transcripts( \@transcript_clusters );
  print STDERR scalar(@lists)." lists returned from link_Transcripts\n";

  # merge each list into a putative transcript
  my @merged_transcripts  = $self->_merge_Transcripts(\@lists);
  print STDERR scalar(@merged_transcripts)." transcripts returned from _merge_Transcripts\n";
	 
  $self->output(@merged_transcripts);
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

  unless (@mytranscripts){
    return;
  }

  # first sort the transcripts
  my @transcripts = sort by_transcript_high @mytranscripts;
  print STDERR "clustering ".scalar(@transcripts)." transcripts\n";

  # create a new cluster 
  my $cluster=Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster->new();
  my $count = 0;
  my @cluster_starts;
  my @cluster_ends;
  my @clusters;
  
  # put the first transcript into these cluster
  $cluster->put_Transcripts( $transcripts[0] );
  

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
  if ( $tran->start_Exon->strand == 1){
    $low = $tran->start_Exon->start;
  }
  else{
    $low = $tran->end_Exon->start;
  }
  return $low;
}


############################################################
#
# CLUSTERING OF TRANSCRIPTS BY CONSECUTIVE EXON OVERLAP
#
############################################################

=head2 link_Transcripts

 description: reads all the transcripts per cluster and merges those that
              are redundant with each other, producing brand new transcripts from that,
              see below the description of the algorithm

 Algorithm:

    pre-processing: cluster transcripts by genomic extent overlap (rough clustering), 
                    done in _cluster_Transcripts_by_genomic_range

    for each cluster of transcripts C

      pre-processing: sort the transcripts in C in ascending order by their start coordinate. If this coincides
                      we sort desc by the end coordinate.

      L = { all the lists that we are going to create }

      for each transcript t_i in C {

        start list list(i)(0) with t_i
	list(i)(a) represents every list that begins with t_i
 
        for each t_j in C, j>i {
  
          for each list(i)(a), a = 0,...,N {
  
            compare t_j with every element t_k in list(i)(a)  
	  }
	  there are 3 possibilities: 
	      
	  1. if t_j merges with at least one element in list(i)(a) and
	     with those that do not merge, it does even not overlap,
	     ==> add t_j to the end of list(i)(a)
	         go to the next t_j
		
	  2. if t_j has the properties as in 1) but with a proper sub-list of list(i)(a)
             ==> create a new list list(i)(a_max + 1) = sublist of list(i)(a) and add t_j at the end
                 go to the next t_j
  
          3. if not 1. and not 2. ==> go to the next list list(i)(a+1)
	}

        put the lists list(i)(a) a=0,...,N in L	
        go to the next t_i
      }

      remove lists embedded in longer lists (e.g. 3->4 is embedded in 1->3->4 )	       
      sort the lists in L descending by the number of elements
      accept the first list list_0
      for each list_m in L m>0 {   
        if list_m is not included in (or equal to) list_n for n = 0,...,m-1 ==> accept list_m
      }
   
=cut

############################################################

sub link_Transcripts{
  my ($self,$transcript_clusters) = @_;

  my @final_lists;  

  # look in each cluster
 CLUSTER:
  foreach my $cluster ( @{ $transcript_clusters} ){

     
    my %overlap_matrix;
    $self->matrix(\%overlap_matrix);

    # keep track of all the lists
    my @lists;
    
    # keep track of every list starting in each transcript:
    my %lists;

    # get the transcripts in this cluster
    my @transcripts = @{$cluster->get_Transcripts};
    
    # sort the transcripts by the left most coordinate in ascending order
    # for equal value, order by the right most coordinate in descending order
    @transcripts = sort { my $result = ( $self->transcript_low($a) <=> $self->transcript_low($b) );
	  		  if ($result){
		 	    return $result;
			  }
			  else{
			    return ( $self->transcript_high($b) <=> $self->transcript_high($a) );
			  }
			} @transcripts;
    
    print STDERR "Cluster: (sorted transcripts)\n";
    foreach my $t ( @transcripts ){
      Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($t);
    }

    # for each transcript
  TRAN1:
    for (my $i=0; $i<scalar(@transcripts); $i++){

      # store this one in the first list
      my @first_list = ();
      push (@first_list, $transcripts[$i]);

      # store all the lists created started in this transcript:
      my @current_lists = ();
      push (@current_lists, \@first_list );
            
      # go over the rest
    TRAN2:
      for (my $j=$i+1; $j<scalar(@transcripts); $j++){
	
	# loop over each of the lists{i}
      LIST:
	foreach my $list ( @current_lists ){
	  
	  # check whether this trans can be linked to this list or to a proper sublist
	  # or maybe to none of the above
	  my $new_list = $self->_test_for_link( $list, $transcripts[$j] );
	  
	  unless ( $new_list ){
	    next LIST;
	  }
	  # if it returns the same list
	  if ( $new_list == $list ){
	    
	    # add the it to the current list:
	    push( @$list, $transcripts[$j] );
	    next TRAN2;
	  }
	  # it could return a proper sub list
	  elsif( $new_list != $list ){
	    
	    # add it to the sublist:
	    push ( @$new_list , $transcripts[$j] );
	    
	    # add this new list to the $lists{$i}
	    push ( @current_lists, $new_list );
	  }
	  
	} # end of LIST
      }   # end of TRAN2
      
      #put @current_lists in @lists
      push ( @lists, @current_lists );

      #print STDERR "current lists:\n";
      #foreach my $list ( @current_lists ){
      # foreach my $t (@$list){
      #  Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($t);
      # }
      #}
      
    }  
    
    # remove lists embedded in longer lists (e.g. 3->4 is embedded in 1->3->4 )	       
    
    # sort the lists in descending by the number of elements
    my @sorted_lists = map  { $_->[1] } sort { $b->[0] <=> $a->[0] } map  { [ scalar( @{$_} ), $_] } @lists;
    
    my @accepted_lists;

    # accept the longest list
    push ( @accepted_lists, shift @sorted_lists );
    
    # check the rest
  LIST:
    while ( @sorted_lists ){
      my $list = shift @sorted_lists;
      my $found_embedding = 0;
      
    ACCEPTED_LIST:
      foreach my $accepted_list ( @accepted_lists ){
	if ( $self->_check_embedding( $list, $accepted_list ) ){
	  next LIST;
	}
      }
      
      # if we get to this point, it means that this list is genuine
      push( @accepted_lists, $list );
    }
    
    # store the lists for this cluster into the big final list:
    push (@final_lists, @accepted_lists);
    
  } # end of CLUSTER
		     
	     
  # @final_lists contain a list of listrefs, 
  # each one containing the transcripts that can merge with each other		      


  $self->sub_clusters( @final_lists );

  print STDERR "final lists:\n";
  my $count = 0;
  foreach my $list (@final_lists){
    $count++;
    print STDERR "list $count\n";
    foreach my $t ( @$list ){
      Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript( $t );
    }
  }



		       
  return @final_lists;	      
}

############################################################

=head2 _test_for_link

description: this method computes the largest proper sublist in $list 
             starting in the same element as $list
             to which $transcript can be linked
             
     Arg[1]: a listref with the list elements
     Arg[2]: a transcript object
     Return: it can return $list, a sub_list, or undef

=cut

sub _test_for_link{
  my ($self, $list, $transcript) = @_;
  
  my $comparator = Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptComparator->new( -comparison_level         => 3,
										      -exon_match               => 0,
										      -splice_mismatch          => 0,
										      -intron_mismatch          => 0,
										    );
  my $sublist;
  my $can_merge = 0;
  
  my %overlap_matrix = %{$self->matrix};

 LINK:
  foreach my $trans_link ( @$list ){
    my ($merge,$overlaps);
    
    if ( defined( $overlap_matrix{$transcript}{$trans_link} ) ){
      ($merge,$overlaps) = @{$overlap_matrix{$transcript}{$trans_link}};
      print STDERR "using cached matrix[ $transcript ][ $trans_link ] = ( $merge,$overlaps )\n";
    }
    else{
      ($merge,$overlaps) = $comparator->compare($trans_link, $transcript, $self->_merge_type, $self->_mismatch_allowed);
      $overlap_matrix{$transcript}{$trans_link} = [$merge,$overlaps];
      print STDERR "calculating matrix[ $transcript ][ $trans_link ] = ( $merge,$overlaps )\n";
    }
    
    if ($merge == 1 ){
      $can_merge = 1;
    }
    # the transcript can be appended to this list if
    # it merges to at least one member of the list
    # and it does not overlap with those member with which
    # it cannot merge
    if ( $merge == 1 || $overlaps == 0 ){
      
      # it links to this one
      push ( @{$sublist}, $trans_link );
    }
    else{
      # stop searching as soon as we find an element to which it does not link
      last LINK;
    }
  }
  
  # update the overlap matrix
  $self->matrix(\%overlap_matrix);
 
  if ( $can_merge == 0 || scalar( @$sublist ) == 0 ){
    return undef;
  }
  elsif ( scalar( @$sublist ) > 0 &&  scalar( @$sublist ) < scalar( @$list ) ){
    return $sublist;
  }
  elsif ( scalar( @$sublist ) == scalar( @$list ) ){
    return $list;
  }
  else{
    $self->throw("I didn't expect to be here!!");
  }
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

=head2 _check_embedding

     Arg[1]: a listref of transcripts
     Arg[2]: a listref of transcripts which potentially contains Arg[1]
description: this method checks whether Arg[1] is included in Arg[2], in the following sense
             if list2 = ( 1,3,4,6 ) and list1 = ( 3,4) ==> list2 includes list1. 
             It mimics the idea of the Boyer-Moore algorithm

=cut

sub _check_embedding{
  my ( $self, $list, $bigger_list ) = @_;
  
  my @list = @$list;
  my @bigger_list = @$bigger_list;
  unless ( scalar(@list) <= scalar( @bigger_list ) ){
      $self->throw("problem with the sorting, please check!");
  }
  
  # with a bit of preprocessing this can be made pretty quick
  
  # last element in list:
  my $last_in_list = $list[$#list];

  # store the positions where this element occurs in bigger_list starting from the right
  my @positions;
  for (my $k= $#bigger_list; $k>=0; $k-- ){
      if ( $bigger_list[$k] == $last_in_list ){
	  push( @positions, $k);
	  last;
      }
  }

  # if it does not occur in bigger_list, return false
  unless (@positions){
      return 0;
  }

  # else, start matching from that position:
  my $i = $#list;  
  my $j = shift @positions;

  while ( $j >=0 ){
            
      while ( $list[$i] == $bigger_list[$j] && $i >=0 && $j >= 0){
	  $i--;
	  $j--;
      }

      if ( $i == -1 ){
	  # list is embedded in bigger_list
	  return 1;
      }
      elsif ( @positions ){
	  # start looking from the next occurrence of $last_in_list in @bigger_list
	  $j = shift @positions;
	  $i = $#list;
      }
      else{
	  # no more occurrences
	  return 0;
      }
  }
  
  # sorry, we got to the end of bigger_list, and no match was found
  return 0;
}

############################################################


=head2 _test_for_Merge
 Function: this function is called from link_Transcripts and actually checks whether two transcripts
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
#
# MERGE TRANSCRIPTS
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
	
    # $list is an arrayref of the ests/cdnas that we can merge
    my @merged_transcripts;

  LIST:
    foreach my $list ( @$lists ){
      
      my @allexons;
      my %exon2transcript;			
      my %is_first;
      my %is_last;			
      
      # collect all exons
      foreach my $tran (@{ $list }){
	
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
      
      # cluster the exons
      my $first_cluster_list = $self->_cluster_Exons( @allexons );
      
      # set start and end of the clusters (using the info collected above)
      my $cluster_list = $self->_set_splice_Ends($first_cluster_list,\%exon2transcript,\%is_first,\%is_last);
      
      # we turn each cluster into an exon and create a new transcript with these exons
      my $transcript    = Bio::EnsEMBL::Transcript->new();
      my @exon_clusters = $cluster_list->sub_SeqFeature;
      
      foreach my $exon_cluster (@exon_clusters){
	
	my $new_exon = Bio::EnsEMBL::Exon->new();
	$new_exon->start ($exon_cluster->start );
	$new_exon->end   ($exon_cluster->end   );
	
	# set the strand to be the same as the exon_cluster
	$new_exon->strand($exon_cluster->strand);
	
	my %evidence_hash;
	my %evidence_obj;
	foreach my $exon ( $exon_cluster->sub_SeqFeature ){
	  $self->_transfer_Supporting_Evidence($exon,$new_exon);
	}
	
	$transcript->add_Exon($new_exon);
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
   my $fused  = 0;
   
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
   
   # if it is not fused, simply collect the exons
   if ( $fused != 1 ){
     push( @exon_list, $cluster->sub_SeqFeature);
   }    
   
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
    
    $new_end = shift( @ends );
    $max_end = $end{ $new_end };
    
    # if we have too little exons to obtain the end, take the original value
    if ( !defined $max_end ){
      print STDERR "In last position, cluster end wins!\n";
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
	print STDERR "Reset: new_start: $new_start\t new_end: $new_end\n";
      }
      else{
	## ok, you tried to change, but you got nothing, what can we do about it?
	print STDERR "Could not reset the start and end coordinates\n";
	if ( $other_end <= $other_start ){
	  print STDERR "Sorry will have to put the end = end of cluster\n";
	  $new_end  = $cluster->end;
	  if ( $new_start >= $new_end ){
	    print STDERR "Last resort: I'm afraid we'll also have to put start = start of cluster, good luck\n";
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
#
# GET/SET METHODS
#
############################################################

sub _merge_type{
  my ($self, $type ) = @_;
  if ($type){
    $self->{_merge_type} = $type;
  }
  return $self->{_merge_type};
}

############################################################

sub _mismatch_allowed{
  my ($self, $mismatch) = @_;
  if ($mismatch){
    $self->{_mismatch_allowed} = $mismatch;
  }
  return $self->{_mismatch_allowed};
}

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


1;
