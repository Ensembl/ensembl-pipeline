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

Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptComparator

=head1 SYNOPSIS
    
    my $level; # one number between 1 and 4
    my $comparator = Bio::EnsEMBL::Pipeline::::GeneComparison::TranscriptComparator->new(
                                                       -comparison_level         => $level,
                                                       -exon_match               => 0,
                                                       -splice_mismatch          => 1,
                                                       -intron_mismatch          => 1,
                                                                                        );


    my ($merge,$overlaps) = $comparator->compare($transcript1,$transcript2);

    $merge = 1 if the two transcripts are equivalent in the sense specified by the parameters above
             0 otherwise
    $overlaps is the number of exon overlaps found by the comparison
    

there are four parameters that we can pass to a comparison method:

exon_match = BOOLEAN ------------> TRUE if we want both transcripts to match 1-to-1 all their exons 
                                   Appropriate for comparison levels 1, 2 and 3

splice_mismatch = INT -----------> Maximum number of bases mismatch that we allow in the internal splice sites
                                   Alignment programs are sometimes not able to resolve some ambiguities in
                                   the splice sites, this might help to identify truly equivalent splice sites.

intron_mismatch = INT -----------> Maximum number of bases that we consider for an intron to be non-real.
                                   Any intron of this size or smaller will be allowed to be merged with exons covering it.
                                   The reason for this is that we do not expect very very small intron to be
                                   real

comparison_level = INT ----------> There are currently 4 comparison levels:
 
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



=head1 CONTACT

ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptComparator;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Pipeline::GeneComparison::GeneCluster;
use Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster;

@ISA = qw(Bio::EnsEMBL::Root);

######################################################################

sub new{
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  my ( $comparison_level, $exon_match, $splice_mismatch, $intron_mismatch ) = 
    $self->_rearrange([qw(COMPARISON_LEVEL
			  EXON_MATCH
			  SPLICE_MISMATCH
			  INTRON_MISMATCH
			 )],
		      @args);
  
  if (defined $comparison_level){
    $self->comparison_level($comparison_level);
  }
  else{
    $self->throw("you must define a comparison_level. See documentation for more info");
  }

  if ( defined $exon_match ){
    $self->exon_match($exon_match);
  }

  if (defined $splice_mismatch){
    $self->splice_mismatch($splice_mismatch);
  }

  if (defined $intron_mismatch){
    $self->intron_mismatch($intron_mismatch);
  }
  
  return $self;  
}

############################################################

sub comparison_level{
  my ($self, $level) = @_;
  if ( defined $level ){
     $self->{_comparison_level} = $level;
  }
  return $self->{_comparison_level};
}

sub exon_match{
  my ($self, $boolean) = @_;
  if ( defined $boolean ){
     $self->{_exon_match} = $boolean;
  }
  return $self->{_exon_match};
}

sub splice_mismatch{
  my ($self, $int) = @_;
  if ( defined $int ){ 
    $self->{_splice_mismatch} = $int;
  }
  return $self->{_splice_mismatch};
}

sub intron_mismatch{
  my ($self, $int) = @_;
  if ( defined $int ){
     $self->{_intron_mismatch} = $int;
  }
  return $self->{_intron_mismatch};
}

############################################################

=head2 compare
  Arg[1] and Arg[2] : 2 transcript objects to compare
  Arg[3]: the mode ( a string of text). Possible modes:

=cut 

sub compare{
  my ($self, $tran1, $tran2) = @_;
  
  my ($merge, $overlaps);

  # switch on comparison level
  if (   $self->comparison_level == 3 ){
    ($merge, $overlaps) = $self->_test_for_fuzzy_semiexact_Merge( $tran1, $tran2 );
  }
  elsif( $self->comparison_level == 2 ){
    ($merge, $overlaps) = $self->_test_for_semiexact_Merge(  $tran1, $tran2 );
  }
  elsif( $self->comparison_level == 4 ){
    ($merge, $overlaps) = $self->_test_for_Merge_allow_gaps( $tran1, $tran2 );
  }
  elsif( $self->comparison_level == 1 ){
    ($merge, $overlaps) = $self->_test_for_strict_merge(     $tran1, $tran2 );
  }
  return ($merge,$overlaps);
}
  
#########################################################################
# this function checks whether two transcripts merge
# with exact exon matches, except for
# possible mismatches in the extremal exons

sub _test_for_semiexact_Merge{
  my ($self,$est_tran,$ens_tran) = @_;

  # allowed_exterior_mismatch is the number of bases that we allow the first or last exon
  # of a transcript to extend beyond an overlapping exon in the other transcript
  # which is an internal exon. We default it to zero:
  my $allowed_exterior_mismatch;
  unless ($allowed_exterior_mismatch){
    $allowed_exterior_mismatch = 0;
  }

  my @exons1 = @{$est_tran->get_all_Exons};
  my @exons2 = @{$ens_tran->get_all_Exons};	
  
  @exons1 = sort {$a->start <=> $b->start} @exons1;
  @exons2 = sort {$a->start <=> $b->start} @exons2;

  my $foundlink = 0; # flag that gets set when starting to link exons
  my $start     = 0; # start looking at the first one
  my $merge     = 0; # =1 if they merge
  my $overlaps  = 0;
  
 EXON1:
  for (my $j=0; $j<=$#exons1; $j++) {
      
    EXON2:
      for (my $k=$start; $k<=$#exons2; $k++){
	#print STDERR "comparing j = $j : ".$exons1[$j]->start."-".$exons1[$j]->end." and k = $k : ".$exons2[$k]->start."-".$exons2[$k]->end."\n";
	  
	  # we allow some mismatches at the extremities
	  #                        ____     ____     ___   
	  #              exons1   |____|---|____|---|___|  $j
	  #                         ___     ____     ____     __
	  #              exons2    |___|---|____|---|____|---|__|  $k
	  
	  # if there is no overlap, go to the next EXON2
	  if ( $foundlink == 0 && !( $exons1[$j]->overlaps($exons2[$k]) ) ){
	      #print STDERR "foundlink = 0 and no overlap --> go to next EXON2\n";
	      next EXON2;
	  }
	  # if there is no overlap and we had found a link, there is no merge
	  if ( $foundlink == 1 && !($exons1[$j]->overlaps($exons2[$k]) ) ){
	      #print STDERR "foundlink = 1 and no overlap --> leaving\n";
	      $merge = 0;
	      last EXON1;
	  }	
	  
	  # the first exon can have a mismatch in the start 
	  if ( ( ($k == 0 && ( $exons1[$j]->start - $exons2[$k]->start ) <= $allowed_exterior_mismatch ) || 
		 ($j == 0 && ( $exons2[$k]->start - $exons1[$j]->start ) <= $allowed_exterior_mismatch )   )
	       && 
	       $exons1[$j]->end == $exons2[$k]->end 
	     ){
	    
	      # but if it is also the last exon
	    if ( ( ( $k == 0 && $k == $#exons2 )   || 
		   ( $j == 0 && $j == $#exons1 ) ) ){
	      
	      # we force it to match the start
	      if ( $exons1[$j]->start == $exons2[$k]->start ){
		$foundlink  = 1;
		$merge      = 1;
		$overlaps++;
		#print STDERR "merged single exon transcript\n";
		last EXON1;
	      }
	      # we call it a non-merge
	      else{
		$foundlink = 0;
		$merge     = 0;
		#print STDERR "non-merged single exon transcript\n";
		last EXON1;
	      }
	    }
	    else{
	      #else, we have a link
	      $foundlink = 1;
	      $start = $k+1;
	      $overlaps++;
	      #print STDERR "found a link\n";
	      next EXON1;
	    }
	  }
	  # the last one can have a mismatch on the end
	  elsif ( ( $k == $#exons2 && ( $exons2[$k]->end - $exons1[$j]->end ) <= $allowed_exterior_mismatch ) || 
		  ( $j == $#exons1 && ( $exons1[$j]->end - $exons2[$k]->end ) <= $allowed_exterior_mismatch ) 
		  &&
		  ( $foundlink == 1 )                  
		  &&
		  ( $exons1[$j]->start == $exons2[$k]->start ) 
		){
	    #print STDERR "link completed, merged transcripts\n";
	    $overlaps++;
	    $merge = 1;
	    last EXON1;
	  }
	  # the middle one must have exact matches
	  elsif ( ($k != 0 && $k != $#exons2) && 
		  ($j != 0 && $j != $#exons1) &&
		  ( $foundlink == 1)          &&
		  ( $exons1[$j]->start == $exons2[$k]->start ) &&
		  ( $exons1[$j]->end   == $exons2[$k]->end   )
		){
	    $overlaps++;
	    $start = $k+1;
	    #print STDERR "continue link\n";
	    next EXON1;
	  }
	  
	} # end of EXON2 
      
      if ($foundlink == 0){
	$start = 0;
      }
      
    }   # end of EXON1      
   if ( $self->exon_match ){
    if ( $merge == 1 && scalar(@exons1) == scalar(@exons2) && scalar(@exons1) == $overlaps ){
      return ( $merge, $overlaps );
    }
    else{
      return ( 0, $overlaps );
    }
  }

  


  return ($merge, $overlaps);
}

#########################################################################
# this function checks whether two transcripts merge
# with fuzzy exon matches: there is consecutive exon overlap 
# but there are mismatches of $self->splice_mismatch bases allowed at the edges of any exon pair

sub _test_for_fuzzy_semiexact_Merge{
  my ($self,$est_tran,$ens_tran) = @_;
  my $allowed_mismatch = 0;
  if ( defined $self->splice_mismatch ){
    $allowed_mismatch =  $self->splice_mismatch;
  }
  
  #print STDERR "=========== comparing ================\n";
  #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript( $est_tran );
  #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript( $ens_tran );
  
  my @exons1 = @{$est_tran->get_all_Exons};
  my @exons2 = @{$ens_tran->get_all_Exons};	
  
  @exons1 = sort {$a->start <=> $b->start} @exons1;
  @exons2 = sort {$a->start <=> $b->start} @exons2;

  my $foundlink = 0; # flag that gets set when starting to link exons
  my $start     = 0; # start looking at the first one
  my $merge     = 0; # =1 if they merge
  my $overlaps  = 0; # number of exon overlaps
  
 EXON1:
  for (my $j=0; $j<=$#exons1; $j++) {
    
  EXON2:
    for (my $k=$start; $k<=$#exons2; $k++){
      #print STDERR "comparing j = $j : ".$exons1[$j]->start."-".$exons1[$j]->end.
      #  " and k = $k : ".$exons2[$k]->start."-".$exons2[$k]->end."\n";
      
      # we allow some mismatches at the extremities
      #                        ____     ____     ___   
      #              exons1   |____|---|____|---|___|  $j
      #                         ___     ____     ____     ____
      #              exons2    |___|---|____|---|____|---|____|  $k
      
      # if there is no overlap, go to the next EXON2
      if ( $foundlink == 0 && !($exons1[$j]->overlaps($exons2[$k]) ) ){
	#print STDERR "foundlink = 0 and no overlap --> go to next EXON2\n";
	next EXON2;
      }
      # if there is no overlap and we had found a link, there is no merge
      elsif ( $foundlink == 1 && !($exons1[$j]->overlaps($exons2[$k]) ) ){
	#print STDERR "foundlink = 1 and no overlap --> leaving\n";
	$merge = 0;
	last EXON1;
      }	
      # the first exon can have a mismatch ( any number of bases ) in the start
      # and $allowed_mismatch bases mismatch at the end
      elsif ( ($k == 0 || $j == 0) ){
	
	# if one of them is a single-exon transcript
	# we allow any mismatch
	if( ( $k == 0 && $k == $#exons2 ) || 
	    ( $j == 0 && $j == $#exons1 )    ){
	  $foundlink = 1;
	  $overlaps++;
	  $merge = 1;
	  #print STDERR "found link\n";
	  #print STDERR "merged single exon transcript\n";
	  last EXON1;
	}
	elsif ( abs($exons1[$j]->end - $exons2[$k]->end)<= $allowed_mismatch ){
	  # else we force it to match the end (with a mismatch of $allowed_mismatch bases allowed)
	  
	  $foundlink = 1;
	  $overlaps++;
	  $start = $k+1;
	  #print STDERR "found link\n";
	  next EXON1;
	}
      }
      # the last one can have any mismatch on the end
      # but must have a match at the start (with $allowed_mismatch mismatches allowed)
      elsif ( ( $k == $#exons2 || $j == $#exons1 ) &&
	      ( $foundlink == 1 )                  &&
	      ( abs($exons1[$j]->start - $exons2[$k]->start)<= $allowed_mismatch ) 
	    ){
	#print STDERR "link completed, merged transcripts\n";
	$merge = 1;
	$overlaps++;
	last EXON1;
      }
      # the middle one must have exact matches
      # (up to an $allowed_mismatch mismatch)
      elsif ( ($k != 0 && $k != $#exons2) && 
	      ($j != 0 && $j != $#exons1) &&
	      ( $foundlink == 1)          &&
	      abs( $exons1[$j]->start - $exons2[$k]->start )<= $allowed_mismatch &&
	      abs( $exons1[$j]->end   - $exons2[$k]->end   )<= $allowed_mismatch
	    ){
	$overlaps++;
	$start = $k+1;
	#print STDERR "continue link\n";
	next EXON1;
      }
      
    } # end of EXON2 
    
    if ($foundlink == 0){
      $start = 0;
    }
    
  }   # end of EXON1      
  
  # check whether we only want the same number of exons:
  if ( $self->exon_match ){
    if ( $merge == 1 && scalar(@exons1) == scalar(@exons2) && scalar(@exons1) == $overlaps ){
      return ( $merge, $overlaps );
    }
    else{
      return ( 0, $overlaps );
    }
  }

  return ($merge, $overlaps);
}

#########################################################################
# this function checks whether two transcripts merge
# according to consecutive exon overlap (just overlap, without looking at the 
# exon positions) and it only considers 1-to-1 matches, so things like
#                        ____     ____        
#              exons1 --|____|---|____|------ etc... $j
#                        ____________  
#              exons2 --|____________|------ etc...  $k
#
# are considered a mismatch

sub _test_for_Simple_Merge{
  my ($self,$tran1,$tran2) = @_;
  my @exons1 = @{$tran1->get_all_Exons};
  my @exons2 = @{$tran2->get_all_Exons};	
 
  my $foundlink = 0; # flag that gets set when starting to link exons
  my $start     = 0; # start looking at the first one
  my $overlaps  = 0; # independently if they merge or not, we compute the number of exon overlaps
  my $merge     = 0; # =1 if they merge

  my $one2one_overlap = 0;
  my $one2two_overlap = 0;
  my $two2one_overlap = 0;
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
	  $two2one_overlap++;
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
	  $one2two_overlap++;
	  $overlaps++;
          $addition++;
	}      
	$start = $k+1+$addition;
	next EXON1;
      }    
      
    } # end of EXON2 
    
    # if you haven't found any match for this exon1, start again from the first exon2:
    if ($foundlink == 0){
      $start = 0;
    }
 
  }   # end of EXON1      

  # we only make them merge if $merge = 1 and the 2-to-1 and 1-to-2 overlaps are zero;
  if ( $merge == 1 && $one2two_overlap == 0 && $two2one_overlap == 0 ){
    return ( 1, $overlaps );
  }
  else{
    return ( 0, $overlaps);
  }
}


#########################################################################
# this function checks whether two transcripts merge
# according to consecutive exon overlap
# this time, matches like this:
#                        ____     ____        
#              exons1 --|____|---|____|------ etc... $j
#                        ____________  
#              exons2 --|____________|------ etc...  $k
#
# are considered a match


=head2 _test_for_Merge_allow_gaps
 Function: this function is called from link_Transcripts and actually checks whether two transcripts
           inputs merge.
 Returns : It returns two numbers ($merge,$overlaps), where
           $merge = 1 (0) when they do (do not) merge,
           and $overlaps is the number of exon-overlaps.

=cut

sub _test_for_Merge_allow_gaps{
  my ($self,$tran1,$tran2) = @_;
  
  my @exons1 = sort{ $a->start <=> $b->start } @{$tran1->get_all_Exons};
  my @exons2 = sort{ $a->start <=> $b->start } @{$tran2->get_all_Exons};
  
  my $foundlink = 0; # flag that gets set when starting to link exons
  my $start     = 0; # start looking at the first one
  my $overlaps  = 0; # independently if they merge or not, we compute the number of exon overlaps
  my $merge     = 0; # =1 if they merge

  #print STDERR "comparing ".$tran1->dbID." ($tran1)  and ".$tran2->dbID." ($tran2)\n";


  my $splice_mismatch = 0;
  if ( defined $self->splice_mismatch ){
    $splice_mismatch =  $self->splice_mismatch;
  }
  my $intron_mismatch = 0;
  if ( defined $self->intron_mismatch ){
    $intron_mismatch =  $self->intron_mismatch;
  }

EXON1:
  for (my $j=0; $j<=$#exons1; $j++) {
    
  EXON2:
    for (my $k=$start; $k<=$#exons2; $k++){
      print STDERR "comparing ".($j+1)." and ".($k+1)."\n";
	    
      # if exon 1 is not the first, check first whether it matches the previous exon2 as well, i.e.
      #                        ____     ____        
      #              exons1 --|____|---|____|------ etc... $j
      #                        ____________  
      #              exons2 --|____________|------ etc...  $k
      #
      # we allow to go over an intron if it is <= $intron_mismatch
      # and can have a mismatch at the other end of <= $splice_mismatch
      if ($foundlink == 1 && $j != 0){
	if ( $k != 0 
	     && $exons1[$j]->overlaps($exons2[$k-1])
	     && ( $exons1[$j]->start - $exons1[$j-1] -1 ) <= $intron_mismatch ){
	  print STDERR ($j+1)." <--> ".($k)."\n";
	  $overlaps++;
          next EXON1;
	}
      }
      
      # if exons1[$j] and exons2[$k] overlap go to the next exon1 and next $exon2
      if ( $exons1[$j]->overlaps($exons2[$k]) ){
	print STDERR ($j+1)." <--> ".($k+1)."\n";
        $overlaps++;
	
        # in order to merge the link always start at the first exon of one of the transcripts
        # we allow some mismatch at the start, and a mismatch of 
	# $splice_mismatch at the end
	if ( $j == 0 || $k == 0 
	     && abs($exons1[$j]->end - $exons2[$k]->end)<= $splice_mismatch){
	  
	  # but if it is also the last exon
	  if ( ( ( $k == 0 && $k == $#exons2 )   || 
		 ( $j == 0 && $j == $#exons1 ) ) ){
	    
	    # we force it to match the start as well (with a mismatch of $allowed_mismatch bases allowed)
	    # this might happen when a single exon est matches one exon in the other est
	    if ( abs($exons1[$j]->start - $exons2[$k]->start)<= $splice_mismatch ){
	      $foundlink = 1;
	      $merge     = 1;
	      $overlaps++;
	      last EXON1;
	    }
	    # else, it is non merge
	    else{
	      $foundlink = 0;
	      $merge     = 0;
	      #print STDERR "non-merged single exon transcript\n";
	      last EXON1;
	    }
	  }
	  else { 
	    #else, we have a link
	    $foundlink = 1;
	    $overlaps++;
	    #print STDERR "found a link\n";
	  }
	}
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
	while ( $k+1+$addition < scalar(@exons2) 
		&& $exons1[$j]->overlaps($exons2[$k+1+$addition])
		&& ( $exons2[$k+1+$addition]->start - $exons2[$k+$addition]->end - 1 ) <= $intron_mismatch
		&& abs( $exons1[$j]->end - $exons2[$k]->end ) <= $splice_mismatch
	      ){
	  print STDERR ($j+1)." <--> ".($k+2+$addition)."\n";
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

  return ($merge,$overlaps);
}



#########################################################################
# this function checks whether two transcripts merge
# according to consecutive exon overlap
# this time, matches like this:
#                        ____     ____        
#              exons1 --|____|---|____|------ etc... $j
#                        ____________  
#              exons2 --|____________|------ etc...  $k
#
# are checked, it won't be considered a merge, but it will count how many of those occur

sub _test_for_Merge_with_gaps{
  my ($self,$tran1,$tran2) = @_;
  my @exons1 = @{$tran1->get_all_Exons};
  my @exons2 = @{$tran2->get_all_Exons};	
 
  my $foundlink = 0; # flag that gets set when starting to link exons
  my $start     = 0; # start looking at the first one
  my $overlaps  = 0; # independently if they merge or not, we compute the number of exon overlaps
  my $merge     = 0; # =1 if they merge

  my $one2one_overlap = 0;
  my $one2two_overlap = 0;
  my $two2one_overlap = 0;
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
	  $two2one_overlap++;
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
	  $one2two_overlap++;
	  $overlaps++;
          $addition++;
	}      
	$start = $k+1+$addition;
	next EXON1;
      }    
      
    } # end of EXON2 
    
    # if you haven't found any match for this exon1, start again from the first exon2:
    if ($foundlink == 0){
      $start = 0;
    }
 
  }   # end of EXON1      

  # we only make them merge if $merge = 1 and the 2-to-1 and 1-to-2 overlaps are zero;
  if ( $merge == 1 ){
    return ( 1, $overlaps );
  }
  else{
    return ( 0, $overlaps);
  }
}
  

#########################################################################
   
# this compares both transcripts and calculate the number of overlapping exons and
# the length of the overlap

sub _compare_Transcripts {         
  my ($tran1, $tran2) = @_;
  my @exons1   = @{$tran1->get_all_Exons};
  my @exons2   = @{$tran2->get_all_Exons};
  my $overlaps = 0;
  my $overlap_length = 0;
  foreach my $exon1 (@exons1){
    foreach my $exon2 (@exons2){
      if ( ($exon1->overlaps($exon2)) && ($exon1->strand == $exon2->strand) ){
	$overlaps++;
	
	# calculate the extent of the overlap
	if ( $exon1->start > $exon2->start && $exon1->start <= $exon2->end ){
	  if ( $exon1->end < $exon2->end ){
	    $overlap_length += ( $exon1->end - $exon1->start + 1);
	  }
	  elsif ( $exon1->end >= $exon2->end ){
	    $overlap_length += ( $exon2->end - $exon1->start + 1);
	  }
	}
	elsif( $exon1->start <= $exon2->start && $exon2->start <= $exon1->end ){
	  if ( $exon1->end < $exon2->end ){
	    $overlap_length += ( $exon1->end - $exon2->start + 1);
	  }
	  elsif ( $exon1->end >= $exon2->end ){
	    $overlap_length += ( $exon2->end - $exon2->start + 1);
	  }
	}
      }
    }
  }
  
  return ($overlaps,$overlap_length);
}    

#########################################################################

	
########################################################################

1;
