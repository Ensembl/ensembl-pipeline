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


=head1 AUTHOR - Eduardo Eyras

This module is part of the Ensembl project http://www.ensembl.org

=head1 CONTACT

ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptComparator;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Pipeline::GeneComparison::GeneCluster;
use Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster;
use Bio::EnsEMBL::Pipeline::GeneComparison::ObjectMap;

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
    #print STDERR "############### comparison level $comparison_level ###################\n";
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
    ($merge, $overlaps) = $self->_test_for_Merge_allow_small_introns( $tran1, $tran2 );
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

  my $verbose = 1;
  
  print STDERR "=========== comparing ================\n";
  Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript( $est_tran ) if $verbose;
  Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript( $ens_tran ) if $verbose;
  
  my @exons1 = @{$est_tran->get_all_Exons};
  my @exons2 = @{$ens_tran->get_all_Exons};	
  
  # the simplest check is to see whether they are in the same genomic region:
  if ( $exons1[0]->start > $exons2[$#exons2]->end 
       ||
       $exons2[0]->start > $exons1[$#exons1]->end 
     ){
    print STDERR "transcript genomic regions do not overlap\n" if $verbose;
    return (0,0);
  }

  # the simplest case is with two single-exon transcripts:
  if ( scalar(@exons1) == 1 && scalar(@exons2) == 1 ){
    if ( $exons1[0]->overlaps($exons2[0] )){
      return (1,1);
    }
  }

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
      print STDERR "comparing j = $j : ".$exons1[$j]->start."-".$exons1[$j]->end.
        " and k = $k : ".$exons2[$k]->start."-".$exons2[$k]->end."\n" if $verbose;
      
      # we allow some mismatches at the extremities
      #                        ____     ____     ___   
      #              exons1   |____|---|____|---|___|  $j
      #                         ___     ____     ____     ____
      #              exons2    |___|---|____|---|____|---|____|  $k
      
      # if there is no overlap, go to the next EXON2
      if ( $foundlink == 0 && !($exons1[$j]->overlaps($exons2[$k]) ) ){
	print STDERR "foundlink = 0 and no overlap --> go to next EXON2\n" if $verbose;
	next EXON2;
      }
      # if there is no overlap and we had found a link, there is no merge
      elsif ( $foundlink == 1 && !($exons1[$j]->overlaps($exons2[$k]) ) ){
	print STDERR "foundlink = 1 and no overlap --> leaving\n" if $verbose;
	$merge = 0;
	last EXON1;
      }	
      # the first exon can have a mismatch ( any number of bases ) in the start
      # and $allowed_mismatch bases mismatch at the end
      elsif ( ($k == 0 || $j == 0) ){
	
	# if one of them is a single-exon transcript...
	if( ( $k == 0 && $k == $#exons2 ) || ( $j == 0 && $j == $#exons1 )    ){
	  if ( $k == 0 && $k == $#exons2 
	       && (
		   ( $j>0  && $exons2[$k]->overlaps($exons1[$j-1]) )
		   ||
		   ( $j<$#exons1 && $exons2[$k]->overlaps($exons1[$j+1]) )
		  )
	     ){
	    print STDERR "single exon transcript overlapping internally more than one exon. Not merging\n" if $verbose;
	    $merge = 0;
	    last EXON1;
	  }
	  elsif ( $j == 0 && $j == $#exons1 
		  && (
		      ( $k>0  && $exons1[$j]->overlaps($exons2[$k-1]) )
		      ||
		      ( $k<$#exons2 && $exons1[$j]->overlaps($exons2[$k+1]) )
		     )
		){
	    print STDERR "single exon transcript overlapping internally more than one exon. Not merging\n" if $verbose;
	    $merge = 0;
	    last EXON1;
	  }
	  elsif( ( $j==0 && $j== $#exons1 
		   && ( $exons2[$k]->start - $exons1[$j]->start <=  $allowed_mismatch )
		   && ( $exons1[$j]->end   - $exons2[$k]->end   <=  $allowed_mismatch )
		 )
		 ||
		 ( $k==0 && $k== $#exons2 
		   && ( $exons1[$j]->start - $exons2[$k]->start <=  $allowed_mismatch )
		   && ( $exons2[$k]->end   - $exons1[$j]->end   <=  $allowed_mismatch )
		 ) ){
	    $foundlink = 1;
	    $overlaps++;
	    $merge = 1;
	    print STDERR "merged single exon transcript\n" if $verbose;
	    last EXON1;
	  }
	  else{
	    print STDERR "single exon transcript overlapping beyond the allowed splic-site mismatches\n" if $verbose;
	    $merge = 0;
	    last EXON1;
	  }
	}
	
	# if the first overlaps with the last, we allow any overlap
	if ( ( $k==0 && $j == $#exons1 ) || ( $j==0 && $k == $#exons2 ) ){
	  $foundlink = 1;
	  $overlaps++;
	  $merge = 1;
	  print STDERR "found link --> merged\n" if $verbose;
	  last EXON1;
	}
	
	if ( ( $j != $#exons1 && $exons2[$k]->overlaps( $exons1[$j+1] ) )
	     ||
	     ( $k != $#exons2 && $exons1[$j]->overlaps( $exons2[$k+1] ) )
	   ){
	  print STDERR "Exon overlaps two exons, not merging\n";
	  $merge = 0;
	  last EXON1;
	  
	}
	
	if( abs($exons1[$j]->end - $exons2[$k]->end)<= $allowed_mismatch ){
	  
	  # else we force it to match the end (with a mismatch of $allowed_mismatch bases allowed)
	  
	  $foundlink = 1;
	  $overlaps++;
	  $start = $k+1;
	  print STDERR "found link\n" if $verbose;
	  next EXON1;
	}
      }
      # the last one can have any mismatch on the end
      # but must have a match at the start (with $allowed_mismatch mismatches allowed)
      if ( ( $k == $#exons2 || $j == $#exons1 ) ){
	if ( ( $k != $#exons2 && $exons1[$j]->overlaps( $exons2[$k+1] ) )
	     ||
	     ( $j != $#exons1 && $exons2[$k]->overlaps( $exons1[$j+1] ) )
	   ){
	  print STDERR "Exon overlaps two exons, not merging\n";
	  $merge = 0;
	  last EXON1;
	}
	elsif( $foundlink == 1 
	       &&
	       ( abs($exons1[$j]->start - $exons2[$k]->start)<= $allowed_mismatch ) 
	     ){
	  print STDERR "link completed, merged transcripts\n" if $verbose;
	  $merge = 1;
	  $overlaps++;
	  last EXON1;
	}
      }
      # the middle one must have exact matches
      # (up to an $allowed_mismatch mismatch)
      if ( ($k != 0 && $k != $#exons2) && 
	   ($j != 0 && $j != $#exons1) &&
	   ( $foundlink == 1)          &&
	   abs( $exons1[$j]->start - $exons2[$k]->start )<= $allowed_mismatch &&
	   abs( $exons1[$j]->end   - $exons2[$k]->end   )<= $allowed_mismatch
	 ){
	$overlaps++;
	$start = $k+1;
	print STDERR "continue link\n" if $verbose;
	next EXON1;
      }
      else{
	$merge = 0;
	last EXON1;
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
  if ($merge ){
    print STDERR $est_tran->dbID."V".$ens_tran->dbID." MERGE ". $ens_tran->dbID."V".$est_tran->dbID." MERGE\n" if $verbose;
  }
  else{
    print STDERR $est_tran->dbID."V".$ens_tran->dbID." NO MERGE ". $ens_tran->dbID."V".$est_tran->dbID." NO MERGE\n" if $verbose;
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
  my @exons1 = sort { $a->start <=> $b->start } @{$tran1->get_all_Exons};
  my @exons2 = sort { $a->start <=> $b->start } @{$tran2->get_all_Exons};	
 
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
      #print STDERR "merge, overlaps = $overlaps\n";
      return ( 1, $overlaps );
  }
  else{
      #print STDERR "no merge, overlaps = $overlaps\n";
    return ( 0, $overlaps);
  }
}


############################################################

=head2 _test_for_Merge_allow_small_introns
 Function: this function is called at the level 4 comparison
           it will first bridge over small introns (found in _disfuse_small_introns() )
           and then call the comparison at level 3, which is _test_for_fuzzy_semiexact_Merge().
 Returns: Like the other comparison methods it returns the values $merge = BOOLEAN (whether they merge or not)
            and $overlaps = INT (Number of exon overlaps found.

=cut

sub _test_for_Merge_allow_small_introns{
  my ($self,$tran1,$tran2) = @_;
  
  $tran1 = $self->_difuse_small_introns( $tran1 );
  $tran2 = $self->_difuse_small_introns( $tran2 );
  #return $self->_test_for_fuzzy_semiexact_Merge( $new_tran1, $new_tran2 );
  return $self->_test_for_merge(  $tran1, $tran2 );
}


############################################################

=head2 _difuse_small_introns
 Function: this function is called at the level 4 comparison
           In order to simplfy things, we difuse small introns, according to the
           value of 'intron_mismatch', since LEVEL 4 merges introns of this size.
           inputs merge.
 WARNING: It does not preserve translations. This was made with ests and cdnas in mind.
 Returns : a Bio::EnsEMBL::Transcript object. A new one if the introns have been difused or
           the same one we pass in if 'intron_mismatch' is not defined.

=cut

sub _difuse_small_introns{
  my ($self,$tran) = @_;
  my $modified = 0;
  if ( $self->intron_mismatch ){
    
    my $newtran = Bio::EnsEMBL::Transcript->new();
    if ( $tran->dbID ){
      $newtran->dbID($tran->dbID);
    }
    my @exons = sort{ $a->start <=> $b->start } @{$tran->get_all_Exons};
    my $exon_count = 0;
    my $current_exon;
    for (my $i=0; $i<=$#exons; $i++){
      my $exon = Bio::EnsEMBL::Pipeline::Tools::ExonUtils->_clone_Exon( $exons[$i] );
      if ( $i>0 ){
	if ( $exon->start - $current_exon->end - 1 <= $self->intron_mismatch ){
	  $current_exon->end( $exon->end );
	  $modified++;
	}
	else{
	  $current_exon = $exon;
	  $newtran->add_Exon( $current_exon );
	}
      }
      else{
	$current_exon = $exon;
	$newtran->add_Exon( $current_exon );
      }
    }
    if ($modified){
      print STDERR "difused transcript:\n";
      print STDERR "before:\n";
      Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript($tran);
      Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript($newtran);
      return $newtran;
    }
    else{
      return $tran;
    }
  }
  else{
    $self->("LEVEL 4 invoked but no intron_mismatch value defined. Doing level 3 instead");
    return $tran;
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

 ############################################################
 NOTE: STILL NOT WORKING PROPERLY. 
 ############################################################

=cut

sub _test_for_Merge_allow_gaps{
  my ($self,$tran1,$tran2) = @_;
  
  my @exons1 = sort{ $a->start <=> $b->start } @{$tran1->get_all_Exons};
  my @exons2 = sort{ $a->start <=> $b->start } @{$tran2->get_all_Exons};
  
  my $foundlink = 0; # flag that gets set when starting to link exons
  my $start     = 0; # start looking at the first one
  my $overlaps  = 0; # independently if they merge or not, we compute the number of exon overlaps
  my $merge     = 0; # =1 if they merge

  print STDERR "comparing:\n";
 Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript($tran1);
 Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript($tran2);

  my $splice_mismatch = 0;
  if ( defined $self->splice_mismatch ){
    $splice_mismatch =  $self->splice_mismatch;
  }
  my $intron_mismatch = 0;
  if ( defined $self->intron_mismatch ){
    $intron_mismatch =  $self->intron_mismatch;
  }

  # check first single exon genes:
  if ( scalar(@exons1) == 1 && scalar(@exons2) == 1 ){
    if ( $exons1[0]->overlaps( $exons2[0] ) ){
      return (1,1);
    }
  }

 EXON1:
  for (my $j=0; $j<=$#exons1; $j++) {
    
  EXON2:
    for (my $k=$start; $k<=$#exons2; $k++){
      print STDERR "comparing ".($j+1).": ".
	$exons1[$j]->start."-".$exons1[$j]->end.
	  " and ".($k+1).": ".
	    $exons2[$k]->start."-".$exons2[$k]->end."\n";
      
      # if exon 1 is not the first, check first whether it matches the previous exon2 as well, i.e.
      #                        ____     ____        
      #              exons1 --|____|---|____|------ etc... $j
      #                        ____________  
      #              exons2 --|____________|------ etc...  $k
      #
      # we allow to go over an intron if it is <= $intron_mismatch
      # and can have a mismatch at the other end of <= $splice_mismatch
      if ($foundlink == 1 && $j != 0){
	if ( $k != 0 && $exons1[$j]->overlaps($exons2[$k-1]) ){
	  print STDERR "checking 2-to-1 overlap:\n";
	  if ( $exons1[$j]->start - $exons1[$j-1]->end - 1 <= $intron_mismatch 
	       && ( $exons2[$k-1]->end - $exons1[$j]->end >=0
		    ||
		    $exons1[$j]->end - $exons2[$k-1]->end <= $splice_mismatch  
		  )
	     ){
	    print STDERR ($j+1)." <--> ".($k)."\n";
	    $overlaps++;
	    next EXON1;
	  }
	  else{
	    print STDERR "bad overlap, not merging\n";
	    merge = 0;
	    last EXON1;
	  }
	}
      }
      
      if ( $foundlink == 0 && !($exons1[$j]->overlaps($exons2[$k]) ) ){
	print STDERR "foundlink = 0 and no overlap --> go to next EXON2\n";
	next EXON2;
      }
      # if there is no overlap and we had found a link, there is no merge
      elsif ( $foundlink == 1 && !($exons1[$j]->overlaps($exons2[$k]) ) ){
	print STDERR "foundlink = 1 and no overlap --> leaving\n";
	$merge = 0;
	last EXON1;
      }	
      
      if ( $exons1[$j]->overlaps( $exons2[$k] )){
	
	if (  ($k == 0 && $k == $#exons2 ) || ( $j ==0 && $j == $#exons1 ) ){
	  if ( $self->compare_exon( $exons1[$j], $exons2[$k] ) ){
	    $merge = 1;
	    $foundlink =1;
	    print STDERR ($j+1)." <--> ".($k+1)."\n";
	    print STDERR "merged single exon transcript\n";
	    $overlaps++;
	    last EXON1;
	  }
	  elsif( ( $j ==0 && $j == $#exons1 ) 
		 && $exons1[$j]->overlaps( $exons2[$k+1] ) 
		 && $self->compare_left_exon( $exons1[$j], $exons2[$k] )
	       ){
	    $foundlink =1;
	    print STDERR ($j+1)." <--> ".($k+1)."\n";
	  }
	}	
	elsif ( ($k == 0 || $j == 0) ){    
	  if ( $self->compare_right_exon( $exons1[$j], $exons2[$k] ) ){
	    $foundlink = 1;
	    print STDERR ($j+1)." <--> ".($k+1)."\n";
	    $overlaps++;
	  }
	  elsif( $k < $#exons2 && $exons1[$j]->overlaps( $exons2[$k+1] ) ){
	    $foundlink = 1;
	    print STDERR ($j+1)." <--> ".($k+1)."\n";
	    $overlaps++;
	  }
	}
	elsif( ( $k == $#exons2 || $j == $#exons1 ) ){
	  if ( $self->compare_left_exon( $exons1[$j], $exons2[$k] ) ){
	    $foundlink = 1;
	    print STDERR  ($j+1)." <--> ".($k+1)."\n";
	    $overlaps++;
	  }
	  elsif( $k < $#exons2 && $exons1[$j]->overlaps( $exons2[$k+1] ) ){
	    $foundlink = 1;
	    print STDERR  ($j+1)." <--> ".($k+1)."\n";
	    $overlaps++;
	  }
	}
	# check whether in exons2 there are further exons overlapping exon1, i.e.
	#                       ____________        
	#             exons1 --|____________|------ etc...
	#                       ____     ___  
	#             exons2 --|____|---|___|------ etc...
	# 
	# we keep on linking until it does not overlap, or exon2 end falls beyond $splice_mismatch bases
	# from the exon1
	if ( $foundlink == 1 && $exons1[$j]->overlaps($exons2[$k+1]) ){
	  my $addition = 0;
	FORWARD:
	  while ( $k+1+$addition < scalar(@exons2) ){
	    if ( $exons1[$j]->overlaps($exons2[$k+1+$addition]) ){
	      print STDERR "checking 1-to-2 overlap:\n";
	      if ( $exons2[$k+1+$addition]->start - $exons2[$k+$addition]->end - 1  <= $intron_mismatch
		   && ( $exons1[$j]->end - $exons2[$k+1+$addition]->end >= 0 
			|| 
			$exons2[$k+1+$addition]->end - $exons1[$j] <= $splice_mismatch 
		      )
		 ){
		print STDERR ($j+1)." <--> ".($k+2+$addition)."\n";
		$foundlink = 1;
		$overlaps++;
		$addition++;
		
	      }
	      else{
		print STDERR "Bad overlap found, not merging\n";
		$merge = 0;
		last EXON1;
	      }
	    }
	    else{
	      last FORWARD;
	    }
	  }
	  $start = $k+1+$addition;
	  if ( $foundlink && ( $j == $#exons1 || $start > $#exons2 ) ){
	    print STDERR "link completed, merged transcripts\n";
	    last EXON1;
	  }    
	  next EXON1;
	}
	elsif( $foundlink == 1 ){
	  $start = $k+1;
	  next EXON2;
	}
      }
      # if they don't overlap, try the next one
      else{
	next EXON2;
      }
    }  # end of EXON2      

    if ($foundlink == 0){
      $start = 0;
    }

  }    # end of EXON1
  return ($merge,$overlaps);
}

############################################################

sub compare_left_exon{
  my ($self, $exon1, $exon2 ) = @_;
  my $splice_mismatch = 0;
  if ( defined $self->splice_mismatch ){
    $splice_mismatch =  $self->splice_mismatch;
  }
  #my $intron_mismatch = 0;
  #if ( defined $self->intron_mismatch ){
  #  $intron_mismatch =  $self->intron_mismatch;
  #}
  
  unless ( $exon1->overlaps( $exon2 ) ){
    return 0;
  }

  if ( abs( $exon1->start - $exon2->start ) <= $splice_mismatch ){
    return 1;
  }
  else{
    return 0;
  }
}

############################################################

sub compare_right_exon{
  my ($self, $exon1, $exon2 ) = @_;
  my $splice_mismatch = 0;
  if ( defined $self->splice_mismatch ){
    $splice_mismatch =  $self->splice_mismatch;
  }
  #my $intron_mismatch = 0;
  #if ( defined $self->intron_mismatch ){
  #  $intron_mismatch =  $self->intron_mismatch;
  #}

  unless ( $exon1->overlaps( $exon2 ) ){
    return 0;
  }

  if ( abs( $exon1->end - $exon2->end ) <= $splice_mismatch ){
    return 1;
  }
  else{
    return 0;
  }
}

############################################################

sub compare_exon{
  my ($self, $exon1, $exon2 ) = @_;
  my $splice_mismatch = 0;
  if ( defined $self->splice_mismatch ){
      $splice_mismatch =  $self->splice_mismatch;
  }
  #my $intron_mismatch = 0;
  #if ( defined $self->intron_mismatch ){
  #  $intron_mismatch =  $self->intron_mismatch;
  #}
  
  unless ( $exon1->overlaps( $exon2 ) ){
      return 0;
  }
  
  if ( abs( $exon1->start - $exon2->start ) <= $splice_mismatch 
       &&
       abs( $exon1->end  -  $exon2->end   ) <= $splice_mismatch
     ){
    return 1;
  }
  else{
    return 0;
  }
}
 
############################################################

sub _max{
    my ($self,$a,$b) = @_;
    if ( $a > $b ){
	return $a;
    }
    else{
	return $b;
    }
}

############################################################

sub _min{
    my ($self,$a,$b) = @_;
    if ( $a < $b ){
	return $a;
    }
    else{
	return $b;
    }
}

############################################################


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
  my ($self, $tran1, $tran2) = @_;
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


sub _test_for_merge{
  my ($self,$tran1,$tran2) = @_;

  my $verbose   = 1;

  if ($verbose){
      print STDERR "comparing ".
	  $tran1->dbID."-".$tran2->dbID." ( ".$tran2->dbID."-".$tran1->dbID." )\n";
    Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript($tran1);
    Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript($tran2);
  }
  my @exons1 = sort { $a->start <=> $b->start } @{$tran1->get_all_Exons};
  my @exons2 = sort { $a->start <=> $b->start } @{$tran2->get_all_Exons};	

  if ( $exons1[0]->start > $exons2[$#exons2]->end 
       ||
       $exons2[0]->start > $exons1[$#exons1]->end 
       ){
      print STDERR "transcript genomic regions do not overlap\n" if $verbose;
      return (0,0);
  }
  
  # the simplest case is with two single-exon transcripts:
  if ( scalar(@exons1) == 1 && scalar(@exons2) == 1 ){
      if ( $exons1[0]->overlaps($exons2[0] )){
	  print STDERR "--- single-exon transcripts --- merge ---\n" if $verbose;
	  return (1,1);
      }
  }
  
  my $object_map = Bio::EnsEMBL::Pipeline::GeneComparison::ObjectMap->new();

  my $splice_mismatch = $self->splice_mismatch;
  my $foundlink = 0; # flag that gets set when starting to link exons
  my $start     = 0; # start looking at this one
  my $overlaps  = 0; # independently if they merge or not, we compute the number of exon overlaps
  my $merge     = 0; # =1 if they merge
  
  my %is_first;
  my %is_last;

  # we follow a greedy aproach to try to match all the exons
  # we jump out as soon as we find a problem with the matching
  
 EXON1:
  for (my $j=0; $j<=$#exons1; $j++) {
  
      # index the first and last position
      if ( $j==0 ){
	  $is_first{ $exons1[$j] } = 1;
      }
      else{
	  $is_first{ $exons1[$j] } = 0;
      }
      if ( $j == $#exons1 ){
	  $is_last{ $exons1[$j] } = 1;
      }
      else{
	  $is_last{ $exons1[$j] } = 0;
      }
    
  EXON2:
      for (my $k=$start; $k<=$#exons2; $k++){
	  
	  # index the first and last position
	  if ( $k==0 ){
	      $is_first{ $exons2[$k] } = 1;
	  }
	  else{
	      $is_first{ $exons2[$k] } = 0;
	  }
	  if ( $k == $#exons2 ){
	      $is_last{ $exons2[$k] } = 1;
	  }
	  else{
	      $is_last{ $exons2[$k] } = 0;
	  }
      
	  if ( $foundlink == 0 && !($exons1[$j]->overlaps($exons2[$k])) ){
	      print STDERR "go to next exon2\n";
	      $foundlink = 0;
	      next EXON2;
	  }
	  elsif ( $foundlink == 1 && !($exons1[$j]->overlaps($exons2[$k])) ){
	      print STDERR "link is broken, not merging\n";
	      $foundlink = 0;
	      $merge = 0;
	      last EXON1;
	  }
	  elsif ( $exons1[$j]->overlaps($exons2[$k]) ){
	      if ( $j == $#exons1 || $k == $#exons2 ){
		  #if ( $j != $#exons1 && $exons2[$k]->overlaps( $exons1[$j+1] )
		  #     ||
		  #     $k != $#exons2 && $exons1[$j]->overlaps( $exons2[$k+1] )
		  #     ){
		  #    print STDERR "terminal exon overlaps two exons. Not merging\n" if $verbose;
		  #    $merge = 0;
		  #    last EXON1;
		  #}
		  #else{
		  print STDERR ($j+1)." <--> ".($k+1)."\n" if $verbose;
		  $object_map->match( $exons1[$j], $exons2[$k] );
		  print STDERR "end of transcripts - there is potential merge|\n";
		  $merge = 1;
		  $overlaps++;
		  $foundlink = 1;
		  last EXON1;
		  #}
	      }
	      else{
		  print STDERR ($j+1)." <--> ".($k+1)."\n" if $verbose;
		  $object_map->match( $exons1[$j], $exons2[$k] );
		  $overlaps++;
		  $foundlink = 1;
		  $start++;
		  next EXON1;
	      }
	  }
	  
      } # end of EXON2 
      
      # if you haven't found any match for this exon1, start again from the first exon2:
      if ($foundlink == 0){
	  $start = 0;
      }
      
  }   # end of EXON1      
  
  unless ( $merge ){
      print STDERR "No merge\n" if $verbose;
      return ( 0, $overlaps );
  }
  
  my @list1 = $object_map->list1();
  my @list2 = $object_map->list2();
  
  print STDERR scalar(@list1)." elements in list 1\n";
  print STDERR scalar(@list2)." elements in list 2\n";
  
  ############################################################
  # the simplest case: when they match over one exon only:
  ############################################################
  if ( scalar(@list1)==1 && scalar(@list2)==1 ){
      
      ############################################################
      # if it is a single-exon transcript overlap: if it matches 
      # to the first or the last, we leave the open end unconstrained
      ############################################################
      if( (  $is_first{ $list1[0] } && $is_last{ $list1[0] } 
	     &&
	     $is_first{ $list2[0] }
	     &&
	     $list1[0]->end - $list2[0]->end <=$splice_mismatch 
	     )
	  ||
	  (  $is_first{ $list1[0] } && $is_last{ $list1[0] } 
	     &&
	     $is_last{ $list2[0] }
	     &&
	     $list2[0]->start - $list1[0]->start <=$splice_mismatch 
	     )
	  ||
	  (  $is_first{ $list2[0] } && $is_last{ $list2[0] } 
	     &&
	     $is_first{ $list1[0] }
	     &&
	     $list2[0]->end - $list1[0]->end <=$splice_mismatch 
	     )
	  ||
	  (  $is_first{ $list2[0] } && $is_last{ $list2[0] } 
	     &&
	     $is_last{ $list1[0] }
	     &&
	     $list1[0]->start - $list2[0]->start <=$splice_mismatch 
	     )
	  ){
	  print STDERR "here 1 --- merge ---\n" if $verbose;
	  return (1,1);
      }
      ############################################################
      # else, it is a single-exon overlap to an 'internal' exon
      ############################################################
      elsif( (  $is_first{ $list1[0] } && $is_last{ $list1[0] } 
		&&
		$list1[0]->end   - $list2[0]->end   <= $splice_mismatch 
		&&
		$list2[0]->start - $list1[0]->start <=$splice_mismatch 
		)
	     ||
	     (  $is_first{ $list2[0] } && $is_last{ $list2[0] } 
		&&
		$list2[0]->end   - $list1[0]->end   <= $splice_mismatch 
		&&
		$list1[0]->start - $list2[0]->start <=$splice_mismatch 
		)
	     ){
	  print STDERR "here 2 --- merge ---\n" if $verbose;
	  return (1,1);
      }
      ############################################################
      # if the first overlaps with the last:
      ############################################################
      elsif ( ( $is_first{ $list1[0] } && !$is_last{ $list1[0] } 
		&&
		$is_last{ $list2[0] }  && !$is_first{$list2[0] }
		)
	      ||
	      ( $is_first{ $list2[0] } && !$is_last{ $list2[0] }  
		&& 
		$is_last{ $list1[0] }  && !$is_first{$list1[0] }
		)
	      ){
	  print STDERR "here 3 --- merge ---\n" if $verbose;
	  return (1,1);
      }
      ############################################################
      # we have already dealt with single-exon against single-exon overlap
      ############################################################
      else{
	  print STDERR "No merge\n" if $verbose;
	  return (0,1);
      }
  }
  
  
  ############################################################
  # go over each pair stored in the object map
  ############################################################
 PAIR:
  foreach my $exon1 ( @list1 ){
      my @partners = $object_map->partners( $exon1 );
      if ( scalar( @partners ) > 1 ){
	  $self->warn("One exon has been matched to two exons");
      }
      my $exon2 = shift ( @partners );
  
      ############################################################
      # exon1 and exon2 are a pair, they overlap, we need to check that
      # they actually overlap as we want
      ############################################################
      
      ############################################################
      # both of them could be the first one
      ############################################################
      if ( $is_first{ $exon1} && $is_first{ $exon2 }
	   &&
	   abs( $exon2->end - $exon1->end ) <= $splice_mismatch
	   ){
	  next PAIR;
      }
      ############################################################
      # one of them could be the first one
      ############################################################
      elsif ( ( $is_first{ $exon1 } 
	     &&
	     $exon2->start - $exon1->start <= $splice_mismatch
	     &&
	     abs( $exon2->end - $exon1->end ) <= $splice_mismatch
	     )
	   ||
	   ( $is_first{ $exon2 }
	     &&
	     $exon1->start - $exon2->start <= $splice_mismatch
	     &&
	     abs( $exon1->end - $exon2->end ) <= $splice_mismatch
	     )
	   ){
	  next PAIR;
      }
      ############################################################
      # both could be the last one
      ############################################################
      elsif ( $is_last{ $exon1} && $is_last{ $exon2 }
	   &&
	   abs( $exon2->start - $exon1->start ) <= $splice_mismatch
	   ){
	  next PAIR;
      }
      ############################################################
      # one of them could be the last one
      ############################################################
      elsif ( ( $is_last{ $exon1 } 
	     &&
	     $exon1->end - $exon2->end <= $splice_mismatch
	     &&
	     abs( $exon2->start - $exon1->start ) <= $splice_mismatch
	     )
	   ||
	   ( $is_last{ $exon2 }
	     &&
	     $exon2->end - $exon1->end <= $splice_mismatch
	     &&
	     abs( $exon1->start - $exon2->start ) <= $splice_mismatch
	     )
	   ){
	  next PAIR;
      }
      ############################################################
      # we have already covered the case: first overlaps last
      ############################################################
      elsif( abs( $exon1->start - $exon2->start ) <= $splice_mismatch
	     &&
	     abs( $exon1->end - $exon2->end ) <= $splice_mismatch
	     ){
	  next PAIR;
      }
      else{
	  print STDERR "Failed to find a proper match. Not merging\n" if $verbose;
	  return (0,$overlaps);
      }
  } # end of PAIR
  
  print STDERR "--- merge ---\n" if $verbose;
  return (1,scalar(@list1));
}

########################################################################

1;
