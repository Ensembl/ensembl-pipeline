#
# Cared for by EnsEMBL 
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
    
    my $comparator = Bio::EnsEMBL::Pipeline::::GeneComparison::TranscriptComparator->new();
);

so far there are two methods that check whether two transcripts have consecutive exon overlap:

my ($merge,$overlaps,$exact) = $comparator->test_for_Merge($transcript1,$transcript2)

this second allows two exons to merge into one at the other transcript
my ($merge,$overlaps,$exact) = $comparator->test_for_Merge_with_gaps($transcript1,$transcript2)
    
    $merge =1 if the exons overlap consecutively
    $overlaps is the number of exon overlaps
    $exact = 1 if the exon matches are exact


=head1 CONTACT

ensembl-dev@ebi.ac.uk

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptComparator;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Pipeline::GeneComparison::GeneCluster;
use Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster;
use Bio::EnsEMBL::Pipeline::GeneCombinerConf;

# config file; parameters searched for here if not passed in as @args

@ISA = qw(Bio::EnsEMBL::Root);

######################################################################

sub new{
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  return $self;  
}

#########################################################################
# this function checks whether two transcripts merge
# with exact exon matches, except for
# possible mismatches in the extremal exons

sub _test_for_semiexact_Merge{
  my ($self,$est_tran,$ens_tran) = @_;
  
  my @exons1 = @{$est_tran->get_all_Exons};
  my @exons2 = @{$ens_tran->get_all_Exons};	
  
  @exons1 = sort {$a->start <=> $b->start} @exons1;
  @exons2 = sort {$a->start <=> $b->start} @exons2;

  my $foundlink = 0; # flag that gets set when starting to link exons
  my $start     = 0; # start looking at the first one
  my $merge     = 0; # =1 if they merge
  
 EXON1:
  for (my $j=0; $j<=$#exons1; $j++) {
      
    EXON2:
      for (my $k=$start; $k<=$#exons2; $k++){
	#print STDERR "comparing j = $j : ".$exons1[$j]->start."-".$exons1[$j]->end." and k = $k : ".$exons2[$k]->start."-".$exons2[$k]->end."\n";
	  
	  # we allow some mismatches at the extremities
	  #                        ____     ____     ___   
	  #              exons1   |____|---|____|---|___|  $j
	  #                         ___     ____     ____  
	  #              exons2    |___|---|____|---|____|  $k
	  
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
	  if ( ($k == 0 || $j == 0) && $exons1[$j]->end == $exons2[$k]->end ){
	      
	      # but if it is also the last exon
	      if ( ( ( $k == 0 && $k == $#exons2 )   || 
		     ( $j == 0 && $j == $#exons1 ) ) ){
		  
		  # we force it to match the start
		  if ( $exons1[$j]->start == $exons2[$k]->start ){
		      $foundlink  = 1;
		      $merge      = 1;
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
		  #print STDERR "found a link\n";
		  next EXON1;
	      }
	  }
	  # the last one can have a mismatch on the end
	  elsif ( ( $k == $#exons2 || $j == $#exons1 ) &&
		  ( $foundlink == 1 )                  &&
		  ( $exons1[$j]->start == $exons2[$k]->start ) 
		  ){
	      #print STDERR "link completed, merged transcripts\n";
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
	      $start = $k+1;
	      #print STDERR "continue link\n";
	      next EXON1;
	  }

      } # end of EXON2 
    
      if ($foundlink == 0){
	  $start = 0;
      }
      
  }   # end of EXON1      
  
  return $merge;
}

#########################################################################
# this function checks whether two transcripts merge
# with fuzzy exon matches: there is consecutive exon overlap 
# but there are mismatches of 2 base allowed at the edges of any exon pair
#
# Why 2 bases: 2 bases is perhaps not meaningful enough to be considered
# a biological difference, and it is possibly an artifact of any of the
# analysis previously run: genomewise, est2genome,... it is more likely to
# happen 

sub _test_for_fuzzy_semiexact_Merge{
  my ($self,$est_tran,$ens_tran) = @_;
  
  my @exons1 = @{$est_tran->get_all_Exons};
  my @exons2 = @{$ens_tran->get_all_Exons};	
  
  @exons1 = sort {$a->start <=> $b->start} @exons1;
  @exons2 = sort {$a->start <=> $b->start} @exons2;

  my $foundlink = 0; # flag that gets set when starting to link exons
  my $start     = 0; # start looking at the first one
  my $merge     = 0; # =1 if they merge
  
 EXON1:
  for (my $j=0; $j<=$#exons1; $j++) {
      
    EXON2:
      for (my $k=$start; $k<=$#exons2; $k++){
	#print STDERR "comparing j = $j : ".$exons1[$j]->start."-".$exons1[$j]->end.
	#  " and k = $k : ".$exons2[$k]->start."-".$exons2[$k]->end."\n";
	
	  # we allow some mismatches at the extremities
	  #                        ____     ____     ___   
	  #              exons1   |____|---|____|---|___|  $j
	  #                         ___     ____     ____  
	  #              exons2    |___|---|____|---|____|  $k
	  
	  # if there is no overlap, go to the next EXON2
	  if ( $foundlink == 0 && !($exons1[$j]->overlaps($exons2[$k]) ) ){
	    #print STDERR "foundlink = 0 and no overlap --> go to next EXON2\n";
	    next EXON2;
	  }
	  # if there is no overlap and we had found a link, there is no merge
	  if ( $foundlink == 1 && !($exons1[$j]->overlaps($exons2[$k]) ) ){
	      #print STDERR "foundlink = 1 and no overlap --> leaving\n";
	      $merge = 0;
	      last EXON1;
	  }	
	  
	  # the first exon can have a mismatch ( any number of bases) in the start
	  # and a 2base mismatch at the end
	  if ( ($k == 0 || $j == 0) && abs($exons1[$j]->end - $exons2[$k]->end)<3 ){
	      
	      # but if it is also the last exon
	      if ( ( ( $k == 0 && $k == $#exons2 )   || 
		     ( $j == 0 && $j == $#exons1 ) ) ){
		  
		  # we force it to match the start (with a mismatch of 2bases allowed)
		  if ( abs($exons1[$j]->start - $exons2[$k]->start)< 3 ){
		      $foundlink  = 1;
		      $merge      = 1;
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
		  #print STDERR "found a link\n";
		  next EXON1;
	      }
	  }
	  # the last one can have any mismatch on the end
	  # but must have a match at the start (wiht 2bases mismatch allowed)
	  elsif ( ( $k == $#exons2 || $j == $#exons1 ) &&
		  ( $foundlink == 1 )                  &&
		  ( abs($exons1[$j]->start - $exons2[$k]->start)<3 ) 
		  ){
	      #print STDERR "link completed, merged transcripts\n";
	      $merge = 1;
	      last EXON1;
	  }
	# the middle one must have exact matches
	# (up to a 2base mismatch)
	elsif ( ($k != 0 && $k != $#exons2) && 
		($j != 0 && $j != $#exons1) &&
		( $foundlink == 1)          &&
		abs( $exons1[$j]->start - $exons2[$k]->start )<3 &&
		abs( $exons1[$j]->end   - $exons2[$k]->end   )<3
		  ){
	      $start = $k+1;
	      #print STDERR "continue link\n";
	      next EXON1;
	  }

      } # end of EXON2 
    
      if ($foundlink == 0){
	  $start = 0;
      }
      
  }   # end of EXON1      
  
  return $merge;
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

sub _test_for_Merge{
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
