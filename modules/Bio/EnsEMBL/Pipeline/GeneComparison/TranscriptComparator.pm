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

restrict_external_splice_site ---> TRUE if we want to restrict how much can exceed
                                   an external exon an internal splice site:
    
                                                    |--d--|
                                    ######-------##########
                                   #######-------####-------------#######
    
                                   If TRUE 'd' must be <= splice_mismatch, else 'd' can be anything
                                   The difference 'd' could be due to a true splice site, but with
                                   ESTs it is most likely to be an artifact
    

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
use Bio::EnsEMBL::Pipeline::Tools::ExonUtils;

@ISA = qw(Bio::EnsEMBL::Root);

######################################################################

sub new{
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  my ( $comparison_level, $exon_match, $splice_mismatch, $intron_mismatch, $restrict_external_splice_site ) = 
    $self->_rearrange([qw(COMPARISON_LEVEL
			  EXON_MATCH
			  SPLICE_MISMATCH
			  INTRON_MISMATCH
			  RESTRICT_EXTERNAL_SPLICE_SITE
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
  
  if( defined $restrict_external_splice_site ){
      $self->restrict_external_splice($restrict_external_splice_site);
  }

  $self->verbose(0);

  return $self;  
}

############################################################
#
# PARAMETRIZATION
#
############################################################

sub  verbose{
  my ($self,$boolean) = @_;
  if (defined $boolean){
    $self->{_verbose} = $boolean;
  }
  return $self->{_verbose};
}

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

sub restrict_external_splice_site{
    my ($self, $boolean) = @_;
    if ( defined $boolean ){ 
	$self->{_restrict_external_splice_site} = $boolean;
    }
    return $self->{_restrict_external_splice_site};
}

sub intron_mismatch{
  my ($self, $int) = @_;
  if ( defined $int ){
     $self->{_intron_mismatch} = $int;
  }
  return $self->{_intron_mismatch};
}

############################################################
#
# PUBLIC COMPARE METHOD
#
############################################################

=head2 compare
  Arg[1] and Arg[2] : 2 transcript objects to compare
  Arg[3]: the mode ( a string of text). Possible modes:

=cut 

sub compare{
  my ($self, $tran1, $tran2) = @_;
  
  my ($merge, $overlaps);

  # switch on comparison level
  if ( $self->comparison_level == 3 ){
    ($merge, $overlaps) = $self->_test_for_merge( $tran1, $tran2 );
  }
  elsif( $self->comparison_level == 2 ){
      $self->splice_mismatch(0);
      ($merge, $overlaps) = $self->_test_for_merge(  $tran1, $tran2 );
    }
  elsif( $self->comparison_level == 4 ){
      # this calls _test_for_merge but after eliminating small introns
      ($merge, $overlaps) = $self->_test_for_Merge_allow_small_introns( $tran1, $tran2 );
      
  }
  elsif( $self->comparison_level == 1 ){
    ($merge, $overlaps) = $self->_test_for_strict_merge( $tran1, $tran2 );
  }
  elsif( $self->comparison_level == 5 ){
    # this will not check any differences in the splice sites
    if ( defined $self->intron_mismatch ){
      ($merge, $overlaps) = $self->_test_for_Merge_allow_small_introns( $tran1, $tran2 );
    }
    else{
      ($merge, $overlaps) = $self->_test_for_merge(     $tran1, $tran2 );
    }
  }
  
  return ($merge,$overlaps);
}
  
############################################################
#
# METHODS CALLED BY 'COMPARE' - DOING THE ACTUAL COMPARISON, ETC...
#
############################################################



sub _test_for_strict_merge{
  my ( $self, $tran1, $tran2 ) = @_;
  my @exons1 = sort { $a->start <=> $b->start } @{$tran1->get_all_Exons};
  my @exons2 = sort { $a->start <=> $b->start } @{$tran2->get_all_Exons};	
  
  unless ( scalar(@exons1) == scalar(@exons2) ){
    return 0;
  }

  for ( my $i=0; $i<=$#exons1; $i++ ){
    unless ( $exons1[$i]->start == $exons2[$i]->start && $exons1[$i]->end == $exons2[$i]->end ){
      return 0;
    }
  }
  return 1;
}
      


=head2 _test_for_Merge_allow_small_introns
 Function: this function is called at the level 4 comparison
           it will first bridge over small introns 
           and then call the function _test_for_merge
 Returns: Like the other comparison methods it returns the values $merge = BOOLEAN (whether they merge or not)
            and $overlaps = INT (Number of exon overlaps found)

=cut

sub _test_for_Merge_allow_small_introns{
  my ($self,$tran1,$tran2) = @_;
  $tran1 = $self->_difuse_small_introns( $tran1 );
  $tran2 = $self->_difuse_small_introns( $tran2 );
  
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
      #print STDERR "difused transcript:\n";
      #print STDERR "before:\n";
      #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript($tran);
      #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript($newtran);
      return $newtran;
    }
    else{
      return $tran;
    }
  }
  else{
    #$self->warn("LEVEL 4 invoked but no intron_mismatch value defined. Doing level 3 instead");
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

  my $verbose   = $self->verbose;

  if ($verbose){
      my $id1 = $tran1->dbID || "no-id";
      my $id2 = $tran2->dbID || "no-id";
      print STDERR "comparing ".$id1."-".$id2." ( ".$id2."-".$id1." )\n";
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
    else{
      print STDERR "--- single-exon transcripts --- No merge ---\n" if $verbose;
      return (0,0);
    }
  }

  ############################################################
  # do first a loose comparison in a greedy manner
  ############################################################
  my ($object_map,$merge,$overlaps,$is_first,$is_last) = $self->_fast_compare($tran1,$tran2);
  
  unless ( $merge ){
      print STDERR "No merge\n" if $verbose;
      return ( 0, $overlaps );
  }

  ############################################################
  # if we are in loose mode, we can return already as we have checked
  # for consecutive exon overlap, regardless of the splice mismatches
  ############################################################
  if ( $self->comparison_level == 5 && $merge ){
    print STDERR "Loose mode merge\n" if $verbose;
    return (1,$overlaps);
  }
  
  ############################################################
  # else, check splice site mismatches
  ############################################################
  return ( $self->_process_comparison( $object_map, $tran1, $tran2, $merge, $is_first, $is_last ), $overlaps);
  
}

########################################################################

=head2 discrete_compare
    
    This method compares two transcripts (tran1,tran2) and returns a discrete value which is
    the relation of tran1 with respect to tran2:
    'extension' if tran1 extends tran2
    'inclusion' if tran1 is included in tran2
    'clash'     if tran1 and tran2 have incompatible splicing structure
    'no-overap' if tran1 and tran2 do not overlap
    This method is only used by the algorithm ClusterMerge.
    It assumes that the transcripts are ordered as specified in that module.

=cut

sub discrete_compare{
  my ($self, $tran1, $tran2 ) = @_;
  
  my $verbose   =  $self->verbose;
  
  if ( ( 4 == $self->comparison_level || 5 == $self->comparison_level )
       && $self->intron_mismatch > 0 ){
    $tran1 = $self->_difuse_small_introns( $tran1 );
    $tran2 = $self->_difuse_small_introns( $tran2 );
  }
  
  print STDERR "comparing ".$tran1->dbID." and ".$tran2->dbID."\n";
  if ($verbose){
    #print STDERR "comparing ".
    #$tran1->dbID."-".$tran2->dbID." ( ".$tran2->dbID."-".$tran1->dbID." )\n";
    print STDERR "comparing\n";
    print STDERR $tran1->dbID.": ";
    Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript($tran1);
    print STDERR $tran2->dbID.": ";
    Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript($tran2);
  }
  my @exons1 = sort { $a->start <=> $b->start } @{$tran1->get_all_Exons};
  my @exons2 = sort { $a->start <=> $b->start } @{$tran2->get_all_Exons};	
  
  #######################
  # check for no-overlap
  #######################
  if ( $exons1[0]->start > $exons2[$#exons2]->end 
       ||
       $exons2[0]->start > $exons1[$#exons1]->end 
     ){
    print STDERR "transcript genomic regions do not overlap\n" if $verbose;
    return 'no-overlap';
  }
  
  ############################################################
  # the simplest case is with two single-exon transcripts:
  ############################################################
  if ( scalar(@exons1) == 1 && scalar(@exons2) == 1 ){
    if ( $exons1[0]->overlaps($exons2[0] )){
      
      if ( $exons1[0]->end > $exons2[0]->end ){
	print STDERR "EXTENSION\n" if $verbose;
	return 'extension';
      }
      elsif( $exons1[0]->end <= $exons2[0]->end ){
	print STDERR "INCLUSION\n" if $verbose;
	return 'inclusion';
      }
    }
    else{
      print STDERR "NO_OVERLAP\n" if $verbose;
      return 'no-overlap';
    }
  }
    
  ############################################################
  # do first a loose comparison in a greedy manner
  ############################################################
  my ($object_map,$merge,$overlaps,$is_first,$is_last) = $self->_fast_compare($tran1,$tran2);

  my %is_first = %$is_first;
  my %is_last  = %$is_last;
  
  unless ( $merge ){
    print STDERR "here1: CLASH\n" if $verbose;
    return 'clash';
  }
  
  ############################################################
  # if we are in loose mode, we can return already as we have checked
  # for consecutive exon overlap, regardless of the splice mismatches
  ############################################################
  if ( $self->comparison_level == 5 && $merge ){
    print STDERR "Loose mode merge\n" if $verbose;
    if ( $exons1[$#exons1]->end > $exons2[$#exons2]->end ){
	  print STDERR "EXTENSION\n" if $verbose;  
	  return 'extension';
	}
    elsif( $exons1[$#exons1]->end <= $exons2[$#exons2]->end ){
      print STDERR "INCLUSION\n" if $verbose;
      return 'inclusion';
    }
  }
  
  ############################################################
  # else, check splice site mismatches
  ############################################################
  if ( $self->_process_comparison($object_map, $tran1, $tran2, $merge, $is_first, $is_last) ){
    if ( $exons1[$#exons1]->end > $exons2[$#exons2]->end ){
      print STDERR "EXTENSION\n" if $verbose;
      return 'extension';
    }
    elsif( $exons1[$#exons1]->end <= $exons2[$#exons2]->end ){
      print STDERR "INCLUSION\n" if $verbose;
      return 'inclusion';
    }
    }
  else{
    print STDERR "here2: CLASH\n" if $verbose;
    return 'clash';
  }
}

############################################################

sub _fast_compare{
  my ($self, $tran1, $tran2) = @_;
  my @exons1 = sort { $a->start <=> $b->start } @{$tran1->get_all_Exons};
  my @exons2 = sort { $a->start <=> $b->start } @{$tran2->get_all_Exons};	
  my %is_first;      # flag the first exon
  my %is_last;       # flag the last exon
  my $foundlink = 0; # flag that gets set when starting to link exons
  my $start     = 0; # start looking at this one
  my $overlaps  = 0; # independently if they merge or not, we compute the number of exon overlaps
  my $merge     = 0; # =1 if they merge

  # an ObjectMap holds an application between any two sets of objects
  my $object_map = Bio::EnsEMBL::Pipeline::GeneComparison::ObjectMap->new();

  my $verbose = $self->verbose;

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
      
      print STDERR "foundlink = $foundlink\n" if $verbose;
      print STDERR "comparing ".$exons1[$j]->start."-".$exons1[$j]->end." and ".$exons2[$k]->start."-".$exons2[$k]->end."\n" if $verbose;
      
      if ( $foundlink == 0 && !($exons1[$j]->overlaps($exons2[$k])) ){
	print STDERR "go to next exon2\n" if $verbose;
	$foundlink = 0;
	next EXON2;
      }
      elsif ( $foundlink == 1 && !($exons1[$j]->overlaps($exons2[$k])) ){
	print STDERR "link is broken, not merging\n" if $verbose;
	$foundlink = 0;
	$merge = 0;
	last EXON1;
      }
      elsif ( $exons1[$j]->overlaps($exons2[$k]) ){
	  
	  ############################################################
	  # check there is no 2-to-1 overlap
	  if ( $j != $#exons1 && $exons2[$k]->overlaps( $exons1[$j+1] )
	       ||
	       $k != $#exons2 && $exons1[$j]->overlaps( $exons2[$k+1] )
	       ){
	      print STDERR "exon overlaps two exons. Not merging\n" if $verbose;
	      $merge = 0;
	      last EXON1;
	  }	
	  elsif ( $foundlink == 0 && $j != 0 && $k != 0 ){
	      print STDERR "there is overlap - but it is not consecutive - no merge\n" if $verbose;
	      $merge = 0;
	      last EXON1;
	  }
	  elsif ( $j == $#exons1 || $k == $#exons2 ){
	      print STDERR ($j+1)." <--> ".($k+1)."\n" if $verbose;
	      $object_map->match( $exons1[$j], $exons2[$k] );
	      print STDERR "end of transcripts - there is potential merge!\n" if $verbose;
	      $merge = 1;
	      $overlaps++;
	      $foundlink = 1;
	      last EXON1;
	  }
	  else{
	      print STDERR ($j+1)." <--> ".($k+1)."\n" if $verbose;
	      $object_map->match( $exons1[$j], $exons2[$k] );
	      $overlaps++;
	      $foundlink = 1;
	      $start = $k + 1;
	      next EXON1;
	  }
      }
      
    } # end of EXON2 
    
    # if you haven't found any match for this exon1, start again from the first exon2:
    if ($foundlink == 0){
      $start = 0;
    }
    
  }   # end of EXON1      
 
  return ($object_map,$merge, $overlaps, \%is_first, \%is_last);
 
}



############################################################

sub _process_comparison{
  my ($self, $object_map, $tran1, $tran2, $merge, $is_first, $is_last) = @_;
  
  my $verbose = $self->verbose;

  my $splice_mismatch = $self->splice_mismatch;
  
  my %is_first = %$is_first;
  my %is_last  = %$is_last;

  # exons mapped from tran1:
  my @list1 = $object_map->list1();
  
  # exons mapped from tran2:
  my @list2 = $object_map->list2();
  
  #print STDERR scalar(@list1)." elements in list 1\n";
  #print STDERR scalar(@list2)." elements in list 2\n";
  
  ############################################################
  # the simplest case: when they match over one exon only:
  ############################################################
  if ( scalar(@list1)==1 && scalar(@list2)==1 ){
    
    my $merge = 0;
    ############################################################
    # if it is a single-exon transcript overlap: if it matches 
    # to the first or the last, we leave the open end unconstrained
    ############################################################
    if( (  $is_first{ $list1[0] } && $is_last{ $list1[0] } 
	   &&
	   $is_first{ $list2[0] }
	   &&
	   $self->_check_high_site( $list1[0], $list2[0], $tran2 )
	)
	||
	(  $is_first{ $list1[0] } && $is_last{ $list1[0] } 
	   &&
	   $is_last{ $list2[0] }
	   &&
	   $self->_check_low_site(  $list1[0], $list2[0], $tran2 )
	)
	||
	(  $is_first{ $list2[0] } && $is_last{ $list2[0] } 
	   &&
	   $is_first{ $list1[0] }
	   &&
	   $self->_check_high_site( $list2[0], $list1[0], $tran1 )
	)
	||
	(  $is_first{ $list2[0] } && $is_last{ $list2[0] } 
	   &&
	   $is_last{ $list1[0] }
	   &&
	   $self->_check_low_site(  $list2[0], $list1[0], $tran1 )
	)
      ){
      print STDERR "here 1 --- merge ---\n" if $verbose;
      $merge = 1;
    }
    ############################################################
    # else, it is a single-exon overlap to an 'internal' exon
    ############################################################
    elsif( (  $is_first{ $list1[0] } && $is_last{ $list1[0] } 
	      &&
	      $self->_check_high_site( $list1[0], $list2[0], $tran2 )
	      &&
	      $self->_check_low_site(  $list1[0], $list2[0], $tran2 )
	   )
	   ||
	   (  $is_first{ $list2[0] } && $is_last{ $list2[0] } 
	      &&
	      $self->_check_high_site( $list2[0], $list1[2], $tran1 )
	      &&
	      $self->_check_low_site(  $list2[0], $list1[0], $tran1 )
	   )
	 ){
      print STDERR "list2[0]: ".$list2[0]->start."-".$list2[0]->end." is_first: ".$is_first{ $list2[0] }." is_last: ".$is_last{ $list2[0] }."\n";
      print STDERR "list1[0]: ".$list1[0]->start."-".$list1[0]->end." is_first: ".$is_first{ $list1[0] }." is_last: ".$is_last{ $list1[0] }."\n";
      
      print STDERR "here 2 --- merge ---\n" if $verbose;
      $merge = 1 ;
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
      $merge = 1;
    }
    ############################################################
    # we have already dealt with single-exon against single-exon overlap
    ############################################################
    else{
      print STDERR "No merge\n" if $verbose;
      $merge = 0;
    }
    return $merge;  
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
    elsif ( $is_first{ $exon1 } ||   $is_first{ $exon2 } ){
      
      if ( ( $is_first{ $exon1 } 
	     &&
	     $self->_check_low_site( $exon1, $exon2, $tran2 )
	     &&
	     abs( $exon2->end - $exon1->end ) <= $splice_mismatch
	   )
	   ||
	   ( $is_first{ $exon2 }
	     &&
	     $self->_check_low_site( $exon2, $exon1, $tran1 )
	     &&
	     abs( $exon1->end - $exon2->end ) <= $splice_mismatch
	   )
	 ){
	print STDERR "yes1\n";
	next PAIR;
      }
      else{
	print STDERR "CLASH\n" if $verbose;
	return 0;
      }
    
    }
    ############################################################
    # both could be the last one
    ############################################################
    elsif ( $is_last{ $exon1} && $is_last{ $exon2 } 
	    &&
	    abs( $exon2->start - $exon1->start ) <= $splice_mismatch
	  ){
      print STDERR "yes2\n";
      next PAIR;
    }
    ############################################################
    # one of them could be the last one
    ############################################################
    elsif (  $is_last{ $exon1 } || $is_last{ $exon2 } ){ 
      
      if ( ( $is_last{ $exon1 } 
	     &&
	     $self->_check_high_site( $exon1, $exon2, $tran2 )
	     &&
	     abs( $exon2->start - $exon1->start ) <= $splice_mismatch
	   )
	   ||
	   ( $is_last{ $exon2 }
	     &&
	     $self->_check_high_site( $exon2, $exon1, $tran1 )
	     &&
	     abs( $exon1->start - $exon2->start ) <= $splice_mismatch
	   )
	 ){
	print STDERR "yes3\n";
	next PAIR;
      }
      else{
        print STDERR "CLASH\n" if $verbose;
	return 0;
      }
    
    }
    ############################################################
    # we have already covered the case: first overlaps last
    ############################################################
    elsif( abs( $exon1->start - $exon2->start ) <= $splice_mismatch
	   &&
	   abs( $exon1->end - $exon2->end ) <= $splice_mismatch
	 ){
      print STDERR "yes4\n";
      next PAIR;
    }
    else{	    
      print STDERR "here3: CLASH\n" if $verbose;
      return 0;
    }
  } # end of PAIR
  
  return $merge;
}

############################################################

# this method checks whether the high end of an external exon overlaps two exons:
#
#                       last_exon
#          #######-----################
#         ########-----#########------#########
#                      middle_exon

sub _check_high_site{
  my ($self, $last_exon, $middle_exon, $tran ) = @_;
  
  # middle_exon is from tran, and last_exon is from the other transcript (not used here)
  #print STDERR "check_high_site(): checking exon ".$last_exon->start."-".$last_exon->end."\n";
  #print STDERR "against:\n";
  #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript($tran);


  ############################################################
  # allow any excess, as long as it does not overlap the next exon
  ############################################################
  my @exons = sort { $a->start <=> $b->start } @{$tran->get_all_Exons};
  #print STDERR "Exons: @exons\n";
  #print STDERR "middle_exon = $middle_exon : ".$middle_exon->start."-".$middle_exon->end."\n";
  for ( my $i=0; $i<= $#exons; $i++){
    if ( $middle_exon == $exons[$i] ){
      if ( $i+1 <= $#exons && $last_exon->overlaps($exons[$i+1]) ){
	print STDERR "overlaps with next exon\n"; 
	return 0;
      }
      else{
	
	############################################################
	# if we restrict external splice sites, 
	# we simply check that they do not exceed the splice_mismatch
	############################################################
	if ( $self->restrict_external_splice_site ){
	  if ($last_exon->end - $middle_exon->end <= $self->splice_mismatch ){
	    return 1;
	  }
	  else{
	    return 0;
	  }
	}
	return 1;
	
      }
    }
  }
  return 1;
}

############################################################

# this method checks whether the low end of an external exon overlaps two exons:
#
#                    first_exon
#                ###############------########
#         ########-----#########------#########
#                      middle_exon

sub _check_low_site{
  my ($self, $first_exon, $middle_exon, $tran ) = @_;
  #print STDERR "check_low_site(): checking exon ".$first_exon->start."-".$first_exon->end."\n";
  #print STDERR "against:\n";
  #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript($tran);
  

  # middle_exon is from tran, and first_exon is from the other transcript (not used here)
  
  ############################################################
  # allow any excess, as long as it does not overlap the previous exon
  ############################################################
  my @exons = sort { $a->start <=> $b->start } @{$tran->get_all_Exons};
  #print STDERR "Exons: @exons\n";
  #print STDERR "middle_exon = $middle_exon : ".$middle_exon->start."-".$middle_exon->end."\n";
  for ( my $i=0; $i<= $#exons; $i++){
    if ( $middle_exon == $exons[$i] ){
      if ( $i-1 >= 0 && $first_exon->overlaps($exons[$i-1]) ){
	print STDERR "overlaps with previous exon\n"; 
	return 0;
      }
      else{
	
	############################################################
	# if we restrict external splice sites, 
	# we simply check that they do not exceed the splice_mismatch
	############################################################
	if ( $self->restrict_external_splice_site ){
	  if ($middle_exon->start - $first_exon->start <= $self->splice_mismatch ){
	    return 1;
	  }
	  else{
	    return 0;
	  }
	}
	return 1;
      }
    }
  }
  
  return 1;
}

############################################################



1;
