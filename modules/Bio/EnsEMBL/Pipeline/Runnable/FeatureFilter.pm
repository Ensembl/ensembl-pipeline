
#
# Ensembl module for Bio::EnsEMBL::Pipeline::Runnable::SearchFilter
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::SearchFilter - Filters a search runnable

=head1 SYNOPSIS

    $search = Bio::EnsEMBL::Pipeline::Runnable::SearchFilter->new( -coverage  => 5,
								   -minscore  => 100,
								   -maxevalue => 0.001,
								   -prune     => 1
								 );
    

   my @filteredfeatures = $search->run(@features);

=head1 DESCRIPTION

Filters search results, such as Blast, on a number of criteria. The
most important ones are minscore, maxevalue, coverage. Coverage means
that only XX number of completely containing higher scores will be
permitted for this feature. The option prune performs a clustering of 
the features according to overlap with respect to the genomic sequence
and take only a maximum number of features per cluster; this number being
specified by coverage.

=head1 CONTACT

Ensembl - ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Pipeline::Runnable::FeatureFilter;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::EnsEMBL::Pipeline::RunnableI;

use Bio::EnsEMBL::Pipeline::RunnableI;


@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);  

  my($minscore,$maxevalue,$coverage,$prune) = $self->_rearrange([qw(MINSCORE
								    MAXEVALUE
								    COVERAGE
								    PRUNE
							     )],
							 @args);



  $minscore  = -100000 unless $minscore;
  $maxevalue = 0.1     unless $maxevalue;
  $coverage  = 10      unless $coverage;
  $prune     = 0       unless $prune;

  $self->minscore($minscore);
  $self->maxevalue($maxevalue);
  $self->coverage($coverage);
  $self->prune   ($prune);
  
  $self->{'_output'} = [];

  return $self;
}


=head2 run

 Title   : run
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub run{
  my ($self,@input) = @_;
  
  my ($minscore,$maxevalue,$coverage);
  
  $minscore = $self->minscore;
  $maxevalue= $self->maxevalue;
  $coverage = $self->coverage;
  
  my %validhit;
  my %hitarray;
  
  
  # first- scan across all features, considering
  # valid to be > minscore < maxevalue
  
  my $maxend   = 0;
  
  # all featurepairs have a score. 
  # some may have an evalue.
  
  # valid hits are stored in a hash of arrays
  # we sort by score to know that the first score for a hseqname is its best  
  @input = sort { $b->score <=> $a->score } @input;
  
  foreach my $f ( @input ) {
    
    if( $f->score > $minscore ) {
      
      unless ( $validhit{$f->hseqname} ){
	$validhit{$f->hseqname} = 0;
      }
      
      if( $f->can('evalue') && defined $f->evalue ) {
	if( $f->evalue < $maxevalue ) {
	  
	  if( $validhit{$f->hseqname} < $f->score ) {
	    $validhit{$f->hseqname} = $f->score;
	  }
	  if( $f->end > $maxend ) {
	    $maxend = $f->end;
	  }
	  	  
	}
      }

      else {
	if( $validhit{$f->hseqname} < $f->score ) {
	  $validhit{$f->hseqname} = $f->score;
	}
	if( $f->end > $maxend ) {
	  $maxend = $f->end;
	}
	
      }
    }
    
    # irregardless of score, take if this hseqname is valid
    if( exists $validhit{$f->hseqname} == 1 ) {
      if( ! exists $hitarray{$f->hseqname} ) {
	$hitarray{$f->hseqname} = [];
      }
      push(@{$hitarray{$f->hseqname}},$f);
    }
    
  }
  
  # empty input array - saves on memory!
  
  @input = ();
  
  # perl will automatically extend this array 
  my @list;
  $list[$maxend] = 0;
  
  # sort the list by highest score
  my @inputids = sort { $validhit{$b} <=> $validhit{$a} } keys %validhit; 
  
  # this now holds the accepted hids ( a much smaller array )
  my @accepted_hids;
  
  # we accept all feature pairs which are valid and meet coverage criteria
  FEATURE :
    foreach my $hseqname ( @inputids ) {
      
      my $hole = 0;
      
      foreach my $f ( @{$hitarray{$hseqname}} ) {
	# only mark if this feature is valid
	if( $f->score > $minscore || ($f->can('evalue') && defined $f->evalue && $f->evalue<$maxevalue ) ) {
	  for my $i ( $f->start .. $f->end ) {
	    unless( $list[$i] ){
	      $list[$i] = 0;
	    }

	    if( $list[$i] < $coverage ) {
	      # accept!
	      $hole = 1;
	      last;
	    }
	  }
	}
      }
      
      if( $hole == 0 ) {
	# completely covered 
	next;
      }
      
      push ( @accepted_hids, $hseqname );
      foreach my $f ( @{$hitarray{$hseqname}} ) {
	for my $i ( $f->start .. $f->end ) {
	  $list[$i]++; 
	}
      }
    }
  
  # drop this huge array to save memory
  @list = ();
  
  if ($self->prune) {
    my @new;
     
    # prune the features per hid (per hseqname)
    foreach my $hseqname ( @accepted_hids ) {
      my @tmp = $self->prune_features(@{$hitarray{$hseqname}});
      push(@new,@tmp);
      #push(@{$self->{'_output'}},@tmp);
    }
    @accepted_hids = ();
    return @new;
  
  } 
  else {
    my @accepted_features;
    foreach my $hid ( @accepted_hids ){
      push(@accepted_features, @{$hitarray{$hid}} );
    }
    @accepted_hids = ();
    #push(@{$self->{'_output'}},@accepted_features);
    return @accepted_features;
  }
}

sub prune {
  my ($self,$arg) = @_;

  if (defined($arg)) {
    $self->{_prune} = $arg;
  }
  return $self->{_prune};
}




sub prune_features {
  my ($self,@input) = @_;

  # define the depth of the coverage
  my $depth = $self->coverage;

  # here we store the created clusters
  my @clusters;
  my @cluster_starts;
  my @cluster_ends;

  #print STDERR "Before:" . scalar(@input) . "features\n";

  # sort the features by start coordinates, this is crucial
  @input = sort {$a->start <=> $b->start} @input;

  ## this is handy to compare the two clustering methods quickly
  #my $old_method = 0;
  #my $new_method = 1;

#  if ( $new_method == 1 ){

  # create the first cluster
  my $count = 0;
  my $cluster = [];
  
  # start it off with the first feature
  my $first_feat = shift( @input );
  push (@$cluster, $first_feat);
  $cluster_starts[$count] = $first_feat->start;
  $cluster_ends[$count]   = $first_feat->end;
  
  # store the list of clusters
  push(@clusters,$cluster);
  
  # loop over the rest of the features
  
 FEATURE:
  foreach my $f ( @input ){
    #print STDERR "trying to place feature: $f ".$f->start."-".$f->end."\n";
    
    # add $f to the current cluster if overlaps and strand matched
    #print STDERR "comparing with cluster $count : "
    #  .$cluster_starts[$count]."-".$cluster_ends[$count]."\n";
    
    if (!($f->end < $cluster_starts[$count] || $f->start > $cluster_ends[$count])) {
      
      push(@$cluster,$f);
      
      # re-adjust size of cluster
      if ($f->start < $cluster_starts[$count]) {
	$cluster_starts[$count] = $f->start;
      }
      if ($f->end  > $cluster_ends[$count]) {
	$cluster_ends[$count] = $f->end;
      }
      
    }
    else{
      # else, start create a new cluster with this feature
      $count++;
      $cluster = [];
      push (@$cluster, $f);
      $cluster_starts[$count] = $f->start;
      $cluster_ends[$count]   = $f->end;
          
      # store it in the list of clusters
      push(@clusters,$cluster);
    }
  }

#}  

#### the previous clustering method ######

#if ( $old_method == 1 ){    
    
#  FEAT: 
#    foreach my $f (@input) {
#      #print STDERR "Processing feature " . $f->gffstring . "\n";
#      #print STDERR "trying to place feature: ".$f->hseqname." ".$f->start."-".$f->end."\n";
#      my $found = 0;
      
#      my $count = 0;
#    CLUS: foreach my $clus (@clusters) {
		
#	#print STDERR "comparing with cluster $count : "
#	#  .$cluster_starts[$count]."-".$cluster_ends[$count]."\n";
#	foreach my $f2 ( @$clus) {
#	  print STDERR "       ".$f2->start."-".$f2->end."\n";
#	}

#	if ($f->end < $cluster_starts[$count] || $f->start > $cluster_ends[$count]) {
#	  #print STDERR "No, go to next one\n";
	  
#	  #next CLUS;
#	}
#	if (!($f->end < $cluster_starts[$count] || $f->start > $cluster_ends[$count])) {
	  
#	  #print STDERR "Yes, we put it in this one\n";
#	  $found = 1;
#	  push(@$clus,$f);
	  
#	  if ($f->start < $cluster_starts[$count]) {
#	    $cluster_starts[$count] = $f->start;
#	  }
#	  if ($f->end   > $cluster_ends[$count]) {
#	    $cluster_ends[$count] = $f->end;
#	  }
	  
#	  next FEAT;
#	}
#	$count++;
#      }
#      if ($found == 0) {
#	#print STDERR "found new cluster\n";
#	my $newclus = [];
#	push (@$newclus,$f);
#	push(@clusters,$newclus);
	
#	$cluster_starts[$count] = $f->start;
#	$cluster_ends[$count]   = $f->end;
	
#      }
#    }
#  }
    
  # put here the final result
  my @new;

  # we take up to a maximum of $depth features per cluster
  foreach my $clus (@clusters) {
    my $count = 0;
    my @tmp = @$clus;
    @tmp = sort {$b->score <=> $a->score} @tmp;
    
    ## test
    #print STDERR "cluster:\n";
    #foreach my $f (@$clus){
    #   print STDERR "feature ".$f->hseqname." ".$f->start."-".$f->end." score(".$f->score.")\n";
    # }
      
    while ($count < $depth && $#tmp >= 0) {
      my $f = shift( @tmp );
      #print STDERR "accepting ".$f->start."-".$f->end." score(".$f->score.")\n";
      push(@new, $f);
      $count++;
    }
    #if ($#tmp >= 0) {
      #foreach my $f ( @tmp ){
      #print STDERR "rejecting ".$f->start."-".$f->end." score(".$f->score.")\n";
      #}
      #print STDERR "Removing " . scalar(@tmp) . " features\n";
    #}
  }
  # drop array to save memory
  @clusters = ();

  return @new;
}

  


=head2 output

 Title   : output
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub output{
   my ($self,@args) = @_;
   unless ( $self->{'_output'} ){
     $self->{'_output'} = [];
   }
   if ( @args ){
     push( @{ $self->{'_output'} }, @args);
   } 
   return @{$self->{'_output'}};
}




=head2 minscore

 Title   : minscore
 Usage   : $obj->minscore($newval)
 Function: 
 Returns : value of minscore
 Args    : newvalue (optional)


=cut

sub minscore{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'minscore'} = $value;
    }
    return $obj->{'minscore'};

}

=head2 maxevalue

 Title   : maxevalue
 Usage   : $obj->maxevalue($newval)
 Function: 
 Returns : value of maxevalue
 Args    : newvalue (optional)


=cut

sub maxevalue{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'maxevalue'} = $value;
    }
    return $obj->{'maxevalue'};

}

=head2 coverage

 Title   : coverage
 Usage   : $obj->coverage($newval)
 Function: 
 Returns : value of coverage
 Args    : newvalue (optional)


=cut

sub coverage{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'coverage'} = $value;
    }
    return $obj->{'coverage'};

}
1;
