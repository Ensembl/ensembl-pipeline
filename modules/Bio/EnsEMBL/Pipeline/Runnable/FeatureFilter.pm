
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

    $search = Bio::EnsEMBL::Pipeline::Runnable::SearchFilter->new( -coverage => 5,
								  -minscore => 100,
								  -maxevalue => 0.001);
    

   my @filteredfeatures = $search->run(@features);

=head1 DESCRIPTION

Filters search results, such as Blast, on a number of criteria. The
most important ones are minscore, maxevalue, coverage. Coverage means
that only XX number of completely containing higher scores will be
permitted for this feature

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



  $minscore  = -100000 if !defined $minscore;
  $maxevalue = 0.1     if !defined $maxevalue;
  $coverage  = 10      if !defined $coverage;
  $prune     = 0       if !defined $prune;

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
   # sort by score to know that the first score for a hseqname is its best
   
   foreach my $f ( @input ) {

       if( $f->score > $minscore ) {
	   if( $f->can('evalue') && defined $f->evalue ) {
	       if( $f->evalue < $maxevalue ) {

		   if( $validhit{$f->hseqname} < $f->score ) {
		       $validhit{$f->hseqname} = $f->score;
		   }
		   if( $f->end > $maxend ) {
		       $maxend = $f->end;
		   }


	       }
	   } else {
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
   my @accepted;

   # we accept all feature pairs which are valid and meet coverage criteria
   FEATURE :
   foreach my $hseqname ( @inputids ) {

       my $hole = 0;

       foreach my $f ( @{$hitarray{$hseqname}} ) {
	   # only mark if this feature is valid
	   if( $f->score > $minscore || ($f->can('evalue') && defined $f->evalue) ) {
	       for my $i ( $f->start .. $f->end ) {
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

       foreach my $f ( @{$hitarray{$hseqname}} ) {
	 push(@accepted,$f);
	 for my $i ( $f->start .. $f->end ) {
	   $list[$i]++; 
	 }
       }
     }


   if ($self->prune) {

     my %hit;

     foreach my $f (@accepted) {
       if (!defined($hit{$f->hseqname})) {
	 $hit{$f->hseqname} = [];
       }
       push(@{$hit{$f->hseqname}},$f);
     }
     
     my @new;

     foreach my $hseqname (keys %hit) {
#       print STDERR "Pruning $hseqname\n";
       my @tmp = $self->prune_features(@{$hit{$hseqname}});
       push(@new,@tmp);
       push(@{$self->{'_output'}},@tmp);
     }
     return @new;
   } else {
     push(@{$self->{'_output'}},@accepted);
     return @accepted;
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

  my @new;
  my $depth = 5;

  my @clusters;
  my @cluster_starts;
  my @cluster_ends;

#  print STDERR "Before " . scalar(@input) . "\n";

  @input = sort {$a->start <=> $b->start} @input;

  FEAT: foreach my $f (@input) {
      #print STDERR "Processing feature " . $f->gffstring . "\n";
    my $found = 0;

    my $count = 0;
    CLUS: foreach my $clus (@clusters) {
      my @clusf = @$clus;

      if ($f->end < $cluster_starts[$count] || $f->start > $cluster_ends[$count]) {
#	print STDERR "Skipping cluster\n";

	next CLUS;
      }
      if (!($f->end < $cluster_starts[$count] || $f->start > $cluster_ends[$count])) {
	
	#foreach my $f2 (@clusf) {
	#if ($f->overlaps($f2)) {
#	print STDERR "Found existing cluster\n";
	$found = 1;
	push(@$clus,$f);
	
	if ($f->start < $cluster_starts[$count]) {
	  $cluster_starts[$count] = $f->start;
	}
	if ($f->end   < $cluster_ends[$count]) {
	  $cluster_ends[$count] = $f->end;
	}

	next FEAT;
      }
      #      }
      $count++;
    }
    if ($found == 0) {
#      print STDERR "found new cluster\n";
      my $newclus = [];
      push (@$newclus,$f);
      push(@clusters,$newclus);

      $cluster_starts[$count] = $f->start;
      $cluster_ends[$count] = $f->end;

    }
  }


  foreach my $clus (@clusters) {
    my $count = 0;

    my @tmp = @$clus;

    @tmp = sort {$b->score <=> $a->score} @tmp;

    while ($count < 5 && $#tmp >= 0) {
      push(@new,shift @tmp);
      $count++;
    }
    if ($#tmp >= 0) {
#      print STDERR "Removing " . scalar(@tmp) . " features\n";
    }
  }

#  print STDERR "After " . scalar(@new) . "\n";
  
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
