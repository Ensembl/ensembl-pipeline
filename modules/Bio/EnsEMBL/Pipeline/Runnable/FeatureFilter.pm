#
# Ensembl module for Bio::EnsEMBL::Pipeline::Runnable::FeatureFilter
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::FeatureFilter - Filters a search runnable

=head1 SYNOPSIS

   $search = Bio::EnsEMBL::Pipeline::Runnable::FeatureFilter->new(
                                                    -coverage  => 5,
					            -minscore  => 100,
					            -maxevalue => 0.001,
					            -prune     => 1
				                                 );

   my @filteredfeatures = $search->run(@features);

=head1 DESCRIPTION

Filters search results, such as Blast, on several criteria. The most
important ones are minscore, maxevalue, coverage. Crudely, coverage
reduces redundant data (e.g., almost-identical ESTs with different
accession numbers) and prune reduces overlapping features (e.g., hits
to a repetitive sequence).

In detail, coverage filtering acts as follows:

  sort hit-sequence-accessions in decreasing order of
  maximum feature score;

  for each hit-sequence-accession in turn:

    if all parts of all features for hit-sequence-accession are
    already covered by other features to a depth of <coverage>

      remove all features for this hit-sequence-accession;

Within the set of features for a given hit-sequence-accession, features
are considered in decreasing order of score. Where two or more
hit-sequence-accessions have equal maximum feature score, secondary
sorting is in decreasing order of total score for the
hit-sequence-accession's features, followed by alphabetical order for
hit-sequence accession number.

The option prune allows only a maximum number of features per
strand per genomic base per hit sequence accession, this number also
being specified by the coverage parameter. Prune works on a
per-hit-sequence-accession basis and removes features (not entire
hit-sequence-accessions) until the criterion is met for each
hit-sequence-accession. Prune filtering occurs after coverage
filtering.

=head1 CONTACT

Ensembl - ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

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
  
  # print "XXX start filter with ".scalar(@input)." features \n";

  # first- scan across all features, considering
  # valid to be > minscore < maxevalue
  
  my $maxend   = 0;
  my %totalscore;     # total score per hid
  
  # all featurepairs have a score. 
  # some may have an evalue.
  
  # valid hits are stored in a hash of arrays
  # we sort by score to know that the first score for a hseqname is its best  
  @input = sort { $b->score <=> $a->score } @input;
  foreach my $f ( @input ) {

    if( $f->score > $minscore ) {
      
      unless ( $validhit{$f->hseqname} ){
	$validhit{$f->hseqname} = 0;
	$totalscore{$f->hseqname} = 0;
      }

      $totalscore{$f->hseqname} =+ $f->score;
      
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
  my @feature_counts;
  foreach my $hid(keys(%hitarray)){
    push(@feature_counts, @{$hitarray{$hid}});
  }

  # print "have ".scalar(@feature_counts)." after filtering by score\n";

  @input = ();
  
  # perl will automatically extend this array 
  my @list;
  $list[$maxend] = 0;
  
  # sort the list by highest score, then by total score for ties, and
  # alphabetically as a last resort
  my @inputids = sort {    $validhit{$b}   <=> $validhit{$a}
                        or $totalscore{$b} <=> $totalscore{$a}
			or $a cmp $b } keys %validhit;
  
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
  @feature_counts = ();

  foreach my $hid(@accepted_hids){
    push(@feature_counts, @{$hitarray{$hid}});
  }

  # drop this huge array to save memory
  @list = ();
  
  if ($self->prune) {
    my @new;

    my @all_features;
    
    # collect all the features
    foreach my $hseqname ( @accepted_hids ){
      my @tmp = $self->prune_features(@{$hitarray{$hseqname}});
      push(@new,@tmp);
      push ( @all_features, @{$hitarray{$hseqname}} );
    }
    # and prune all together taking the first '$self->coverage' according to score 
    #@new = $self->prune_features( @all_features );

    ## prune the features per hid (per hseqname)
    #foreach my $hseqname ( @accepted_hids ) {
    #  my @tmp = $self->prune_features(@{$hitarray{$hseqname}});
    #  push(@new,@tmp);
    #  #push(@{$self->{'_output'}},@tmp);
    #}
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

=head2 prune_features

 Title   : prune_features
 Usage   : @pruned_fs = $self->prune_features(@f_array);
 Function: reduce coverage of each base to a maximum of
           $self->coverage features on each strand, by
	   removing low-scoring features as necessary.
 Returns : array of features
 Args    : array of features

=cut

sub prune_features {
  my ($self, @input) = @_;
  $self->throw('interface fault') if @_ < 1;	# @input optional

  $self->warn('experimental new prune_features method!');

    my @plus_strand_fs = $self->_prune_features_by_strand(+1, @input);
    my @minus_strand_fs = $self->_prune_features_by_strand(-1, @input);
    @input = ();
    push @input, @plus_strand_fs;
    push @input, @minus_strand_fs;
    return @input;
}

=head2 _prune_features_by_strand

 Title   : _prune_features_by_strand
 Usage   : @pruned_neg = $self->prune_features(-1, @f_array);
 Function: reduce coverage of each genomic base to a maximum of
           $self->coverage features on the specified strand;
	   all features not on the specified strand are discarded
 Returns : array of features
 Args    : strand, array of features

=cut

sub _prune_features_by_strand {
   my ($self, $strand, @in) = @_;
   $self->throw('interface fault') if @_ < 2;	# @in optional

   my @input_for_strand = ();
   foreach my $f (@in) {
     push @input_for_strand, $f if $f->strand eq $strand;
   }

   return () if !@input_for_strand;

   # get the genomic first and last bases covered by any features
   my @sorted_fs = sort{ $a->start <=> $b->start } @input_for_strand;
   my $first_base = $sorted_fs[0]->start;
   my @sorted_fs = sort{ $a->end <=> $b->end } @input_for_strand;
   my $last_base = $sorted_fs[$#sorted_fs]->end;

   # fs_per_base: set element i to the number of features covering base i
   my @fs_per_base = ();
   foreach  my $base ($first_base..$last_base) {
     $fs_per_base[$base] = 0;	# initialise
   }
   foreach my $f (@input_for_strand) {
     foreach my $covered_base ($f->start..$f->end) {
       $fs_per_base[$covered_base]++;
     }
   }

   # put the worst features first, so they get removed with priority
   # we use score, which assumes they're from the same database search
   @sorted_fs = sort { $a->score <=> $b->score } @input_for_strand;

   # over_covered_bases: list of base numbers where coverage must be
   # reduced, listed worst-case-first
   my $max_coverage = $self->coverage;
   my @over_covered_bases = ();
   foreach my $base ($first_base..$last_base) {
     my $excess_fs = $fs_per_base[$base] - $max_coverage;
     if ($excess_fs > 0) {
       push @over_covered_bases, $base;
     }
   }
   @over_covered_bases = reverse sort { $fs_per_base[$a] <=> $fs_per_base[$b] }
     @over_covered_bases;

   foreach my $base (@over_covered_bases) {
     my $f_no = 0;
     while ($fs_per_base[$base] > $max_coverage) {
       my $start = $sorted_fs[$f_no]->start;
       my $end = $sorted_fs[$f_no]->end;
       if ($start <= $base and $end >= $base) {	# cut this feature
         splice @sorted_fs, $f_no, 1;	# same index will give next feature
         foreach my $was_covered ($start..$end) {
           $fs_per_base[$was_covered]--;
         }
       } else {	# didn't overlap this base, move on to next feature
         $f_no++;
       }
     }
   }
   return @sorted_fs;
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

