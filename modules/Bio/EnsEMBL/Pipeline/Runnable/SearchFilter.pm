
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

    $search = Bio::EnsEMBL::Pipeline::Runnable::SearchFilter->new(
					     -runnable => $blastrunnable);
                                             -coverage => 5,
                                             -minscore => 100,
                                             -maxevalue => 0.001);
    


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


package Bio::EnsEMBL::Pipeline::Runnable::SearchFilter;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::EnsEMBL::Pipeline::RunnableI;

use Bio::EnsEMBL::Pipeline::RunnableI;


@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

sub new {
  my($class,@args) = @_;

  my $self = {};
  bless $self,$class;
  

  my($runnable,$minscore,$maxevalue,$coverage) = $self->_rearrange(
								   [qw( RUNNABLE
									MINSCORE
									MAXEVALUE
									COVERAGE
									)]
								   ,@args);



  $minscore  = -100000 if !defined $minscore;
  $maxevalue = 0.1     if !defined $maxevalue;
  $coverage  = 10      if !defined $coverage;

  if( !defined $runnable ) {
      $self->throw("Must have a runnable for search filter");
  }



  $self->runnable($runnable);
  $self->minscore($minscore);
  $self->maxevalue($maxevalue);
  $self->coverage($coverage);

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
   my ($self,$dir) = @_;

   $self->runnable->run($dir);
   
   my ($minscore,$maxevalue,$coverage);
   my (@accepted);

   $minscore = $self->minscore;
   $maxevalue= $self->maxevalue;
   $coverage = $self->coverage;


   my %validhit;
   my %hitarray;
   my @input = $self->runnable->output;

   # first- scan across all features, considering
   # valid to be > minscore < maxevalue

   my $maxend   = 0;

   # all featurepairs have a score. 
   # some may have an evalue.

   # valid hits are stored in a hash of arrays
   # sort by score to know that the first score for a hseqname is its best
   
   @input = sort { $a->score <=> $b->score } @input; 

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

   my @inputids = sort { $validhit{$a} <=> $validhit{$b} } keys %validhit; 


   # we accept all feature pairs which are valid and meet coverage criteria

   my @accepted;
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

   push(@{$self->{'_output'}},@accepted);


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




=head2 runnable

 Title   : runnable
 Usage   : $obj->runnable($newval)
 Function: 
 Returns : value of runnable
 Args    : newvalue (optional)


=cut

sub runnable{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'runnable'} = $value;
    }
    return $obj->{'runnable'};

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


