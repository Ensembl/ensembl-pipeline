#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::ECR

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::Runnable::ECR->new
    (
    '-features'       => $features,
    '-union'          => 1 # 0 => intersection
    '-merge'          => 10,
    
    $obj->run

    my @newfeatures = $obj->output;


=head1 DESCRIPTION

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::Runnable::ECR;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::SimpleFeature;


@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

sub new {
  my ($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  
  my ($features, $union, $merge) = $self->_rearrange([qw(FEATURES
                                                         UNION
                                                         MERGE)],
                                                         @args);

  $union = 0 if not defined $union;
  $merge = 0 if not defined $merge;

  $self->throw("You must supply a reference to an array of features\n") if not defined $features;

  $self->merge($merge);
  $self->union($union);
  $self->features($features);
  
  return $self;
}





=head2 run

  Title   : run
  Usage   : $self->run()
  Function: 
  Returns : none
  Args    : 

=cut

sub run {
  my ($self) = @_;

  # convert the features into a more malleable form
  my (@feat_lists, $min_start, $max_end, $seqname);
  foreach my $flist (@{$self->features}) {
    my @this_list;
    foreach my $f (@$flist) {
      if (not defined $min_start or $f->start < $min_start) {
        $min_start = $f->start;
      }
      if (not defined $max_end or $f->end > $max_end) {
        $max_end = $f->end;
      }
      if (defined $seqname) {
        if ($seqname ne $f->seqname) {
          $self->warn("ECR supplied with featuers from different sequences. ".
                      "Results will be meaningless");
        } 
      } else {
        $seqname = $f->seqname;
      }

      push @this_list, { start => $f->start,
                         end   => $f->end,
                       };

    }
    @this_list = sort {$a->{start} <=> $b->{start}} @this_list;
    push @feat_lists, \@this_list;
  }

  my $result;
  if ($self->union) {
    $result = $self->_project_features($min_start, $max_end, @feat_lists);
  } else {
    $result = $self->_merge_features($min_start, $max_end, @feat_lists);
  }

  # create the simple features here
  my @output_features;
  foreach my $f (@$result) {
    my $sf = Bio::EnsEMBL::SimpleFeature->new();
    $sf->seqname($seqname);
    $sf->start($f->{start});
    $sf->end($f->{end});
    $sf->strand(0);

    push @output_features, $sf;
  }
  
  $self->output(@output_features);

  return 1;
}



# _merge_features
#   takes two more lists of features and returns a list of features 
#   corresponding to the regions that are covered in all of the input lists

sub _merge_features {
  my ($self, $min_start, $max_end, @lists) = @_;

  my $range = $max_end - $min_start + 1;

  my $bitstring = "1" x $range; 
  my @retsegs;
  
  # return $_[0] if scalar(@_) == 1; 

  foreach my $list (@lists) {
    my $thisbitstring = "0" x $range;
    
    if ($list) {
      foreach my $seg (sort {$a->{start} <=> $b->{start}} @$list) {
        my ($st, $en) = ($seg->{start} - $min_start, $seg->{end} - $min_start);
        substr($thisbitstring, $st, $en - $st + 1) = 
            "1" x ($en - $st + 1);
      }
      $bitstring &= $thisbitstring;
    }
  }
  
  # finally, locate the strings of ones
  while( $bitstring =~ /(1+)/g ) {
    my $end = pos($bitstring) + $min_start - 1;
    my $match = $1;
    
    my $start = $end - length($match) + 1;
    
    if (@retsegs and ($start - $retsegs[-1]->{end} - 1 <= $self->merge)) {
      $retsegs[-1]->{end} = $end;
    }
    else {
      push @retsegs, { start => $start, end => $end  };
    }
    
  }
  
  return \@retsegs;
}



# _project_features
#   Takes two or more lists of non-overlapping segments and
#   returns a new list of segments comprising regions for
#   which all the bases occur in at least one list 

sub _project_features {
  my ($self, $min_start, $max_end, @lists) = @_;
  
  my $range = $max_end - $min_start + 1;

  my $bitstring = "0" x $range; 
  my @retsegs;
  
  foreach my $list (@lists) {
    if ($list) {
      foreach my $seg (sort {$a->{start} <=> $b->{start}} @$list) {
        my ($st, $en) = ($seg->{start} - $min_start, $seg->{end} - $min_start);
        substr($bitstring, $st, $en - $st + 1) = 
            "1" x ($en - $st + 1);
      }
    }
  }
  
  # finally, locate the strings of ones
  while( $bitstring =~ /(1+)/g ) {
    my $end = pos($bitstring) + $min_start - 1;
    my $match = $1;
    
    my $start = $end - length($match) + 1;
    
    if (@retsegs and ($start - $retsegs[-1]->{end} - 1 <= $self->merge)) {
      $retsegs[-1]->{end} = $end;
    }
    else {
      push @retsegs, { start => $start, end => $end  };
    }
  }
  
  return \@retsegs;
}




=head2 features

    Title   :   features
    Usage   :   $self->features($features)
    Function:   ref. to the array fo features to be projected
    Returns :   
    Args    :   

=cut

sub features {
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    $self->{'_features'} = $arg;
  }
  
  if (!exists($self->{'_features'})) {
    $self->{'_features'} = [];
  }    

  return $self->{'_features'};
}



=head2 union

    Title   :   union
    Usage   :   $self->union($union)
    Function:   Binary flag that determines whether the given feature set is 
       projected as a union (if true) ir an intersection (if false).
    Returns :   
    Args    :   

=cut

sub union {
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    $self->{'_union'} = $arg;
  }
  
  if (!exists($self->{'_union'})) {
    $self->{'_union'} = 0;
  }    

  return $self->{'_union'};
}


=head2 merge

    Title   :   merge
    Usage   :   $self->merge($merge)
    Function:   After the projection, resulting features that are within $merge 
      of each other will be merged into one
    Returns :   
    Args    :   

=cut

sub merge {
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    $self->{'_merge'} = $arg;
  }
  
  if (!exists($self->{'_merge'})) {
    $self->{'_merge'} = 0;
  }    

  return $self->{'_merge'};
}


=head2 output

    Title   :   merge
    Usage   :   $self->merge($output)
    Function:
    Returns :   
    Args    :   

=cut

sub output {
  my ($self,@feats) = @_;
  
  if (@feats) {
    $self->{'_output'} = [@feats];
  }
  
  if (!exists($self->{'_output'})) {
    $self->{'_output'} = [];
  }    

  return @{$self->{'_output'}};
}



1;
