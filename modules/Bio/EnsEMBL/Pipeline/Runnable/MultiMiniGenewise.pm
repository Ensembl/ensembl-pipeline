#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::MultiMiniGenewise

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::Runnable::MultiMiniGenewise->new(-genomic  => $genseq,
								  -features => $features)

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

package Bio::EnsEMBL::Pipeline::Runnable::MultiMiniGenewise;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Pipeline::Runnable::MiniGenewise;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI );

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@_);    
           
    my( $genomic, $features,$seqfetcher, $endbias, $minimum_intron, $exon_padding, $minimum_feature_length) = 
      $self->_rearrange([qw(GENOMIC
                            FEATURES
                            SEQFETCHER
                            ENDBIAS  
                            MINIMUM_INTRON
                            EXON_PADDING
                            MINIMUM_FEATURE_LENGTH
			   )],
			@args);


    $self->throw("No genomic sequence input")                     unless defined($genomic);
    $self->throw("No seqfetcher provided")                        unless defined($seqfetcher);
    $self->throw("No features input")                             unless defined($features);
    
    $self->throw("[$genomic] is not a Bio::PrimarySeqI")          unless $genomic->isa("Bio::PrimarySeqI");
    $self->throw("[$seqfetcher] is not a Bio::DB::RandomAccessI") unless $seqfetcher->isa("Bio::DB::RandomAccessI");
    
    $self->genomic_sequence($genomic)           if defined($genomic);
    $self->seqfetcher($seqfetcher)              if defined($seqfetcher);
    $self->endbias($endbias)                    if defined($endbias);
    $self->features($features)                  if defined($features);
    $self->_minimum_intron($features)           if defined($minimum_intron);
    $self->_exon_padding($features)             if defined($exon_padding);
    $self->_minimum_feature_length($features)   if defined($minimum_feature_length);
    #print STDERR @$features." have be passed into MultiMiniGenewise\n";
    return $self;
  }


=head2 features

  Arg [1]   : arrayref 
  Function  : sets varible to arrayref
  Returntype: arrayref
  Exceptions: throws if not given an arrayref or if elements of array aren't featurepairs'
  Caller    : $self
  Example   : my $features = $self->features();

=cut




sub features {
  my ($self,$features) = @_;
  
  if (!defined($self->{_features})) {
    $self->{_features} = [];
  }
  if (defined($features)) {
    if (ref($features) eq "ARRAY") {
      foreach my $f (@$features) {
	if ($f->isa("Bio::EnsEMBL::FeaturePair")) {
	  push(@{$self->{_features}},$f);
	} else {
	  $self->throw("Object [$f] is not a Bio::EnsEMBL::FeaturePair");
	}
      }
    } else {
      $self->throw("[$features] is not an array ref.");
    }
  }
  return $self->{_features};
}



=head2 object accessors

  Arg [1]   : Object of correct type
  Function  : get/set object
  Returntype: object
  Exceptions: throws if object not of correct type
  Caller    : $self
  Example   : my $genomic = $self->genomic_sequence;

=cut


  
sub genomic_sequence {
    my( $self, $value ) = @_;    
    if ($value) {
      #need to check if passed sequence is Bio::Seq object
      $value->isa("Bio::PrimarySeqI") || $self->throw("Input isn't a Bio::PrimarySeqI");
      $self->{'_genomic_sequence'} = $value;
    }
    return $self->{'_genomic_sequence'};
  }


sub seqfetcher {

  my( $self, $value ) = @_;    
  
  if (defined($value)) {
    #need to check if we are being passed a Bio::DB::RandomAccessI object
    $self->throw("[$value] is not a Bio::DB::RandomAccessI") unless $value->isa("Bio::DB::RandomAccessI");
    $self->{'_seqfetcher'} = $value;
  }
  return $self->{'_seqfetcher'};
}



=head2 get_all_feature_by_id

  Arg [1]   : none
  Function  : arranges all feature into hash keyed by hseqname, each element containing an anonymous array of features with that name
              also produces an hash of key hseqname and value of score which is used to sort an array of hseqname 
  Returntype: hasfref and array ref
  Exceptions: warns and skips if a feature doesn't have a hseqname'
  Caller    : $self
  Example   : my ($idhash, $idarray) = $self->get_all_features_by_id;

=cut



sub get_all_features_by_id {
    my( $self) = @_;
    
    my  %idhash;
    my  %scorehash;
    my $feature_count = 0;
  FEAT: foreach my $f (@{$self->features}) {
      $feature_count++;
      if (!$f->hseqname) {
        $self->warn("No hit name for " . $f->seqname . "\n");
        next FEAT;
      } 
      if ($idhash{$f->hseqname}) {
        push(@{$idhash{$f->hseqname}},$f);
      } else {
        #print STDERR "Dealing with ".$f->hseqname."\n";
        $idhash{$f->hseqname} = [];
        push(@{$idhash{$f->hseqname}},$f);
      }
      if ($scorehash{$f->hseqname}) {
        if ($f->score > $scorehash{$f->hseqname}) {
          $scorehash{$f->hseqname} = $f->score;
        }
      } else {
        $scorehash{$f->hseqname} = $f->score;
      }
        }
      
      my @sorted_ids = keys %idhash;
      #print STDERR "there are ".$feature_count." features being sorted\n";
      @sorted_ids = sort {$scorehash{$b} <=> $scorehash{$a}} @sorted_ids;
      
      return (\%idhash,\@sorted_ids);
    }
    


=head2 get_Sequence

  Arg [1]   : sequence id/accession
  Function  : calls seqfetcher to get sequence of id passed in
  Returntype: Bio::PrimarySeq
  Exceptions: throws if fetch sequence failed
  Caller    : $self
  Example   : $seq = $self->get_Sequence('AB10323');

=cut




sub get_Sequence {
  my ($self,$id) = @_;
  
  if (defined($self->{'_seq_cache'}{$id})) {
    return $self->{'_seq_cache'}{$id};
  } 
  
  my $seqfetcher = $self->seqfetcher;    
  
  my $seq;
  #print STDERR "Fetching ".$id." sequence\n";
  eval {
    $seq = $seqfetcher->get_Seq_by_acc($id);
  };
  
  if ($@) {
    $self->throw("Problem fetching sequence for [$id]: [$@]\n");
  }
  if(!$seq){
    print STDERR "have had problems fetching sequence for ".$id."\n";
  }
  return $seq;
}



 
=head2 run

  Arg [1]   : none 
  Function  : creates and runs genewise
  Returntype: 1
  Exceptions: throws if features don't' have hseqnames
  Caller    : $runnableDB
  Example   : $runnable->run;

=cut


sub run {
  my ($self) = @_;

  my ($fhash,$ids) = $self->get_all_features_by_id;
  #print STDERR "have ".@$ids." ids\n";
  my $failed_count = 0;
  foreach my $id (@$ids) {
    #print STDERR "$id\n";
    
    my @features = @{$fhash->{$id}};
    #print STDERR "have ".@features." features\n";
   
  
    my $pepseq = $self->get_Sequence($features[0]->hseqname);
      
    my @forward;
    my @reverse;
    if ($pepseq) {
      foreach my $f (@features){
        if($f->strand == 1){
          push(@forward, $f);
        }elsif($f->strand == -1){
          push(@reverse, $f);
        }else{
          $self->throw("unstranded feature not much use for gene building\n") 
        }
      }
      if(@forward){
	my @extras = $self->_find_extras(@forward);
	#print STDERR "Number of features       = ".scalar(@forward)."\n";
	#print STDERR "Number of extra features = ".scalar(@extras)   ."\n";
	if(@extras){
	  my $runnable  = new Bio::EnsEMBL::Pipeline::Runnable::MiniGenewise(
									     -genomic => $self->genomic_sequence,
									     -protein => $pepseq,
									     -features=> \@forward,
									     -endbias => $self->endbias);
	  
	  $runnable->run;
	  ##print STDERR "MiniGenewise output " . $runnable->output . "\n";
	
	  push(@{$self->{_output}},$runnable->output);
	}else{
	  print STDERR $id." had no extra features on the forward strand\n";
    
	}
      }

      if(@reverse){
	my @extras = $self->_find_extras(@reverse);
	#print STDERR "Number of features       = ".scalar(@reverse)."\n";
	#print STDERR "Number of extra features = ".scalar(@extras)   ."\n";
	if(@extras){
	  my $runnable  = new Bio::EnsEMBL::Pipeline::Runnable::MiniGenewise(
									     -genomic => $self->genomic_sequence,
									     -protein => $pepseq,
									     -features=> \@reverse,
									     -endbias => $self->endbias);
	
	  $runnable->run;
	  #print STDERR "MiniGenewise output " . $runnable->output . "\n";
	
	  push(@{$self->{_output}},$runnable->output);
	}else{
	  print STDERR $id." had no extra features on the reverse strand\n";
	 
	}
      }
      
    } else {
      $self->warn("Can't fetch sequence for " 
                  . $features[0]->hseqname . "\n");
      $failed_count++;
    }
  
  }
  if($failed_count == @$ids){
    $self->throw("Can't find any sequences for the ids which match ".
                 $self->genomic_sequence->name); 
  }
  return 1;
  
}




=head2 _find_extras

  Arg [1]   : array of FeaturePairs
  Function  : checks for overlaps of feature with genewise output
  Returntype: array of FeaturePair
  Exceptions: none
  Caller    : $self
  Example   : @new_features = $self->_find_extras(@features);

=cut



#
# checks feature length is long enough then that feature doesn't overlapped with any of output evidence from genewise 
# before adding to the array which will be used to make the miniseq
#

sub _find_extras {
  my ($self,@features) = @_;
  
  my @output = $self->output;
  my @new;
  
 FEAT: 
  foreach my $f (@features) {
    
    $f->slice($self->genomic_sequence);
    $f->seqname($f->slice->name);
    my $found = 0;

    # does this need to be hardcoded?
    # smaller for virgin genomes?
    if (($f->end - $f->start) < $self->_minimum_feature_length) {
      next FEAT;
    }      
    
    foreach my $out (@output) {
      foreach my $sf ($out->sub_SeqFeature) {
        $sf->slice($self->genomic_sequence);
        $sf->seqname($sf->slice->name);
        #print STDERR "Comparing ".$sf." to ".$f."\n";
        if ($f->overlaps($sf)) {
          $found = 1;
        }
      }
    }
    
    if ($found == 0) {
      push(@new,$f);
    }
  }
  return @new;
}





=head2 accessors this is pod for the simple acessors

  Arg [1]   : int (value to be set)
  Function  : set/get a int value
  Returntype: int
  Exceptions: non
  Caller    : $self
  Example   : $self->_minimum_intron;

=cut





sub _minimum_intron {
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    $self->{'_minimum_intron'} = $arg;
  }

  return $self->{'_minimum_intron'};
}


sub _exon_padding {
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    $self->{'_padding'} = $arg;
  }

  return $self->{'_padding'};
  
}


sub _minimum_feature_length {
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    $self->{'_feature_length'} = $arg;
  }
  #does the default want to be 50?
  return $self->{'_feature_length'} || 50;
  
}

sub endbias {
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    $self->{'_endbias'} = $arg;
  }
  if (!defined($self->{'_endbias'})) {
    $self->{'_endbias'} = 0;
  }
  return $self->{'_endbias'};
}


1;

