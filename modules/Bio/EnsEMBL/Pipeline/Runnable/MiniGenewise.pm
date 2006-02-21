#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::MiniGenewise

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::Runnable::MiniGenewise->new(
      -genomic  => $genseq,
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

package Bio::EnsEMBL::Pipeline::Runnable::MiniGenewise;

use vars qw(@ISA);
use strict;


use Bio::EnsEMBL::Pipeline::Runnable::Genewise;
use Bio::EnsEMBL::Pipeline::MiniSeq;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::DnaPepAlignFeature;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI );

sub new {
    my ($class,@args) = @_;

    my $self = $class->SUPER::new(@args);    
           
    $self->{'_features'} = [];
    
    my( $genomic, $protein, $features,$endbias, $gap, $extension, $matrix, 
        $minimum_intron, $terminal_padding, $exon_padding, $max_iterate_dist)
      = $self->_rearrange([qw(GENOMIC
			      PROTEIN
			      FEATURES
			      ENDBIAS
			      GAP
			      EXTENSION
			      MATRIX
			      MINIMUM_INTRON
			      TERMINAL_PADDING
			      EXON_PADDING
                              MAX_SPLIT_ITERATE_DIST
			     )],
			  @args);
    
    $self->throw("No genomic sequence input")                     unless defined($genomic);
    $self->throw("No protein sequence input")                     unless defined($protein);
    $self->throw("No features input")                             unless defined($features);
    
    $self->throw("[$genomic] is not a Bio::PrimarySeqI")          unless $genomic->isa("Bio::PrimarySeqI");
    $self->throw("[$protein] is not a Bio::PrimarySeqI")          unless $genomic->isa("Bio::PrimarySeqI");

    $max_iterate_dist = 0 if not defined $max_iterate_dist;
    
    $self->genomic_sequence($genomic)             if defined($genomic);
    $self->protein_sequence($protein)             if defined($protein);
    $self->endbias($endbias)                      if defined($endbias);
    $self->gap($gap)                              if defined($gap);
    $self->extension($extension)                  if defined($extension);
    $self->matrix($matrix)                        if defined($matrix);
    $self->features($features)                    if defined($features);
    $self->_minimum_intron($minimum_intron)       if defined($minimum_intron); 
    $self->_exon_padding($exon_padding)           if defined($exon_padding);
    $self->_terminal_padding($terminal_padding)     if defined($terminal_padding);
    $self->_max_split_iterate_distance($max_iterate_dist);

    return $self;
  }



sub make_MiniSeq {
    my ($self, @features) = @_;

    my $pairaln  = new Bio::EnsEMBL::Analysis::PairAlign;
    
    my @genomic_features;
		
    my $prevend     = 0;
    my $count       = 0;
    my $mingap      = $self->_minimum_intron;
    my $seqname     = $features[0]->seqname;

    @features = sort {$a->start <=> $b->start} @features;
    
  FEAT: foreach my $f (@features) {
      my $start = $f->start - $self->_exon_padding;
      my $end   = $f->end   + $self->_exon_padding;

      if ($start < 1) { 
        $start = 1;
      }
      
      if ($end   > $self->genomic_sequence->length) {
        $end = $self->genomic_sequence->length;
      }
      
      my $gap     =    ($start - $prevend);
      
      # Extend the region is the gap between features is
      # below a certain size - otherwise start a new region
      
      if ($count > 0 && ($gap < $mingap)) {
	
	if ($end < $prevend) { 
	  $end = $prevend;
	}
	
	$genomic_features[$#genomic_features]->end($end);
	$prevend     = $end;
	
      } else {
	
	my $newfeature = new Bio::EnsEMBL::Feature;
	
        $newfeature->seqname   ($f->seqname);
        $newfeature->start     ($start);
        $newfeature->end       ($end);
        $newfeature->strand    (1);
        $newfeature->slice($self->genomic_sequence);
	
	push(@genomic_features,$newfeature);
	
	$prevend     = $end;
      }
      $count++;
    }

    # make a forward strand sequence, but tell genewise to run reversed if the 
    # features are on the reverse strand - handled by _is_reversed

    @genomic_features = sort {$a->start <=> $b->start } @genomic_features;

    # pad the termini of the sequence with 20K (configurable?) of genomic 
    # sequence - to catch small terminal exons that are missed by blast
    my $adjusted_start = $genomic_features[0]->start;
    $adjusted_start -= $self->_terminal_padding;
    if($adjusted_start < 1 ) {
      $adjusted_start = 1;
    }
    $genomic_features[0]->start($adjusted_start);

    my $adjusted_end = $genomic_features[$#genomic_features]->end;
    $adjusted_end += $self->_terminal_padding;
    if($adjusted_end > $self->genomic_sequence->length) {
        $adjusted_end = $self->genomic_sequence->length;
      }
    $genomic_features[$#genomic_features]->end($adjusted_end);

    my $current_coord = 1;

    foreach my $f (@genomic_features) {

      my $cdna_start = $current_coord;
      my $cdna_end   = $current_coord + ($f->end - $f->start);
      
      my $fp  = new Bio::EnsEMBL::FeaturePair();
      $fp->start($f->start);
      $fp->end($f->end);
      $fp->seqname($f->seqname);
      $fp->slice($self->genomic_sequence);
      $fp->strand($f->strand);
      $fp->hstart($cdna_start);
      $fp->hend($cdna_end);
      $fp->hstrand(1);
      $fp->hseqname($f->seqname."cdna");
      
      $pairaln->addFeaturePair($fp);
      
      $current_coord = $cdna_end+1;
    }
    
    my $miniseq = new Bio::EnsEMBL::Pipeline::MiniSeq(-id        => 'test',
						      -pairalign => $pairaln);
    
    return $miniseq;
    
  }

=head2 run

  Title   : run
  Usage   : $self->run()
  Function: Runs genewise on a MiniSeq
  Returns : none
  Args    :

=cut

sub run {
  my ($self) = @_;


  my @features = @{$self->features};
  my ($miniseq, $minigen, @gw_output, $must_try_again);


  do {
    $must_try_again = 0;

    $miniseq = $self->make_MiniSeq(@features);
    $minigen = $miniseq->get_cDNA_sequence;

    my $gw = new Bio::EnsEMBL::Pipeline::Runnable::Genewise(  -slice    => $self->genomic_sequence,
                                                              -genomic  => $minigen,
                                                              -protein  => $self->protein_sequence,
                                                              -reverse  => $self->_is_reversed,
                                                              -gap      => $self->gap,
                                                              -extension=> $self->extension,
                                                              -matrix   => $self->matrix,
                                                              -endbias  => $self->endbias);
    
    eval{
      $gw->run;
    };

    # if this failed, no point in going on!
    if($@){
      $self->warn("Genewise run failed - getting out of here\n");
      return;
    }

    # output is an array of genes in miniseq coordinates
    @gw_output = $gw->output;

    if ($self->_max_split_iterate_distance) {
      # check that we have no split feature
      my (@bridge_features);
      
      foreach my $g (@gw_output) {
        foreach my $e (@{$g->get_all_Exons}) {        
          my @pieces = $miniseq->convert_SeqFeature($e);
          
          if (@pieces > 1) {
            
            @pieces = sort {$a->start <=> $b->start} @pieces;
            
            my $can_merge = 1;
            for(my $i=1; $i < @pieces; $i++) {
              my $dist = $pieces[$i]->start - $pieces[$i-1]->end - 1;
              if ($dist > $self->_max_split_iterate_distance) {
                $can_merge = 0;
                last;
              }
            }
            
            if ($can_merge) {
              push @bridge_features, Bio::EnsEMBL::Feature->
                  new(-start => $pieces[0]->start,
                      -end   => $pieces[-1]->end);
            }                                                            
          }
        }
      }
      
      if (@bridge_features) {
        push @features, @bridge_features;
        $must_try_again = 1;
        print STDERR ("Warning : gene has feature that splits when mapped; trying again\n");
      }
    }

  } while ($must_try_again);

  my $strand = $self->_is_reversed ? -1 : 1;
  
  GENE:
  foreach my $gene (@gw_output) {
    my $converted_gene = new Bio::EnsEMBL::Gene;
    $converted_gene->type($gene->type);
    
    foreach my $transcript(@{$gene->get_all_Transcripts}){
      my @converted_exons;
      my @all_supp_features;
        
      # need to convert all the exons and all the supporting features
      my @exons = @{$transcript->get_all_Exons};
      my $nexon = scalar(@exons);

      TEXON:
      for(my $i=0; $i < $nexon; $i++) {
        my $exon = $exons[$i];
        $exon->strand($strand);
        
        # need to convert whole exon back to genomic coordinates
        my @genomics = $miniseq->convert_SeqFeature($exon);
                
        if(!@genomics){
          print STDERR "Don't have anything returned by MiniSeq\n";
          next GENE;
        }
        if ($#genomics > 0) {
          # all hell will break loose as the sub alignments will probably not map cheerfully
          # and we may start introducing in frame stops ...
          # for now, ignore this feature.
          print STDERR ("Warning : feature converts into " . scalar(@genomics), " features; ");
          if ($i == $nexon - 1) {
            print STDERR "Skipping last exon $i\n";
            next TEXON;
          } elsif ($exon->phase == $exons[$i+1]->phase) {
            print STDERR "Skipping compatible exon $i\n";
            next TEXON;
          } else {
            print STDERR "Cannot skip exon without adjusting phase; skipping gene\n";
            next GENE;
          }
        }  else {
          my $genomic_exon = new Bio::EnsEMBL::Exon;
          $genomic_exon->start    ($genomics[0]->start);
          $genomic_exon->end      ($genomics[0]->end);
          $genomic_exon->strand   ($strand);
          $genomic_exon->seqname  ($self->genomic_sequence->seq_region_name);
          $genomic_exon->phase    ($exon->phase);
          $genomic_exon->end_phase($exon->end_phase);
          $genomic_exon->slice    ($self->genomic_sequence);
          
          # also need to convert each of the sub alignments back to genomic coordinates
          foreach my $sf (@{$exon->get_all_supporting_features}) {
            my @features;
            my @ungapped = $sf->ungapped_features;
            foreach my $aln (@ungapped) {
              my @alns = $miniseq->convert_PepFeaturePair($aln);
              if ($#alns > 0) {
                print STDERR "Warning : sub_align feature converts into > 1 features " . scalar(@alns) . "\n";
              }
              
              my $align = new Bio::EnsEMBL::DnaPepAlignFeature(-features => \@alns);
              $align->seqname($self->genomic_sequence->seq_region_name);
              $align->hseqname($self->protein_sequence->id);
              $align->slice($self->genomic_sequence);
              # needs fix
              $align->score(100);
              push @features,$align;
            }
            push @all_supp_features,@features;
            my $gapped = new Bio::EnsEMBL::DnaPepAlignFeature(-features => \@features);
            $genomic_exon->add_supporting_features($gapped);
          }
          
          push(@converted_exons,$genomic_exon);
        }
      }
      
      # make a new transcript from @converted_exons
      my $converted_transcript  = new Bio::EnsEMBL::Transcript;
      my $converted_translation = new Bio::EnsEMBL::Translation;
      $converted_transcript->translation($converted_translation);
      
      if (scalar(@all_supp_features)) {
        my $daf = Bio::EnsEMBL::DnaPepAlignFeature->new(-features => \@all_supp_features);
        $converted_transcript->add_supporting_features($daf);
      }
      
      if ($#converted_exons < 0) {
        print STDERR "Odd.  No exons found\n";
        return undef;
        
      } else {
        
        if ($converted_exons[0]->strand == -1) {
          @converted_exons = sort {$b->start <=> $a->start} @converted_exons;
	} else {
	  @converted_exons = sort {$a->start <=> $b->start} @converted_exons;
	}
	
	$converted_translation->start_Exon($converted_exons[0]);
	$converted_translation->end_Exon  ($converted_exons[$#converted_exons]);
        
	# phase is relative to the 5' end of the transcript (start translation)
	if ($converted_exons[0]->phase == 0) {
	  $converted_translation->start(1);
	} elsif ($converted_exons[0]->phase == 1) {
	  $converted_translation->start(3);
	} elsif ($converted_exons[0]->phase == 2) {
	  $converted_translation->start(2);
	}

	$converted_translation->end  ($converted_exons[$#converted_exons]->end -
				      $converted_exons[$#converted_exons]->start + 1);
	foreach my $exon(@converted_exons){
	  $converted_transcript->add_Exon($exon);
	}
      }

      $converted_transcript->slice($self->genomic_sequence);
      $converted_gene->add_Transcript($converted_transcript);
      push(@{$self->{_output}}, $converted_gene);
    }
  }
}

sub _is_reversed {
  my ($self) = @_;
  
  if (!defined($self->{_reverse})) {
    
    my $strand = 0;
    my $fcount = 0;
    my $rcount = 0;
    
    foreach my $f (@{$self->features}) {
      if ($f->strand == 1) {
	$fcount++;
      } elsif ($f->strand == -1) {
	$rcount++;
      }
    }
    
    if ($fcount > $rcount) {
      $self->{_reverse} = 0;
    } else {
      $self->{_reverse} = 1;
    }
  } 
  return $self->{_reverse};
}

sub features {
  my ($self,$arg) = @_;
  
  if (!defined($self->{_features})) {
    $self->{_features} = [];
  }
  
  if (defined($arg)) {
    if (ref($arg) eq "ARRAY") {
      my @f = @$arg;
      
      foreach my $f (@f) {  
	
	$f->isa("Bio::EnsEMBL::FeaturePair") || $self->throw("Input [$f] isn't a Bio::EnsEMBL::FeaturePair");
	push(@{$self->{'_features'}},$f);
      }
    } else {
      $arg->isa("Bio::EnsEMBL::FeaturePair") || $self->throw("Input [$arg] isn't a Bio::EnsEMBL::FeaturePair");
      push(@{$self->{'_features'}},$arg);
    }
  }
  return $self->{_features};
}    


sub _minimum_intron {
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    $self->{'_minimum_intron'} = $arg;
  }
  # does this need to be hardcoded?
  return $self->{'_minimum_intron'} || 1000;
}


sub _exon_padding {
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    $self->{'_padding'} = $arg;
  }
  # does this need to be hardcoded
  return $self->{'_padding'} || 200;
  
}

sub _terminal_padding {
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    $self->{'_terminal_padding'} = $arg;
  }
  # does this need to be hardcoded
  return $self->{'_terminal_padding'} || 20000;
  
}


sub _max_split_iterate_distance {
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    $self->{'_iterate'} = $arg;
  }

  return $self->{'_iterate'};
  
}



=head2 genomic_sequence

Title   :   genomic_sequence
  Usage   :   $self->genomic_sequence($seq)
 Function:   Get/set method for genomic sequence
    Returns :   Bio::Seq object
  Args    :   Bio::Seq object

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
=head2 genomic_sequence

    Title   :   genomic_sequence
    Usage   :   $self->genomic_sequence($seq)
    Function:   Get/set method for genomic sequence
    Returns :   Bio::Seq object
    Args    :   Bio::Seq object

=cut

sub protein_sequence {
    my( $self, $value ) = @_;    
    if ($value) {
        #need to check if passed sequence is Bio::Seq object
        $value->isa("Bio::PrimarySeqI") || $self->throw("Input isn't a Bio::PrimarySeqI");
        $self->{'_protein_sequence'} = $value;
    }
    return $self->{'_protein_sequence'};
}

=head2 endbias

    Title   :   endbias
    Usage   :   $self->endbias($endbias)
    Function:   Get/set method for endbias
    Returns :   
    Args    :   

=cut

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

=head2 gap

    Title   :   gap
    Usage   :   $self->gap($gap)
    Function:   Get/set method for gap
    Returns :   
    Args    :   

=cut

sub gap {
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    $self->{'_gap'} = $arg;
  }
  if (!defined($self->{'_gap'})) {
    $self->{'_gap'} = 0;
  }
  return $self->{'_gap'};
  }

=head2 extension

    Title   :   extension
    Usage   :   $self->extension($extension)
    Function:   Get/set method for extension
    Returns :   
    Args    :   

=cut

sub extension {
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    $self->{'_extension'} = $arg;
  }
  if (!defined($self->{'_extension'})) {
    $self->{'_extension'} = 0;
  }
  return $self->{'_extension'};
  }

=head2 matrix

    Title   :   matrix
    Usage   :   $self->matrix($matrix)
    Function:   Get/set method for matrix
    Returns :   
    Args    :   

=cut

sub matrix {
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    $self->{'_matrix'} = $arg;
  }
  if (!defined($self->{'_matrix'})) {
    $self->{'_matrix'} = 0;
  }
  return $self->{'_matrix'};
  }

1;
