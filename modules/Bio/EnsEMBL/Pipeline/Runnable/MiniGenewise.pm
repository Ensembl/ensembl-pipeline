#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::MiniGenewise

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::Runnable::MiniGenewise->new(-genomic  => $genseq,
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
use Bio::EnsEMBL::SeqFeature;

use Bio::PrimarySeqI;
use Bio::SeqIO;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI );

sub new {
    my ($class,@args) = @_;

    my $self = $class->SUPER::new(@args);    
           
    $self->{'_features'} = [];
    
    my( $genomic, $protein, $features,$endbias)
      = $self->_rearrange([qw(GENOMIC
			      PROTEIN
			      FEATURES
			      ENDBIAS
			     )],
			  @args);
    
    
    $self->throw("No genomic sequence input")                     unless defined($genomic);
    $self->throw("No protein sequence input")                     unless defined($protein);
    $self->throw("No features input")                             unless defined($features);
    
    $self->throw("[$genomic] is not a Bio::PrimarySeqI")          unless $genomic->isa("Bio::PrimarySeqI");
    $self->throw("[$protein] is not a Bio::PrimarySeqI")          unless $genomic->isa("Bio::PrimarySeqI");
    
    $self->genomic_sequence($genomic) if defined($genomic);
    $self->protein_sequence($protein) if defined($protein);
    $self->endbias($endbias)          if defined($endbias);
    $self->features($features)        if defined($features);
    
    return $self;
  }



sub make_miniseq {
    my ($self) = @_;

    my $pairaln  = new Bio::EnsEMBL::Analysis::PairAlign;
    
    my @genomic_features;
		
    my $prevend     = 0;
    my $count       = 0;
    my $mingap      = $self->minimum_intron;
    my @features    = $self->features;
    my $seqname     = $features[0]->seqname;

    @features = sort {$a->start <=> $b->start} @features;
    
  FEAT: foreach my $f (@features) {

      my $start = $f->start - $self->exon_padding;
      my $end   = $f->end   + $self->exon_padding;

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
	
	my $newfeature = new Bio::EnsEMBL::SeqFeature;
	
        $newfeature->seqname   ($f->hseqname);
        $newfeature->start     ($start);
	$newfeature->end       ($end);
	$newfeature->strand    (1);
	$newfeature->attach_seq($self->genomic_sequence);
	
	push(@genomic_features,$newfeature);
	
	$prevend     = $end;
      }
      $count++;
    }
    
    my $current_coord = 1;
    
    # make a forward strand sequence, but tell genewise to run reversed if the 
    # features are on the reverse strand - handled by is_reversed

    @genomic_features = sort {$a->start <=> $b->start } @genomic_features;

    foreach my $f (@genomic_features) {
      
      my $cdna_start = $current_coord;
      my $cdna_end   = $current_coord + ($f->end - $f->start);
      
      my $tmp = new Bio::EnsEMBL::SeqFeature( #-seqname => $f->seqname.'.cDNA',
					     -start => $cdna_start,
					     -end   => $cdna_end,
					     -strand => 1);
      
      my $fp  = new Bio::EnsEMBL::FeaturePair(-feature1 => $f,
					      -feature2 => $tmp);
      
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

  my $miniseq = $self->make_miniseq;
	my $minigen = $miniseq->get_cDNA_sequence;

  my $gw = new Bio::EnsEMBL::Pipeline::Runnable::Genewise(  -genomic => $minigen,
							    -protein => $self->protein_sequence,
							    -reverse => $self->is_reversed,
							    -endbias => $self->endbias);
  
  
  $gw->run;
  
  # output is a list of Features (one per exon) with subseqfeatures 
  # representing the ungapped sub alignments for each exon
  
  my @f = $gw->output;
  
  my @newf;
  
  my $strand = 1;
  
  if ($self->is_reversed == 1) {
    $strand = -1;
  }
  
  my $ec = 0;
  
 FEAT: 
  foreach my $f (@f) {
    $ec++;
    
    $f->strand($strand); 
    
    # need to convert whole exon back to genomic coordinates
    my @genomics = $miniseq->convert_SeqFeature($f);         
    my $gf;
    
    if ($#genomics > 0) {
      
      # all hell will break loose as the sub alignments will probably not map cheerfully 
      # and we may start introducing in frame stops ...
      # for now, ignore this feature.
      
      print STDERR ("Warning : feature converts into > 1 features " . scalar(@genomics) . " Ignoring exon $ec\n");
      
      next FEAT;
    }  else {
      $gf = $genomics[0];
    }
    
    $gf->phase    ($f->phase);
    $gf->end_phase($f->end_phase);
    $gf->strand   ($strand);
    $gf->seqname  ($self->genomic_sequence->id);
    $gf->score    (100);
    
    # also need to convert each of the sub alignments back to genomic coordinates
    
    foreach my $aln ($f->sub_SeqFeature) {
      my @alns = $miniseq->convert_PepFeaturePair($aln);
      
      if ($#alns > 0) {
	print STDERR "Warning : sub_align feature converts into > 1 features " . scalar(@alns) . "\n";
      }
      
      foreach my $a(@alns) {
	$a->strand($strand); 
	$a->hstrand(1);      
	$a->seqname($self->genomic_sequence->id);
	$a->hseqname($self->protein_sequence->id);
	
	# Maybe put a check in that this really is a sub feature
	
	$gf->add_sub_SeqFeature($a,'');
      }
    }
    
    push(@newf,$gf);
    
  }
  
  
  # $fset holds a list of (genomic) SeqFeatures (one fset per gene) plus their constituent exons and
  # sub_SeqFeatures representing ungapped alignments making up the exon:protein alignment
  
  my $fset = new Bio::EnsEMBL::SeqFeature();
  
  foreach my $nf (@newf) {
    $fset->add_sub_SeqFeature($nf,'EXPAND');
    $fset->seqname($nf->seqname);
  }
  
  # Hmm - only add to the output if > 1 exon???
  
  if(scalar($fset->sub_SeqFeature) > 0){
    push(@{$self->{'_output'}},$fset);
  } else { 
    print STDERR $fset." won't be outputed has ".$fset->sub_SeqFeature." exons\n";
  }
  
}

sub is_reversed {
  my ($self) = @_;
  
  if (!defined($self->{_reverse})) {
    
    my $strand = 0;
    my $fcount = 0;
    my $rcount = 0;
    
    foreach my $f ($self->features) {
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
  return @{$self->{_features}};
}    


sub minimum_intron {
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    $self->{'_minimum_intron'} = $arg;
  }
  return $self->{'_minimum_intron'} || 1000;
}


sub exon_padding {
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    $self->{'_padding'} = $arg;
  }
  
  return $self->{'_padding'} || 200;
  
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


1;

