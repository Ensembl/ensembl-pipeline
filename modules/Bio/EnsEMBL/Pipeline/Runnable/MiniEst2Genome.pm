#!/usr/local/bin/perl

#
#
# Cared for by EnsEMBL  <ensembl-dev@ebi.ac.uk>
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::MiniEst2Genome

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::Runnable::MiniEst2Genome->new(-genomic  => $genseq,
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

package Bio::EnsEMBL::Pipeline::Runnable::MiniEst2Genome;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::RootI;
use Bio::EnsEMBL::Pipeline::Runnable::Est2Genome;
use Bio::EnsEMBL::Pipeline::MiniSeq;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Pipeline::SeqFetcher;

#compile time check for executable
use Bio::EnsEMBL::Analysis::Programs qw(pfetch efetch); 
use Bio::PrimarySeqI;
use Bio::SeqIO;

use Data::Dumper;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI Bio::Root::RootI);

sub new {
  my ($class,@args) = @_;
  my $self = bless {}, $class;
  
  $self->{'_fplist'} = []; #create key to an array of feature pairs
  
  my( $genomic, $features ) = $self->_rearrange(['GENOMIC',
						 'FEATURES'], @args);
       
  $self->throw("No genomic sequence input")            unless defined($genomic);
  $self->throw("[$genomic] is not a Bio::PrimarySeqI") unless $genomic->isa("Bio::PrimarySeqI");
  
  $self->genomic_sequence($genomic) if defined($genomic);
  
  if (defined($features)) {
    if (ref($features) eq "ARRAY") {
      my @f = @$features;
      
      foreach my $f (@f) {
	$self->addFeature($f);
      }
    } else {
      $self->throw("[$features] is not an array ref.");
    }
  }
  
  return $self; # success - we hope!
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

=head2 addFeature 

    Title   :   addFeature
    Usage   :   $self->addFeature($f)
    Function:   Adds a feature to the object for realigning
    Returns :   Bio::EnsEMBL::FeaturePair
    Args    :   Bio::EnsEMBL::FeaturePair

=cut

sub addFeature {
    my( $self, $value ) = @_;
    
    if(!defined($self->{_features})) {
	$self->{_features} = [];
    }

    if ($value) {
        $value->isa("Bio::EnsEMBL::FeaturePair") || $self->throw("Input isn't a Bio::EnsEMBL::FeaturePair");
	push(@{$self->{_features}},$value);
    }
}


=head2 get_all_FeaturesbyId

    Title   :   get_all_FeaturesById
    Usage   :   $hash = $self->get_all_FeaturesById;
    Function:   Returns a ref to a hash of features.
                The keys to the hash are distinct feature ids
    Returns :   ref to hash of Bio::EnsEMBL::FeaturePair
    Args    :   none

=cut

sub get_all_FeaturesById {
    my( $self) = @_;
    
    my  %idhash;

    FEAT: foreach my $f ($self->get_all_Features) {
	print STDERR ("Feature is $f " . $f->seqname . "\t" . $f->hseqname ."\n");
    if (!(defined($f->hseqname))) {
	$self->warn("No hit name for " . $f->seqname . "\n");
	    next FEAT;
	} 
	if (defined($idhash{$f->hseqname})) {
	    push(@{$idhash{$f->hseqname}},$f);
	} else {
	    $idhash{$f->hseqname} = [];
	    push(@{$idhash{$f->hseqname}},$f);
	}

    }

    return (\%idhash);
}


=head2 get_all_Features

    Title   :   get_all_Features
    Usage   :   @f = $self->get_all_Features;
    Function:   Returns the array of features
    Returns :   @Bio::EnsEMBL::FeaturePair
    Args    :   none

=cut


sub get_all_Features {
    my( $self, $value ) = @_;
    
    return (@{$self->{_features}});
}


=head2 get_all_FeatureIds

  Title   : get_all_FeatureIds
  Usage   : my @ids = get_all_FeatureIds
  Function: Returns an array of all distinct feature hids 
  Returns : @string
  Args    : none

=cut

sub get_all_FeatureIds {
    my ($self) = @_;

    my %idhash;

    foreach my $f ($self->get_all_Features) {
	if (defined($f->hseqname)) {
	    $idhash{$f->hseqname} = 1;
	} else {
	    $self->warn("No sequence name defined for feature. " . $f->seqname . "\n");
	}
    }

    return keys %idhash;
}

=head2 make_miniseq

  Title   : make_miniseq
  Usage   : 
  Function: makes a mini genomic from the genomic sequence and features list
  Returns : 
  Args    : 

=cut

sub make_miniseq {
    my ($self,@features) = @_;

    my $seqname = $features[0]->seqname;

    print STDERR "seqname $seqname\n";

    @features = sort {$a->start <=> $b->start} @features;
    my $count  = 0;
    my $mingap = $self->minimum_intron;
    
    my $pairaln  = new Bio::EnsEMBL::Analysis::PairAlign;

    my @genomic_features;

    my $prevend     = 0;
    my $prevcdnaend = 0;
    
  FEAT: foreach my $f (@features) {
      print STDERR "Found feature - " . $f->hseqname . "\t" . $f->start . "\t" . $f->end . "\t" . $f->strand . "\n"; 

      my $start = $f->start;
      my $end   = $f->end;
      
      $start = $f->start - $self->exon_padding;
      $end   = $f->end   + $self->exon_padding;

      if ($start < 1) { $start = 1;}
      if ($end   > $self->genomic_sequence->length) {$end = $self->genomic_sequence->length;}

      my $gap     =    ($start - $prevend);

      print STDERR "Feature hstart is " . $f->hstart . "\t" . $prevcdnaend . "\n";
      print STDERR "Padding feature - new start end are $start $end\n";

      print STDERR "Count is $count : $mingap " . $gap  . "\n";

      if ($count > 0 && ($gap < $mingap)) {
	# STRANDS!!!!!
	  if ($end < $prevend) { $end = $prevend;}
	  print(STDERR "Merging exons in " . $f->hseqname . " - resetting end to $end\n");
	    
	  $genomic_features[$#genomic_features]->end($end);
	  $prevend     = $end;
	  $prevcdnaend = $f->hend;
	  print STDERR "Merged start end are " . $genomic_features[$#genomic_features]->start . "\t" .  $genomic_features[$#genomic_features]->end . "\n";
      } else {
	
	    my $newfeature = new Bio::EnsEMBL::SeqFeature;

        $newfeature->seqname ($f->hseqname);
        $newfeature->start     ($start);
	    $newfeature->end       ($end);
	    $newfeature->strand    (1);
# ???	    $newfeature->strand    ($strand);
	    $newfeature->attach_seq($self->genomic_sequence);

	    push(@genomic_features,$newfeature);
	    
	    print(STDERR "Added feature $count: " . $newfeature->start  . "\t"  . 
		  $newfeature->end    . "\t " . 
		  $newfeature->strand . "\n");

	    $prevend = $end;
	    $prevcdnaend = $f->hend; 
	    print STDERR "New end is " . $f->hend . "\n";

	}
	$count++;
    }

    # Now we make the cDNA features
    # but presumably only if we actually HAVE any ... 
    return unless scalar(@genomic_features);

    my $current_coord = 1;

    # make a forward strand sequence, est2genome runs -reverse 
    @genomic_features = sort {$a->start <=> $b->start } @genomic_features;

    foreach my $f (@genomic_features) {
	$f->strand(1);
	my $cdna_start = $current_coord;
	my $cdna_end   = $current_coord + ($f->end - $f->start);
	
	my $tmp = new Bio::EnsEMBL::SeqFeature(
					       -seqname => $f->seqname.'.cDNA',
					       -start => $cdna_start,
					       -end   => $cdna_end,
					       -strand => 1);
	
	my $fp  = new Bio::EnsEMBL::FeaturePair(-feature1 => $f,
						-feature2 => $tmp);
	
	$pairaln->addFeaturePair($fp);
	
	$self->print_FeaturePair($fp);

	$current_coord = $cdna_end+1;
    }
	
    #changed id from 'Genomic' to seqname
    my $miniseq = new Bio::EnsEMBL::Pipeline::MiniSeq(-id        => $seqname,
						      -pairalign => $pairaln);

    my $newgenomic = $miniseq->get_cDNA_sequence->seq;
    $newgenomic =~ s/(.{72})/$1\n/g;
#    print ("New genomic sequence is " . $newgenomic. "\n");
    return $miniseq;

}


=head2 minimum_introm

  Title   : minimum_intron
  Usage   : 
  Function: Defines minimum intron size for miniseq
  Returns : 
  Args    : 

=cut

sub minimum_intron {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_minimum_intron} = $arg;
    }

    return $self->{_minimum_intron} || 1000;
}

=head2 exon_padding

  Title   : exon_padding
  Usage   : 
  Function: Defines exon padding extent for miniseq
  Returns : 
  Args    : 

=cut
   
sub exon_padding {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_padding} = $arg;
    }

    return $self->{_padding} || 100;
#    return $self->{_padding} || 1000;

}

=head2 print_FeaturePair

  Title   : print_FeaturePair
  Usage   : 
  Function: for debugging
  Returns : 
  Args    : 

=cut

sub print_FeaturePair {
    my ($self,$nf) = @_;
    #changed $nf->id to $nf->seqname
    print(STDERR "FeaturePair is " . $nf->seqname    . "\t" . 
	  $nf->start . "\t" . 
	  $nf->end   . "\t(" . 
	  $nf->strand . ")\t" .
	  $nf->hseqname  . "\t" . 
	  $nf->hstart   . "\t" . 
	  $nf->hend     . "\t(" .
	  $nf->hstrand  . ")\n");
}


=head2 get_Sequence

  Title   : get_Sequence
  Usage   : my $seq = get_Sequence($id)
  Function: Fetches sequences with id $id
  Returns : Bio::PrimarySeq
  Args    : none

=cut

sub get_Sequence {
    my ($self,$id) = @_;
    my $seq;
    my $seqfetcher = new Bio::EnsEMBL::Pipeline::SeqFetcher;

    if (defined($self->{_seq_cache}{$id})) {
      return $self->{_seq_cache}{$id};
    } 
    
    $seq = $seqfetcher->run_pfetch($id);
    
    if (!defined($seq)) {
      # try efetch
      $seq = $seqfetcher->run_efetch($id);
    }
    
    if (!defined($seq)) {
      $self->throw("Couldn't find sequence for [$id]");
    }
    
    print (STDERR "Found sequence for $id [" . $seq->length() . "]\n");
    
    return $seq;

}

=head2 get_all_Sequences

  Title   : get_all_Sequences
  Usage   : my $seq = get_all_Sequences(@id)
  Function: Fetches sequences with ids in @id
  Returns : nothing, but $self->{_seq_cache}{$id} has a Bio::PrimarySeq for each $id in @id
  Args    : array of ids

=cut

sub get_all_Sequences {
  my ($self,@id) = @_;
  
 SEQ: foreach my $id (@id) {
    my $seq = $self->get_Sequence($id);
    if(defined $seq) {
      $self->{_seq_cache}{$id} = $seq;
    }
  }
}

=head2 run

  Title   : run
  Usage   : $self->run()
  Function: Runs est2genome on MiniSeq representation of genomic sequence for each EST
  Returns : none
  Args    : 

=cut

sub run {
  my ($self) = @_;
  
  my ($esthash) = $self->get_all_FeaturesById;
  
  my @ests    = keys %$esthash;
  
  $self->get_all_Sequences(@ests);

  my $analysis_obj    = new Bio::EnsEMBL::Analysis
    (-db              => undef,
     -db_version      => undef,
     -program         => "est2genome",
     -program_version => 1,
     -gff_source      => 'est2genome',
     -gff_feature     => 'similarity',);
  
 ID: foreach my $est (@ests) {
    
    my $features = $esthash->{$est};
    my @exons;
    
    print(STDERR "Processing $est\n");
    next ID unless (ref($features) eq "ARRAY");
    
    print(STDERR "Features = " . scalar(@$features) . "\n");

    # why > not >= 1?
    next ID unless (scalar(@$features) >= 1);
    
    eval {
      $self->run_blaste2g($est, $features, $analysis_obj);
    };

    if ($@) {
      print STDERR "Error running blaste2g on" . $features->[0]->hseqname . " [$@]\n";
    }

  }

#  # check the output
#  print STDERR "MiniE2G:\n";
#  foreach my $gene($self->output){
#    print STDERR "Gene: " . $gene->gffstring . "\n";
#    foreach my $exon($gene->sub_SeqFeature) {
#      print STDERR "Exon: " . $exon->gffstring . "\n";
#      foreach my $segment($exon->sub_SeqFeature){
#	print STDERR "Segment: " . $segment->gffstring . "\n";
#      }
#    }
#    print STDERR "\n";
#  }

}

=head2 run_blaste2g

  Title   : run_blaste2g
  Usage   : $self->run_blaste2g()
  Function: Runs est2genome on a MiniSeq
  Returns : none
  Args    : 

=cut

sub run_blaste2g {
  my ($self,$est,$features,$analysis_obj) = @_;
  
  my @extras  = $self->find_extras (@$features);
  print STDERR "Number of extra features = " . scalar(@extras) . "\n";
  return unless (scalar(@extras) >= 1);
  
  my $miniseq = $self->make_miniseq(@$features);
  my $hseq    = $self->get_Sequence($est);
  
  if (!defined($hseq)) {
    $self->throw("Can't fetch sequence for id [$est]\n");
  }
  
  my $eg = new Bio::EnsEMBL::Pipeline::Runnable::Est2Genome(  -genomic => $miniseq->get_cDNA_sequence,
							      -est     => $hseq);
  
  $eg->run;
  
  # each element of @exons represents an exon, with sub_seqfeatures representing the 
  # individual ungapped alignments from est_genome
  my @exons = $self->convert_output($eg, $miniseq);

  # make a SeqFeature representing a gene to hold the exons
  my $gene = new Bio::EnsEMBL::SeqFeature();
  foreach my $ex (@exons) {
    $gene->add_sub_SeqFeature($ex,'EXPAND');
    $gene->seqname($ex->seqname);
    $gene->analysis($analysis_obj);
  }
  
  push(@{$self->{_output}},$gene);
  
}

=head2 convert_output

  Title   : convert_output
  Usage   : $self->convert_output($runnable, $miniseq)
  Function: Converts 
  Returns : none
  Args    : 

=cut

sub convert_output {
  my ($self, $runnable, $miniseq) = @_;

  my @genes = $runnable->output;
  print STDERR "number of genes: " . scalar(@genes)  . "\n";
  
  if ( scalar(@genes) >1 ) {
    $self->throw("more than one gene predicted fropm est_genome - I'm outta here!\n");
  }
  
  my @genomic_exons;
  my $excount = 0;
  my $strand;
  
  GENE: foreach my $gene(@genes) {
      my @converted = $self->_convert_exons($gene, $miniseq);
      push (@genomic_exons, @converted);
    }

  return (@genomic_exons);
}

=head2 _convert_exons

  Title   : _convert_exons
  Usage   : $self->_convert_exons($gene, $miniseq)
  Function: Converts the exons of a SeqFeature representing a gene from cDNA to genomic coordinates
  Returns : Array of SeqFeatures representing exons
  Args    : Bio::SeqFeature, Bio::EnsEMBL::Pipeline::MiniSeq

=cut

sub _convert_exons {
  my ($self, $gene, $miniseq) = @_;
  my @exons = $gene->sub_SeqFeature;
  my $excount = 0;
  my @genomic_exons;

  # we need to keep track of strand information. Exons should not switch strand!
  foreach my $sf($exons[0]->sub_SeqFeature){
    if($sf->strand != $sf->hstrand){
      $sf->strand(-1);
      $sf->hstrand(1);
      $exons[0]->strand(-1);
    }
  }
  my $strand = $exons[0]->strand;
  
  # convert each exon and its component ungapped alignments back to VC coordinates
  EXON :    foreach my $exon(@exons){
    $excount++;
    foreach my $sf($exon->sub_SeqFeature){
      if($sf->strand != $sf->hstrand){
	$sf->strand(-1);
	$sf->hstrand(1);
	$exon->strand(-1);
      }
    }
    $self->throw("Trying to switch strands mid gene!\n") unless $exon->strand == $strand;
    
    # est_genome has no concept of phase, but remapping will fail if this is unset
    $exon->phase(0);
    
    my @converted_exons = $miniseq->pairAlign->convert_cDNA_feature($exon);
    if ($#converted_exons > 0) {
      # all hell will break loose as the sub alignments will probably not map cheerfully 
      # for now, ignore this feature.
      print STDERR "Warning : feature converts into > 1 features " . scalar(@converted_exons) . " ignoring exon $excount\n";
      next EXON;
    }
  
    # now deal with sub alignments
    $self->_convert_segments($exon, $miniseq, $converted_exons[0]);
  
    foreach my $ce (@converted_exons) {
      $ce->{_phase} = $exon->phase; # doesn't give us phase
      $ce->strand($strand);
      #BUGFIX: This should probably be fixed in Bio::EnsEMBL::Analysis
      $ce->seqname($exon->seqname); # urrmmmm?
      $ce->score($exon->score); # probably not
      #end BUGFIX
    }
    
    
    push (@genomic_exons, @converted_exons);
    
  }

  return (@genomic_exons);
}


=head2 _convert_segments

  Title   : _convert_segment
  Usage   : $self->_convert_segment($exon, $miniseq, $converted_exon)
  Function: Converts the subSeqFeatures of a SeqFeature representing an exon from cDNA to genomic coordinates
  Returns : 
  Args    : Bio::SeqFeature, Bio::EnsEMBL::Pipeline::MiniSeq

=cut

sub _convert_segments {
  my ($self, $exon, $miniseq, $converted_exon) = @_;

  my @f = $exon->sub_SeqFeature;
  print STDERR "***_convert_segments: " . scalar(@f). "\n";

  # each gapped exon has a set of subSeqFeatures representing ungapped component alignments to the est.
  SEGMENT: foreach my $segment($exon->sub_SeqFeature){
    # double check strands 
    if($segment->strand != $segment->hstrand){
      $segment->strand(-1);
      $segment->hstrand(1);
      $self->throw("Segment strands don't match exon strand\n") unless $exon->strand == -1;
    }
    
    # convert to genomic coords
    my @converted_segments = $miniseq->convert_FeaturePair($segment);
    if ($#converted_segments > 0) {
      # we're in for fun
      print STDERR "Warning : sub_align feature converts into > 1 features " . scalar(@converted_segments) . "\n";
    }
    
    foreach my $s(@converted_segments) {
      if($s->strand != $s->hstrand){
	$s->strand(-1);
	$s->hstrand(1);
      }

      if( ($s->strand != $segment->strand) || ($s->hstrand != $segment->hstrand) ) {
	$self->print_FeaturePair($s);
	$self->print_FeaturePair($segment);
	$self->throw("Converted segment strands don't match unconverted segment strands\n");
      }

      $s->seqname($segment->seqname);
      $s->hseqname($segment->hseqname);
      $s->score($segment->score);
      
      # try to add to the appropriate converted_exon
      # shouldn't need to expand $converted_exon when adding $s ... unless something's gone horribly wrong
      if($s->start >= $converted_exon->start && $s->end <= $converted_exon->end){
	  $converted_exon->add_sub_SeqFeature($s,'');
	}
      else {
	$self->warn("Sub align feature could not be added ...\n");
      }
      
    }
  }
}

=head2 find_extras

  Title   : find_extras
  Usage   : $self->find_extras(@features)
  Function: 
  Returns : 
  Args    : 

=cut

sub find_extras {
  my ($self,@features) = @_;
  
  my @output = $self->output;
  my @new;
  
 FEAT: foreach my $f (@features) {
    my $found = 0;
    if (($f->end - $f->start) < 50) {
      next FEAT;
    }
    #	print ("New feature\n");
    
    #$self->print_FeaturePair($f);
    foreach my $out (@output) {
      foreach my $sf ($out->sub_SeqFeature) {
	
	if (!($f->end < $out->start || $f->start >$out->end)) {
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

=head2 output

  Title   : output
  Usage   : $self->output
  Function: Returns results of est2genome as array of FeaturePair
  Returns : An array of Bio::EnsEMBL::FeaturePair
  Args    : none

=cut

sub output {
    my ($self) = @_;
    if (!defined($self->{_output})) {
	$self->{_output} = [];
    }
    return @{$self->{'_output'}};
}

1;


