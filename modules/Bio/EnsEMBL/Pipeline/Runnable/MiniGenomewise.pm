#
#
# Cared for by EnsEMBL <ensembl-dev@ebi.ac.uk>
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::MiniGenomewise

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::Runnable::MiniGenomewise->new(-genomic     =
> $genseq,
								    -transcripts => $transcripts)

    $obj->run

    my @newtranscripts = $obj->output;


=head1 DESCRIPTION

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::Runnable::MiniGenomewise;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::RootI;
use Bio::EnsEMBL::Pipeline::Runnable::Genomewise;
use Bio::EnsEMBL::Pipeline::MiniSeq;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::Analysis;
use Bio::PrimarySeqI;
use Bio::SeqIO;
use Bio::Root::RootI;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI );

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@_);    
           
    my( $genomic, $transcripts, $analysis ) = $self->_rearrange([qw(GENOMIC
							 TRANSCRIPTS
							 ANALYSIS)],
						     @args);
    
    $self->throw("No genomic sequence input")                     unless defined($genomic);
    $self->throw("[$genomic] is not a Bio::PrimarySeqI")          unless $genomic->isa("Bio::PrimarySeqI");

    $self->genomic_sequence($genomic) if defined($genomic);
    $self->analysis($analysis);

    if (defined($transcripts)) {
	if (ref($transcripts) eq "ARRAY") {
	    my @t = @$transcripts;
	    
	    foreach my $t (@t) {
		$self->add_Transcript($t);
	    }
	} else {
	    $self->throw("[$transcripts] is not an array ref.");
	}
    }
    
    return $self;
}

sub analysis{
  my ($self,$analysis) = @_;
  if ( $analysis ){
    $self->{analysis} = $analysis;
  }
  return $self->{analysis};
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

=head2 add_Transcript

    Title   :   add_Transcript
    Usage   :   $self->add_Transcript($t)
    Function:   Adds a transcript to the object for realigning
    Returns :   Bio::EnsEMBL::Transcript
    Args    :   Bio::EnsEMBL::Transcript

=cut

sub add_Transcript {
    my( $self, $value ) = @_;
    
    if(!defined($self->{'_transcripts'})) {
	$self->{'_transcripts'} = [];
    }

    if ($value) {
        $value->isa("Bio::EnsEMBL::Transcript") || $self->throw("Input isn't a Bio::EnsEMBL::Transcript");
	push(@{$self->{'_transcripts'}},$value);
    }
}


=head2 get_all_Transcripts

    Title   :   get_all_Transcripts
    Usage   :   @t = $self->get_all_Transcripts
    Function:   Returns the array oftranscripts
    Returns :   @Bio::EnsEMBL::Transcript
    Args    :   none

=cut


sub get_all_Transcripts {
    my( $self, $value ) = @_;
    
    return (@{$self->{'_transcripts'}});
}

=head2 miniseq

    Title   :   miniseq
    Usage   :   $self->miniseq($ms)
    Function:   get/set for miniseq
    Returns :   Bio::EnsEMBL::Pipeline::MiniSeq
    Args    :   Bio::EnsEMBL::Pipeline::MiniSeq

=cut

sub miniseq {
    my( $self, $miniseq ) = @_;
    
    if ($miniseq) {
        $miniseq->isa("Bio::EnsEMBL::Pipeline::MiniSeq") || $self->throw("[$miniseq] isn't a Bio::EnsEMBL::Pipeline::MiniSeq");
	$self->{'_miniseq'} = $miniseq;
    }

    return $self->{'_miniseq'};
}

=head2 make_miniseq

    Title   :   make_miniseq
    Usage   :   $self->make_miniseq
    Function:   produces a miniseq 
    Returns :   nothing
    Args    :   none

=cut

sub make_miniseq {
    my ($self) = @_;

    my @exons;
    my $strand;

    foreach my $transcript($self->get_all_Transcripts){
    EXON:
      foreach my $exon($transcript->get_all_Exons){
	# we expect to be given a set of transcripts that are all on the same 
	# strand - genomewise can't deal with mixed strand predictions and we 
	# really dont; want to be doing with it here
	if(!defined $strand){
	  # this is the first exon
	  $strand = $exon->strand;
	}
	
	if($exon->strand != $strand){
	  $self->throw("Mixed strands in input exons! I'm outta here!\n");
	}

	if ( $exon->start > $exon->end){
	  my $start = $exon->start;
	  my $end   = $exon->end;
	  $self->throw("Exon $exon has start $start larger than end $end");
	}
	push(@exons, $exon);
      }
    }



    my $count  = 0;
    my $mingap = $self->minimum_intron;
    my $pairaln  = new Bio::EnsEMBL::Analysis::PairAlign;

    my @genomic_features;

    my $prevend     = 0;
    
  FEAT: foreach my $exon (@exons) {

      my $start = $exon->start;
      my $end   = $exon->end;

      #print STDERR "In make_miniseq(), exon -> start: $start, end: $end\n";

      $start = $exon->start - $self->exon_padding;
      $end   = $exon->end   + $self->exon_padding;
      
      if ($start < 1) { $start = 1;}
      if ($end   > $self->genomic_sequence->length) {$end = $self->genomic_sequence->length;}

      #print STDERR "In make_miniseq(), feature-> start: $start, end: $end\n";

      my $gap = ($start - $prevend);

      if ($count > 0 && ($gap < $mingap)) {
	  if ($end < $prevend) { $end = $prevend;}
	  $genomic_features[$#genomic_features]->end($end);
	  $prevend = $end;
      } else {
	    my $newfeature = new Bio::EnsEMBL::SeqFeature;
	    
	    $newfeature->seqname   ('genomic');
	    $newfeature->start     ($start);
	    $newfeature->end       ($end);
	    $newfeature->strand    (1);
	    $newfeature->attach_seq($self->genomic_sequence);

	    push(@genomic_features,$newfeature);
	    
	    $prevend = $end;
	}
	$count++;
    }

    # Now we make the cDNA features
    # but presumably only if we actually HAVE any ... 
    return unless scalar(@genomic_features);

    my $current_coord = 1;
    
    # make a forward strand sequence - genomewise can't handle reversed strand
    # strand flipping has to be dealt with outside Genomewise/MiniGenomewise 
    # eg see EST_GeneBuilder
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
	$current_coord = $cdna_end+1;
    }
	
    #changed id from 'Genomic' to seqname
    my $miniseq = new Bio::EnsEMBL::Pipeline::MiniSeq(-id        => 'miniseq',
						      -pairalign => $pairaln);

    $self->miniseq($miniseq);

}

=head2 minimum_intron

    Title   :   minimum_intron
    Usage   :   my $minimum_intron = $self->minimum_intron
    Function:   get/set for minimum intron length
    Returns :   integer
    Args    :   optional integer

=cut


sub minimum_intron {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{'_minimum_intron'} = $arg;
    }
    return $self->{'_minimum_intron'} || 1000;
}

=head2 exon_padding

    Title   :   exon_padding
    Usage   :   my $exon_padding = $self->exon_padding
    Function:   get/set for exon_padding length
    Returns :   integer
    Args    :   optional integer

=cut
    
sub exon_padding {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{'_padding'} = $arg;
    }

    return $self->{'_padding'} || 1000;

}

=head2 run

  Title   : run
  Usage   : $self->run()
  Function: Runs genomewise on MiniSeq representation of genomic sequence, plus the array of transcripts
  Returns : none
  Args    : 

=cut

sub convert_transcript_to_miniseq{
  my ($self, $transcript) = @_;
  foreach my $exon($transcript->get_all_Exons){
    my $start = $self->miniseq->pairAlign->genomic2cDNA($exon->start);
    my $end   = $self->miniseq->pairAlign->genomic2cDNA($exon->end);

    $exon->start($start);
    $exon->end($end);
  }

  return $transcript;
}

=head2 run

  Title   : run
  Usage   : $self->run()
  Function: Runs genomewise on MiniSeq representation of genomic sequence, plus the array of transcripts
  Returns : none
  Args    : 

=cut

sub run {
  my ($self) = @_;
  

  # make miniseq representation to cover all transcripts
  $self->make_miniseq;

  my $genomewise = new Bio::EnsEMBL::Pipeline::Runnable::Genomewise;
  $genomewise->seq($self->miniseq->get_cDNA_sequence);
  $genomewise->analysis($self->analysis);

#  $genomewise->seq($self->genomic_sequence);
  foreach my $t($self->get_all_Transcripts){
    my $converted = $self->convert_transcript_to_miniseq($t);
    $genomewise->add_Transcript($converted);
    #$genomewise->add_Transcript($t);
  }
  
  $genomewise->run;

  $self->convert_output($genomewise->output);
}


=head2 convert_output

  Title   : convert_output
  Usage   : $self->convert_output(@transcripts)
  Function: Converts exons predicted by genomewsie back into virtual contig coordinates
  Returns : nothing
  Args    : @transcripts

=cut

sub convert_output{ 
  my ($self, @transcripts) = @_;
  my $analysis_obj    = new Bio::EnsEMBL::Analysis
    (-db              => undef,
     -db_version      => undef,
     -program         => "genomewise",
     -program_version => 1,
     -gff_source      => 'genomewise',
     -gff_feature     => 'similarity');

  #print STDERR "number of transcripts: " . scalar(@transcripts)  . "\n";
  
 TRANSCRIPT: 
  foreach my $transcript (@transcripts) {
    #print STDERR "In MiniGenomewise: transcript is a ".ref($transcript)."\n";
    my @newexons;

    # test
#    print STDERR "\nIn MiniGenomewise.convert_output\n";
#    print STDERR " Transcript        : ".$transcript."\n";
#    print STDERR " Translation       : ".$transcript->translation."\n";
#    print STDERR " translation starts: ".$transcript->translation->start."\n";
#    print STDERR " translation ends  : ".$transcript->translation->end."\n";
#    print STDERR " start exon: ".$transcript->translation->start_exon
#      ." starts: ".$transcript->translation->start_exon->start
#      ." ends: ".$transcript->translation->start_exon->end."\n";
#    print STDERR " end  exon : ".$transcript->translation->end_exon
#      ." starts: ".$transcript->translation->end_exon->start
#      ." ends: ".$transcript->translation->end_exon->end."\n";
    
    # get the translation from the transcript
    my $translation = $transcript->translation;
    my $start_exon  = $translation->start_exon;
    my $end_exon    = $translation->end_exon;
    my ($ss,$se)    = ( $start_exon->start, $start_exon->end );
    my ($es,$ee)    = ( $end_exon->start, $end_exon->end );
    

    # we need not create a new translation for the new coordinate system
    #my $new_translation = new Bio::EnsEMBL::Translation;    
 
    # redefine the one we had
    $translation->start($translation->start);
    $translation->end  ($translation->end);
    
    # convert coordinates exon by exon
    my $ec = 0;
  EXON:
    foreach my $exon ($transcript->get_all_Exons) {
      $ec++;
   
      my $end_phase = $exon->end_phase;
      my $phase  = $exon->phase;
      my $strand = $exon->strand;

      # get the supporting evidence
      my @evidence = $exon->each_Supporting_Feature;

      my @genomics = $self->miniseq->convert_SeqFeature($exon);         
      if ($#genomics > 0) {
	# for now, ignore this exon.
	print STDERR "Warning : miniseq exon converts into > 1 genomic exon " 
	             . scalar(@genomics) . " ignoring exon $ec\n";
	
	# this can screw up the translation as the phases may not be consistent anymore
	# there is a check for this below
	next EXON;
      }
      
      foreach my $f (@genomics) {
	my $new_exon = new Bio::EnsEMBL::Exon;
	$new_exon->start($f->start);
	$new_exon->end($f->end);
	$new_exon->phase($phase);
	$new_exon->strand($strand);
	$new_exon->end_phase($end_phase);
	
	# transfer the supporting evidence!!!
	foreach my $evi ( @evidence ){
	  $new_exon->add_Supporting_Feature( $evi );
	}
	
	#BUGFIX: This should probably be fixed in Bio::EnsEMBL::Analysis
	$new_exon->seqname($exon->seqname);
	$new_exon->analysis($analysis_obj);
	#end BUGFIX
	push(@newexons, $new_exon);    
	
	# check whether $exon is the start or end exon
	if ( $exon->start == $ss && $exon->end == $se ) {
	  #print STDERR " >> start_exon found, converting ".$exon." into ".$new_exon."\n";
	  $translation->start_exon($new_exon);
	}
	if ( $exon->start == $es && $exon->end == $ee ) {
	  #print STDERR " >> end_exon found, converting   ".$exon." into ".$new_exon."\n";
	  $translation->end_exon($new_exon);
	}
      }
      
      # flush out old exons from transcript and replace them with newly remapped exons
      $transcript->flush_Exon;
      foreach my $exon(@newexons){
	$transcript->add_Exon($exon);
      }
    
      # check the consistency of the phases:
      my @exons = $transcript->get_all_Exons;
      for (my $i = 1; $i <= $#exons; $i++) {
	if ( $exons[$i-1]->end_phase != $exons[$i]->phase  ){
	  print STDERR "transcript has phase inconsistency, skipping it...\n";
	  $self->_print_Transcript($transcript);
	  next TRANSCRIPT;
	}
      }

    }
    $transcript->sort;
    
    # include the modified translation
    $transcript->translation($translation);
    
    $self->output($transcript);  
  
  } # end of TRANSCRIPT
}

sub _print_Transcript{
  my ($self,$tran) = @_;
  $tran->sort;
  my @exons = $tran->get_all_Exons;
  foreach my $exon (@exons){
    print STDERR $exon->start."-".$exon->end." phase: ".$exon->phase." end_phase: ".$exon->end_phase." strand: ".$exon->strand."\n";
  }
}


=head2 output

  Title   : output
  Usage   : $self->output
  Function: Get/set for output array
  Returns : An array of Bio::EnsEMBL::Transcript
  Args    : optional Bio::EnsEMBL::Transcript

=cut

sub output {
    my ($self, $transcript) = @_;

    if (!defined($self->{'_output'})) {
	$self->{'_output'} = [];
    }

    if(defined $transcript){
      if($transcript->isa("Bio::EnsEMBL::Transcript")){
	push(@{$self->{'_output'}},$transcript);
      }
      else{
	$self->warn("[$transcript] is not a Bio::EnsEMBL::Transcript - not adding to output array\n");
      }
    }

    return @{$self->{'_output'}};
}

1;

