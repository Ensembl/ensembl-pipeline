#
# BioPerl module for PredictionGeneBuilder
#
# Cared for by EnsEMBL <ensembl-dev@ebi.ac.uk>
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

  Bio::EnsEMBL::Pipeline::Runnable::PredictionGeneBuilder

=head1 SYNOPSIS

given 

@predictions = a set of transcripts, normally ab initio predictions in the form of PredictionTranscripts

and

@features = features, any set of features, normally from the precompute-pipeline

@annotations (OPTIONAL) = is a set of genes you want to give preference to over the predictions. 

my $genecooker = new Bio::EnsEMBL::Pipeline::Runnable::PredictionGeneBuilder(
									     -predictions => \@predictions,
									     -features    => \@features,
									     -annotations => \@annotations,
									     );

Note: features and predictions should be in the same slice coordinate system.

my @supported_predictions = $genecooker->run;

=head1 DESCRIPTION

This module is a revampped version of a piece of history within 
Ensembl automatic annotation. Not only brings you back to the days where
genome annotation was being devised, but also shows how versatile and powerful
can be the Ensembl code to manipulate biological information.

This module reads features and predictions (genscan, genefinder, etc ).
Predictions are read as PredictionTranscript objects.
It will take only exons which are confirmed by the features and it will make Exon pairs
accordign to consecutive feature overlap. An exon pair is constructed as follows :

  ---------          --------    prediction exons (genscan, genefinder,... any ab initio prediction )
    -------          ----->      feature that spans an intron
    1     10        11    22        

For an exon pair to make it into a gene there must be at least 2 blast
hits (features) that span across an intron. This is called the
coverage of the exon pair. After all exon pairs have been generated for 
all the genscan exons there is a recursive routine (_recurseTranscript)
that looks for all exons that are the start of an exon pair with no preceding exons
It will then recursively link the
Exon pairs into full transcripts. The recursion will also generate alternative splicing.
The run method returns the generated transcripts.

=head1 CONTACT


=head1 CONTACT

ensembl-dev@ebi.ac.uk

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::Runnable::PredictionGeneBuilder;

use Bio::EnsEMBL::Pipeline::ExonPair;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::GeneConf qw (
					 TRANSCRIPT_ID_SUBSCRIPT
					 GB_MIN_GENSCAN_EXONS
					 GB_GENSCAN_MAX_INTRON
					 GB_TARGETTED_GW_GENETYPE
					 GB_SIMILARITY_GENETYPE
					 GB_COMBINED_GENETYPE
					 GB_MIN_FEATURE_SCORE
					 GB_MIN_FEATURE_LENGTH
					 GB_INPUTID_REGEX
					 GB_ABINITIO_TYPE
					 GB_ABINITIO_SUPPORTED_TYPE
					);
use vars qw(@ISA);
use strict;

@ISA = qw(Bio::EnsEMBL::Root);

############################################################

sub new {
  my ($class,@args) = @_;
  
  my $self = $class->SUPER::new(@args);
  
  my ($predictions,$features,$annotations) = $self->_rearrange([qw(PREDICTIONS FEATURES ANNOTATIONS)],
							       @args);
  
  unless ( @$predictions && @$features ){
    $self->throw("Must input predictions and features to PredictionGeneBuilder");
  }
  
  $self->predictions(@{$predictions});
  $self->features(@{$features});
  $self->{_annotations} = [];
  if ( @$annotations ){
    $self->annotations( @$annotations );
  }
  return $self;
}

############################################################

sub run{
  my ($self) = @_;
  
  # get all exons from the PredictionTranscripts, take only exons with overlapping similarity features, which
  # are incorporated as supporting evidence
  print STDERR "making exons...\n";
  my @supported_exons = $self->make_Exons;

  print STDERR "Number of exons supported out of predictions ".scalar( @supported_exons )."\n";
  
  # pair up the exons according to consecutive overlapping supporting evidence
  # pairs can be retrieved using 'get_all_ExonPairs'
  print STDERR "making exonpairs...\n";
  $self->make_ExonPairs(@supported_exons);
  print STDERR scalar($self->get_all_ExonPairs)." Exon pairs created\n";

  # link exons recursively according to the exon pairs with shared evidence to form transcripts
  my @linked_predictions = $self->link_ExonPairs(@supported_exons);
  print STDERR scalar(@linked_predictions)." linked transcripts generated\n";

  # check the generated transcripts:
  my @checked_predictions;
  foreach my $prediction ( @linked_predictions ){
    next unless ( $self->_check_Transcript($prediction) && $self->_check_Translation($prediction) );
    push( @checked_predictions, $prediction );
  }

  print STDERR scalar(@checked_predictions)." checked predictions generated\n";
  return @checked_predictions;
}  

############################################################
#
# METHODS DEALING WITH PREDICTION TRANSCRIPTS AND FEATURES
#
############################################################ 

=head2 make_Exons

Example : my @exons = $self->make_Exons;
Function: Turns features into exons with the help of the ab initio predictions

=cut
  
sub make_Exons {
  my ($self) = @_;
  
  my @supported_exons;
  my @features    = sort { $a->start <=> $b->start } $self->features;
  my @predictions = $self->predictions;
  my $gscount  = 1;
  
  my $ignored_exons = 0;
  
 PREDICTION:  
  foreach my $prediction (@predictions) {
    my $excount    = 1;
    unless ( @{$prediction->get_all_Exons} ){
      next PREDICTION;
    }
    
  EXON: 
    foreach my $prediction_exon (@{$prediction->get_all_Exons}) {
      
      # Don't include any genscans that are inside a genewise/combined transcript
      if ($self->annotations){
      OTHER_GENES:
	foreach my $gene ($self->annotations) {
	  my @exons    = sort {$a->start <=> $b->start} @{$gene->get_all_Exons};
	  my $g_start  = $exons[0]->start;
	  my $g_end    = $exons[$#exons]->end;
	  my $g_strand = $exons[0]->strand;
	  
	  if (!(($g_end < $prediction_exon->start) || $g_start > $prediction_exon->end)) {
	    if ($g_strand == $prediction_exon->strand) {
	      $ignored_exons++;
	      next EXON;
	    }
	  }
	}
      }
      
      my $newexon = $self->_make_Exon($prediction_exon,$excount,"genscan." . $gscount . "." . $excount );
      $newexon->find_supporting_evidence(\@features);
      
      # take only the exons that get supporting evidence
      if ( @{ $newexon->get_all_supporting_features } ){
	push(@supported_exons,$newexon);
	$excount++;
      }
    }
    
    $gscount++;
  }
  
  if ( $self->annotations ){
    print "\nIgnoring $ignored_exons genscan exons due to overlaps with genewise genes\n";
  }
  $self->supported_exons(@supported_exons);
  return @supported_exons;
}

############################################################

sub _make_Exon { 
  my ($self,$exon_prediction,$stub) = @_;
  
  my $exon       = new Bio::EnsEMBL::Exon;
  
  $exon->seqname   ($exon_prediction->seqname);
  $exon->contig    ($exon_prediction->contig);
  $exon->start     ($exon_prediction->start);
  $exon->end       ($exon_prediction->end  );
  $exon->strand    ($exon_prediction->strand);
  $exon->phase     ($exon_prediction->phase);
  $exon->end_phase ( ( $exon->end - $exon->start + 1 + $exon->phase )%3 );

  return $exon;
}

############################################################
=head2 make_ExonPairs

 Description: Links exons with supporting evidence into ExonPairs

=cut
  
sub  make_ExonPairs {
  my ($self,@exons) = @_;
  
  my $gap = 5;  
  my %pairhash;
  my @forward;
  my @reverse;
  
 EXON: 
  for (my $i = 0; $i < scalar(@exons)-1; $i++) {
    my %idhash;
    my $exon1 = $exons[$i];
    my $jstart = $i - 2;  
    
    if ($jstart < 0) {
      $jstart = 0;
    }
    
    my $jend   = $i + 2;  
    if ( $jend >= scalar(@exons)) {
      $jend    = scalar(@exons) - 1;
    }
    
  J: 
    for (my $j = $jstart ; $j <= $jend; $j++) {
      next J if ($i == $j);
      next J if ($exons[$i]->strand != $exons[$j]->strand);
      next J if ($exons[$i]->{'temporary_id'}  eq $exons[$j]->{'temporary_id'});
      
      my $exon2 = $exons[$j];
      my %doneidhash;
      
      # For the two exons we compare all of their supporting features.
      # If any of the supporting features of the two exons
      # span across an intron a pair is made.
      my @f1 = sort {$b->score <=> $a->score} @{$exon1->get_all_supporting_features};
      
    F1: 
      foreach my $f1 (@f1) {
	next F1 if (!$f1->isa("Bio::EnsEMBL::FeaturePair"));
	my @f = sort {$b->score <=> $a->score} @{$exon2->get_all_supporting_features};
       	
      F2: 
	foreach my $f2 (@f) {
	  next F2 if (!$f2->isa("Bio::EnsEMBL::FeaturePair"));
	  
	  my @pairs = $self->get_all_ExonPairs;		
	  
	  # Do we have hits from the same sequence
	  # n.b. We only allow each database hit to span once
	  # across the intron (%idhash) and once the pair coverage between
	  # the two exons reaches $minimum_coverage we 
	  # stop finding evidence. (%pairhash)
	  
	  if ($f1->hseqname eq $f2->hseqname &&
	      $f1->strand   == $f2->strand   &&
	      !(defined($idhash{$f1->hseqname})) &&
	      !(defined($pairhash{$exon1}{$exon2}))) {
	    
	    my $ispair = 0;
	    my $thresh = $self->threshold;
	    
	    if ($f1->strand == 1) {
	      if (abs($f2->hstart - $f1->hend) < $gap) {
		
		if (!(defined($doneidhash{$f1->hseqname}))) {
		  $ispair = 1;
		}
	      }
	    } 
	    elsif ($f1->strand == -1) {
	      if (abs($f1->hend - $f2->hstart) < $gap) {
		if (!(defined($doneidhash{$f1->hseqname}))) {
		  $ispair = 1;
		}
	      }
	    }
	    
	    # This checks if the coordinates are consistent if the 
	    # exons are on the same contig
	    if ($ispair == 1) {
	      if ($exon1->contig->name eq $exon2->contig->name) {
		if ($f1->strand == 1) {
		  if ($f1->end >  $f2->start) {
		    $ispair = 0;
		  }
		} 
		else {
		  if ($f2->end >  $f1->start) {
		    $ispair = 0;
		  }
		}
	      }
	    }
	    
	    # We finally get to make a pair
	    if ($ispair == 1) {
	      eval {
		my $check = $self->check_link($exon1,$exon2,$f1,$f2);
		next J unless $check;
		
		# we make an exon pair
		my $pair = $self->makePair($exon1,$exon2,"ABUTTING");
		
		if ( $pair) {
		  $idhash    {$f1->hseqname} = 1;
		  $doneidhash{$f1->hseqname} = 1;
		  
		  $pair->add_Evidence($f1);
		  $pair->add_Evidence($f2);
		  
		  if ($pair->is_Covered == 1) {
		    $pairhash{$exon1}{$exon2}  = 1;
		  }
		  next EXON;
		}
		
		
	      };
	      if ($@) {
		warn("Error making ExonPair from [" . $exon1->{'temporary_id'} . "][" .$exon2->{'temporary_id'} ."] $@");
	      }
	    }
	  }
	}
      }
    }
  }
  return $self->get_all_ExonPairs;
}

############################################################

=head2 makePair

  Example    : my $pair = $self->makePair($exon1,$exon2)
  Description: it creates a ExonPair and checks whether this has been included yet or not.
  Returns    : Bio::EnsEMBL::Pipeline::ExonPair
  Args       : Bio::EnsEMBL::Exon,Bio::EnsEMBL::Exon

=cut

sub makePair {
  my ($self,$exon1,$exon2,$type) = @_;
  
  if (!defined($exon1) || !defined($exon2)) {
    $self->throw("Wrong number of arguments [$exon1][$exon2] to makePair");
  }
  
  $self->throw("[$exon1] is not a Bio::EnsEMBL::Exon") unless $exon1->isa("Bio::EnsEMBL::Exon");
  $self->throw("[$exon2] is not a Bio::EnsEMBL::Exon") unless $exon2->isa("Bio::EnsEMBL::Exon");
  
  # create a new pair
  my $tmppair = new Bio::EnsEMBL::Pipeline::ExonPair(-exon1 => $exon1,
						     -exon2 => $exon2,
						     -type  => $type,
						    );
  
  my $found = 0;
  foreach my $pair ($self->get_all_ExonPairs) {
    if ($pair->compare($tmppair) == 1) {
      $pair->add_coverage;
      $tmppair = $pair;
      $found = 1;
    }
  }
  
  if ($found == 0 && $self->check_ExonPair($tmppair)) {
    $self->add_ExonPair($tmppair);
    return $tmppair;
  }
  else{
    return 0;
  }
}

############################################################

############################################################

=head2 prune_features

 Description: prunes out duplicated features
 Returntype : array of Bio::EnsEMBL::SeqFeature
 Args       : array of Bio::EnsEMBL::SeqFeature

=cut

sub prune_features {
  my ($self,@features)  = @_;
    
  my @pruned;

  @features = sort {$a->start <=> $b->start} @features;

  my $prev = -1;

  F: 
  foreach  my $f (@features) {
    if ($prev != -1 && $f->hseqname eq $prev->hseqname &&
	$f->start   == $prev->start &&
	$f->end     == $prev->end   &&
	$f->hstart  == $prev->hstart &&
	$f->hend    == $prev->hend   &&
	$f->strand  == $prev->strand &&
	$f->hstrand == $prev->hstrand) {
    } 
    else {
      push(@pruned,$f);
      $prev = $f;
    }
  }
  return @pruned;
}

############################################################

=head2 check_link

 Example    : $self->check_link($exon1, $exon2, $feature1, $feature2)
 Description: checks to see whether the 2 exons can be linked by the 2 features, i.e.
              if they haven not been linked yet by other feature or independently in another pair.
 Returns    : 1 if exons can be linked, otherwise 0
 Args       : two Bio::EnsEMBL::Exon, two Bio::EnsEMBL::FeaturePair

=cut

sub check_link {
  my ($self,$exon1,$exon2,$f1,$f2) = @_;

  my @pairs = $self->get_all_ExonPairs;
  
  # are these 2 exons already linked in another pair
  foreach my $pair (@pairs) {
    
    if ($exon1->strand == 1) {
      if ($exon1 == $pair->exon1) {
	my @linked_features = @{$pair->get_all_Evidence};
	
	foreach my $f (@linked_features) {
	  
	  if ($f->hseqname eq $f2->hseqname && $f->hstrand == $f2->hstrand) {
	    return 0;
	  }
	}
      }
    } 
    else {
      if ($exon2 == $pair->exon2) {
	my @linked_features = @{$pair->get_all_Evidence};
	
	foreach my $f (@linked_features) {
	  
	  if ($f->hseqname eq $f2->hseqname && $f->hstrand == $f2->hstrand) {
	    return 0;
	  }
	}
      }
    }
    
    # if we're still here, are these 2 exons already part of a pair but linked by different evidence?
    if(($exon1 == $pair->exon1 && $exon2 == $pair->exon2) || 
       ($exon1 == $pair->exon2 && $exon2 == $pair->exon1)    ){
      
      # add in new evidence
      $pair->add_Evidence($f1);
      $pair->add_Evidence($f2);
      
      return 0;
    }
  }
  
  # exons are not linked
  return 1;
}

############################################################

=head2 link_ExonPairs
  
Usage   : my @transcripts = $self->make_ExonPairs(@exons);
Function: It takes the @supported_exons and links them according to the ExonPairs 
          created previously., Links ExonPairs into Transcripts, validates transcripts, 
          rejects any with exons < $GB_MIN_GENSCAN_EXONS
Returns : Array of Bio::EnsEMBL::Pipeline::ExonPair
  Args    : Array of Bio::EnsEMBL::Exon, 

=cut

sub link_ExonPairs {
  my ($self,@exons) = @_;
  

  my @tmpexons;
  
 EXON: 
  foreach my $exon (@exons) {
    $self->throw("[$exon] is not a Bio::EnsEMBL::Exon") unless $exon->isa("Bio::EnsEMBL::Exon");
    
    if ($self->isHead($exon) == 1) {
      
      my $transcript = new Bio::EnsEMBL::Transcript;
      $transcript->type($GB_ABINITIO_SUPPORTED_TYPE);
      
      $self      ->linked_predictions($transcript);
      $transcript->add_Exon       ($exon);
      
      $self->_recurseTranscript($exon,$transcript);
    }
  }
  my $count = 1;
  
  # the recursion is ready, now check the results:
  my @t = $self->linked_predictions;
  
  # create a translation
  foreach my $tran ( @t ) {
    $self->make_Translation($tran);
  }
  
  # flush transcripts & re-add valid ones.
  $self->flush_linked_predictions;  
  
  # validate the transcripts
  foreach my $transcript (@t){
    my @valid = $self->validate_transcript($transcript);
    foreach my $vt(@valid){
      @tmpexons = @{$vt->get_all_Exons};
      if(scalar (@tmpexons) >= $GB_MIN_GENSCAN_EXONS){
	$self->linked_predictions($vt);
      }
    }
  }
  return $self->linked_predictions;
}

############################################################


=head2 _recurseTranscript

 Usage   : $self->_recurseTranscript($exon,$transcript)
 Function: Follows ExonPairs recursively to form a new transcript
 Args    : Bio::EnsEMBL::Exon Bio::EnsEMBL::Transcript

=cut
  
  
sub _recurseTranscript {
  my ($self,$exon,$tran) = @_;
  
  if (defined($exon) && defined($tran)) {
    $self->throw("[$exon] is not a Bio::EnsEMBL::Exon")       unless $exon->isa("Bio::EnsEMBL::Exon");
    $self->throw("[$tran] is not a Bio::EnsEMBL::Transcript") unless $tran->isa("Bio::EnsEMBL::Transcript");
  } 
  else {
    $self->throw("Wrong number of arguments [$exon][$tran] to _recurseTranscript");
  }
  
  # Checks for circular genes here.
  my %exonhash;
  
  foreach my $exon (@{$tran->get_all_Exons}) {
    $exonhash{$exon->{'temporary_id'}}++;
  }
  
  foreach my $exon (keys %exonhash) {
    if ($exonhash{$exon} > 1) {
      $self->warn("Eeeek! Found exon " . $exon . " more than once in the same gene. Bailing out");
      $tran = undef;
      return;
    }
  }
  
  # First copy all the exons into a new transcript
  my $tmptran = new Bio::EnsEMBL::Transcript;
  $tmptran->type($GB_ABINITIO_SUPPORTED_TYPE);
  
  foreach my $ex (@{$tran->get_all_Exons}) {
    $tmptran->add_Exon($ex);
  }
  
  my $count = 0;
  my @pairs = $self->_getPairs($exon);
  
  #    print STDERR "Pairs are @pairs\n";
  
  my @exons = @{$tran->get_all_Exons};
  
  if ($exons[0]->strand == 1) {
    @exons = sort {$a->start <=> $b->start} @exons;
    
  } 
  else {
    @exons = sort {$b->start <=> $a->start} @exons;
  }
  
  
 PAIR: 
  foreach my $pair (@pairs) {
    next PAIR if ($exons[$#exons]->end_phase != $pair->exon2->phase);
    
    if ($count > 0) {
      my $newtran = new Bio::EnsEMBL::Transcript;
      $newtran->type($GB_ABINITIO_SUPPORTED_TYPE);
      $self->add_Transcript($newtran);
      
      foreach my $tmpex (@{$tmptran->get_all_Exons}) {
	$newtran->add_Exon($tmpex);
      }
      
      $newtran->add_Exon($pair->exon2);
      $self->_recurseTranscript($pair->exon2,$newtran);
    } 
    else {
      $tran->add_Exon($pair->exon2);
      $self->_recurseTranscript($pair->exon2,$tran);
    }
    $count++;
  }
}

############################################################



=head2 _getPairs

 Usage   : my @pairs = $self->_getPairs($exon)
 Function: Returns an array of all the ExonPairs 
           in which this exon is exon1
 Returns : @Bio::EnsEMBL::Pipeline::ExonPair
 Args    : Bio::EnsEMBL::Exon

=cut

sub _getPairs {
  my ($self,$exon) = @_;
  
  my $minimum_coverage = 1;
  my @pairs;
  
  $self->throw("No exon input") unless defined($exon);
  $self->throw("Input must be Bio::EnsEMBL::Exon") unless $exon->isa("Bio::EnsEMBL::Exon");
  
  foreach my $pair ($self->get_all_ExonPairs) {
    if (($pair->exon1->{'temporary_id'} eq $exon->{'temporary_id'}) && ($pair->is_Covered == 1)) {
      push(@pairs,$pair);
    }
  } 
  @pairs = sort { $a->exon2->start <=> $b->exon2->start} @pairs;
  return @pairs;
}

############################################################	
	
=head2 isHead

 Usage   : my $foundhead = $self->isHead($exon)
 Function: checks through all ExonPairs to see whether this
           exon is connected to a preceeding exon in a pair, i.e. it returns FALSE if this exon is
           pair->exon2 in a pair. Returns TRUE otherwise.
           pair, where head is the 3prime exon or exon2 in the pair
 Returns : BOOLEAN
 Args    : Bio::EnsEMBL::Exon

=cut
  
sub isHead {
  my ($self,$exon) = @_;
  my $minimum_coverage = 1;
  foreach my $pair ($self->get_all_ExonPairs) {
    my $exon2 = $pair->exon2;
    if (($exon == $exon2 && ($pair->is_Covered == 1) ) ) {
      return 0;
    }
  }
  return 1;
}

############################################################

=head2 isTail

 Title   : isTail
 Usage   : my $foundtail = $self->isTail($exon)
 Function: checks through all ExonPairs to see whether this
           exon is connected to a following exon, i.e. it returns FALSE if this
	   exon is pair->exon1 in a pair. Returns TRUE otherwise.
 Returns : BOOLEAN
 Args    : Bio::EnsEMBL::Exon

=cut

sub isTail {
  my ($self,$exon) = @_;
  
  my $minimum_coverage = 1;
  foreach my $pair ($self->get_all_ExonPairs) {
    my $exon1 = $pair->exon1;
    if ($exon == $exon1 && $pair->is_Covered == 1) {
      return 0;
    }
  }
  return 1;
}

############################################################

=head2 make_Translation

 Function: builds a translation for a Transcript object
 Returns : Bio::EnsEMBL::Translation
 Args    : $transcript - Bio::EnsEMBL::Transcript object
           
=cut

sub make_Translation{
    my ($self,$transcript) = @_;

    my $translation = new Bio::EnsEMBL::Translation;    
    my @exons       = @{$transcript->get_all_Exons};
        
    if ($exons[0]->strand == 1) {
      @exons = sort {$a->start <=> $b->start} @exons;
    } 
    else {
      @exons = sort {$b->start <=> $a->start} @exons;
    }
    
    # we assume here that translation starts in the first exon and ends in the last
    $translation->start_Exon($exons[0]);
    
    # set translation start according to exon_phase
    if( $exons[0]->phase == 0 ) {
      $translation->start(1);
    } 
    elsif ( $exons[0]->phase == 1 ) {
      $translation->start(3);
    } 
    elsif ( $exons[0]->phase == 2 ) {
      $translation->start(2);
    } 
    else {
      $self->throw("Nasty exon phase ".$exons[0]->phase);
    }
    
    $translation->end_Exon  ($exons[$#exons]);
    $translation->end($exons[$#exons]->end - $exons[$#exons]->start + 1 - $exons[$#exons]->end_phase );
    $transcript->translation($translation);
}   


############################################################

=head2 check_ExonPair

  example    : unless ( $self->check_ExonPair($pair) ){ ... };
 Description : builds a translation for a Transcript object
               At the moment it is merely informative and 
               only print outs if there is any problem, and always returns 1.
 Returns     : a BOOLEAN
 Args        : a Bio::EnsEMBL::Pipeline::ExonPair

=cut
  
sub check_ExonPair {
  my ($self,$pair) = @_;
    
  my $exon1  = $pair->exon1;
  my $exon2  = $pair->exon2;
  
  my $trans1 = $pair->exon1->translate();
  my $trans2 = $pair->exon2->translate();
  
  my $splice1;
  my $splice2;
  
  my $spliceseq;
  
  
  # check the splice site sequence
  if ($pair->exon1->strand == 1) {
    $splice1 = $exon1->seq->subseq($exon1->end+1,$exon1->end+2);
    $splice2 = $exon2->seq->subseq($exon2->start-2,$exon2->start-1);
    $spliceseq = new Bio::Seq('-id' => "splice",
			      '-seq' => "$splice1$splice2");
  } 
  else {
    $splice1 = $exon1->seq->subseq($exon1->start-2,$exon1->start-1);
    $splice2 = $exon2->seq->subseq($exon2->end+1,  $exon2->end+2  );
    $spliceseq = new Bio::Seq('-id' => "splice",
			      '-seq' => "$splice2$splice1");
    $spliceseq = $spliceseq->revcom;
  }
  
  $pair->splice_seq($spliceseq);
  
  if ($spliceseq ->seq ne "GTAG"){
    print STDERR "splice sites are not GT-AG for this pair\n";
  }
  
  if ($pair->exon1->end_phase == $pair->exon2->phase){
    print STDERR "exon pair has a phase inconsistency\n";
  }
  return 1;
}

############################################################

sub each_ExonFrame {
    my ($self,$exon) = @_;

    return $self->{'_framehash'}{$exon};
}

############################################################

sub add_ExonPhase {
    my ($self,$exon) = @_;

    if (defined($self->{'_exonphase'}{$exon})) {
#	print STDERR "Already defined phase : old phase " . $self->{'_exonphase'}{$exon} . " new " . $exon->phase . "\n";
	if ($self->{'_exonphase'}{$exon} != $exon->phase) {
	    return 0;
	}
    } else {
	$self->{'_exonphase'}{$exon} = $exon->phase;
	return 1;
    }


}
############################################################


=head2 split_transcript

 Title   : split_transcript 
 Usage   : my @splits = $self->split_transcript($transcript)
 Function: splits a transcript into multiple transcripts at long introns. Rejects single exon 
           transcripts that result. 
 Returns : @Bio::EnsEMBL::Transcript
 Args    : Bio::EnsEMBL::Transcript

=cut


sub split_transcript{
  my ($self, $transcript) = @_;
  $transcript->sort;
  my @split_transcripts   = ();

  if(!($transcript->isa("Bio::EnsEMBL::Transcript"))){
    $self->warn("[$transcript] is not a Bio::EnsEMBL::Transcript - cannot split");
    return (); # empty array
  }
  
  my $prev_exon;
  my $exon_added = 0;
  my $curr_transcript = new Bio::EnsEMBL::Transcript;
  my $translation     = new Bio::EnsEMBL::Translation;
  $curr_transcript->type($transcript->type);
  $curr_transcript->translation($translation);


EXON:   
  foreach my $exon( @{$transcript->get_all_Exons} ){
    
    $exon_added = 0;
      # is this the very first exon?
    if($exon == $transcript->start_Exon){
      $prev_exon = $exon;      
      # set $curr_transcript->translation start and start_exon
      $curr_transcript->add_Exon($exon);
      $exon_added = 1;
      $curr_transcript->translation->start_Exon($exon);
      $curr_transcript->translation->start($transcript->translation->start);
      push(@split_transcripts, $curr_transcript);
      next EXON;
    }
    
    if ($exon->strand != $prev_exon->strand){
      return (); # empty array
    }

    # We need to start a new transcript if the intron size between $exon and $prev_exon is too large
    my $intron = 0;
    if ($exon->strand == 1) {
      $intron = abs($exon->start - $prev_exon->end + 1);
    } else {
      $intron = abs($exon->end   - $prev_exon->start + 1);
    }
    
    if ($intron > $GB_GENSCAN_MAX_INTRON) {
      $curr_transcript->translation->end_Exon($prev_exon);
      # need to account for end_phase of $prev_exon when setting translation->end
      $curr_transcript->translation->end($prev_exon->end - $prev_exon->start + 1 - $prev_exon->end_phase);
      
      # start a new transcript 
      my $t  = new Bio::EnsEMBL::Transcript;
      my $tr = new Bio::EnsEMBL::Translation;
      $t->type($transcript->type);
      $t->translation($tr);

      # add exon unless already added, and set translation start and start_exon
      $t->add_Exon($exon) unless $exon_added;
      $exon_added = 1;

      $t->translation->start_Exon($exon);

      if ($exon->phase == 0) {
	$t->translation->start(1);
      } elsif ($exon->phase == 1) {
	$t->translation->start(3);
      } elsif ($exon->phase == 2) {
	$t->translation->start(2);
      }

      # start exon always has phase 0
      $exon->phase(0);      

      # this new transcript becomes the current transcript
      $curr_transcript = $t;

      push(@split_transcripts, $curr_transcript);
    }

    if($exon == $transcript->end_Exon){
      # add it unless already added
      $curr_transcript->add_Exon($exon) unless $exon_added;
      $exon_added = 1;

      # set $curr_transcript end_exon and end
      $curr_transcript->translation->end_Exon($exon);
      $curr_transcript->translation->end($transcript->translation->end);
    }

    else{
      # just add the exon
      $curr_transcript->add_Exon($exon) unless $exon_added;
    }
    
    # this exon becomes $prev_exon for the next one
    $prev_exon = $exon;

  }

  # discard any single exon transcripts
  my @t = ();
  my $count = 1;
  
  foreach my $st(@split_transcripts){
    $st->sort;
    my @ex = @{$st->get_all_Exons};
    if(scalar(@ex) > 1){
      $st->{'temporary_id'} = $transcript->dbID . "." . $count;
      $count++;
      push(@t, $st);
    }
  }

  return @t;

}

############################################################

=head2 validate_transcript

 Title   : validate_transcript 
 Usage   : my @valid = $self->validate_transcript($transcript)
 Function: Validates a transcript - rejects if mixed strands, splits if long introns
 Returns : @Bio::EnsEMBL::Transcript
 Args    : Bio::EnsEMBL::Transcript

=cut

sub validate_transcript{
  my ($self, $transcript) = @_;
  my @valid_transcripts;

  my $valid = 1;
  my $split = 0;

  # check exon phases:
  my @exons = @{$transcript->get_all_Exons};
  $transcript->sort;
  for (my $i=0; $i<(scalar(@exons-1)); $i++){
    my $end_phase = $exons[$i]->end_phase;
    my $phase    = $exons[$i+1]->phase;
    if ( $phase != $end_phase ){
      $self->warn("rejecting transcript with inconsistent phases( $phase <-> $end_phase) ");
      return undef;
    }
  }
  

  my $previous_exon;
  foreach my $exon (@exons){
    if (defined($previous_exon)) {
      my $intron;
      
      if ($exon->strand == 1) {
	$intron = abs($exon->start - $previous_exon->end + 1);
      } else {
	$intron = abs($exon->end   - $previous_exon->start + 1);
      }
      
      if ($intron > $GB_GENSCAN_MAX_INTRON) {
	print STDERR "Intron too long $intron  for transcript " . $transcript->{'temporary_id'} . "\n";
	$split = 1;
	$valid = 0;
      }
      
      if ($exon->strand != $previous_exon->strand) {
	print STDERR "Mixed strands for gene " . $transcript->{'temporary_id'} . "\n";
	$valid = 0;
	return;
      }
    }
    $previous_exon = $exon;
  }
  
       if ($valid) {
	 push(@valid_transcripts,$transcript);
       }
       elsif ($split){
	 # split the transcript up.
	 my @split_transcripts = $self->split_transcript($transcript);
    push(@valid_transcripts, @split_transcripts);
  }
       return @valid_transcripts;
}

############################################################


sub _check_Transcript{
  my ($self,$transcript) = @_;
  my $slice = $self->slice;
  
  my $id = $self->transcript_id( $transcript );
  
  my $valid = 1;

  # check that transcripts are not completely outside the slice
  if ( $transcript->start > $slice->length || $transcript->end < 1 ){
    print STDERR "transcript $id outside the slice\n";
    $valid = 0;
  }
  # allow transcripts that fall partially off the slice only at one end, the 'higher' end of the slice
  elsif ( $transcript->start < 1 && $transcript->end > 1 ){
      print STDERR "transcript $id falls off the slice by its lower end\n";
    $valid = 0;
  }
  
  # sort the exons 
  $transcript->sort;
  my @exons = @{$transcript->get_all_Exons};
  
  if ($#exons > 0) {
    for (my $i = 1; $i <= $#exons; $i++) {
      
      # check phase consistency:
      if ( $exons[$i-1]->end_phase != $exons[$i]->phase  ){
	print STDERR "transcript $id has phase inconsistency\n";
	$valid = 0;
	last;
      }
      
      # check for folded transcripts
      if ($exons[0]->strand == 1) {
	if ($exons[$i]->start < $exons[$i-1]->end) {
	  print STDERR "transcript $id folds back on itself\n";
	  $valid = 0;
	} 
      } 
      elsif ($exons[0]->strand == -1) {
	if ($exons[$i]->end > $exons[$i-1]->start) {
	  print STDERR "transcript $id folds back on itself\n";
	  $valid = 0;
	} 
      }
    }
  }
  if ($valid == 0 ){
    $self->_print_Transcript($transcript);
  }
  return $valid;
}


############################################################

sub _check_Translation{
  my ($self,$transcript) = @_;
  
  my $id = $self->transcript_id( $transcript );
  
  my $valid = 1;
  
  # check that they have a translation
  my $translation = $transcript->translation;
  my $sequence;
  eval{
    $sequence = $transcript->translate;
  };
  unless ( $sequence ){
    print STDERR "transcript $id has no translation\n";
    return 0;
  }
  if ( $sequence ){
    my $peptide = $sequence->seq;
    if ( $peptide =~ /\*/ ){
      print STDERR "translation of transcript $id has STOP codons\n";
      $valid = 0;
    }
  }
  if ($valid == 0 ){
    $self->_print_Transcript($transcript);
  }
  return $valid;
}

############################################################
#
# GETSET METHODS
#
############################################################

=head2 linked_predictions

 Description: get/set for the transcripts built from linking the exon-pairs
              taken from the prediction transcripts according to consecutive
              feature overlap
=cut

sub linked_predictions {
  my ($self,@linked_predictions) = @_;

  if(!$self->{_linked_predictions}){
    $self->{_linked_predictions} = [];
  }
  if ( @linked_predictions ) {
     push(@{$self->{_linked_predictions}},@linked_predictions);
  }
  return @{$self->{_linked_predictions}};
}

############################################################

=head2 flush_linked_predictions

=cut

sub flush_linked_predictions{
  my ($self) = @_;
  $self->{_linked_predictions} = [];
  return;
}

###################################################

=head2 predictions

 Description: get/set for the PredictionTranscripts. It is  set in new()

=cut

sub predictions {
  my ($self,@predictions) = @_;

  if ( @predictions ) {
     push(@{$self->{_predictions}},@predictions);
  }
  return @{$self->{_predictions}};
}

############################################################

=head2 predictions

 Description: get/set for the similarity features

=cut

sub features {
  my ($self,@features) = @_;
  
  if (!defined($self->{_feature})) {
    $self->{_feature} = [];
  }
  if ( @features ) {
    push(@{$self->{_feature}},@features);
  }
  return @{$self->{_feature}};
}

############################################################

sub annotations{
  my ($self,@annotations) = @_;
  if (@annotations){
    push( @{$self->{_annotations}}, @annotations);
  }
  return @{$self->{_annotations}};
}


############################################################

sub supported_exons {
    my ($self,@exons) = @_;

    if (!defined($self->{'_supported_exons'})) {
	$self->{'_supported_exons'} = [];
    }
    if (@exons ) {
	push(@{$self->{'_supported_exons'}},@exons);
    }
    return @{$self->{'_supported_exons'}};
}

############################################################

sub get_all_ExonPairs {
    my ($self) = @_;

    if (!defined($self->{'_exon_pairs'})) {
	$self->{'_exon_pairs'} = [];
    }
    return @{$self->{'_exon_pairs'}};
}

############################################################

sub add_ExonPair {
    my ($self,$arg) = @_;


    if (!defined($self->{'_exon_pairs'})) {
	$self->{'_exon_pairs'} = [];
    }

    if (defined($arg) && $arg->isa("Bio::EnsEMBL::Pipeline::ExonPair")) {
	push(@{$self->{'_exon_pairs'}},$arg);
#        print STDERR "Adding exon pair $arg\n";
    } else {
	$self->throw("[$arg] is not a Bio::EnsEMBL::Pipeline::ExonPair");
    }
}

#############################################################
sub threshold {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{'_threshold'} = $arg;
    }

    return $self->{'_threshold'} || 100;
}


############################################################


sub _print_Transcript{
  my ($self,$transcript) = @_;
  my @exons = @{$transcript->get_all_Exons};
  my $id;
  if ( $transcript->dbID ){
    $id = $transcript->dbID;
  }
  else{
    $id = "no id";
  }
  print STDERR "transcript id: ".$id."\n";
  foreach my $exon ( @exons){
    $self->print_Exon($exon);
  }
  #print STDERR "Translation : ".$transcript->translation."\n";
  if ( $transcript->translation ){
    print STDERR "translation start exon: ".
      $transcript->translation->start_Exon->start."-".$transcript->translation->start_Exon->end.
	" start: ".$transcript->translation->start."\n";
    print STDERR "translation end exon: ".
      $transcript->translation->end_Exon->start."-".$transcript->translation->end_Exon->end.
	  " end: ".$transcript->translation->end."\n";
    print STDERR "\n";
  }
}

1;
