#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils - 

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 CONTACT

ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Pipeline::Runnable::Protein::Seg;
use Bio::EnsEMBL::DnaPepAlignFeature;

@ISA = qw(Bio::EnsEMBL::Root);

 



###########################################################
#
# METHODS DOING CHECKS
#
###########################################################

# parameter slice is optional. It makes sense to use it when working on fixed length slices.
# If it is not used, the method can still be used to check consistency of the transcript
# although always on chromosomal/slice coordinates, never in rawcontig coordinates.

sub _check_Transcript{
    my ($self,$transcript, $slice) = @_;
    
    # hardcoded stuff, to go in a config file
    my $MAX_EXON_LENGTH   = 20000;
    my $UNWANTED_EVIDENCE = "NG_";
    my $MAX_INTRON_LENGTH = 2000000;
    
    my $id = $self->transcript_id( $transcript );
    my $valid = 1;
    
    my $strand =  $transcript->start_Exon->strand;

    $transcript->sort;
    
    ############################################################
    # check that transcripts are not completely outside the slice
    # allow transcripts that fall partially off the slice only at 
    # one end, the 'higher' end of the slice
    ############################################################
    if ( $slice ){
	if ( $transcript->start > $slice->length || $transcript->end < 1 ){
	    print STDERR "transcript $id outside the slice\n";
	    $valid = 0;
	}
	elsif ( $transcript->start < 1 && $transcript->end > 1 ){
	    print STDERR "transcript $id falls off the slice by its lower end\n";
	    $valid = 0;
	}
    }
    

    my @exons = @{$transcript->get_all_Exons};
    
    if (scalar(@exons) > 1 ) {

      EXON:
	for (my $i = 0; $i <= $#exons; $i++) {

	    ##############################
	    # check exon length
	    ##############################
	    my $length = $exons[i]->end - $exons[i] + 1;
	    if ( $length > $MAX_EXON_LENGTH ){
		print STDERR "exon too long: length = $length >  MAX_EXON_ENGTH = $MAX_EXON_LENGTH\n";
		$valid = 0;
		last EXON;
	    }
	    
	    
	    if ( $i>0 ){
		##############################
		# check phase consistency:
		##############################
		if ( $exons[$i-1]->end_phase != $exons[$i]->phase  ){
		    print STDERR "transcript $id has phase inconsistency\n";
		    $valid = 0;
		    last EXON;
		}
	    

	    
		##############################
		# check intron length
		##############################
		if ( $strand == 1 ){
		    my $intron_length = $exons[i]->start - $exons[i-1]->end -1;
		    if ( $intron_length > $MAX_INTRON_LENGTH ){
			print STDERR "intron too long: length = $intron_length >  MAX_INTRON_ENGTH = $MAX_INTRON_LENGTH\n";
			$valid = 0;
			last EXON;
		    }
		}
		elsif( $strand = -1 ){
		    my $intron_length = $exons[i-1]->start - $exons[i]->end -1;
		    if ( $intron_length > $MAX_INTRON_LENGTH ){
			print STDERR "intron too long: length = $intron_length >  MAX_INTRON_ENGTH = $MAX_INTRON_LENGTH\n";
			$valid = 0;
			last EXON;
		    }
		}
		
		##############################
		# check for folded transcripts
		##############################
		if ($exons[0]->strand == 1) {
		    if ($exons[$i]->start < $exons[$i-1]->end) {
			print STDERR "transcript $id folds back on itself\n";
			$valid = 0;
			last EXON;
		    } 
		} 
		elsif ($exons[0]->strand == -1) {
		    if ($exons[$i]->end > $exons[$i-1]->start) {
			print STDERR "transcript $id folds back on itself\n";
			$valid = 0;
			last EXON;
		    } 
		}
	    }
	    ############################################################
	    # we don't want the NG_ entries going through, they are evil
	    ############################################################
	    foreach my $evidence (@{$exons[i]->get_all_supporting_features}){
		if ( $evidence->hseqname=~/$UNWANTED_EVIDENCE/ ){
		    print STDERR "transcript with evil evidence: ".$evidence->hseqname." skippping\n";
		    $valid = 0;
		    next EXON;
		}
	    }
	}
	
    }
    elsif( scalar(@exons) == 1 ){
	my $length =  $exons[0]->end - $exons[0]->start + 1;
	if ( $length >  $MAX_EXON_LENGTH ){
	    print STDERR "single exon transcript is too long: length = $length >  MAX_EXON_LENGTH = $MAX_EXON_LENGTH\n";
	    $valid = 0;
	}
    }
    else{
	print STDERR "transcript with no exons\n";
	$valid = 0;
    }
    if ($valid == 0 ){
	$self->_print_Transcript($transcript);
    }
    return $valid;
}

############################################################

=head2 _check_Translation

Description : it returns TRUE if a transcript has a translation, and this has 
              no  stop codons. It returns FALSE otherwise. 
  IMPORTANT : we want to check translation independently from other
              properties of the transcripts. Basically because
              we may have some transcripts whih are valid but
              for which we haven not assigned a translation yet.
ReturnType  : a BOOLEAN.
=cut

sub _check_Translation{
  my ($self,$transcript) = @_;
  
  my $id = $self->transcript_id( $transcript );
  
  my $valid = 1;
  
  my $translation = $transcript->translation;

  # double check sane translation start & end
  if( $translation->start < 1){
      print STDERR "dodgy translation start: ".$translation->start."\n";
      $valid = 0;
  }
  
  if( $translation->end < 1 || $translation->end > $translation->end_Exon->length ){
     print STDERR "dodgy translation end " . $translation->end . "\n";
     $valid = 0;
  }
    
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
    #print STDERR "peptide: $peptide\n";
    # check only terminal stops
    if ( $peptide =~ /\*./ ){
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

sub transcript_id {
  my ( $self, $t ) = @_;
  my $id;
  if ( $t->stable_id ){
    $id = $t->stable_id;
  }
  elsif( $t->dbID ){
    $id = $t->dbID;
  }
  elsif( $t->temporary_id ){
    $id = $t->temporary_id;
  }
  else{
    $id = 'no-id';
  }
  
  if ($t->type){
      $id .= " ".$t->type."\n";
  }
  return $id;
}

############################################################

=head2 split_transcript

  Function: splits a transcript into multiple transcripts at long introns. Rejects single exon 
    transcripts that result. 
    Returns : Ref to @Bio::EnsEMBL::Transcript
    Args    : Bio::EnsEMBL::Transcript

=cut
    
sub split_Transcript{
  my ($transcript,$max_intron) = @_;
      
  $transcript->sort;
      
      my @split_transcripts   = ();
      
      if(!($transcript->isa("Bio::EnsEMBL::Transcript"))){
	#$self->throw("[$transcript] is not a Bio::EnsEMBL::Transcript - cannot split");
      }
      
      my $prev_exon;
      my $exon_added = 0;
      
      my $curr_transcript = new Bio::EnsEMBL::Transcript;
      my $translation     = new Bio::EnsEMBL::Translation;
      
      $curr_transcript->translation($translation);
      
    EXON: foreach my $exon (@{$transcript->get_all_Exons}){
	
	$exon_added = 0;
	
	# Start a new transcript if we are just starting out
	
	if($exon == $transcript->start_Exon){
	  
	  $prev_exon = $exon;
	  
	  $curr_transcript->add_Exon($exon);
	  $exon_added = 1;
	  $curr_transcript->translation->start_Exon($exon);
	  $curr_transcript->translation->start($transcript->translation->start);
	  
	  push(@split_transcripts, $curr_transcript);
	  next EXON;
	}
	
	# We need to start a new transcript if the intron size between $exon and $prev_exon is too large
	my $intron = 0;
	
	if ($exon->strand == 1) {
	  $intron = abs($exon->start - $prev_exon->end - 1);
	} else {
	  $intron = abs($prev_exon->start - $exon->end - 1);
	}
	
	if ($intron > $max_intron) {
	  $curr_transcript->translation->end_Exon($prev_exon);
	  $curr_transcript->translation->end($prev_exon->end - $prev_exon->start + 1 - $prev_exon->end_phase);
	  
	  my $t  = new Bio::EnsEMBL::Transcript;
	  my $tr = new Bio::EnsEMBL::Translation;
	  
	  $t->translation($tr);
	  
	  # add exon unless already added, and set translation start and start_Exon
	  # But the exon will nev er have been added ?
	  
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
	  
	  $exon->phase(0);
	  
	  $curr_transcript = $t;
	  
	  push(@split_transcripts, $curr_transcript);
	}
	
	if ($exon == $transcript->end_Exon){
	  $curr_transcript->add_Exon($exon) unless $exon_added;
	  $exon_added = 1;
	  
	  $curr_transcript->translation->end_Exon($exon);
	  $curr_transcript->translation->end($transcript->translation->end);
	} else {
	  $curr_transcript->add_Exon($exon) unless $exon_added;
	}
	
	foreach my $sf(@{$exon->get_all_supporting_features}){
	  $sf->seqname($exon->contig_id);
	}
	
	$prev_exon = $exon;
	
      }
      
      # discard any single exon transcripts
      my @final_transcripts = ();
      my $count = 1;
      
      foreach my $st (@split_transcripts){
	$st->sort;
	
	my @ex = @{$st->get_all_Exons};

	if(scalar(@ex) > 1){
	  $st->{'temporary_id'} = $transcript->dbID . "." . $count;
	  $count++;
	  push(@final_transcripts, $st);

	}
  }

 return \@final_transcripts;
      
}
############################################################
#
# METHODS DOING THE PRINTING
#
############################################################

sub _print_SimpleTranscript{
    my ($self,$transcript) = @_;
    my @exons = @{$transcript->get_all_Exons};
    my $id;
    if ($transcript->stable_id){
	$id = $transcript->stable_id;
    }
    elsif ( $transcript->dbID ){
	$id = $transcript->dbID;
    }
    else{
	$id = "no id";
    }
    if ( defined( $transcript->type ) ){
	$id .= " ".$transcript->type;
    }
    print STDERR "transcript: ".$id.": ";
    foreach my $exon ( @exons){
	print STDERR $exon->start."-".$exon->end." ";
    }
    print STDERR "\n";
}

############################################################

sub _print_Transcript{
  my ($self,$transcript) = @_;
  my @exons = @{$transcript->get_all_Exons};
  my $id;
  if ($transcript->stable_id){
    $id = $transcript->stable_id;
  }
  elsif ( $transcript->dbID ){
    $id = $transcript->dbID;
  }
  else{
    $id = "no id";
  }
  if ( defined( $transcript->type ) ){
    $id .= " ".$transcript->type;
  }
  print STDERR "transcript: ".$id."\n";
  foreach my $exon ( @exons){
    print STDERR $exon->gffstring."\n";
  }
  if ( $transcript->can('translation') && $transcript->translation){
    $self->_print_Translation($transcript->translation);
  }
}

############################################################

sub _print_Translation{
  my ($self,$translation) = @_;
  
  
  if ( $translation->start_Exon ){
      print STDERR "translation start exon: ".
	  $translation->start_Exon->start."-".$translation->start_Exon->end.
	      " start: ".$translation->start."\t phase: ".$translation->start_Exon->phase.
		  " end_phase: ".$translation->start_Exon->end_phase."\n";
  }
  else{
      print STDERR "translation->start_Exon does not exist\n";
  }

  if ( $translation->end_Exon ){
      print STDERR "translation end exon: ".
	  $translation->end_Exon->start."-".$translation->end_Exon->end.
	      " end: ".$translation->end."\t phase: ".$translation->end_Exon->phase.
		  " end_phase: ".$translation->end_Exon->end_phase."\n";
  }
  else{
      print STDERR "translation->end_Exon does not exist\n";
  }
  
}

############################################################

sub _print_Evidence{
  my ($self,$transcript) = @_;
  my @exons = @{$transcript->get_all_Exons};
  my $id;
  if ($transcript->stable_id){
    $id = $transcript->stable_id;
  }
  elsif ( $transcript->dbID ){
    $id = $transcript->dbID;
  }
  else{
    $id = "no id";
  }
  if ( defined( $transcript->type ) ){
    $id .= " ".$transcript->type;
  }
  print STDERR "transcript: ".$id."\n";
  my $count = 0;
  foreach my $exon ( @exons){
    $count++;
    my $exon_id;
    if ($exon->stable_id){
      $exon_id = $exon->stable_id;
    }
    elsif ( $exon->dbID ){
      $exon_id = $exon->dbID;
    }
    else{
      $exon_id = "no id";
    }
    print STDERR "Exon $exon_id: ".$exon->gffstring."\n";
    my @evidence = @{$exon->get_all_supporting_features};
    if (@evidence){
      foreach my $evi ( @evidence ){
	print STDERR "Evidence: ".$evi->gffstring."\n";
      }
    }	
    else{
      print STDERR "no evidence for exon ".$count."\n";
    }
  }
}

############################################################

sub _print_TranscriptEvidence{
    my ($self,$transcript) = @_;
    my @exons = @{$transcript->get_all_Exons};
    my %evidence;
    my %score;
    my %percent_id;
    foreach my $exon ( @exons){
	my @evidences = @{$exon->get_all_supporting_features};
	if (@evidences){
	    foreach my $evi ( @evidences ){
		$evidence{$evi->hseqname} = 1;
		unless( $score{$evi->hseqname} ){
		    $score{$evi->hseqname} = $evi->score;
		}
		if ( $score{$evi->hseqname} < $evi->score ){
		    $score{$evi->hseqname} = $evi->score;
		}
		unless( $percent_id{$evi->hseqname} ){
		    $percent_id{$evi->hseqname} = $evi->percent_id;
		}
		if ( $percent_id{$evi->hseqname} < $evi->percent_id ){
		    $percent_id{$evi->hseqname} = $evi->percent_id;
		}
	    }
	}
    }
    foreach my $evidence ( keys %evidence ){
	print STDERR "hit_name: ".$evidence." score: ".$score{$evidence}." percent_id: ".$percent_id{$evidence}."\n";
    }
}

############################################################

sub _print_Peptide{
  my ($self, $transcript) = @_;
  
  my $seqout = new Bio::SeqIO->new(-fh => \*STDERR);
  my $translation;
  
  eval {
    $translation = $transcript->translate;
    print "translation is a $translation\n";
  };  
  if ($@) {
    print STDERR "Couldn't translate transcript\n";
  }
  else{
    #unless ( $translation->display_id ){
    #  $translation->display_id($self->transcript_id($transcript));
    #}
    $seqout->write_seq($translation);
  }
}


############################################################

#sub _print_Exon{
#  my ($self,$exon) = @_;

#  print STDERR $exon->contig->
#  if ($exon->isa('Bio::EnsEMBL::Sticky'){

#  }
1;
