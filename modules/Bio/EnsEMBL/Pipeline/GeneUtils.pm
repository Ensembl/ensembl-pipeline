#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::GeneUtils - 

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 CONTACT

ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::Pipeline::GeneUtils;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Pipeline::Runnable::Protein::Seg;
use Bio::EnsEMBL::DnaPepAlignFeature;

@ISA = qw(Bio::EnsEMBL::Root);

=head2 validate_transcript

Title   : validate_transcript 
  Usage   : my @valid = $self->validate_transcript($transcript)
 Function: Validates a transcript - rejects if mixed strands, 
  rejects if low coverage, 
  rejects if stops in translation
  splits if long introns and insufficient coverage of parental protein
  Returns : Ref to @Bio::EnsEMBL::Transcript
  Args    : Bio::EnsEMBL::Transcript

=cut
  
  sub validate_transcript {
    my ($transcript,$coverage,$low_complexity,$maxintron,$seqfetchers ) = @_;
    
    my @valid_transcripts;
    
    my $valid = 1;
    my $split = 0;
    
    Bio::EnsEMBL::Pipeline::GeneUtils::check_coverage      ($transcript,$coverage,$seqfetchers);
    Bio::EnsEMBL::Pipeline::GeneUtils::check_translation   ($transcript);
    Bio::EnsEMBL::Pipeline::GeneUtils::check_low_complexity($transcript,$low_complexity);
    Bio::EnsEMBL::Pipeline::GeneUtils::check_strand       ($transcript);
    
    my @tran = Bio::EnsEMBL::Pipeline::GeneUtils::check_introns ($transcript,$maxintron);
    
    
    if ($valid) {
      # Do we really need to do this? Can we just take out the dbID adaptor stuff
      # make a new transcript that's a copy of all the important parts of the old one
      # but without all the db specific gubbins
      
      my $newtranscript  = new Bio::EnsEMBL::Transcript;
      my $newtranslation = new Bio::EnsEMBL::Translation;
      
      $newtranscript->translation($newtranslation);
      $newtranscript->translation->start_Exon($transcript->translation->start_Exon);
      $newtranscript->translation->end_Exon($transcript->translation->end_Exon);
      $newtranscript->translation->start($transcript->translation->start);
      $newtranscript->translation->end($transcript->translation->end);
      
      foreach my $exon(@{$transcript->get_all_Exons}){
	$newtranscript->add_Exon($exon);
	foreach my $sf(@{$exon->get_all_supporting_features}){
	  $sf->seqname($exon->contig_id);
	}
      }
      
      push(@valid_transcripts,$newtranscript);
    }
    
    return \@valid_transcripts;
  }
  
  sub check_introns {
    my ($transcript,$maxintron) = @_;
    
    my $previous_exon;
    my $split = 0;
    
    foreach my $exon (@{$transcript->get_all_Exons}){
      
      if (defined($previous_exon)) {
	my $intron;
	
	if ($exon->strand == 1) {
	  $intron = abs($exon->start - $previous_exon->end - 1);
	} else {
	  $intron = abs($previous_exon->start - $exon->end - 1);
	}
	
	if ( $intron > $maxintron ) {
	  print STDERR "Intron too long $intron  for transcript " . $transcript->dbID . "\n";
	  $split = 1;
	}
	
      }
      $previous_exon = $exon;
    }
    
    if ($split) {
      return @{Bio::EnsEMBL::Pipeline::GeneUtils::split_Transcript($transcript,$maxintron)};
    } else {
      return ($transcript);
    }
  }

=head2 split_transcript

  Title   : split_transcript 
    Usage   : my @splits = $self->split_transcript($transcript)
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
    
    sub check_strand {
      my ($transcript) = @_;
      
      my $previous_exon;
      
      foreach my $exon (@{$transcript->get_all_Exons}){
	
	if (defined($previous_exon)) {
	  my $intron;
	  
	  if ($exon->strand == 1) {
	    $intron = abs($exon->start - $previous_exon->end - 1);
	  } else {
	    $intron = abs($previous_exon->start - $exon->end - 1);
	  }
	  
	  if ($exon->strand != $previous_exon->strand) {
	    print STDERR "Mixed strands for gene " . $transcript->{'temporary_id'} . "\n";
	    return 0;
	  }
	}
	$previous_exon = $exon;
      }
      return 1;
    }

=head2 check_translation

    Title   : check_translation
      Usage   :
    Function: 
    Example :
      Returns : 1 if transcript translates with no stops, otherwise 0
      Args    :
      

=cut
      
      sub check_translation {
	my ($transcript) = @_;
	
	my $tseq;
	
	eval{
	  $tseq = $transcript->translate;
	};
	
	if((!defined $tseq) || ($@)){
	  my $msg = "problem translating :\n$@\n";
	  #    $self->warn($msg);
	  return 0;
	}
	
	if ($tseq->seq =~ /\*/ ) {
	  #   $self->warn("discarding transcript - translation has stop codons\n");
	  return 0;
	}
	else{
	  return 1;
	}
      }
      

=head2 check_coverage

      Title   : check_coverage
	Usage   :
      Function: returns how much of the parent protein is covered by the genewise prediction
	Example :
	Returns : percentage
	Args    :
	

=cut
	
	sub check_coverage {
	  my ( $transcript, $coverage, $seqfetchers) = @_;
	  
	  my $matches = 0;
	  my $pstart  = 0;
	  my $pend    = 0;
	  my $protname;
	  my $plength;
	  
	  
	  foreach my $exon(@{$transcript->get_all_Exons}) {
	    $pstart = 0;
	    $pend   = 0;
	    
	    foreach my $f(@{$exon->get_all_supporting_features}){
	      
	      if (!defined($protname)){
		$protname = $f->hseqname;
	      }
	      if($protname ne $f->hseqname){
		warn("$protname ne " . $f->hseqname . "\n");
	      }
	      
	      if((!$pstart) || $pstart > $f->hstart){
		$pstart = $f->hstart;
	      }
	      
	      if((!$pend) || $pend < $f->hend){
		$pend= $f->hend;
	      }
	    }
	    $matches += ($pend - $pstart + 1);
	  }
	  
	  my $seq; 
	  
	SEQFETCHER:
	  foreach my $seqfetcher (@$seqfetchers){
	    print STDERR "FPC_BlastMiniGenewise: getting sequence for $protname\n";
	    eval{
	      $seq = $seqfetcher->get_Seq_by_acc($protname);
	    };
	    if ($@) {
	      print("FPC_BMG:Error fetching sequence for [$protname] - trying next seqfetcher:[$@]\n");
	    }
	    
	    if (defined $seq) {
	      last SEQFETCHER;
	    }
	    
	  }
	  
	  if(!defined $seq){
	    print("FPC_BMG:No sequence fetched for [$protname] - can't check coverage, letting gene through\n");
	    return 100;
	  }
	  
	  $plength = $seq->length;
	  
	  if(!defined($plength) || $plength == 0){
	    warn("no sensible length for $protname - can't get coverage\n");
	    return 0;
	  }
	  
	  my $realcoverage = 100 * $matches/$plength;
	  
	  if ($realcoverage < $coverage){
	    print(" rejecting transcript for low coverage: $coverage\n");
	  }
	  
	  return $realcoverage;
	  
	}

=head2 check_low_complexity

	Title   : check_complexity
	  Usage   :
	Function: uses seg to find low complexity regions in transcript->translate. 
	  Calculates overall %low complexity of the translation
	  Example :
	  Returns : percentage low complexity sequence
	  Args    :
	  

=cut
	  
	  sub check_low_complexity{
	    my ($transcript,$complexity) = @_;
	    
	    my $low_complexity;
	    
	    eval{
	      
	      my $protseq = $transcript->translate;
	      
	      # Ugh! 
	      my $analysis = Bio::EnsEMBL::Analysis->new(
							 -db           => 'low_complexity',
							 -program      => '/usr/local/ensembl/bin/seg',
							 -program_file => '/usr/local/ensembl/bin/seg',
							 -gff_source   => 'Seg',
							 -gff_feature  => 'annot',
							 -module       => 'Seg',
							 -logic_name   => 'Seg'
							 
							);
	      
	      my $seg = new  Bio::EnsEMBL::Pipeline::Runnable::Protein::Seg(    
									    -query    => $protseq,
									    -analysis => $analysis,
									   );
	      
	      $seg->run;
	      
	      
	      if($seg->get_low_complexity_length > $complexity){
		print("discarding transcript - translation has $low_complexity% low complexity sequence\n");
		
	      }
	      
	      my $compl = $seg->get_low_complexity_length;
	      print "Comple $compl\n";
	      
	      return $compl;
	      
	      
	    };
	    
	    if($@){
	      print STDERR "problem running seg: \n[$@]\n";
	      return 0;		# let transcript through
	    }
	    
	  }

=head S

	  Title   : make_transcript
	    Usage   : $self->make_transcript($gene, $contig, $genetype, $analysis_obj)
	  Function: makes a Bio::EnsEMBL::Transcript from a SeqFeature representing a gene, 
	    with sub_SeqFeatures representing exons.
	    Example :
	    Returns : Bio::EnsEMBL::Transcript with Bio::EnsEMBL:Exons(with supporting feature 
								       data), and a Bio::EnsEMBL::translation
	    Args    : $gene: Bio::EnsEMBL::SeqFeatureI, $contig: Bio::EnsEMBL::RawContig,
	    $genetype: string, $analysis_obj: Bio::EnsEMBL::Analysis
	    

=cut
	    
	    sub SeqFeature_to_Transcript {
	      my ($gene, $contig, $analysis_obj,$db,$phase) = @_;
	      
	      unless ($gene->isa ("Bio::EnsEMBL::SeqFeatureI")){
		print "$gene must be Bio::EnsEMBL::SeqFeatureI\n";
	      }
	      
	      my $transcript   = new Bio::EnsEMBL::Transcript;
	      my $translation  = new Bio::EnsEMBL::Translation;    
	      
	      $transcript->translation($translation);
	      
	      my $excount = 1;
	      my @exons;
	      
	      my $ea;
	      
	      if ($db ne "" && defined($db)) {
		$ea = $db->get_ExonAdaptor;
	      }
	      
	      my $curr_phase = 0;
	      
	      my @pred = $gene->sub_SeqFeature;
	      
	      if ($pred[0]->strand ==1 ) {
		@pred = sort {$a->start <=> $b->start} @pred;
	      } else {
		@pred = sort {$b->start <=> $a->start} @pred;
	      }
	      foreach my $exon_pred (@pred) {
		my $exon = new Bio::EnsEMBL::Exon;
		
		$exon->start    ($exon_pred->start);
		$exon->end      ($exon_pred->end);
		$exon->strand   ($exon_pred->strand);
		
		if ($phase) {
		  $exon->phase($curr_phase);
		  my $exon_length =  $exon->end - $exon->start + 1;
		  my $end_phase   = ( $exon_length + $exon->phase ) %3;
		  
		  $exon->end_phase($end_phase);
		  
		  $curr_phase = $end_phase;
		} else {
		  $exon->phase    ($exon_pred->phase);
		  $exon->end_phase($exon_pred->end_phase);
		}
		
		
		$exon->contig   ($contig); 
		$exon->adaptor  ($ea);
		
		
		# sort out supporting evidence for this exon prediction
		# Now these should be made in the correct type to start with
		
		my @sf = $exon_pred->sub_SeqFeature;
		my $prot_adp;
		
		if ($db ne "" && defined($db)) {
		  $prot_adp = $db->get_ProteinAlignFeatureAdaptor;    
		}
		
		if(@sf){
		  my $align = new Bio::EnsEMBL::DnaPepAlignFeature(-features => \@sf); 
		  
		  $align->seqname($contig->dbID);
		  $align->contig($contig);
		  $align->adaptor($prot_adp);
		  $align->score(100); # Hmm!!!!
		  $align->analysis($analysis_obj);
		  
		  $exon->add_supporting_features($align);
		  
		}
		
		push(@exons,$exon);
		
		$excount++;
	      }
	      
	      if ($#exons < 0) {
		print STDERR "Odd.  No exons found\n";
		return;
	      } else {
		
		if ($exons[0]->strand == -1) {
		  @exons = sort {$b->start <=> $a->start} @exons;
		} else {
		  @exons = sort {$a->start <=> $b->start} @exons;
		}
		
		foreach my $exon (@exons) {
		  $transcript->add_Exon($exon);
    }
    
    $translation->start_Exon($exons[0]);
    $translation->end_Exon  ($exons[$#exons]);
    
    if ($exons[0]->phase == 0) {
      $translation->start(1);
    } elsif ($exons[0]->phase == 1) {
      $translation->start(3);
    } elsif ($exons[0]->phase == 2) {
      $translation->start(2);
    }
    
    $translation->end  ($exons[$#exons]->end - $exons[$#exons]->start + 1);
  }
  
  return $transcript;

}

1;

