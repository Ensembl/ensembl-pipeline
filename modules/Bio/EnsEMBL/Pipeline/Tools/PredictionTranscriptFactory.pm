# Ensembl module for Bio::EnsEMBL::PredictionTranscriptFactory
#
# Cared for by EnsEMBL (www.ensembl.org)
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code



=head1 NAME

Bio::EnsEMBL::Pipeline::Tools::PredictionTranscriptFactory - Module having the fset2transcript*
subroutines

=head1 SYNOPSIS

    use Bio::EnsEMBL::Pipeline::Tools::PredictionTranscriptFactory;

    &Bio::EnsEMBL::Pipeline::Tools::PredictionTranscriptFactory::fset2transcript($fset_id);

=head1 DESCRIPTION

Module containing the subroutines fset2transcript*, 
which create transcripts from features (formally housed in
Bio::EnsEMBL::DBSQL::Utils).

=head1 CONTACT

Ensembl - ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut



package Bio::EnsEMBL::Pipeline::Tools::PredictionTranscriptFactory;



use strict;
use Bio::EnsEMBL::PredictionTranscript;
use Bio::EnsEMBL::PredictionExon;


sub fset2transcript {
    my ($genscan,$contig)=@_;

  
    unless ($genscan->isa ("Bio::EnsEMBL::SeqFeatureI"))
    {print "$genscan must be Bio::EnsEMBL::SeqFeatureI\n";}
     
    my $transcript = new Bio::EnsEMBL::PredictionTranscript; 
    my @exons;
    my $count= 1;
    
    foreach my $f ($genscan->sub_SeqFeature) {
  
	my $exon  = new Bio::EnsEMBL::Exon;
	$transcript->add_Exon($exon);
        $exon->slice   ($contig);
	$exon->start    ($f->start);
	$exon->end      ($f->end  );
	$exon->strand   ($f->strand);
	$exon->phase    ($f->phase);
	$exon->end_phase();
	$exon->score($f->score);
	$exon->p_value($f->p_value);
	$exon->slice($contig->primary_seq);
	
	push(@exons,$exon);
	$count++;
	
    }
    
    if( $count == 1 ) {
	$genscan->throw("Got a 0 exon genscan");
    }

    my $translation = new Bio::EnsEMBL::Translation;
    #
    # This code got changed due to Translation convention changing. Should work...
    #
    
    if ($exons[0]->strand == 1) {
	@exons = sort {$a->start <=> $b->start} @exons;
    } else {
	@exons = sort {$b->start <=> $a->start} @exons;
    }
    
    $translation->start(1);
    $translation->end($exons[scalar(@exons)-1]->length);
    
    $translation->start_Exon($exons[0]);
    $translation->end_Exon($exons[$#exons]);

    my $endphase = $exons[0]->end_phase;
    
    foreach my $exon (@exons) {
	
      if ( $exon == $exons[0] ){
	next;
      }
      $exon->phase($endphase);
      $endphase = $exon->end_phase;
    }
    
    $transcript->translation($translation);
    
    return $transcript;
}

sub fset2transcript_guess_phases {
    my ($fset,$contig) = @_;

    my $transcript = new Bio::EnsEMBL::PredictionTranscript;


    my @exons;
    my $count    = 1;
    my $endphase = 0;
    foreach my $f ($fset->sub_SeqFeature) {

      my $exon  = new Bio::EnsEMBL::PredictionExon;
      $exon->slice   ($contig);
      $exon->start    ($f->start);
      $exon->end      ($f->end  );
      $exon->strand   ($f->strand);
      $exon->score($f->score);
      $exon->p_value($f->p_value);
      $exon->slice($contig);
      $exon->phase($f->phase); 
      push(@exons,$exon);
      $count++;
      $exon->phase   ($endphase);
      $transcript->add_Exon($exon);
      $endphase = $exon->end_phase();
    }
	
    my $translation = new Bio::EnsEMBL::Translation;
	
    if ($exons[0]->strand == 1) {
      @exons = sort {$a->start <=> $b->start} @exons;
    } else {
      @exons = sort {$b->start <=> $a->start} @exons;
    }
	
    $translation->start        (1);
    $translation->end          ($exons[$#exons]->end - $exons[$#exons]->start + 1);
    $translation->start_Exon($exons[0]);
    $translation->end_Exon($exons[$#exons]);
    
    
    
    
    foreach my $exon (@exons) {
	
      
	
    }


    if ($transcript->translate->seq !~ /\*/) {
	return $transcript;
    }  	

    $endphase = 1;
    
    foreach my $exon (@exons) {
	$exon->phase($endphase);
	$endphase = $exon->end_phase();
    }

    if ($transcript->translate->seq !~ /\*/) {
	return $transcript;
    }  	

    $endphase = 2;
    
    foreach my $exon (@exons) {
	$exon->phase($endphase);
	$endphase = $exon->end_phase();
    }
    
    if ($transcript->translate->seq !~ /\*/) {
	return $transcript;
    }  	
}

sub fset2transcript_3frame {
  my ($fset,$contig) = @_;
  my @f = $fset->sub_SeqFeature;
  
  if ($f[0]->strand == 1) {
    @f = sort {$a->start <=> $b->start} @f;
  } else {
    @f = sort {$b->start <=> $a->start} @f;
  }

  my @transcripts;

  my $startphase = 0;

  while ($startphase < 3) {
    my $endphase = $startphase;

    my $transcript = new Bio::EnsEMBL::PredictionTranscript;

    push(@transcripts,$transcript);
    my $count    = 1;
    my @exons;

    #print STDERR "Dealing with ".@f." exons\n";
    foreach my $f (@f) {
     
      my $exon  = new Bio::EnsEMBL::PredictionExon;
     
      push(@exons,$exon);
      $exon->seqname($f->seqname);
      $exon->slice   ($contig);
      $exon->start    ($f->start);
      $exon->end      ($f->end  );
      $exon->strand   ($f->strand);
      $exon->slice($contig);
      $exon->phase    ($endphase);
      $exon->end_phase();
      $exon->score    ($f->score);
      $exon->p_value  ($f->p_value);
      $endphase = $exon->end_phase;

      $transcript->add_Exon($exon);
      $count++;

      #print STDERR "exon ".$exon->start." ".$exon->end ." phase ".
      #  $exon->phase."\n";
      #print STDERR "Exon length ".($exon->end - $exon->start +1)."\n";
      #Bio::EnsEMBL::Pipeline::Tools::PredictionTranscriptFactory::display_exon($exon);
    }
       
    my $translation = new Bio::EnsEMBL::Translation;

    my $contig_id = "";
    my $fset_id   = "";

    if ($contig->name) {
       $contig_id = $contig->name;
    }
    if (defined($fset->id)) {
       $fset_id = $fset->id;
    }

    $translation->start        (1);
    $translation->end          ($exons[$#exons]->end - $exons[$#exons]->start + 1);
    $translation->start_Exon($exons[0]);
    $translation->end_Exon  ($exons[$#exons]);
    $transcript->translation($translation);

 #  print STDERR "Phase $startphase " . $transcript->translate->seq . "\n";

    $startphase++;
  }
  #print "finshed  fset2transcript_3frame\n";
  return @transcripts;
}


sub fset2transcript_with_seq {
    my ($genscan,$seq)=@_;

  
    unless ($genscan->isa ("Bio::EnsEMBL::SeqFeatureI"))
    {print "$genscan must be Bio::EnsEMBL::SeqFeatureI\n";}
    unless ($seq->isa ("Bio::PrimarySeqI") || $seq->isa ("Bio::SeqI"))
    {print "$seq must be Bio::SeqI or a Bio::PrimarySeqI\n";}

    #print STDERR "running fset2transcript\n";
    my $transcript = new Bio::EnsEMBL::PredictionTranscript;
    
        
    my @exons;
    my $count= 1;
    
    foreach my $f ($genscan->sub_SeqFeature) {
  
	my $exon  = new Bio::EnsEMBL::PredictionExon;
  $exon->slice   ($seq);
	$exon->start    ($f->start);
	$exon->end      ($f->end  );
	$exon->strand   ($f->strand);
	$exon->phase    ($f->phase);
	$exon->end_phase();
	$exon->score ($f->score);
	#print STDERR "contig is a = ".$seq."\n";
	$exon->slice($seq);
	
	push(@exons,$exon);
	$count++;
	
    }

    foreach my $exon (@exons) {
       	
      $transcript->add_Exon($exon);
	
	
    }
    return $transcript;
   
}


sub display_exon{
  my ($e) = @_;

  #print "exon ".$e->start." ".$e->end." ".$e->strand."\n";
  my $seq = $e->seq;
  my $pep0 = $seq->translate('*', 'X', 0);
  my $pep1 = $seq->translate('*', 'X', 1);
  my $pep2 = $seq->translate('*', 'X', 2);
  #print "exon sequence :\n".$e->seq->seq."\n\n";
  #print $e->seqname." ".$e->start." : ".$e->end." translation in 0 frame\n ".$pep0->seq."\n\n";
  #print $e->seqname." ".$e->start." : ".$e->end." translation in 1 phase\n ".$pep2->seq."\n\n";
  #print $e->seqname." ".$e->start." : ".$e->end." translation in 2 phase\n ".$pep1->seq."\n\n";
  #print "\n\n";
  
}



1;
