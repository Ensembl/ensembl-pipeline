#
# Written by Eduardo Eyras
#
# Copyright GRL/EBI 2002
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod

=head1 NAME

Bio::EnsEMBL::Pipeline::Tools::PseudoGeneTests;


=head1 SYNOPSIS


=head1 DESCRIPTION

For a given protein id, it finds all the genewise modesl that have been built from that protein id.
It will locate cases where there is at least one spliced-case and one or more unspliced 
(or only with frameshifts) transcripts. The latter cases are potential pseudogenes.
For these potential pseudogeens, we find the syntenic region in Mouse and run genewise on that
sequence. If there is no match or very poor alignment, we have a very likely pseudogene.

=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::Tools::PseudoGeneTests;

use strict;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Bio::EnsEMBL::Pipeline::Config::GeneBuild::PseudoGenes;
use Bio::EnsEMBL::Pipeline::Runnable::Blast;
use Bio::EnsEMBL::Pipeline::Runnable::BlastDB;
use Bio::EnsEMBL::Pipeline::Runnable::NewExonerate;
use Bio::EnsEMBL::Pipeline::Runnable::Genewise;
use Bio::EnsEMBL::Pipeline::GeneComparison::ComparativeTools;

use vars qw(@ISA);

@ISA = qw (Bio::EnsEMBL::Pipeline::RunnableI);

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  
  return $self;
}
############################################################
#
# METHOD WITH THE LOGIC OF THE TESTS
#
############################################################

sub pseudogene_test{
  my ( $self, 
       $transcript, 
       $db,
       $compara_db, 
       $focus_db, $focus_species, 
       $target_db, $target_species,
       $target_db2, $target_species2
       $threshold
     ) = @_;
  
  my ( $frameshift, $polyA, $Met, $spliced_elsewhere, $mouse_homology, $rat_homology, $break_synteny_mouse, $break_synteny_rat, $repeat );
    
    ############################################################
    # is it spliced?
    if ( Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->is_spliced($transcript) ){
      print STDERR "Spliced:\tYES\n";
    }
    else{
      print STDERR "Spliced:\tNO\n";
    }

    ############################################################
    # does it have frameshifts?
    if ( !Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->is_spliced($transcript) 
	 && 
	 scalar( @{$transcript->get_all_Exons} > 1) 
       ){
      print STDERR "Frameshifts:\tYES\n";
      $frameshift = 1;
    }
    else{
      print STDERR "Frameshifts:\tNO\n";
      $frameshift = 0;
    }

    ############################################################
    # does it have an in-frame stop codon?
    
    ############################################################
    # has polyA?
    if ( $self->has_polyA_track( $transcript, $db )){
	print STDERR "polyA track:\tYES\n";
	$polyA = 1;
    }
    else{
	print STDERR "polyA track:\tNO\n";
	$polyA = 0;
    }
    
    ############################################################
    # starts with Methionine?
    if ($self->starts_with_Methionine( $transcript ) ){
	print STDERR "starts with Methionine: YES\n";
	$Met = 1;
    }
  else{
	print STDERR "starts with Methionine: NO\n";
	$Met = 0;
    }
  
  ############################################################
  # is the evidence spliced elsewhere in the genome?
  my $spliced_ones = $self->check_for_spliced_real_gene( $transcript, $focus_db );
  if ( $spliced_ones && @$spliced_ones ){
    print STDERR "Evidence is spliced elsewhere: YES\n";
    $spliced_elsewhere = 1;
  }
  else{
    print STDERR "Evidence is spliced elsewhere: NO\n";
    $spliced_elsewhere = 0;
  }
  
    ############################################################
    # if is spliced, has it got canonical splice sites?
    if ( scalar( @{$transcript->get_all_Exons} ) > 1 ){
	if ( Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils
	     ->check_splice_sites( $transcript ) ){
	    print STDERR "Canonical splice sites:\tYES\n";
	}
	else{
	    print STDERR "Canonical splice sites:\tNO\n";
	}
    }
    
  ############################################################
  # has it got homology in mouse?
  my $orthologues = Bio::EnsEMBL::Pipeline::GeneComparison::ComparativeTools
    ->test_for_orthology($transcript, $db, $focus_db, $focus_species, $compara_db, $target_db, $target_species, $threshold );
  if ( $orthologues && @{$orthologues} ){
    $mouse_homology = 1;
  }
  else{
    print STDERR "Could not find any ortholog\n";
    $mouse_homology = 0;
  }
  
  ############################################################
  # has it got homology in rat?
  my $orthologues2 = Bio::EnsEMBL::Pipeline::GeneComparison::ComparativeTools
    ->test_for_orthology($transcript, $db, $focus_db, $focus_species, $compara_db, $target_db2, $target_species2, $threshold );
  if ( $orthologues2 && @{$orthologues2} ){
    $rat_homology = 1;
  }
  else{
    print STDERR "Could not find any ortholog\n";
    $rat_homology = 0;
  }
  
  ############################################################
  # Does it break synteny in mouse?
  my $synteny_breaking = Bio::EnsEMBL::Pipeline::GeneComparison::ComparativeTools
    ->test_for_synteny_breaking($transcript, $db, $focus_db, $focus_species, $compara_db, $target_db, $target_species );
  if ( $synteny_breaking ){
    print STDERR "synteny in mouse is broken\n";
    $break_synteny_mouse = 1;
  }
  else{
    print STDERR "synteny in mouse is preserved\n";
    $break_synteny_mouse = 0;
  }
  
  ############################################################
    # Does it break synteny in rat?
    my $synteny_breaking2 = Bio::EnsEMBL::Pipeline::GeneComparison::ComparativeTools
      ->test_for_synteny_breaking($transcript, $db, $focus_db, $focus_species, $compara_db, $target_db2, $target_species2 );
  if ( $synteny_breaking2 ){
    print STDERR "synteny in rat is broken\n";
    $break_synteny_rat= 1;
  }
  else{
    print STDERR "synteny in rat is preserved\n";
    $break_synteny_rat = 0;
  }
  
  ############################################################
  # Does it overlap with repeats?
  if ( $self->overlap_repeats( $transcript, $focus_db ) ){
    print STDERR "it overlaps with repeats\n";
    $repeat = 1;
  }
  else{
    print STDERR "it does not overlap with repeats\n";
    $repeat = 0;
  }
  
  return ( $frameshift, $polyA, $Met, $spliced_elsewhere, $mouse_homology, $rat_homology, $break_synteny_mouse, $break_synteny_rat, $repeat );
}

############################################################
# 
# METHODS DOING THE ACTUAL WORK
#
############################################################

sub is_spliced{
  my ($self,$transcript) =@_;
  return Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->is_spliced($transcript);
}

############################################################
# method for checking whether a transcript
# has a polyA attached. This polyA can be in general
# outside the annotated transcript, so we must take extra 
# downstream sequence.

sub has_polyA_track{
  my ($self,$transcript,$db) = @_;
  my @exons = sort { $a->start <=> $b->end } @{$transcript->get_all_Exons};
  my $chr_name  =  $exons[0]->contig->chr_name;
  my $end   =  $exons[0]->contig->chr_start + $exons[$#exons]->start - 1;
  
  # we take 20 bases downstream:
  my $seq = $db->get_SliceAdaptor->fetch_by_chr_start_end( $chr_name, $end - 5, $end + 10 )->seq;
  print STDERR "20bp downstream: $seq\n";
  my $length = length($seq);
  my $a_count = $seq =~ tr/Aa//;
  
  if ( ($a_count/$length) >=7/10 ){
    return 1;
  }
  else{
    return 0;
  }
}

############################################################

sub starts_with_Methionine{
  my ($self,$transcript) = @_;
  
  if ( $transcript->translation ){
    
    my $tseq;
    eval{
      $tseq = $transcript->translate;
    };
    unless($tseq){
      return 0;
    }
    
    my $first_aa = substr ( $tseq->seq, 0, 1 );
    print STDERR "first AA: $first_aa\n";
    if ( $first_aa eq 'M' || $first_aa eq 'm' ){
      return 1;
    }
    else{
      return 0;
    }
  }
  else{
    
    my $trans_seq;
    eval{
      $trans_seq = $transcript->seq;
    };
    unless( $trans_seq ){
      return 0;
    }

    my $first_codon = substr ( $trans_seq->seq, 0, 3 );
    print STDERR "first codon $first_codon\n";
    if ( $first_codon eq 'ATG' || $first_codon eq 'atg' ){
      return 1;
    }
    else{
      return 0;
    }
  }
}

############################################################
# this method checks whether the pseudogene overlaps
# with repeats. You need to pass a slice where the pseudogene is
# and a db where to get the repeats from (they should be both on the same assembly)

sub overlap_repeats{
  my ($self,$pseudogene, $db) = @_;

  ############################################################
  # pseudogene is a transcript
  my @exons  = sort { $a->start <=> $b->start} @{$pseudogene->get_all_Exons};
  my $start  = $exons[0]->start;
  my $end    = $exons[$#exons]->end;
  my $strand = $exons[0]->strand;
  my $length = $end - $start + 1;

  ############################################################
  # need to get a slice where pseudogene sits:
  my $chr_name  = $exons[0]->contig->chr_name;
  my $chr_start = $exons[0]->contig->chr_start;
  my $chr_end   = $exons[0]->contig->chr_end;
  
  ############################################################
  # get slice from the same db where our pseudogene is
  my $slice = $db->get_SliceAdaptor->fetch_by_chr_start_end( $chr_name, $chr_start + $start - 1 - 100, $chr_start + $end - 1 + 100);
  
  my $rpfa = $db->get_RepeatFeatureAdaptor();
  my @features = @{$rpfa->fetch_all_by_Slice($slice, 'RepeatMask')};
  
  ############################################################
  # find repeats that overlap at least 50% of the transript
  my @overlapping;
  print STDERR "found ".scalar(@features)." repeat features\n";
  foreach my $feature ( @features ){
    next unless ( $strand == $feature->strand );
    my $overlap = 0;
    my $overlap_length = 0;
    my $feature_start = $slice->chr_start + $feature->start - 1;
    my $feature_end   = $slice->chr_start + $feature->end - 1;
    
    if ( $feature_start <= $start && $feature_end >= $start ){
      $overlap = 1;
      if ( $feature_end <= $end ){
	$overlap_length = $feature_end - $start + 1;
      }
      if ( $feature_end > $end ){
	$overlap_length = $end - $start + 1;
      }
    }
    if ( $feature_start >= $start && $feature_end >= $start ){
      $overlap = 1;
      if ( $feature->end >= $end ){
	$overlap_length = $end - $feature_start + 1;
      }
      if ( $feature_end < $end ){
	$overlap_length = $feature_end - $feature_start + 1;
      }
    }
    if ( $overlap ){
      my $coverage = sprintf "%.2f", 100 * ( $overlap_length ) / $length;
      if ( $coverage > 50 ){
	push(@overlapping, $feature);
      }
    }
  }
  
  if ( @overlapping ){
    return 1;
  }
  else{
    return 0;
  }  
}

############################################################
# this method checks whether a transcript ( $transcript ) 
# is based on evidence which is in fact spliced somethere else 
# in the genome (provided by a database $db )

sub check_for_spliced_real_gene{
  my ($self,$transcript,$db) = @_;
  
  # get the evidence split by dna and protein align features:
  my ($prot_ids, $cdna_ids) = $self->_get_split_evidence($transcript); 
  
  # we take only one spliced transcript per evidence (if there is any)
  my @spliced_transcripts;
  if ( $prot_ids ){
    foreach my $id ( @$prot_ids ){
      my @trans =  Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils
	->find_transcripts_by_protein_evidence($id,$db);
      print STDERR scalar(@trans)." found with protein evidence $id\n";
      if ( @trans ){
	my @sorted = 
	  sort { scalar(@{$b->get_all_Exons}) <=> scalar(@{$a->get_all_Exons}) } @trans;
	
	if ( Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->is_spliced($sorted[0]) && !$self->overlap( $transcript, $sorted[0] ) ) {
	  push ( @spliced_transcripts, $sorted[0] );
	}
      }
    }
  }
  
  if ( $cdna_ids ){
    foreach my $id ( @$cdna_ids ){
      my @trans =  Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils
	->find_transcripts_by_dna_evidence($id,$db);
      print STDERR scalar(@trans)." found with cdna evidence $id\n";

      if ( @trans ){
	my @sorted = 
	  sort { scalar(@{$b->get_all_Exons}) <=> scalar(@{$a->get_all_Exons}) } @trans;
	
	if ( Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->is_spliced($sorted[0]) && !$self->overlap( $transcript, $sorted[0] ) ){
	  push ( @spliced_transcripts, $sorted[0] );
	}
      }
    }
  }
  return \@spliced_transcripts;
}
  
############################################################
# it checks the overlap of the extent of 2 transcripts

sub overlap{
  my ($self, $t1, $t2 ) = @_;
  my @e1 = sort { $a->start <=> $b->start } @{$t1->get_all_Exons};
  my @e2 = sort { $a->start <=> $b->start } @{$t2->get_all_Exons};
  unless ( $e1[0]->contig->chr_name eq $e2[0]->contig->chr_name ){
    return 0;
  }
  my $start1 = $e1[0]->start;
  my $end1   = $e1[$#e1]->end;
  my $start2 = $e2[0]->start;
  my $end2   = $e2[$#e2]->end;

  if ( !( $start1 > $end2 ) && !( $end1< $start2)  ){
    return 1;
  }
  return 0;
}

############################################################

sub _get_split_evidence{
    my ($self,$transcript) = @_;
    my %prot;
    my %cdna;
    foreach my $exon (@{$transcript->get_all_Exons}){
      foreach my $evidence ( @{$exon->get_all_supporting_features} ){
	#print STDERR "evidence is a $evidence\n";
	if ( $evidence->isa('Bio::EnsEMBL::DnaDnaAlignFeature') ){
	  $cdna{$evidence->hseqname}++;
	}
	if ( $evidence->isa('Bio::EnsEMBL::DnaPepAlignFeature') ){
	  $prot{$evidence->hseqname}++;
	}
      }
    }
    my @prot_ids = keys %prot;
    my @cdna_ids = keys %cdna;
    return (\@prot_ids, \@cdna_ids);
  }





############################################################

sub read_probabilities{
  my ($self,$positivo, $negativo, $attributes ) = @_;
  
  my @attributes = @$attributes;

  ############################################################
  # target
  my %positivo;
  foreach my $att ( @attributes ){
    $positivo{$att}     = 0;
  }
  
  my $pos_count = 0;
  
  print STDERR "about to open $positivo\n";
  open( POSITIVO, "<$positivo") || die("cannot open $positivo");
  
  while (<POSITIVO>){
    chomp;
    my @att = split;
    next unless ( $att[1] eq '1' || $att[1] eq '0' );
    
    $positivo{frameshift}     += $att[0];
    $positivo{polyA}          += $att[1];
    $positivo{Met}            += $att[2];
    $positivo{mouse_homology} += $att[3];
    $positivo{break_synteny}  += $att[4];
    $positivo{repeat}         += $att[5];
    
    $pos_count++;
  }
  close(POSITIVO);
  
  ############################################################
  # target
  my %negativo;
  
  foreach my $att ( @attributes ){
    $negativo{$att}     = 0;
  }
  
  my $neg_count = 0;
  
  print STDERR "about to open $negativo\n";
  open( NEGATIVO, "<$negativo") || die("cannot open $negativo");
  
  while (<NEGATIVO>){
    chomp;
    my @att = split;
    next unless ( $att[1] eq '1' || $att[1] eq '0' );
    
    $negativo{frameshift}     += $att[0];
    $negativo{polyA}          += $att[1];
    $negativo{Met}            += $att[2];
    $negativo{mouse_homology} += $att[3];
    $negativo{break_synteny}  += $att[4];
    $negativo{repeat}         += $att[5];
    
    $neg_count++;
  }


  close(NEGATIVO);

  ############################################################
  # calculate the likelihoods
  
  my %Ppos;
  my %Pneg;
  
  foreach my $att ( @attributes ){
    $Ppos{$att}{0} = ( $pos_count - $positivo{$att} )/ $pos_count;
    $Ppos{$att}{1} = ( $positivo{$att} )/ $pos_count;
    $Pneg{$att}{0} = ( $neg_count - $negativo{$att} )/ $neg_count;
    $Pneg{$att}{1} = ( $negativo{$att} )/ $neg_count;
  }
  
  return ( $pos_count, \%Ppos, $neg_count, \%Pneg );
}


############################################################




1;
