=head1 NAME

GenePair

=head1 SYNOPSIS

my $gene_pair = Bio::EnsEMBL::Pipeline::GeneComparison::GenePair->new($gene1,$gene2);

$gene_pair->compare_isoforms();

=head1 DESCRIPTION

Class to compare the isoforms between two orthologous genes.
It carry out two comparisons. First it pairs up the transcript by using blastn
and a stable-marriage optimization procedure. Then it
takes each of the found transcript-pairs and pair up the exons by
using a Needlemann-Wusch method on the exon space, using as score matrix the
tblastx comparison between exons.
Alternatively, the exons can be paired up
using blastz on the genomic extent (plus some external padding)
of the transcripts to be compared. Then exon alignments are calculated from the
blastz alignments. In this way, a possible use of non orthologous
exons can be resolved as they would not overlap in the genomic comparison.

=head1 CONTACT

eae@sanger.ac.uk

=cut

# Let the code begin ...

package Bio::EnsEMBL::Pipeline::GeneComparison::GenePair;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::GeneComparison::ObjectMap;
use Bio::EnsEMBL::Pipeline::GeneComparison::ExonPair;
use Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptPair;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Bio::EnsEMBL::Pipeline::Runnable::Blast;
use Bio::EnsEMBL::Pipeline::Runnable::Blastz;

@ISA = qw(Bio::EnsEMBL::Root);

=head1 METHODS

=cut

#########################################################################

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);
  
    my( $gap_penalty ) = $self->_rearrange([qw(
					       GAP_PENALTY
					       )], 
					   @args);
 
    if ( $gap_penalty ){
	$self->gap_penalty( $gap_penalty );
    }
    else{
	$self->gap_penalty( -100 );
    }
    
    return $self;
}

############################################################

sub gap_penalty{
    my ($self,$value) = @_;
    if ( $value ){
	$self->{_gap_penalty} = $value;
    }
    return $self->{_gap_penalty};
}

############################################################
# this method reads two genes suposedly orthologs (in slice coordinates)
# and compares its transcripts all against all using blastn.
# The best possible transcript-pairs are derived using
# the stable-mariage algorithm. 
# For each one of these selected pairs
# ir runs a global dynamic programming alignment algorithm
# on exon space to match up the exons.
#
# You can pass a flag to run only with coding exons.

sub compare{
  my ($self,$human_gene,$mouse_gene, $coding_exons) = @_;
  
  my @human_transcripts = @{$human_gene->get_all_Transcripts};
  my @mouse_transcripts = @{$mouse_gene->get_all_Transcripts};
  
  ############################################################
  # we make a pair only if the transcripts align with gaps 
  # no longer than the smallest exon
  
  my $object_map = Bio::EnsEMBL::Pipeline::GeneComparison::ObjectMap->new();
  my @transcript_matches;
  foreach my $human_t ( @human_transcripts ){
      foreach my $mouse_t ( @mouse_transcripts ){
	#print STDERR "blasting isoforms\n";
	
	############################################################
	# blast transcripts
	my ($score,$pair) = 
	  Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptPair->blast_isoforms( $human_t, $mouse_t, $coding_exons );
	
	if ( $score && $pair ){
	  $object_map->match($human_t, $mouse_t, $score );
	}
	unless ($score){
	  $score = 0;
	}
	my $id1 = $human_t->stable_id || $human_t->dbID;
	my $id2 = $mouse_t->stable_id || $mouse_t->dbID;
	print STDERR "Pair ( $id1 , $id2 ) score = $score\n";
      }
    }
  
  
  unless ( $object_map->list1 && scalar($object_map->list1) ){
    return 0;
  }
  
  my $best_pairs_object = $object_map->stable_marriage;
  #my $best_pairs_object = $object_map;
  ############################################################
  # pairs created:
  my $pair_count = scalar($best_pairs_object->list1);
  print STDERR "pairs created: ".$pair_count."\n";
  foreach my $element1 ( $best_pairs_object->list1 ){
      foreach my $partner ( $best_pairs_object->partners( $element1 ) ){
      # there should be only one
      my $id1 = $element1->stable_id || $element1->dbID;
      my $id2 = $partner->stable_id || $partner->dbID;
      print STDERR "Pair ( $id1 , $id2 ) with score: ".$best_pairs_object->score( $element1, $partner )."\n";
      #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($element1);
      #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($partner);
    }
  }
  
  ############################################################
  # compare the exons for each pair
  print STDERR "comparing exons\n";
  my $conserved_count        = 0;
  my $skipped_exons_count    = 0;
  my $terminal_missing_count = 0;
  my $all_conserved_count    = 0;

  my $human_id = $human_gene->stable_id || $human_gene->dbID;
  my $mouse_id = $mouse_gene->stable_id || $mouse_gene->dbID;
  
  foreach my $element1 ( $best_pairs_object->list1 ){
      foreach my $partner ( $best_pairs_object->partners( $element1 ) ){
	  
      # there should be only one
	  #print STDERR "Pair with score: ".$best_pairs_object->score( $element1, $partner )."\n";
	  #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($element1);
	  #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($partner);
	  
	my ($missing_terminal_exons, $exon_skipping, $all_exons_conserved) = 
	  $self->compare_Exons( $element1, $partner, $self->gap_penalty, $coding_exons, $human_id, $mouse_id);
	  
	  if ($exon_skipping ){
	      $skipped_exons_count++;
	  }
	  if ($missing_terminal_exons){
	      $terminal_missing_count++;
	  }
	  unless ( $exon_skipping || $missing_terminal_exons ){
	      $conserved_count++;
	  }
	  if ( $all_exons_conserved ){
	      $all_conserved_count++;
	  }
	  
	  #$self->print_alignment( $element1, $partner, \@score_matrix);
	  
      }
  }

  my $human_trans_count = scalar(@human_transcripts);
  my $mouse_trans_count = scalar(@mouse_transcripts);
  print STDERR "GENEPAIR\t".
      "$human_id\thuman_trans_count:$human_trans_count\t".
	  "$mouse_id\tmouse_trans_count:$mouse_trans_count\t".
	    "pairs:$pair_count\t".
		  "conserved:$conserved_count\t".
		      "with_skipped_exons:$skipped_exons_count\t".
			  "with_missing_terminals:$terminal_missing_count\n";
  return 1;
}

############################################################
# This method is similar to the one above, except 
# for the exon comparison step, which uses the genomic extent
# of the transcripts in the pair.
#
# This method reads two genes suposedly orthologs (in slice coordinates)
# and compares its transcripts all against all using blastn.
# The best possible transcript-pairs are derived using
# the stable-mariage algorithm. 
# For each one of these selected pairs
# it uses blastz to compare the genomic extent of the transcripts
# and locate whether the exons align with each other.
#
# You can pass a flag to run only with coding exons.

sub find_exact_matches{
  my ($self,$human_gene,$mouse_gene, $coding_exons) = @_;
  
  my @human_transcripts = @{$human_gene->get_all_Transcripts};
  my @mouse_transcripts = @{$mouse_gene->get_all_Transcripts};

  my $gene_id1 = $human_gene->stable_id || $human_gene->dbID;
  my $gene_id2 = $mouse_gene->stable_id || $mouse_gene->dbID;
  
  ############################################################
  # we make a pair only if the transcripts align with gaps 
  # no longer than the smallest exon
  
  my $object_map = Bio::EnsEMBL::Pipeline::GeneComparison::ObjectMap->new();
  my @transcript_matches;
  foreach my $human_t ( @human_transcripts ){
    foreach my $mouse_t ( @mouse_transcripts ){
      print STDERR "blasting isoforms\n";
      
      ############################################################
      # blast transcripts
      my ($score,$pair) = 
	Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptPair->blast_isoforms( $human_t, $mouse_t, $coding_exons );
      
      if ( $score && $pair ){
	$object_map->match($human_t, $mouse_t, $score );
      }
      unless ($score){
	$score = 0;
      }
      my $id1 = $human_t->stable_id || $human_t->dbID;
      my $id2 = $mouse_t->stable_id || $mouse_t->dbID;
      print STDERR "Pair ( $id1 , $id2 ) score = $score\n";
    }
  }
  unless ( $object_map->list1 && scalar($object_map->list1) ){
    return 0;
  }

  my $best_pairs_object = $object_map->stable_marriage;
  #my $best_pairs_object = $object_map;

  ############################################################
  # pairs created:
  my $pair_count = scalar($best_pairs_object->list1);
  print STDERR "pairs created: ".$pair_count."\n";
  foreach my $element1 ( $best_pairs_object->list1 ){
    foreach my $partner ( $best_pairs_object->partners( $element1 ) ){
      # there should be only one
      my $id1 = $element1->stable_id || $element1->dbID;
      my $id2 = $partner->stable_id || $partner->dbID;
      print STDERR "Pair ( $id1 , $id2 ) with score: ".$best_pairs_object->score( $element1, $partner )."\n";
      #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($element1);
      #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($partner);
    }
  }
  
  ############################################################
  # check with the genomic alignment that the exons are
  # really overlapping each other:
  my $conserved_count        = 0;
  my $skipped_exons_count    = 0;
  my $terminal_missing_count = 0;
  my $all_conserved_count    = 0;

 HUMAN:
  foreach my $human_t ( $best_pairs_object->list1 ){
      
    MOUSE:
      foreach my $mouse_t ( $best_pairs_object->partners( $human_t ) ){
	  
	  ############################################################
	  # blastz genomic extent of transcripts
	  my $id1 = $human_t->stable_id || $human_t->dbID;
	  my $id2 = $mouse_t->stable_id || $mouse_t->dbID;
	  print STDERR "comparing genomic extent of $id1 and $id2:\n";
	  #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($human_t);
	  #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($mouse_t);
	  
	  my ($missing_terminal_exons, $exon_skipping, $all_conserved) = Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptPair->blast_genomic_isoforms( $human_t, $mouse_t, $coding_exons, $gene_id1, $gene_id2);
	  
	  if ($exon_skipping ){
	      $skipped_exons_count++;
	  }
	  if ($missing_terminal_exons){
	      $terminal_missing_count++;
	  }
	  unless ( $exon_skipping || $missing_terminal_exons ){
	      $conserved_count++;
	  }
	  if ( $all_conserved ){
	      $all_conserved++;
	  }
	  
	  #$self->print_alignment( $element1, $partner, \@score_matrix);
	  
      }
  }
  
  my $human_id = $gene_id1;
  my $mouse_id = $gene_id2;
  my $human_trans_count = scalar(@human_transcripts);
  my $mouse_trans_count = scalar(@mouse_transcripts);
  print STDERR "BLASTZ_GENEPAIR\t".
      "$human_id\thuman_trans_count:$human_trans_count\t".
	  "$mouse_id\tmouse_trans_count:$mouse_trans_count\t".
	    "pairs:$pair_count\t".
		  "conserved:$conserved_count\t".
		      "with_skipped_exons:$skipped_exons_count\t".
			  "with_missing_terminals:$terminal_missing_count\n";
  


  

  return 1;
}




############################################################
#
#
############################################################




sub max{
    my ($self, $max, @values ) = @_;
    for (my $i=0; $i<@values; $i++ ){
	$max = $values[$i] if $values[$i]>$max;
    }
    return $max;
}

############################################################

sub min{
    my ($self, $min, @values ) = @_;
    for (my $i=0; $i<@values; $i++ ){
	$min = $values[$i] if $values[$i]<$min;
    }
    return $min;
}

############################################################

sub print_Feature{
    my ($self,$f, $cigar) = @_;
    my $string =
	$f->seqname."\t".
	    $f->start."-".$f->end."\t".
		($f->end - $f->start + 1)."\t".
		    $f->strand."\t".
			$f->hseqname."\t".
			    $f->hstart."-".$f->hend."\t".
				($f->hend - $f->hstart + 1 )."\t".
				    $f->strand."\t".
					"score:".$f->score."\t".
					    "perc_id:".$f->percent_id;
    if ($cigar){
	$string .= "\t".$f->cigar_string;
    }
    print STDERR $string."\n";
}
	    
############################################################
# dynamic programming method to align the exons
# It uses a global alignment algorithm to
# pair up the exons from each transcript

sub compare_Exons{
  my ($self,$human_t, $mouse_t, $gap_penalty, $coding_exons, $human_gene_id, $mouse_gene_id ) = @_;
  

  # get the exons 5' to 3'
  my @human_exons = @{$self->get_Exons($human_t, $coding_exons)};
  my @mouse_exons = @{$self->get_Exons($mouse_t, $coding_exons)};
  
  my @score_matrix;
  my %comparison_score;
  my %features;

  my $human_length = scalar(@human_exons);
  my $mouse_length = scalar(@mouse_exons);
  
  foreach my $i (0..$human_length){
      $score_matrix[$i][0] = $i * $gap_penalty;
  }
  foreach my $j (0..$mouse_length){
      $score_matrix[0][$j] = $j * $gap_penalty;
  }
  
  my $exon_pair = Bio::EnsEMBL::Pipeline::GeneComparison::ExonPair->new();
  foreach my $i ( 1..$human_length ){
    foreach my $j ( 1..$mouse_length ){
      
      my $human_exon = $human_exons[$i-1];
      my $mouse_exon = $mouse_exons[$j-1];

      ############################################################
      # use correct phases when using coding exons to avoid lingering -1
      if ( $coding_exons ){
	if ( $i == 1 && $human_exon->phase == -1 ){
	  $human_exon->phase(0);
	}
	if ( $i == $human_length ){
	    #print STDERR "changing phase of human_exon($i) from ".$human_exon->end_phase." to "
	    #. (($human_exon->phase + $human_exon->length) %3 )."\n";
	  $human_exon->end_phase( ($human_exon->phase + $human_exon->length) %3 );
	}
	if ( $j == 1 && $mouse_exon->phase == -1 ){
	  $mouse_exon->phase(0);
	}
	if ( $j == $mouse_length ){
	  #print STDERR "changing phase of mouse_exon($j) from ".$mouse_exon->end_phase." to "
	    #. (($mouse_exon->phase + $mouse_exon->length) %3 )."\n";
	    $mouse_exon->end_phase( ($mouse_exon->phase + $mouse_exon->length) %3 );
	}
      }

      ($comparison_score{$human_exon}{$mouse_exon},
       $features{$human_exons[$i-1]}{$mouse_exons[$j-1]}) =
	 $exon_pair->blast_Exons( $human_exon, $mouse_exon );
      
      #print STDERR "comparison( ".$human_exons[$i-1]->stable_id."-".$mouse_exons[$j-1]->stable_id." ) = ".
      #$comparison_score{$human_exons[$i-1]}{$mouse_exons[$j-1]} ."\n";
	
      $score_matrix[$i][$j] = 
	$self->max( $score_matrix[$i-1][$j]   + $gap_penalty,
		    $score_matrix[$i][$j-1]   + $gap_penalty,
		    $score_matrix[$i-1][$j-1] + $comparison_score{$human_exons[$i-1]}{$mouse_exons[$j-1]} );
    }
  }
  
  my ($human_list, $mouse_list ) = 
      $self->get_alignment( \@human_exons, \@mouse_exons, \@score_matrix, \%comparison_score, $gap_penalty );
  
  my $human_missing = 0;
  my $mouse_missing = 0;
  my $human_internal_missing = 0;
  my $mouse_internal_missing = 0;
  my $exon_skipping = 0;
  my $conserved     = 0;
  my $same_length   = 0;
  my $same_phases   = 0;

  for ( my $i=0 ; $i< scalar(@$human_list); $i++ ){
    if ( $human_list->[$i] eq 'gap' ){
      $human_missing++;
    }
    if ( $mouse_list->[$i] eq 'gap'){
      $mouse_missing++;
    }
    if ( !( $human_list->[$i] eq 'gap' || $mouse_list->[$i] eq 'gap') ){
      $conserved++;
      my $human_length = $human_list->[$i]->end - $human_list->[$i]->start + 1;
      my $mouse_length = $mouse_list->[$i]->end - $mouse_list->[$i]->start + 1;
      if ( $human_length == $mouse_length ){
	$same_length++;
      }
      if ( $human_list->[$i]->phase == $mouse_list->[$i]->phase 
	   &&
	   $human_list->[$i]->end_phase == $mouse_list->[$i]->end_phase ){
	$same_phases++;
      }
    }
    if ( $i > 0 && $i< ( scalar(@$human_list) - 1 ) ){
	if ( $self->haspredecessor( $human_list, $i) 
	     && 
	     $self->hassuccessor( $human_list, $i )
	     &&
	     $self->haspredecessor( $mouse_list, $i )
	     && 
	     $self->hassuccessor( $mouse_list, $i )
	     ){
	    if ( $human_list->[$i] eq 'gap' ){ 
		$exon_skipping++;
		$human_internal_missing++;
	    }
	    if( $mouse_list->[$i] eq 'gap'){
		$exon_skipping++;
		$mouse_internal_missing++;
	    }
	}
    }
}
  
  my $human_terminal_missing = $human_missing - $human_internal_missing;
  my $mouse_terminal_missing = $mouse_missing - $mouse_internal_missing;
  
  my $human_id = $human_t->stable_id || $human_t->dbID;
  my $mouse_id = $mouse_t->stable_id || $mouse_t->dbID;
  
  print STDERR "TRANPAIR\t".
    "$human_gene_id\t$human_id\thuman_exons:$human_length\thuman_miss_term_exons:$human_terminal_missing\thuman_miss_int_exons:$human_internal_missing\t".
      "conserved_exons:$conserved\twith_same_length:$same_length\twith_same_phase:$same_phases\t".
	"$mouse_gene_id\t$mouse_id\tmouse_exons:$mouse_length\tmouse_miss_term_exons:$mouse_terminal_missing\tmouse_miss_int_exons:$mouse_internal_missing\n";
  
  my $print_report = 1;
  if ( $print_report ){
    for ( my $i=0; $i<scalar(@$human_list); $i++ ){
      my $human_string;
      my $mouse_string;
      my $cigars = '';
      my $score  = 0;
      if ( $human_list->[$i] eq 'gap'){
	$human_string = "             ####GAP####";
      }
      else{
	$human_string = $self->exon_string( $human_list->[$i] );
      }
      if ( $mouse_list->[$i] eq 'gap'){
	$mouse_string = "             ####GAP####";
      }
      else{
	$mouse_string = $self->exon_string( $mouse_list->[$i] );
      }
      if( !($human_string eq "gap" || $mouse_string eq "gap") ){
	#print STDERR "score (".$human_list->[$i]->stable_id.",".$mouse_list->[$i]->stable_id.")=".$comparison_score{$human_list->[$i]}{$mouse_list->[$i]}."\n";
	$score = $comparison_score{$human_list->[$i]}{$mouse_list->[$i]};
	
	foreach my $feat ( @{ $features{$human_list->[$i]}{$mouse_list->[$i]} } ){
	  $cigars = $feat->cigar_string;
	}
      }
      $score = 0 unless $score;
      $score = sprintf "%.2f", $score;
      print STDERR $human_string."\t<---->\t".$mouse_string.
	"\t score= ".$score."\t$cigars\n";
    }
  }
  my $missing_terminal_exons = 0;
  if ( $human_terminal_missing || $mouse_terminal_missing ){
    $missing_terminal_exons = 1;
  }
  my $all_conserved = 0;
  if ( $same_length == scalar( @$human_list ) ){
    $all_conserved = 1;
  }
  return ($missing_terminal_exons, $exon_skipping, $all_conserved);
}

############################################################

sub get_Exons{
  my ( $self, $trans , $coding) = @_;
  my @exons;
  my @newexons;
  my $strand = $trans->start_Exon->strand;

  my $verbose = 0;

  if ( $coding ){
    if ( $strand == 1 ){
      @exons = sort {$a->start <=> $b->start} @{$trans->get_all_translateable_Exons};
    }
    else{
      @exons = sort {$b->start <=> $a->start} @{$trans->get_all_translateable_Exons};
    }
  }
  else{
    if ( $strand == 1 ){
      @exons = sort {$a->start <=> $b->start} @{$trans->get_all_Exons};
    }
    else{
      @exons = sort {$b->start <=> $a->start} @{$trans->get_all_Exons};
    }
  }

  #print STDERR "exons:\n";
  #foreach my $exon (@exons){
  #  print STDERR $exon->contig->name." length: ".($exon->end - $exon->start + 1)
  #    ." exon_seq length: ".$exon->seq->length." seq length ".length($exon->seq->seq)."\n";
  #}

  my $c=0;
  for (my $i=0; $i< scalar(@exons); $i++ ){
    if ( $i>0 && $strand == 1 ){
      if ( $exons[$i]->start - $exons[$i-1]->end - 1 < 10 ){
	
	# adjust_start_end() creates a new exon object!
        my $shift = $exons[$i]->end - $exons[$i-1]->end;
	print STDERR "adding right $shift bp to exon: ".$exons[$i-1]->start."-".$exons[$i-1]->end."\n" if $verbose;
	
	$newexons[$c-1] = $exons[$i-1]->adjust_start_end(0,$shift);
	
	print STDERR "new exon: ".$newexons[$c-1]->start."-".$newexons[$c-1]->end."\n" if $verbose;
	#$exons[$i-1]->end($exons[$i]->end);
	next;
      }
    }
    if ( $i>0 && $strand == -1 ){
      if ( $exons[$i-1]->start - $exons[$i]->end - 1 < 10 ){
	my $shift = $exons[$i-1]->start - $exons[$i]->start;
	print STDERR "adding left $shift bp to exon: ".$exons[$i-1]->start."-".$exons[$i-1]->end."\n" if $verbose;
	# adjust_start_end() creates a new exon object!
	
	$newexons[$c-1] = $exons[$i-1]->adjust_start_end(0,$shift);
	
	print STDERR "new exon: ".$newexons[$c-1]->start."-".$newexons[$c-1]->end."\n" if $verbose;
	#$exons[$i-1]->start($exons[$i]->start);
	next;
      }
    }
    push (@newexons, $exons[$i] );
    $c++;
  }
  
  #print STDERR "New exons:\n";
  #foreach my $exon (@newexons){
  #  print STDERR $exon->contig->name." length: ".($exon->end - $exon->start + 1)
  #    ." exon_seq length: ".$exon->seq->length." seq length ".length($exon->seq->seq)."\n";
  #}
  

  return \@newexons;
}

############################################################

sub haspredecessor{
    my ($self,$list,$i) = @_;
    for( my $j=$i-1; $j>0; $j-- ){
	if ( !( $list->[$j] eq 'gap' ) ){
	    return 1;
	}
    }
    return 0;
}
############################################################

sub hassuccessor{
    my ($self,$list,$i) = @_;
    for( my $j=$i+1; $j< scalar(@$list); $j++ ){
	if ( !( $list->[$j] eq 'gap' ) ){
	    return 1;
	}
    }
    return 0;
}
############################################################
# method to recover the alignment

sub get_alignment{
  my ($self,$human_list, $mouse_list, $matrix, $comparison, $gap_penalty) = @_;
  my @matrix     = @$matrix;  
  my %comparison = %$comparison;
  my @human_list = @$human_list;
  my @mouse_list = @$mouse_list;

  my $human_length = scalar( @human_list );
  my $mouse_length = scalar( @mouse_list );

  unless( $human_length ){
      for ( my $i=1; $i<= $mouse_length; $i++ ){
	  push( @human_list, "gap" );
      }
      return ( \@human_list, \@mouse_list );
  }
  unless( $mouse_length ){
      for ( my $j=1; $j<= $human_length; $j++ ){
	  push( @mouse_list, "gap" );
      }
      return ( \@human_list, \@mouse_list );
  }
  
  my $human_last = $human_list[-1];
  my $mouse_last = $mouse_list[-1];
  
  ############################################################
  # last exons are paried-up in the optimal alignment
  if ( $matrix[$human_length][$mouse_length] 
       == $matrix[$human_length-1][$mouse_length-1] 
       + $comparison{$human_last}{$mouse_last} ){
      pop @human_list;
      pop @mouse_list;
      my ( $human_list2, $mouse_list2) = 
	  $self->get_alignment( \@human_list, \@mouse_list, $matrix, $comparison, $gap_penalty);
      push ( @{$human_list2}, $human_last );
      push ( @{$mouse_list2}, $mouse_last );
      return ( $human_list2, $mouse_list2 );
  }
  ############################################################
  # last exon of the first list is paired-up with a gap
  elsif( $matrix[$human_length][$mouse_length] 
	 == $matrix[$human_length-1][$mouse_length] + $gap_penalty ){
    pop @human_list;
    my ( $human_list2, $mouse_list2) =
      $self->get_alignment( \@human_list, \@mouse_list, $matrix, $comparison, $gap_penalty);
    push ( @{$human_list2}, $human_last );
    push ( @{$mouse_list2}, "gap" );
    return ( $human_list2, $mouse_list2 );
  }
  ############################################################
  # last exons of the second list is paired up with a gap
  else{
    pop @mouse_list;
    my ( $human_list2, $mouse_list2) =
      $self->get_alignment( \@human_list, \@mouse_list, $matrix, $comparison, $gap_penalty);
    push ( @{$human_list2}, "gap" );
    push ( @{$mouse_list2}, $mouse_last );
    return ( $human_list2, $mouse_list2 );
  }
} 

############################################################

sub exon_string{
  my ($self,$exon) = @_;
  #my $id = $exon->stable_id || $exon->dbID;
  my $string = $exon->seqname.":".$exon->start."-".$exon->end.
    " (".($exon->end - $exon->start + 1 ).")".
      " strand:".$exon->strand.
	" phase:".$exon->phase.
	  " endphase:".$exon->end_phase;
  
}    

############################################################

sub print_exons_in_transcript{
    my ($self,$tran) = @_;
    my @exons =  sort { $a->start <=> $b->start } @{$tran->get_all_Exons};
    my $length = 0;
    my $start  = 1;
    my $end;
    foreach my $exon ( @exons ){
	$start += $length;
	$length = $exon->length;
	$end = $start + $length - 1;
	print STDERR "$start-$end ($length) ";
	
    }
    print STDERR "\n";
}


############################################################
# this method will compare the CDSs of two genes
# it will generate transcript pairs according to 
# the protein-protein alignment of the transcripts
# and it will generate output for the transcript pairs and
# a summary output for the gene comparison

sub compare_CDSs{
  my ($self,$human_gene,$mouse_gene) = @_;
  
  my @human_transcripts = @{$human_gene->get_all_Transcripts};
  my @mouse_transcripts = @{$mouse_gene->get_all_Transcripts};
  
  ############################################################
  # we make a pair only if the transcripts align with gaps 
  # no longer than the smallest exon
  
  my $object_map = Bio::EnsEMBL::Pipeline::GeneComparison::ObjectMap->new();
  my @transcript_matches;
  my %target_coverage;
  my %query_coverage;
  my %perc_id;
  foreach my $human_t ( @human_transcripts ){
    foreach my $mouse_t ( @mouse_transcripts ){
      
      ############################################################
      # blast CDSs
      my ($score,$best_features, $target_coverage, $query_coverage, $perc_id) = 
	Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptPair->blast_CDSs( $human_t, $mouse_t);
      
      if ( $score && $best_features ){
	$object_map->match($human_t, $mouse_t, $score );
	$perc_id{$human_t}{$mouse_t} = $perc_id;
	$target_coverage{$human_t}{$mouse_t} = $target_coverage;
	$query_coverage{$human_t}{$mouse_t} = $query_coverage;
      }
      unless ($score){
	$score = 0;
      }
      my $id1 = $human_t->stable_id || $human_t->dbID;
      my $id2 = $mouse_t->stable_id || $mouse_t->dbID;
      print STDERR "Pair ( $id1 , $id2 ) score = $score\n";
    }
  }
  
  
  unless ( $object_map->list1 && scalar($object_map->list1) ){
    return 0;
  }
  
  my $best_pairs_object = $object_map->stable_marriage;
  #my $best_pairs_object = $object_map;
  ############################################################
  # pairs created:
  my $pair_count = scalar($best_pairs_object->list1);
  print STDERR "pairs created: ".$pair_count."\n";
  foreach my $element1 ( $best_pairs_object->list1 ){
    foreach my $partner ( $best_pairs_object->partners( $element1 ) ){
      # there should be only one
      my $id1 = $element1->stable_id || $element1->dbID;
      my $id2 = $partner->stable_id || $partner->dbID;
      print STDERR "Pair ( $id1 , $id2 ) with score: ".$best_pairs_object->score( $element1, $partner )."\n";
      
      ############################################################
      # summary line
      my $gene_id1 = $human_gene->stable_id || $human_gene->dbID;
      my $gene_id2 = $mouse_gene->stable_id || $mouse_gene->dbID;
      my @exons1 = @{Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptPair->get_Exons($element1,1)};
      my @exons2 = @{Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptPair->get_Exons($partner,1)};
      my $exons1   = scalar( @exons1 );
      my $exons2   = scalar( @exons2 );
      
      my @real_exons1 = @{$element1->get_all_translateable_Exons};
      my @real_exons2 = @{$partner->get_all_translateable_Exons};
      
      my ($length1,$length2) = (0,0);
      foreach my $exon ( @exons1 ){
	$length1 += $exon->length;
      }
      foreach my $exon ( @exons2 ){
	$length2 += $exon->length;
      }

      my ($real_length1,$real_length2) = (0,0);
      foreach my $exon ( @real_exons1 ){
	$real_length1 += $exon->length;
      }
      foreach my $exon ( @real_exons2 ){
	$real_length2 += $exon->length;
      }
      
      my $perc_id         = $perc_id{$element1}{$partner};
      my $target_coverage = $target_coverage{$element1}{$partner};
      my $query_coverage  = $query_coverage{$element1}{$partner};
      my $order = 0;
      $order = 1 if ( $exons1 == $exons2 );
      my $length_diff      = $length1 - $length2;
      my $real_length_diff = $real_length1 - $real_length2;
      print STDERR "CDS_PAIR\t".
	"$gene_id1\t$id1\texons1:$exons1\tlength1:$length1\treal_length1:$real_length1\t".
	  "coverage1:$target_coverage\t".
	  "$gene_id2\t$id2\texons2:$exons2\tlength2:$length2\treal_length2:$real_length2\t".
	    "coverage2:$query_coverage\t".
	    "perc_id:$perc_id\tsame_exon_number:$order\t".
	      "length_diff:$length_diff\treal_length_diff:$real_length_diff\n";
      #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($element1);
      #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($partner);
    }
  }
  
  ############################################################
  # compare the exons for each pair
  print STDERR "comparing exons\n";
  my $conserved_count        = 0;
  my $skipped_exons_count    = 0;
  my $terminal_missing_count = 0;
  my $all_conserved_count    = 0;
  my $coding_exons = 1;
  
  my $human_id = $human_gene->stable_id || $human_gene->dbID;
  my $mouse_id = $mouse_gene->stable_id || $mouse_gene->dbID;
  
  foreach my $element1 ( $best_pairs_object->list1 ){
    foreach my $partner ( $best_pairs_object->partners( $element1 ) ){
      
      # there should be only one
      #print STDERR "Pair with score: ".$best_pairs_object->score( $element1, $partner )."\n";
      #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($element1);
      #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($partner);
      
      my ($missing_terminal_exons, $exon_skipping, $all_exons_conserved) = 
	$self->compare_Exons( $element1, $partner, $self->gap_penalty, $coding_exons, $human_id, $mouse_id);
      
      if ($exon_skipping ){
	$skipped_exons_count++;
      }
      if ($missing_terminal_exons){
	$terminal_missing_count++;
      }
      unless ( $exon_skipping || $missing_terminal_exons ){
	$conserved_count++;
      }
      if ( $all_exons_conserved ){
	$all_conserved_count++;
      }
      
      #$self->print_alignment( $element1, $partner, \@score_matrix);
      
    }
  }
  
  my $human_trans_count = scalar(@human_transcripts);
  my $mouse_trans_count = scalar(@mouse_transcripts);
  print STDERR "GENEPAIR\t".
      "$human_id\thuman_trans_count:$human_trans_count\t".
	  "$mouse_id\tmouse_trans_count:$mouse_trans_count\t".
	    "pairs:$pair_count\t".
		  "conserved:$conserved_count\t".
		      "with_skipped_exons:$skipped_exons_count\t".
			  "with_missing_terminals:$terminal_missing_count\n";
  return 1;
}

1;































