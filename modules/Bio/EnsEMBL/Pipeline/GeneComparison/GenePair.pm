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
	  print STDERR "blasting isoforms\n";
	  
	  ############################################################
	  # blast transcripts
	  my ($score,$pair) = $self->blast_isoforms( $human_t, $mouse_t, $coding_exons );
	  
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
  #my $best_pairs_object = $object_map->stable_marriage;
  my $best_pairs_object = $object_map;
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

  foreach my $element1 ( $best_pairs_object->list1 ){
      foreach my $partner ( $best_pairs_object->partners( $element1 ) ){
	  
      # there should be only one
	  #print STDERR "Pair with score: ".$best_pairs_object->score( $element1, $partner )."\n";
	  #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($element1);
	  #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($partner);
	  
	  my ($missing_terminal_exons, $exon_skipping, $all_exons_conserved) = 
	      $self->compare_Exons( $element1, $partner, $self->gap_penalty, $coding_exons);
	  
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

  my $human_id = $human_gene->stable_id || $human_gene->dbID;
  my $mouse_id = $mouse_gene->stable_id || $mouse_gene->dbID;
  my $human_trans_count = scalar(@human_transcripts);
  my $mouse_trans_count = scalar(@mouse_transcripts);
  print STDERR "GENEPAIR\t".
      "$human_id\thuman_trans_count:$human_trans_count\t".
	  "$mouse_id\tmouse_trans_count:$mouse_trans_count\t".
	    "pairs:$pair_count\t".
		  "conserved:$conserved_count\t".
		      "with_skipped_exons:$skipped_exons_count\t".
			  "with_missing_terminals:$terminal_missing_count\n";
}

############################################################

sub find_exact_matches{
  my ($self,$human_gene,$mouse_gene, $coding_exons) = @_;
  
  my @human_transcripts = @{$human_gene->get_all_Transcripts};
  my @mouse_transcripts = @{$mouse_gene->get_all_Transcripts};
  
  ############################################################
  # we make a pair only if the transcripts align with gaps 
  # no longer than the smallest exon
  
  my $object_map = Bio::EnsEMBL::Pipeline::GeneComparison::ObjectMap->new();
  my @transcript_matches;
  my %exon_pairs;
 HUMAN:
  foreach my $human_t ( @human_transcripts ){
    MOUSE:
      foreach my $mouse_t ( @mouse_transcripts ){
	  
	  ############################################################
	  # blastz genomic extent of transcripts
	  my $id1 = $human_t->stable_id || $human_t->dbID;
	  my $id2 = $mouse_t->stable_id || $mouse_t->dbID;
	  print STDERR "comparing genomic extent of $id1 and $id2:\n";
	Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($human_t);
	Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($mouse_t);
	  
	  my ($missing_terminal_exons, $exon_skipping, $all_conserved) = $self->blast_genomic_isoforms( $human_t, $mouse_t, $coding_exons );
	  
	  #if ( $all_exons_match ){
	  #    $object_map->match($human_t, $mouse_t, 1); 
	  #    print STDERR "Pair ( $id1 , $id2 ) all exons match\n";
	  #    $exon_pairs{$human_t} = $exon_pairs;
	  #    next HUMAN;
	  #}
      } 
  }
  my $best_pairs_object = $object_map;
  
  ############################################################
  # compare the exons for each pair
  print STDERR "comparing exons\n";
  my $conserved_count        = 0;
  my $skipped_exons_count    = 0;
  my $terminal_missing_count = 0;
  my $all_conserved_count    = 0;
  
  foreach my $element1 ( $best_pairs_object->list1 ){
      foreach my $partner ( $best_pairs_object->partners( $element1 ) ){
	  
	  # there should be only one
	  #print STDERR "Pair with score: ".$best_pairs_object->score( $element1, $partner )."\n";
	  #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($element1);
	  #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($partner);
	  
	  my ($missing_terminal_exons, $exon_skipping, $all_exons_conserved) = 
	      $self->compare_Exons( $element1, $partner, $self->gap_penalty, $coding_exons, $exon_pairs{$element1});
	  
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
  
  my $human_id = $human_gene->stable_id || $human_gene->dbID;
  my $mouse_id = $mouse_gene->stable_id || $mouse_gene->dbID;
  my $human_trans_count = scalar(@human_transcripts);
  my $mouse_trans_count = scalar(@mouse_transcripts);
  #print STDERR "GENEPAIR\t".
  #    "$human_id\thuman_trans_count:$human_trans_count\t".
#	  "$mouse_id\tmouse_trans_count:$mouse_trans_count\t".
#	      "pairs:$pair_count\t".
#		  "conserved:$conserved_count\t".
#		      "with_skipped_exons:$skipped_exons_count\t".
#			  "with_missing_terminals:$terminal_missing_count\n";
}

############################################################
############################################################


sub blast_isoforms{
    my ( $self,$tran1,$tran2, $coding_exons ) = @_;
    
    # query
    my $id1;
    if ( $tran1->dbID ){
	$id1 = $tran1->stable_id || $tran2->dbID;
    }
    else{
	$id1 = "no id";
    }
    
    #target
    my $id2;
    if ( $tran2->dbID ){
	$id2 = $tran2->stable_id || $tran2->dbID;
    }
    else{
	$id2 = "no id";
    }
    
    print STDERR "\tcomparing $id1 and $id2\n";
    
    my ($seq1, $seq2 );
    if ( $coding_exons ){
	my $string1 = $tran1->translateable_seq;
	$seq1       = Bio::Seq->new(
				    -DISPLAY_ID => $id1,
				    -MOLTYPE    => 'dna',
				  -SEQ        => $string1,
				 );
      
      my $string2 = $tran2->translateable_seq;
      $seq2       = Bio::Seq->new(
				  -DISPLAY_ID => $id2,
				  -MOLTYPE    => 'dna',
				  -SEQ        => $string2,
				 );
    }
    else{
      $seq1    = $tran1->seq;
      unless ( $seq1->display_id ){
	$seq1->display_id($id1);
      }
      $seq2    = $tran2->seq;
      unless ( $seq2->display_id ){
	$seq2->display_id($id2);
      }
    }
    my $length1 = $seq1->length;
    my $length2 = $seq2->length;
    
    ############################################################
    # create database
    my $file = 'seq_'.$$.'.fa';
    my $database = "/tmp/".$file;
    open( DB_SEQ,">$database") || die("Could not open $database $!");
    
    my $seqout = Bio::SeqIO->new('-format' => 'Fasta',
				 '-fh'     => \*DB_SEQ);
    
    $seqout->write_seq($seq2);
    close( DB_SEQ );
    
    system("pressdb $database > /dev/null 2>&1");
    
    ############################################################
    my $options = "-nogap W=5";
    
    # Ian's parameters:
    #my $options = "W=5 M=1 N=-1 Q=3 R=3";
    
    #my $options = "W=5";
    my $blast =  
      Bio::EnsEMBL::Pipeline::Runnable::Blast->new ('-query'          => $seq1,
						    '-program'        => 'wublastn',
						    '-database'       => $database,
						    '-threshold_type' => "PVALUE",
						    #'-threshold'      => 1e-10,
						    '-options'        => $options,
						   );
    
    $blast->add_regex($file,'(\S+)');
    $blast->run();
    
    unlink ( $database );

    my @featurepairs = $blast->output();
    
    ############################################################
    # compute the average score
    my $ave_score = 0;
    foreach my $fp (sort {$a->hstart <=> $b->hstart} @featurepairs) {
      $ave_score += $fp->score;
      #print $fp->gffstring . "\n";
    }
    if ( @featurepairs ){
      $ave_score = $ave_score/scalar( @featurepairs );
    }
    else{
      $ave_score = 0;
    }

    print STDERR "\t$id1 length = $length1\n";
    print STDERR "\t$id2 length = $length2\n";
    
    ############################################################
    # calculate coverage
    my @pos_features = grep { $_->strand == 1 } @featurepairs;  
    my @neg_features = grep { $_->strand == -1} @featurepairs; 
    my @features;
    # hpos - gpos
    @{$features[0]} = grep { $_->hstrand == 1 } @pos_features;
    # hneg - gpos
    @{$features[1]} = grep { $_->hstrand == -1} @pos_features;
    # hpos - gneg
    @{$features[2]} = grep { $_->hstrand == 1 } @neg_features;
    # hneg - hneg
    @{$features[3]} = grep { $_->hstrand == -1} @neg_features;
    
    my $max_score = 0;
    my $pair;
    my %coverage;
    my $spliced = 0;
    for (my $i=0; $i<4; $i++ ){
      unless ( $features[$i] && @{$features[$i]} ){
	next;
      }
      # we use query/target as in feature pairs the target=seqname and query=hseqname
      my ($query_coverage,  $query_spliced)  = 
	$self->process_query( $features[$i], $tran2 , $coding_exons);
      my ($target_coverage, $target_spliced) = 
	$self->process_target( $features[$i], $tran1, $coding_exons);
      my $score = ( $query_coverage + $target_coverage )/2;
      
      foreach my $f ( @{$features[$i]} ){
	$self->print_Feature($f);
      }

      ### this is going to let through also spliced cases - it is just a test
      if ( $score > $max_score ){
      	$max_score = $score;
      	$pair = $features[$i];
      }
  
      if ( $query_spliced || $target_spliced ){
	$spliced = 1;
	#print STDERR "one of them is spliced\n";
      }
      #if ( $score > $max_score && $query_spliced == 0 && $target_spliced == 0 ){
      #	$max_score = $score;
      #	$pair = $features[$i];
      #}
      print STDERR "\tquery:$id1 coverage:$query_coverage spliced:$query_spliced\n";
      print STDERR "\ttarget:$id2 coverage:$target_coverage spliced:$target_spliced\n";
    }
    #return ($max_score,$pair);
    #$ave_score = 0 if $spliced;
    return ( $ave_score, $pair );
  }

############################################################

sub blast_genomic_isoforms{
  my ( $self,$tran1,$tran2, $coding_exons , $padding) = @_;
  
  $padding = 300;
  # query
  my $id1;
  if ( $tran1->dbID ){
    $id1 = $tran1->stable_id || $tran2->dbID;
  }
  else{
    $id1 = "no id";
  }
  
  #target
  my $id2;
  if ( $tran2->dbID ){
    $id2 = $tran2->stable_id || $tran2->dbID;
  }
  else{
    $id2 = "no id";
  }
  
  print STDERR "######### comparing the genomic extent of $id1 and $id2\n";
  
  my @exons;
  my @exons2;
  if ( $coding_exons ){
    @exons  = sort { $a->start <=> $b->start } @{$tran1->get_all_translateable_Exons};
    @exons2 = sort { $a->start <=> $b->start } @{$tran2->get_all_translateable_Exons};
  }
  else{
    @exons  = sort { $a->start <=> $b->start } @{$tran1->get_all_Exons};
    @exons2 = sort { $a->start <=> $b->start } @{$tran2->get_all_Exons};
  }
   
  my $chr_start = $exons[0]->start;
  my $chr_end   = $exons[-1]->end;
  my $seqname   = $exons[0]->seqname;
  $seqname =~/(\S+)\.\d+-\d+/;
  my $chr_name  = $1;
  my $slice_adaptor = $tran1->adaptor->db->get_SliceAdaptor;
  my $seq1;
  if ( $exons[0]->strand == 1 ){
    $seq1 = 
	$slice_adaptor->fetch_by_chr_start_end($chr_name,$chr_start-$padding,$chr_end+$padding);
  }
  else{
    my $seq = 
	$slice_adaptor->fetch_by_chr_start_end($chr_name,$chr_start-$padding,$chr_end+$padding);
    my $id  = "reverse_".$seq->display_id;
    $seq1   = $seq->invert;
    #$seq1->display_id($id);
    #$seq1->name($id);
    #$seq1->desc('');
  }

  my $chr_start2 = $exons2[0]->start;
  my $chr_end2   = $exons2[-1]->end;
  my $seqname2   = $exons2[0]->seqname;
  $seqname2 =~/(\S+)\.\d+-\d+/;
  my $chr_name2  = $1;
  my $slice_adaptor2 = $tran2->adaptor->db->get_SliceAdaptor;
  my $seq2;
  if ( $exons2[0]->strand == 1 ){
    $seq2 = $slice_adaptor2->fetch_by_chr_start_end($chr_name2,$chr_start2,$chr_end2);
  }
  else{
    my $seq = $slice_adaptor2->fetch_by_chr_start_end($chr_name2,$chr_start2,$chr_end2);
    my $id  = "reverse_".$seq->display_id;
    $seq2   = $seq->invert;
    #$seq2->display_id($id);
    #$seq2->name($id);
    #$seq2->desc('');
  }
  my $length1 = $seq1->length;
  my $length2 = $seq2->length;
  ############################################################
  # create database
  my $file = 'seq_'.$$.'.fa';
  my $database = "/tmp/".$file;
  open( DB_SEQ,">$database") || die("Could not open $database $!");
  
  my $seqout = Bio::SeqIO->new('-format' => 'Fasta',
			       '-fh'     => \*DB_SEQ);
  
  $seqout->write_seq($seq2);
  close( DB_SEQ );
  
  system("pressdb $database > /dev/null 2>&1");
  
  
  my $blastz =  Bio::EnsEMBL::Pipeline::Runnable::Blastz->new ('-query'     => $seq1,
							       '-database'  => $database,
							       '-options'   => 'B=0 C=2 K=2200',
							      );
    
  $blastz->run();

  my @featurepairs = $blastz->output();
  
  foreach my $fp (@featurepairs) {
      print STDERR $self->print_Feature($fp);
  }
  print STDERR "##############################\n";
  unlink ( $database );
  
  ############################################################
  # map the exons into the feature coordinates:
  
  my $copy1 = Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_clone_Transcript($tran1);
  my $copy2 = Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_clone_Transcript($tran2);
  my $transcript1 = $self->map_to_slice( $copy1,$seq1 );
  my $transcript2 = $self->map_to_slice( $copy2,$seq2 );
  
  my @exon_pairs;
  my ($exon_object_map, $exon_map) = $self->get_exon_pairs( $transcript1, $transcript2, \@featurepairs ,$coding_exons);
  
  return ($exon_object_map, $exon_map);
}

############################################################

sub map_to_slice{
  my ($self, $tran, $slice ) = @_;

  my $gene = Bio::EnsEMBL::Gene->new();
  $gene->add_Transcript($tran);
  $gene->transform($slice);
  my @transcripts = @{$gene->get_all_Transcripts};
  return $transcripts[0];
}

############################################################

sub get_exon_pairs{
  my ($self, $tran1, $tran2, $features, $coding ) = @_;
    
  my @features = @$features;
  my @exons1 = @{$self->get_Exons( $tran1,$coding)};
  my @exons2 = @{$self->get_Exons( $tran2,$coding)};
  
  my $verbose = 0;
  my %exon_map;

  my $start = 0;
 FEATURE:
  for ( my $k = 0; $k<scalar(@features); $k++){
      my $feat = $features->[$k];
      
      # which exons from tran1 overlap with feat?
    EXON1:
      for ( my $i=0; $i< scalar( @exons1); $i++ ){
	  my $exon1 = $exons1[$i];
	  if ( $exon1->start > $feat->end || $exon1->end < $feat->start ){
	      next EXON1;
	  }
	  if ( $exon1->start <= $feat->end && $exon1->end >= $feat->start ){
	      
	      # which exons from tran2 are in feat?
	    EXON2:
	      for (my $j=$start; $j< scalar( @exons2 ); $j++ ){
		  my $exon2 = $exons2[$j];
		  if ( $exon2->start > $feat->hend || $exon2->end < $feat->hstart ){
		      next EXON2;
		  }
		  if ( $exon2->start <= $feat->hend && $exon2->end >= $feat->hstart ){
		      
		      ############################################################
		      # potential pair - are these overlapping?
		      my $start1 = $exon1->start;
		      my $end1   = $exon1->end;
		      
		      my $start2 = $exon2->start;
		      my $end2   = $exon2->end;
		      
		      print STDERR "potential pair:  exon1($i): $start1-$end1  feature: (".
			  $feat->start."-".$feat->end."):(".
			      $feat->hstart."-".$feat->hend.")  exon2($j): $start2-$end2\n";
		      
		      my $cigar_string = $feat->cigar_string;
		      my @blocks = ( $cigar_string =~ /(\d*[MDI])/g );
		      #print STDERR "parsing string: ",join ( ",", @blocks ),"\n";
		      
		      my $start_1 = $exon1->start - $feat->start + 1;
		      my $end_1   = $exon1->end   - $feat->start + 1;
		      my $start_2 = $exon2->start - $feat->hstart + 1;
		      my $end_2   = $exon2->end   - $feat->hstart + 1;
		      
		      my $pos1 = 0;
		      my $pos2 = 0;
		      my @sub_blocks;
		
		      # flag when we are in an exon-match
		      my $in_exon = 0;
		      
		      ############################################################
		      # coordinates of the blocks on the sequences system:
		      # (s1,e1) is on the seq1 system
		      # (s2,e2) is on the seq2 system
		      my ($s1,$e1,$s2,$e2) = (0,$feat->start - 1,0,$feat->hstart - 1);
		      foreach my $block ( @blocks ){
			  my ($length) = ( $block =~ /^(\d*)/ );
			  $length =1 if $length eq "";
			  
			  if ( $block =~ /M$/ ){
			      $s1 = $e1 + 1;
			      $e1 = $s1 + $length - 1;
			      $s2 = $e2 + 1;
			      $e2 = $s2 + $length - 1 ;
			      print STDERR "block: $block\tfeat1: $s1-$e1\tfeat2: $s2-$e2\n" if $verbose;
			  }
			  if ( $block =~ /I$/ ){
			      $s1 = $e1 + 1;
			      $e1 = $s1 + $length - 1;
			      print STDERR "block: $block\tfeat1: $s1-$e1\n" if $verbose;
			  }
			  if ( $block =~ /D$/ ){
			      $s2 = $e2 + 1;
			      $e2 = $s2 + $length - 1;
			      print STDERR "block: $block\t\t\t\tfeat2: $s2-$e2\n" if $verbose;
			  }
			  
			  ############################################################
			  # match state
			  # simplest case: both exons coincide on a match block:
			  if ( $block =~ /M$/ ){
			      my ($start1,$end1,$start2,$end2);
			      if ( $exon1->start >= $s1 && $exon1->start <= $e1 ){
				  $start1 = $exon1->start - $s1 + 1;
				  print STDERR "exon1 starts at pos $start1\n" if $verbose;
			      }
			      if ( $exon1->end >= $s1 && $exon1->end <= $e1 ){
				  $end1 = $exon1->end - $s1 + 1;
				  print STDERR "exon1 ends at pos $end1\n" if $verbose;
			      }
			      if ( $exon2->start >= $s2 && $exon2->start <= $e2 ){
				  $start2 = $exon2->start - $s2 + 1;
				  print STDERR "exon2 starts at pos $start2\n" if $verbose;
			      }
			      if ( $exon2->end >= $s2 && $exon2->end <= $e2 ){
				  $end2 = $exon2->end - $s2 + 1;
				  print STDERR "exon2 ends at pos $end2\n" if $verbose;
			      }
			      
			      ############################################################
			      # check if they fall entirely on the M state
			      if ( $start1 && $start2 && $end1 && $end2 ){
				  
				  ############################################################
				  # do they overlap?
				  if ( $start1 > $end2 || $end1 < $start2 ){
				      print STDERR "they do not overlap - skipping\n";
				      next EXON2;
				  }
				  
				  
				  ############################################################
				  # exact match?
				  if ( $start1 == $start2 && $end1 == $end2 ){
				      print STDERR "MATCH\n";
				      
				      $exon_map{$i}{$j} = ($end2 - $start2 + 1 )."M";
				      print STDERR "exon_map($i)($j): ".$exon_map{$i}{$j}."\n";
				      $start = $j + 1;
				      next EXON1;
				  }
				  
				  ############################################################
				  # overlap
				  elsif( !( $start1>$end2 || $end1<$start2 ) ){
				      my $overlap = $self->min ($end1, $end2 ) 
					  - $self->max($start1,$start2) + 1; 
				      $exon_map{$i}{$j} = '';
				      my $left = $start2 - $start1;
				      if ($left){
					  $exon_map{$i}{$j} .=      $left."I" if $left>0;
					  $exon_map{$i}{$j} .= abs($left)."D" if $left<0;
				      }
				      $exon_map{$i}{$j} .= $overlap."M";
				      my $right = $end1 - $end2;
				      if ($right){
					  $exon_map{$i}{$j} .=      $right."I" if $right>0;
					  $exon_map{$i}{$j} .= abs($right)."D" if $right<0;
				      }
				      print STDERR "exon_map($i)($j): ".$exon_map{$i}{$j}."\n";
				  }
			      }
			      ############################################################
			      # if one has the start outside the match state
			      # or both are but we are at the first feature:
			      elsif( $end1 && $end2 && ( $start1 || $start2 || $k==0 ) ){ 
				  my $mismatch1 = ($s1 - $exon1->start);
				  my $mismatch2 = ($s2 - $exon2->start);
				  my $left = ($mismatch1 - $mismatch2 );
				  if ($left){
				      my $mismatch = ( $self->max($mismatch1,$mismatch2) -
						       $self->min($mismatch1,$mismatch2) );
				      $exon_map{$i}{$j} .=      $mismatch."I" if $left>0;
				      $exon_map{$i}{$j} .= abs($mismatch)."D" if $left<0;
				  }
				  # potentially alignable bases
				  # not aligned by blastz
				  my $extra_bases = max( 0, $self->min($mismatch1,$mismatch2) );
				  
				  #$exon_map{$i}{$j} .= $self->min($mismatch1,$mismatch2)."m"
				  #    if ( $self->min($mismatch1,$mismatch2) > 0);
				  
				  
				  # same end?
				  if ( $end1 == $end2 ){
				      $exon_map{$i}{$j} .= ($end2 + $extra_bases)."M";
				  }
				  else{
				      $exon_map{$i}{$j} .= 
					  ($self->min($end1,$end2) + $extra_bases)."M";
				      my $right = $end1 - $end2;
				      if ($right){
					  $exon_map{$i}{$j} .=      $right."I" if $right>0;
					  $exon_map{$i}{$j} .= abs($right)."D" if $right<0;
				      }
				  }
				  print STDERR "exon_map($i)($j): ".$exon_map{$i}{$j}."\n";
				  next EXON1;
				  $start = $j + 1;
			      }
			      ############################################################
			      # if both have the start outside the block
			      elsif( $end1 && $end2 && ! ( $start1 || $start2 ) ){
				  # same end?
				  if ( $end1 == $end2 ){
				      $exon_map{$i}{$j} .= $end2."M";
				  }
				  else{
				      $exon_map{$i}{$j} .= $self->min($end1,$end2)."M";
				      my $right = $end1 - $end2;
				      if ($right){
					  $exon_map{$i}{$j} .=      $right."I" if $right>0;
					  $exon_map{$i}{$j} .= abs($right)."D" if $right<0;
				      }
				  }
				  print STDERR "exon_map($i)($j): ".$exon_map{$i}{$j}."\n";
				  next EXON1;
				  $start = $j + 1;
			      }
			      ############################################################
			      # if one has the end outside the feature block
			      # or they are at the last feature
			      elsif ($start1 && $start2 
				     && 
				     ( $end1 || $end2 || $k==$#features ) ){
				  my $mismatch1 = $exon1->end - $e1;
				  my $mismatch2 = $exon2->end - $e2;
				  
				  # alignable bases not aligned by blastz
				  my $extra_bases = max(0, $self->min($mismatch1,$mismatch2) );
				  my $Match = min(
						  min($exon1->length,$exon2->length),
						  min($e1 - $exon1->start + 1,
						      $e2 - $exon2->start + 1
						      )
						  );
				  if ( $start1 == $start2 ){
				      $exon_map{$i}{$j} .= ($Match + $extra_bases)."M";
				  }
				  else{
				      my $left = $start1 - $start2;
				      if ($left){
					  $exon_map{$i}{$j} .=      $left."I" if $left>0;
					  $exon_map{$i}{$j} .= abs($left)."D" if $left<0;
				      }
				      $exon_map{$i}{$j} .= ($Match + $extra_bases)."M";
				  }
				  
				  
				  #$exon_map{$i}{$j} .= $self->min($mismatch1,$mismatch2)."m"
				  #    if ( $self->min($mismatch1,$mismatch2) > 0 );
				  
				  my $right = ($mismatch1 - $mismatch2 );
				  if ($right){
				      my $mismatch = ( $self->max($mismatch1,$mismatch2) -
						       $self->min($mismatch1,$mismatch2) );
				      $exon_map{$i}{$j} .=      $mismatch."I" if $right>0;
				      $exon_map{$i}{$j} .= abs($mismatch)."D" if $right<0;
				  }
				  print STDERR "exon_map($i)($j): ".$exon_map{$i}{$j}."\n";
			      }
			      ############################################################
			      # if both have the end outside the feature block
			      # or they are at the last feature
			      elsif ($start1 && $start2 && !( $end1 || $end2 ) ){
				  $in_exon = 1;
				  my $Match = min(
						  min($exon1->length,$exon2->length),
						  min($e1 - $exon1->start + 1,
						      $e2 - $exon2->start + 1
						      )
						  );
				  if ( $start1 == $start2 ){
				      $exon_map{$i}{$j} .= $Match."M";
				  }
				  else{
				      my $left = $start1 - $start2;
				      if ($left){
					  $exon_map{$i}{$j} .=      $left."I" if $left>0;
					  $exon_map{$i}{$j} .= abs($left)."D" if $left<0;
				      }
				      $exon_map{$i}{$j} .= $Match."M";
				  }
				  
				  print STDERR "exon_map($i)($j): ".$exon_map{$i}{$j}."\n";
			      }
			  }
			  
			  ############################################################
			  # insert state
			  if ( $block =~ /I$/ ){
			      my ($start1,$end1);
			      if ( $exon1->start >= $s1 && $exon1->start <= $e1 ){
				  $start1 = $exon1->start - $s1 + 1;
				  print STDERR "exon1 starts at pos $start1 in I-state\n";
			      }
			      if ( $exon1->end >= $s1 && $exon1->end <= $e1 ){
				  $end1 = $exon1->end - $s1 + 1;
				  print STDERR "exon1 ends at pos $end1 in I-state\n";
			      }
			      if ( !($start1 || $end1 ) && $in_exon ){
				  $exon_map{$i}{$j} .= ($e1 - $s1 + 1)."I";
				  print STDERR "exon_map($i)($j): ".$exon_map{$i}{$j}."\n";
			      }
			      elsif( $start1 && !$end1 ){
				  $exon_map{$i}{$j} .= ($e1 - $exon1->start + 1 )."I";
				  $in_exon = 1;
				  print STDERR "exon_map($i)($j): ".$exon_map{$i}{$j}."\n";
			      }
			      elsif( !$start1 && $end1 && $in_exon ){
				  $exon_map{$i}{$j} .= ($exon1->end - $s1 + 1 )."I";
				  print STDERR "exon_map($i)($j): ".$exon_map{$i}{$j}."\n";
			      }
			  }
			  ############################################################
			  # delete state
			  if ( $block =~ /D$/ ){
			      my ($start2,$end2);
			      if ( $exon2->start >= $s2 && $exon2->start <= $e2 ){
				  $start2 = $exon2->start - $s2 + 1;
				  print STDERR "exon2 starts at pos $start2 in D-state\n";
			      }
			      if ( $exon2->end >= $s2 && $exon2->end <= $e2 ){
				  $end2 = $exon2->end - $s2 + 1;
				  print STDERR "exon2 ends at pos $end2 in D-state\n";
			      }
			      
			      if ( !($start2 || $end2 ) && $in_exon ){
				  $exon_map{$i}{$j} .= ($e2 - $s2 + 1)."D";
				  print STDERR "exon_map($i)($j): ".$exon_map{$i}{$j}."\n";
			      }
			      elsif( $start2 && !$end2 ){
				  $exon_map{$i}{$j} .= ($e2 - $exon2->start + 1 )."D";
				  $in_exon = 1;
				  print STDERR "exon_map($i)($j): ".$exon_map{$i}{$j}."\n";
			      }
			      elsif( !$start2 && $end2 && $in_exon ){
				  $exon_map{$i}{$j} .= ($exon2->end - $s2 + 1 )."D";
				  print STDERR "exon_map($i)($j): ".$exon_map{$i}{$j}."\n";
			      }
			  }
		      }
		  }
	      } # end of EXON2
	  }
      }         # end of EXON1
      
  }             # end of FEATURE
  
  
  my $exon_object_map = Bio::EnsEMBL::Pipeline::GeneComparison::ObjectMap->new();
  print STDERR "pairs found:\n";
  foreach my $i ( keys %exon_map ){
      print STDERR "matches for exon $i\n";
      foreach my $j ( keys %{$exon_map{$i}} ){

	  my $match    = ( $exon_map{$i}{$j} =~ /(\d*M)/g );
	  my $mismatch = ( $exon_map{$i}{$j} =~ /(\d*[DI])/g );
	  my $flag = '';
	  if ( $match && !$mismatch ){
	      $exon_object_map->match( $exons1[$i], $exons2[$j], +1 );
	      $flag = "match";
	  }
	  elsif( $mismatch ){
	      $exon_object_map->match( $exons1[$i], $exons2[$j], -1 );
	      $flag = "mismatch";
	  }
	  print STDERR $self->exon_string($exons1[$i]).
	      "\t<---->\t".
		  $self->exon_string($exons2[$j])."\t".$exon_map{$i}{$j}."\t[$flag]\n";
	  
      }
  }
  my ($human_list,$mouse_list) = $self->get_exon_pair_alignment(\@exons1,\@exons2,\%exon_map);

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
  
  my $human_id = $tran1->stable_id || $tran1->dbID;
  my $mouse_id = $tran2->stable_id || $tran2->dbID;
  
  my $human_count = scalar(@exons1);
  my $mouse_count = scalar(@exons2);

  ############################################################
  # summary line
  print STDERR "TRANS_PAIR\t".
      "$human_id\thuman_exons:$human_count\thuman_miss_term_exons:$human_terminal_missing\thuman_miss_int_exons:$human_internal_missing\t".
	  "conserved_exons:$conserved\twith_same_length:$same_length\twith_same_phase:$same_phases\t".
	      "$mouse_id\tmouse_exons:$mouse_count\tmouse_miss_term_exons:$mouse_terminal_missing\tmouse_miss_int_exons:$mouse_internal_missing\n";
  
  ############################################################
  # summary line for exact matches
  if ( $human_count == $mouse_count && $conserved == $same_length && $human_count == $conserved ){
      my $cigars;
      foreach my $i ( keys %exon_map ){
	  foreach my $j ( keys %{$exon_map{$i}} ){
	      $cigars .= $exon_map{$i}{$j}."\t";
	  }
      }
      print STDERR "TRANS_EXACT_MATCH\t".
	  "$human_id\thuman_exons:$human_count\t".
	      "conserved_exons:$conserved\twith_same_length:$same_length\twith_same_phase:$same_phases\t".
		  "$mouse_id\tmouse_exons:$mouse_count\t".$cigars."\n";
  }
  ############################################################
  # summary line for semi_exact matches
  elsif ( $human_count == $mouse_count && $human_count == $conserved ){
      my $cigars;
      foreach my $i ( keys %exon_map ){
	  foreach my $j ( keys %{$exon_map{$i}} ){
	      $cigars .= $exon_map{$i}{$j}."\t";
	  }
      }
      print STDERR "TRANS_SEMI_MATCH\t".
	  "$human_id\thuman_exons:$human_count\t".
	      "conserved_exons:$conserved\twith_same_length:$same_length\twith_same_phase:$same_phases\t".
		  "$mouse_id\tmouse_exons:$mouse_count\t".$cigars."\n";
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

sub get_exon_pair_alignment{
  my ($self,$human_list, $mouse_list, $exon_map ) = @_;
  my %exon_map = %$exon_map;
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
  # last exons 
  if ( defined( $exon_map{$human_length-1}{$mouse_length-1} ) &&
       $exon_map{$human_length-1}{$mouse_length-1} ){
      pop @human_list;
      pop @mouse_list;
      my ( $human_list2, $mouse_list2) = 
	  $self->get_exon_pair_alignment( \@human_list, \@mouse_list, $exon_map );
      push ( @{$human_list2}, $human_last );
      push ( @{$mouse_list2}, $mouse_last );
      return ( $human_list2, $mouse_list2 );
  }
  ############################################################
  # last exon of the first list is paired-up with a gap
  elsif( !defined( $exon_map{$human_length-1}{$mouse_length-1} ) &&
	 defined(  $exon_map{$human_length-2}{$mouse_length-1} ) &&
	 $exon_map{$human_length-2}{$mouse_length-1} ){
      pop @human_list;
      my ( $human_list2, $mouse_list2) =
	  $self->get_exon_pair_alignment( \@human_list, \@mouse_list, $exon_map );
      push ( @{$human_list2}, $human_last );
      push ( @{$mouse_list2}, "gap" );
      return ( $human_list2, $mouse_list2 );
  }
  ############################################################
  # last exons of the second list is paired up with a gap
  else{
      pop @mouse_list;
      my ( $human_list2, $mouse_list2) =
	  $self->get_exon_pair_alignment( \@human_list, \@mouse_list, $exon_map );
      push ( @{$human_list2}, "gap" );
      push ( @{$mouse_list2}, $mouse_last );
      return ( $human_list2, $mouse_list2 );
  }
} 

############################################################




sub max{
    my ($self, $max, @values ) = @_;
    for (my $i=0; $i<@values; $i++ ){
	$max = $values[$i] if $values[$i]>$max;
    }
    return $max;
}

sub min{
    my ($self, $min, @values ) = @_;
    for (my $i=0; $i<@values; $i++ ){
	$min = $values[$i] if $values[$i]<$min;
    }
    return $min;
}



############################################################
# the target
############################################################

sub process_target{
    my ($self,$feat, $tran, $coding_exons) = @_;
    my $transcript_length;
    my @exons;

    if ( $coding_exons ){
      $transcript_length = length($tran->translateable_seq);
      @exons = sort { $a->length <=> $b->length } @{$tran->get_all_translateable_Exons};
    }
    else{
      $transcript_length = $tran->seq->length;
      @exons = sort { $a->length <=> $b->length } @{$tran->get_all_Exons};
    }
    
    my $min_exon_length = $exons[0]->length;
    
    my $is_spliced;
    
    my @clusters;
    my @cluster_starts;
    my @cluster_ends;
    my @features = sort{ $a->start <=> $b->start} @$feat;
 
    # create the first cluster
    my $count = 0;
    my $cluster = [];
    
    # start it off with the first feature
    my $first_feat = shift( @features );
    push (@$cluster, $first_feat);
    $cluster_starts[$count] = $first_feat->start;
    $cluster_ends[  $count] = $first_feat->end;
 
    # store the list of clusters
    push(@clusters,$cluster);
    
    ############################################################
    # loop over the rest of the features
  FEATURE:
    foreach my $f ( @features ){
	if (!($f->end < $cluster_starts[$count] || $f->start > $cluster_ends[$count])) {      
	    push(@$cluster,$f);
	    
	    # re-adjust size of cluster
	    if ($f->start < $cluster_starts[$count]) {
		$cluster_starts[$count] = $f->start;
	    }
	    if ($f->end  > $cluster_ends[$count]) {
		$cluster_ends[$count]   = $f->end;
	    }
	}
	else{
	    # else, start create a new cluster with this feature
	    $count++;
	    $cluster = [];
	    push (@$cluster, $f);
	    $cluster_starts[$count] = $f->start;
	    $cluster_ends[  $count] = $f->end;
	    
	    # store it in the list of clusters
	    push(@clusters,$cluster);
	}
    }

    ############################################################
    # check whether the transcript has one or more exons unaligned
    if ( scalar( @clusters ) == 1 ){
	$is_spliced = 0;
    }
    else{
      # compute the size of the 'gaps'
      my @gaps;
      $is_spliced = 0;
      for(my $i=0; $i<$#clusters-1; $i++){
	my $gap = $cluster_starts[$i+1] - $cluster_ends[$i] - 1;
	#print STDERR "gap: $gap, min_exon_length = $min_exon_length\n";
	if ( $gap >= $min_exon_length ){
	  $is_spliced = 1;
	  #print STDERR "is spliced\n";
	}
      }
    }
    
    ############################################################
    # calculate the coverage of the transcript
    my $feature_length = 0;
    for(my $i=0; $i<=$#clusters; $i++){
	#print STDERR "target cluster $i: $cluster_starts[$i] - $cluster_ends[$i]\n";
	$feature_length += $cluster_ends[$i] - $cluster_starts[$i] + 1;
    }
    my $coverage = 100*$feature_length/$transcript_length;
    #$self->print_exons_in_transcript($tran);
    return ($coverage,$is_spliced);
   
}

############################################################
# the query 
############################################################
sub process_query{
  my ($self,$feat, $qtran, $coding_exons) = @_;

  my $qtranscript_length;
  my @exons;
  if ( $coding_exons ){
    $qtranscript_length = length($qtran->translateable_seq);
    @exons = sort { $a->length <=> $b->length } @{$qtran->get_all_translateable_Exons};
  }
  else{
    $qtranscript_length = $qtran->seq->length;
    @exons = sort { $a->length <=> $b->length } @{$qtran->get_all_Exons};
  }
  
  my $min_exon_length = $exons[0]->length;
  
  my $is_spliced;
 
 my @clusters;
 my @cluster_hstarts;
 my @cluster_hends;
 my @features = sort{ $a->hstart <=> $b->hstart} @$feat;
 
 # create the first cluster
 my $count = 0;
 my $cluster = [];
  
 # start it off with the first feature
 my $first_feat = shift( @features );
 push (@$cluster, $first_feat);
 $cluster_hstarts[$count] = $first_feat->hstart;
 $cluster_hends[  $count] = $first_feat->hend;
 
 # store the list of clusters
 push(@clusters,$cluster);
 
 ############################################################
 # loop over the rest of the features
 FEATURE:
 foreach my $f ( @features ){
     if (!($f->hend < $cluster_hstarts[$count] || $f->hstart > $cluster_hends[$count])) {      
	 push(@$cluster,$f);
	 
	 # re-adjust size of cluster
	 if ($f->hstart < $cluster_hstarts[$count]) {
	     $cluster_hstarts[$count] = $f->hstart;
	 }
	 if ($f->hend  > $cluster_hends[$count]) {
	     $cluster_hends[$count] = $f->hend;
	 }
    }
     else{
	 # else, start create a new cluster with this feature
	 $count++;
	 $cluster = [];
	 push (@$cluster, $f);
	 $cluster_hstarts[$count] = $f->hstart;
	 $cluster_hends[$count]   = $f->hend;
	 
	 # store it in the list of clusters
	 push(@clusters,$cluster);
     }
 }

 ############################################################
 # check whether the transcript has one or more exons unaligned
 if ( scalar( @clusters ) == 1 ){
     $is_spliced = 0;
 }
 else{
     # compute the size of the 'gaps'
     my @gaps;
     $is_spliced = 0;
     for(my $i=0; $i<$#clusters-1; $i++){
	 my $gap = $cluster_hstarts[$i+1] - $cluster_hends[$i] - 1;
	 #print STDERR "gap: $gap, min_exon_length = $min_exon_length\n";
	 if ( $gap >= $min_exon_length ){
	     $is_spliced = 1;
	     #print STDERR "is spliced\n";
	 }
     }
 }
 
 ############################################################
 # calculate the coverage of the transcript
 my $feature_length = 0;
 for(my $i=0; $i<=$#clusters; $i++){
     #print STDERR "query cluster $i: $cluster_hstarts[$i] - $cluster_hends[$i]\n";
     $feature_length += $cluster_hends[$i] - $cluster_hstarts[$i] + 1;
 }
 my $coverage = sprintf "%.2f", 100*$feature_length/$qtranscript_length;
 #print STDERR "coverage = $feature_length / $qtranscript_length = $coverage\n";

 #$self->print_exons_in_transcript($qtran);
 
 return ($coverage,$is_spliced);
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
  my ($self,$human_t, $mouse_t, $gap_penalty, $coding_exons, $exon_pairs ) = @_;
  
  
  my @exons_pairs = @$exon_pairs if ($exon_pairs);

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
  
  foreach my $i ( 1..$human_length ){
    foreach my $j ( 1..$mouse_length ){
	($comparison_score{$human_exons[$i-1]}{$mouse_exons[$j-1]},
	 $features{$human_exons[$i-1]}{$mouse_exons[$j-1]}) =
	     $self->blast_Exons( $human_exons[$i-1], $mouse_exons[$j-1] );
	
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
    "$human_id\thuman_exons:$human_length\thuman_miss_term_exons:$human_terminal_missing\thuman_miss_int_exons:$human_internal_missing\t".
      "conserved_exons:$conserved\twith_same_length:$same_length\twith_same_phase:$same_phases\t".
	"$mouse_id\tmouse_exons:$mouse_length\tmouse_miss_term_exons:$mouse_terminal_missing\tmouse_miss_int_exons:$mouse_internal_missing\n";
  
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
	      print STDERR "score (".$human_list->[$i]->stable_id.",".$mouse_list->[$i]->stable_id.")=".$comparison_score{$human_list->[$i]}{$mouse_list->[$i]}."\n";
	      $score = $comparison_score{$human_list->[$i]}{$mouse_list->[$i]};
	      
	      foreach my $feat ( @{ $features{$human_list->[$i]}{$mouse_list->[$i]} } ){
		  $cigars = $feat->cigar_string ."|";
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

  for (my $i=0; $i< scalar(@exons); $i++ ){
    if ( $i>0 && $strand == 1 ){
      if ( $exons[$i]->start - $exons[$i-1]->end - 1 < 10 ){
	$exons[$i-1]->end($exons[$i]->end);
	next;
      }
    }
    if ( $i>0 && $strand == -1 ){
      if ( $exons[$i-1]->start - $exons[$i]->end - 1 < 10 ){
	$exons[$i-1]->start($exons[$i]->start);
	next;
      }
    }
    push (@newexons, $exons[$i] );
  }
  
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

sub blast_Exons{
  my ($self,$exon1, $exon2) =@_;
  
  my $id1;
  if ( $exon1->dbID ){
    $id1 = $exon1->stable_id || $exon2->dbID;
  }
  else{
    $id1 = $exon1;
  }
  
  my $id2;
  if ( $exon2->dbID ){
    $id2 = $exon2->stable_id || $exon2->dbID;
  }
  else{
    $id2 = $exon2;
  }
  
  my $seq1    = $exon1->seq;
  my $length1 = $seq1->length;
  unless ( $seq1->display_id ){
    $seq1->display_id($id1);
  }

  my $seq2    = $exon2->seq;
  my $length2 = $seq2->length;
  unless ( $seq2->display_id ){
    $seq2->display_id($id2);
  }
  
  my $min_length = $length1;
  if ( $length2 < $length1 ){
    $min_length = $length2;
  }
  my $word = 5;
  if ( 3*$word >= $min_length ){
    $word = int($min_length/3) - 1;
  }
  if ( $word < 1 ){
    return 0;
  }

  #print STDERR "word: $word\n";
  #print STDERR "comparing $id1: ".$seq1->seq."\n";
  #print STDERR "      and $id2: ".$seq2->seq."\n";


  ############################################################
  # create database
  my $file = 'seq_'.$$.'.fa';
  my $database = "/tmp/".$file;
  open( DB_SEQ,">$database") || die("Could not open $database $!");
  
  my $seqout = Bio::SeqIO->new('-format' => 'Fasta',
			       '-fh'     => \*DB_SEQ);
  
  $seqout->write_seq($seq2);
  close( DB_SEQ );
  
  system("pressdb $database > /dev/null 2>&1");
  
  ############################################################
  # Ian's parameters:
  my $options = "W=5 M=1 N=-1 Q=3 R=3";
  
  #my $options = 'altscore="* any na" altscore="any * na" S2=13';
  #$options .= " V=200 B=200 ";
  #$options .= " -nogap ";
  #$options .= " W=$word ";
  #print STDERR "options: $options\n";
  #my $options = 'V=200 B=200 altscore="* any na" altscore="any * na" W=4 E=0.01 E2=0.01 -nogap';
  #my $options = 'V=200 B=200 W=9 E=0.01 E2=0.01';
  my $blast =  
    Bio::EnsEMBL::Pipeline::Runnable::Blast->new ('-query'          => $seq1,
						  '-program'        => 'wutblastx',
						  #'-program'        => 'wublastn',
						  '-database'       => $database,
						  '-threshold_type' => "PVALUE",
						  #'-threshold'      => 1e-10,
						  '-options'        => $options,
						 );
  
  
  $blast->add_regex($file,'(\S+)');
  $blast->run();
  
  unlink( $database );
  
  my @featurepairs = $blast->output();
  
  if ( @featurepairs ){
      #my @pos_strand = grep { $_->strand == 1} @featurepairs;  
      #my @neg_strand = grep { $_->strand == -1} @featurepairs;  
      #foreach my $fp (sort{ $a->hstart <=> $b->hstart} @pos_strand) {
      #  print $fp->gffstring . "\n";
      #}
      #foreach my $fp (sort{ $a->hstart <=> $b->hstart} @neg_strand) {
      #  print $fp->gffstring . "\n";
      #}
      my $score = 0;
      foreach my $fp ( @featurepairs ){
	  $score += $fp->score;
      }
      $score = $score/scalar(@featurepairs);
      return ($score,\@featurepairs);
      
      #my @feat_by_score = sort { $b->score <=> $a->score } @featurepairs;
      #return $feat_by_score[0]->score;
  }
  else{
      #print STDERR "no hits for:\n";
      #print STDERR "$id1:\t".$seq1->seq."\n";
      #print STDERR "$id2:\t".$seq2->seq."\n";
      return (0);
  }
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

1;































