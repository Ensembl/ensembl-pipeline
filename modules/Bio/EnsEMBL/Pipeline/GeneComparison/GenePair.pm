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
}

############################################################

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
  # check with te genomic alignment that the exons are
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
	  
	  my ($missing_terminal_exons, $exon_skipping, $all_conserved) = $self->blast_genomic_isoforms( $human_t, $mouse_t, $coding_exons, $gene_id1, $gene_id2);
	  
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
  


  
}




############################################################
#
#
############################################################



sub blast_genomic_isoforms{
  my ( $self,$tran1, $tran2, $coding_exons , $gene_id1, $gene_id2) = @_;
  
  my $padding = 500;
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
    $seq1->display_id($id);
    $seq1->name($id);
    $seq1->desc('');
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
    $seq2->display_id($id);
    $seq2->name($id);
    $seq2->desc('');
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
  
  #foreach my $fp (@featurepairs) {
  #    print STDERR $self->print_Feature($fp);
  #}
  #print STDERR "##############################\n";
  unlink ( $database );
  
  ############################################################
  # map the exons into the feature coordinates:
  
  my $copy1 = Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_clone_Transcript($tran1);
  my $copy2 = Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_clone_Transcript($tran2);
  my $transcript1 = $self->map_to_slice( $copy1,$seq1 );
  my $transcript2 = $self->map_to_slice( $copy2,$seq2 );
  
  my @exon_pairs;
  my ($exon_object_map, $exon_map) = $self->get_exon_pairs( $transcript1, $transcript2, \@featurepairs ,$coding_exons, $gene_id1, $gene_id2);
  
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
    my ($self, $tran1, $tran2, $features, $coding , $gene_id1, $gene_id2) = @_;
    
  my @features = @$features;
  my @exons1 = @{$self->get_Exons( $tran1,$coding)};
  my @exons2 = @{$self->get_Exons( $tran2,$coding)};
  
  my $verbose = 1;
  my %exon_map;
  my %exon_pointer_map;

  my $start = 0;
 FEATURE:
  for ( my $k = 0; $k<scalar(@features); $k++){
      my $feat = $features->[$k];
      
      # which exons from tran1 overlap with feat?
    EXON1:
      for ( my $i=0; $i< scalar( @exons1); $i++ ){
	my $exon1 = $exons1[$i];
	if ( $coding ){
	  if ( $i == 0 && $exon1->phase == -1 ){
	    $exon1->phase(0);
	  }
	  if ( $i == $#exons1 && $exon1->end_phase == -1 ){
	    $exon1->end_phase( ($exon1->phase + $exon1->length) %3 );
	  }  
	}
	if ( $exon1->start > $feat->end || $exon1->end < $feat->start ){
	  next EXON1;
	}
	if ( $exon1->start <= $feat->end && $exon1->end >= $feat->start ){
	    	 
	      # which exons from tran2 are in feat?
	    EXON2:
	      for (my $j=$start; $j< scalar( @exons2 ); $j++ ){
		  my $exon2 = $exons2[$j];
		  if ($coding){
		    if ( $j == 0 && $exon2->phase == -1 ){
		      $exon2->phase(0);
		    }
		    if ( $j == $#exons2 && $exon2->end_phase == -1 ){
		      $exon2->end_phase( ($exon2->phase + $exon2->length) %3 );
		    }
		  }
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
			      $feat->hstart."-".$feat->hend.")  exon2($j): $start2-$end2\n" if $verbose;
		      
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
		      
		      my $in_exon1 = 0;
		      my $in_exon2 = 0;
		      my $seen_exon1 = 0;
		      my $seen_exon2 = 0;
		      
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
			

			  if ( ( $seen_exon1 == 1 && $in_exon2 == 0 && $seen_exon2 == 0 )
			       ||
			       ( $seen_exon2 == 2 && $in_exon1 == 0 && $seen_exon1 == 0 )
			       ){
			      print STDERR "exons do not overlap - skipping\n" if $verbose;
			      delete $exon_map{$i}{$j};
			      next EXON2;
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
				#print STDERR "start1: $start1\n" if $start1;
				#print STDERR "start2: $start2\n" if $start2;
				#print STDERR "end1  : $end1\n" if $end1;
				#print STDERR "end2  : $end2\n" if $end2;
				
			      }
			      
			      $seen_exon1 = 1 if ( $end1 );
			      $seen_exon2 = 1 if ( $end2 );


			      ############################################################
			      # do they overlap?
			      if ( ( $start1 && $end2 && $start1 > $end2 )
				   || 
				   ( $start2 && $end1 && $end1 < $start2 )
				   ){
				  print STDERR "they do not overlap - skipping\n" if $verbose;
				  next EXON2;
			      }
			      

			      ############################################################
			      # check if they fall entirely on the M state
			      if ( $start1 && $start2 && $end1 && $end2 ){
				
				  ############################################################
				  # do they overlap?
				  if ( $start1 > $end2 || $end1 < $start2 ){
				      print STDERR "they do not overlap - skipping\n" if $verbose;
				      next EXON2;
				    }
				  
				  
				  ############################################################
				  # exact match?
				  if ( $start1 == $start2 && $end1 == $end2 ){
				      print STDERR "MATCH\n" if $verbose;
				      
				      $exon_map{$i}{$j} = ($end2 - $start2 + 1 )."M";
				      $exon_pointer_map{$exon1}{$exon2} = $exon_map{$i}{$j};
				      print STDERR "exon_map($i)($j): ".$exon_map{$i}{$j}."\n" if $verbose;
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
				      print STDERR "exon_map($i)($j): ".$exon_map{$i}{$j}."\n" if $verbose;
				      $exon_pointer_map{$exon1}{$exon2} = $exon_map{$i}{$j};
				      $start = $j + 1;
				      #next EXON1;
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
				print STDERR "exon_map($i)($j): ".$exon_map{$i}{$j}."\n" if $verbose;
				$exon_pointer_map{$exon1}{$exon2} = $exon_map{$i}{$j};
				$start = $j + 1;
				$in_exon1 = 0;
				$in_exon2 = 0;
				#next EXON1;
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
				  print STDERR "exon_map($i)($j): ".$exon_map{$i}{$j}."\n" if $verbose;
				  $exon_pointer_map{$exon1}{$exon2} = $exon_map{$i}{$j};
				  $start = $j + 1;
				  $in_exon1 = 0;
				  $in_exon2 = 0;
				  #next EXON1;
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
				  
				  
				  if ( $exon_map{$i}{$j} ){
				      $exon_pointer_map{$exon1}{$exon2} = $exon_map{$i}{$j};
				      print STDERR "exon_map($i)($j): ".$exon_map{$i}{$j}."\n" if $verbose;
				  }
				  
				  $in_exon1 = 1 unless ( $end1 );
				  $in_exon2 = 1 unless ( $end2 );
				  $in_exon1 = 0 if ( $end1 );
				  $in_exon2 = 0 if ( $end2 );
				  
			      }
			      ############################################################
			      # if both have the end outside the feature block
			      # or they are at the last feature
			      elsif ($start1 && $start2 && !( $end1 || $end2 ) ){
				$in_exon1 = 1;
				$in_exon2 = 1;
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
				  $exon_pointer_map{$exon1}{$exon2} = $exon_map{$i}{$j};
				  print STDERR "exon_map($i)($j): ".$exon_map{$i}{$j}."\n" if $verbose;
				  $in_exon1 = 1 unless ( $end1 );
				  $in_exon2 = 1 unless ( $end2 );
				}
			      ############################################################
			      # if only one ends here
			      elsif ( ( $end1 && !($start1 || $start2 || $end2 ) )
				      ||
				      ( $end2 && !($start1 || $start2 || $end1 ) )
				    ){
				my $mismatch;
				if ( $end1 ){
				  $mismatch = ($exon1->end - $s1 + 1)."I";
				}
				elsif( $end2 ){
				  $mismatch = ($exon2->end - $s2 + 1)."D";
				}
				$exon_map{$i}{$j} .= $mismatch;
				
				
				if ( $exon_map{$i}{$j} ){
				    $exon_pointer_map{$exon1}{$exon2} = $exon_map{$i}{$j};
				    print STDERR "exon_map($i)($j): ".$exon_map{$i}{$j}."\n" if $verbose;
				}
				$in_exon1 = 0 if ( $end1 );
				$in_exon2 = 0 if ( $end2 );
			      }
			      ############################################################
			      # if only one starts here
			      elsif ( ( $start1 && !($end1 || $start2 || $end2 ) )
				      ||
				      ( $start2 && !($start1 || $end1 || $end2 ) )
				    ){
				my $mismatch;
				if ( $start1 ){
				  $mismatch = ($e1 - $exon1->start + 1)."I";
				}
				elsif( $start2 ){
				  $mismatch = ($e2 - $exon2->start + 1)."D";
				}
				$exon_map{$i}{$j} .= $mismatch;
				$exon_pointer_map{$exon1}{$exon2} = $exon_map{$i}{$j};
				print STDERR "exon_map($i)($j): ".$exon_map{$i}{$j}."\n" if $verbose;
				$in_exon1 = 1 if ( $start1 );
				$in_exon2 = 1 if ( $start2 );
			      }
			      ############################################################
			      # if both starts and end are outside one only feature:
			      elsif( !($start1||$end1||$start2||$end2) && $k==0 && $k==$#features ){
				
				my $mismatch1 = ($s1 - $exon1->start);
				my $mismatch2 = ($s2 - $exon2->start);
				
				# mismatch
			        my $mismatch = $mismatch1 - $mismatch2;
			        if ($mismatch){
				  $exon_map{$i}{$j} .=      $mismatch."I" if $mismatch>0;
				  $exon_map{$i}{$j} .= abs($mismatch)."D" if $mismatch<0;
				}
			      
				# alignable bases not aligned by blastz
				my $extra_bases = $self->min($mismatch1,$mismatch2);
			        
				my $Match = $e1 - $s1 + 1;
			      
				# end alignable bases not aligned by blastz
				my $end_mismatch1 = ( $exon1->end - $e1 );
				my $end_mismatch2 = ( $exon2->end - $e2 );
				my $end_extra_bases = $self->min($end_mismatch1,$end_mismatch2);
				
				$exon_map{$i}{$j} .= ($Match + $extra_bases + $end_extra_bases)."M";
				
				my $end_mismatch = $end_mismatch1 - $end_mismatch2;
			        if ($end_mismatch){
				  $exon_map{$i}{$j} .=      $end_mismatch."I" if $end_mismatch>0;
				  $exon_map{$i}{$j} .= abs($end_mismatch)."D" if $end_mismatch<0;
				}
				print STDERR "exon_map($i)($j): ".$exon_map{$i}{$j}."\n" if $verbose;
				$exon_pointer_map{$exon1}{$exon2} = $exon_map{$i}{$j};
				$in_exon1 = 0;
				$in_exon2 = 0;
				$start = $j + 1;
				#next EXON1;
			      }
			      
			    }
			  
			  ############################################################
			  # insert state
			  if ( $block =~ /I$/ ){
			      my ($start1,$end1);
			      if ( $exon1->start >= $s1 && $exon1->start <= $e1 ){
				  $start1 = $exon1->start - $s1 + 1;
				  print STDERR "exon1 starts at pos $start1 in I-state\n" if $verbose;
			      }
			      if ( $exon1->end >= $s1 && $exon1->end <= $e1 ){
				  $end1 = $exon1->end - $s1 + 1;
				  print STDERR "exon1 ends at pos $end1 in I-state\n" if $verbose;
			      }

			      $seen_exon1 = 1 if $end1;

			      if ( !($start1 || $end1 ) && $in_exon1 ){
				  $exon_map{$i}{$j} .= ($e1 - $s1 + 1)."I";
				  print STDERR "exon_map($i)($j): ".$exon_map{$i}{$j}."\n" if $verbose;
			      }
			      elsif( $start1 && !$end1 ){
				  $exon_map{$i}{$j} .= ($e1 - $exon1->start + 1 )."I" if $verbose;
				  $in_exon1 = 1;
				  print STDERR "exon_map($i)($j): ".$exon_map{$i}{$j}."\n" if $verbose;
			      }
			      elsif( !$start1 && $end1 && $in_exon1 ){
				  $exon_map{$i}{$j} .= ($exon1->end - $s1 + 1 )."I";
				  print STDERR "exon_map($i)($j): ".$exon_map{$i}{$j}."\n" if $verbose;
				  $in_exon1 = 0;
				}
			      if ( $exon_map{$i}{$j} ){
				$exon_pointer_map{$exon1}{$exon2} = $exon_map{$i}{$j};
			      }
			    }
			  ############################################################
			  # delete state
			  if ( $block =~ /D$/ ){
			    my ($start2,$end2);
			    if ( $exon2->start >= $s2 && $exon2->start <= $e2 ){
				$start2 = $exon2->start - $s2 + 1;
			      print STDERR "exon2 starts at pos $start2 in D-state\n" if $verbose;
			  }
			    if ( $exon2->end >= $s2 && $exon2->end <= $e2 ){
			      $end2 = $exon2->end - $s2 + 1;
			      print STDERR "exon2 ends at pos $end2 in D-state\n" if $verbose;
			    }
			    
			    $seen_exon2 = 1 if $end2;

			      if ( !($start2 || $end2 ) && $in_exon2 ){
				$exon_map{$i}{$j} .= ($e2 - $s2 + 1)."D";
				print STDERR "exon_map($i)($j): ".$exon_map{$i}{$j}."\n" if $verbose;
			      }
			      elsif( $start2 && !$end2 ){
				$exon_map{$i}{$j} .= ($e2 - $exon2->start + 1 )."D";
				$in_exon2 = 1;
				print STDERR "exon_map($i)($j): ".$exon_map{$i}{$j}."\n" if $verbose;
			      }
			      elsif( !$start2 && $end2 && $in_exon2 ){
				$exon_map{$i}{$j} .= ($exon2->end - $s2 + 1 )."D";
				print STDERR "exon_map($i)($j): ".$exon_map{$i}{$j}."\n" if $verbose;
				$in_exon2 = 0;
			      }
			    if ( $exon_map{$i}{$j} ){
			      $exon_pointer_map{$exon1}{$exon2} = $exon_map{$i}{$j};
			    }
			  }
			}
		      
		      print STDERR "finished looking at i=$i j=$j\n";
		      if ( $exon_map{$i}{$j} ){
			my $s = $exon_map{$i}{$j};
			print STDERR "checking exon_map($i)($j) for matches: ".$s."\n";
			my $matches    = ( $s =~ /M/g );
			print STDERR "$matches matches found\n";
			unless ( $matches ){
			  print STDERR "rejecting exon_map($i)($j): ".$exon_map{$i}{$j}."\n" if $verbose;
			  delete $exon_map{$i}{$j};
			  
			  }
			
			# if there are previous exon_maps we reject this one:
		      }
		      if ( $exon_map{$i}{$j} ){
			for (my $m=0; $m< $i; $m++ ){
			  if ( $exon_map{$m}{$j} ){
			    print STDERR "exon($j) already taken - rejecting exon_map($i)($j): ".$exon_map{$i}{$j}."\n" if $verbose;
			    delete $exon_map{$i}{$j};
			    last;
			  }
			}
		      }
		      if ( $exon_map{$i}{$j} ){
			for (my $n=0; $n< $j; $n++ ){
			  if ( $exon_map{$i}{$n} ){
			    print STDERR "exon($i) already taken - rejecting exon_map($i)($j): ".$exon_map{$i}{$j}."\n" if $verbose;
			    delete $exon_map{$i}{$j};
			    last;
			  }
			}
		      }
		      
		      if ( $exon_map{$i}{$j} ){
			  $exon_pointer_map{$exon1}{$exon2} = $exon_map{$i}{$j};
		      }
		  }
	      } # end of EXON2
	  }
    }         # end of EXON1
      
  }             # end of FEATURE
    
    
    my $exon_object_map = Bio::EnsEMBL::Pipeline::GeneComparison::ObjectMap->new();
    #print STDERR "pairs found:\n";
  foreach my $i ( keys %exon_map ){
    #print STDERR "matches for exon $i\n";
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
      print STDERR "exon($i) ".$self->exon_string($exons1[$i]).
	"\t<---->\t".
	  "exon($j) " .$self->exon_string($exons2[$j])."\t".$exon_map{$i}{$j}."\t[$flag]\n";
    }
  }

  ############################################################
  # get alignment of exons:
  my ($human_list,$mouse_list) = $self->get_exon_pair_alignment(\@exons1,\@exons2,\%exon_map,\%exon_pointer_map);
  
  #my $alignment_object_map = Bio::EnsEMBL::Pipeline::GeneComparison::ObjectMap->new();
  my @alignment;

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
      my $cigar = $mouse_list->[$i]->length."D";
      push ( @alignment, $cigar );
      $human_missing++;
    }
    if ( $mouse_list->[$i] eq 'gap'){
      my $cigar = $human_list->[$i]->length."I";
      push( @alignment, $cigar );
      $mouse_missing++;
    }
    if ( !( $human_list->[$i] eq 'gap' || $mouse_list->[$i] eq 'gap') ){
      $conserved++;
      
      my $human_length = $human_list->[$i]->end - $human_list->[$i]->start + 1;
      my $mouse_length = $mouse_list->[$i]->end - $mouse_list->[$i]->start + 1;
      
      push( @alignment, $exon_pointer_map{$human_list->[$i]}{$mouse_list->[$i]} );
      
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

  my $alignment = join ' ',@alignment;


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
      print STDERR $human_string."\t<---->\t".$mouse_string."\t".$alignment[$i]."\n";
    }
  }

  ############################################################
  # summary line
  print STDERR "TRANS_PAIR\t".
      "$gene_id1\t$human_id\thuman_exons:$human_count\thuman_miss_term_exons:$human_terminal_missing\thuman_miss_int_exons:$human_internal_missing\t".
	  "conserved_exons:$conserved\twith_same_length:$same_length\twith_same_phase:$same_phases\t".
	      "$gene_id2\t$mouse_id\tmouse_exons:$mouse_count\tmouse_miss_term_exons:$mouse_terminal_missing\tmouse_miss_int_exons:$mouse_internal_missing $alignment\n";
  
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
	  "$gene_id1\t$human_id\thuman_exons:$human_count\t".
	      "conserved_exons:$conserved\twith_same_length:$same_length\twith_same_phase:$same_phases\t".
		  "$gene_id2\t$mouse_id\tmouse_exons:$mouse_count\t".$alignment."\n";
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
	  "$gene_id1\t$human_id\thuman_exons:$human_count\t".
	      "conserved_exons:$conserved\twith_same_length:$same_length\twith_same_phase:$same_phases\t".
		  "$gene_id2\t$mouse_id\tmouse_exons:$mouse_count\t".$alignment."\n";
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
  my ($self,$human_list, $mouse_list, $exon_map, $exon_pointer_map ) = @_;

  my $verbose = 1;
  my %exon_map = %$exon_map;
  my @human_list = @$human_list;
  my @mouse_list = @$mouse_list;
  
  print STDERR "human_list: ".$self->list_string(\@human_list)."\n" if $verbose;
  print STDERR "mouse_list: ".$self->list_string(\@mouse_list)."\n" if $verbose;
  
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
  
  print STDERR "before checking: exon_map[$human_length-1][$mouse_length-1] = ".$exon_map{$human_length-1}{$mouse_length-1}."\n" if $verbose;
  print STDERR "before checking: ".keys( %{ $exon_map{$human_length-1} } )."\n" if $verbose;
  
  ############################################################
  # last exons 
  if ( defined( $exon_map{$human_length-1}{$mouse_length-1} ) &&
       $exon_map{$human_length-1}{$mouse_length-1} ){
   
    print STDERR "exon_map[$human_length-1][$mouse_length-1] = ".$exon_map{$human_length-1}{$mouse_length-1}."\n" if $verbose;
    pop @human_list;
    pop @mouse_list;
    my ( $human_list2, $mouse_list2) = 
      $self->get_exon_pair_alignment( \@human_list, \@mouse_list, $exon_map );
    print STDERR "case1: human_list2: ".$self->list_string($human_list2).
      " mouse_list2: ".$self->list_string($mouse_list2)."\n" if $verbose;
    push ( @{$human_list2}, $human_last );
    push ( @{$mouse_list2}, $mouse_last );
    return ( $human_list2, $mouse_list2 );
  }
  ############################################################
  # last exon of the first list is paired-up with a gap
  elsif( !defined( $exon_map{$human_length-1}{$mouse_length-1} ) &&
	 !( keys( %{ $exon_map{$human_length-1} } ) )
       ){
    
    if ( defined $exon_map{$human_length-2}{$mouse_length-1} ){
      print STDERR "exon_map[$human_length-2][$mouse_length-1] = ".$exon_map{$human_length-2}{$mouse_length-1}."\n" if $verbose;
    }
    pop @human_list;
    my ( $human_list2, $mouse_list2) =
      $self->get_exon_pair_alignment( \@human_list, \@mouse_list, $exon_map );
    print STDERR "case2: human_list2: ".$self->list_string($human_list2).
	" mouse_list2: ".$self->list_string($mouse_list2)."\n" if $verbose;
    push ( @{$human_list2}, $human_last );
    push ( @{$mouse_list2}, "gap" );
    return ( $human_list2, $mouse_list2 );
  }
  ############################################################
  # last exons of the second list is paired up with a gap
  elsif( !defined( $exon_map{$human_length-1}{$mouse_length-1} ) &&
	 !( $self->mouse_exon_has_human_match( $exon_map, $mouse_length-1 ) )
       ){
	
    if ( defined $exon_map{$human_length-1}{$mouse_length-2} ){
      print STDERR "exon_map[$human_length-1][$mouse_length-2] = ".$exon_map{$human_length-1}{$mouse_length-2}."\n";
    }
    pop @mouse_list;
    my ( $human_list2, $mouse_list2) =
      $self->get_exon_pair_alignment( \@human_list, \@mouse_list, $exon_map );
    print STDERR "case3: human_list2: ".$self->list_string($human_list2).
	" mouse_list2: ".$self->list_string($mouse_list2)."\n" if $verbose;
    push ( @{$human_list2}, "gap" );
    push ( @{$mouse_list2}, $mouse_last );
    return ( $human_list2, $mouse_list2 );
  }
} 

############################################################

sub list_string{
    my ($self,$list) = @_;
    my $string;
    foreach my $l ( @$list ){
	if ( $l->isa('Bio::EnsEMBL::Exon' ) ){
	    $string .= $l->length."-";
	}
	else{
	    $string .= "gap-";
	}
    }
    return $string;
}

sub mouse_exon_has_human_match{
  my ($self,$map,$mouse_pos) = @_;
  my %exon_map = %$map;
  foreach my $i ( keys %exon_map ){
    foreach my $j ( keys %{$exon_map{$i}} ){
	print STDERR "mouse_pos=$mouse_pos, j=$j\n"; 
	if ( $j == $mouse_pos && defined( $exon_map{$i}{$mouse_pos} ) && $exon_map{$i}{$mouse_pos} ){
	    print STDERR "exon_map($i)($mouse_pos) = ".$exon_map{$i}{$mouse_pos}."\n";
	    return 1;
	}
    }
}
  return 0;
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
	if ( $i == $human_length && $human_exon->end_phase == -1 ){
	  $human_exon->end_phase( ($human_exon->phase + $human_exon->length) %3 );
	}
	if ( $j == 1 && $mouse_exon->phase == -1 ){
	  $mouse_exon->phase(0);
	}
	if ( $j == $mouse_length && $mouse_exon->end_phase == -1 ){
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
	#print STDERR "adding right $shift bp to exon: ".$exons[$i-1]->start."-".$exons[$i-1]->end."\n";
	
	$newexons[$c-1] = $exons[$i-1]->adjust_start_end(0,$shift);
	
	#print STDERR "new exon: ".$newexons[$c-1]->start."-".$newexons[$c-1]->end."\n";
	#$exons[$i-1]->end($exons[$i]->end);
	next;
      }
    }
    if ( $i>0 && $strand == -1 ){
      if ( $exons[$i-1]->start - $exons[$i]->end - 1 < 10 ){
	my $shift = $exons[$i-1]->start - $exons[$i]->start;
	#print STDERR "adding left $shift bp to exon: ".$exons[$i-1]->start."-".$exons[$i-1]->end."\n";
	# adjust_start_end() creates a new exon object!
	
	$newexons[$c-1] = $exons[$i-1]->adjust_start_end(0,$shift);
	
	#print STDERR "new exon: ".$newexons[$c-1]->start."-".$newexons[$c-1]->end."\n";
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

1;































