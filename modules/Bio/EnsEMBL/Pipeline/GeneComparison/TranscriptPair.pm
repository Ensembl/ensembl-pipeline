=head1 NAME

TranscriptPair

=head1 SYNOPSIS

Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptPair->blast_isoforms( $transcript1, $transcript2 );


=head1 DESCRIPTION


=head1 CONTACT

eae@sanger.ac.uk

=cut

# Let the code begin ...

package Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptPair;

use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Root;
use Bio::Seq;
use Bio::SeqIO;
use Bio::EnsEMBL::Pipeline::Runnable::Blast;
use Bio::EnsEMBL::Pipeline::Runnable::Blastz;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;

@ISA = qw(Bio::EnsEMBL::Root);

=head1 METHODS

=cut

#########################################################################


=head2 new()

=cut

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
#
# compare transcripts with blastn
#
############################################################


sub blast_isoforms{
    my ( $self,$tran1,$tran2, $coding_exons ) = @_;
    
    my $verbose = 1;

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
    #my $options = "-nogap W=5";
    
    # Ian's parameters:
    #my $options = "W=5 M=1 N=-1 Q=3 R=3";
    
    my $options = "W=5";
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

    my @featurepairs = $blast->output();

    
    #    my $blastz =  Bio::EnsEMBL::Pipeline::Runnable::Blastz->new ('-query'     => $seq1,
    #								  '-database'  => $database,
    #								  '-options'   => 'B=0 C=2 K=2200',
    #								 );
    
#    $blastz->run();
    
    #my @featurepairs = $blastz->output();
    


    #foreach my $fp (sort {$a->hstart <=> $b->hstart} @featurepairs) {
    #  print $fp->gffstring . "\n";
    #}


    unlink ( $database );

        
    #print STDERR "\t$id1 length = $length1\n";
    #print STDERR "\t$id2 length = $length2\n";
    
    ############################################################
    # separate by strands

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
    
    my $best_score = 0;
    my $best_features;
    for (my $i=0; $i<4; $i++ ){
      unless ( $features[$i] && @{$features[$i]} ){
	next;
      }

      ############################################################
      # compute the score
      my $score = 0;
      foreach my $fp (sort {$a->hstart <=> $b->hstart} @{$features[$i]}) {
        $score += $fp->score;
        #print $fp->gffstring . "\n";
      }
      if ( $score > $best_score ){
        $best_score    = $score;
        $best_features = $features[$i];
      }
     }
    
     my $coverage = 0;
     my $spliced = 0;

    
    if ( $best_features ){
      
      ############################################################
      # calculate coverage
      # we use query/target as in feature pairs the target=seqname and query=hseqname
      my ($query_coverage,  $query_spliced)  = 
	$self->process_query( $best_features, $tran2 , $coding_exons);
      my ($target_coverage, $target_spliced) = 
	$self->process_target( $best_features, $tran1, $coding_exons);
      $coverage = ( $query_coverage + $target_coverage )/2;
      
      if ($verbose){
	foreach my $f ( sort { $a->start <=> $b->start } @{$best_features} ){
	  $self->print_Feature($f);
	}
      }
      if ( $query_spliced || $target_spliced ){
	$spliced = 1;
      }
    
      ############################################################
      # calculate the perc id
      my $perc_id = 0;
      foreach my $f ( @$best_features ){
	$perc_id += $f->percent_id;
      }
      $perc_id = sprintf "%.2f", ( $perc_id/scalar(@$best_features) );
      
      print STDERR "\tquery:$id1 coverage:$query_coverage spliced:$query_spliced\n";
      print STDERR "\ttarget:$id2 coverage:$target_coverage spliced:$target_spliced\n";
      print STDERR "\taveraged percent id: $perc_id\n";
      $best_score = 0 if ( $target_coverage + $query_coverage < 100 );
    }
    
    $best_score = 0 if $spliced;
    return ( $best_score, $best_features );
  }


############################################################
# the target
############################################################

sub process_target{
    my ($self,$feat, $tran, $coding_exons, $protein) = @_;
    my $transcript_length;
    my @exons;

    if ( $coding_exons || $protein ){
      $transcript_length = length($tran->translateable_seq);
      @exons = sort { $a->length <=> $b->length } @{$tran->get_all_translateable_Exons};
    }
    else{
      $transcript_length = $tran->seq->length;
      @exons = sort { $a->length <=> $b->length } @{$tran->get_all_Exons};
    }
    my $min_exon_length = $exons[0]->length;
    if ( $protein ){
      if ( $tran->translateable_seq=~/TAA$|TGA$|TAG$/i ){
	$transcript_length = ($transcript_length - 3)/3;
      }
      else{
	$transcript_length /= 3;
      }
      $min_exon_length /= 3;
    }  
        
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
      for(my $i=0; $i<$#clusters; $i++){
	my $gap = $cluster_starts[$i+1] - $cluster_ends[$i] - 1;
	print STDERR "gap: $gap, min_exon_length = $min_exon_length\n";
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
    #print STDERR "feature length   : $feature_length\n";
    #print STDERR "transcript length: $transcript_length\n";
    my $coverage = sprintf "%.2f", 100*$feature_length/$transcript_length;
    #$self->print_exons_in_transcript($tran);
    return ($coverage,$is_spliced);
   
}

############################################################
# the query 
############################################################
sub process_query{
  my ($self,$feat, $qtran, $coding_exons, $protein) = @_;

  my $qtranscript_length;
  my @exons;
  if ( $coding_exons || $protein){
    $qtranscript_length = length($qtran->translateable_seq);
    @exons = sort { $a->length <=> $b->length } @{$qtran->get_all_translateable_Exons};
  }
  else{
    $qtranscript_length = $qtran->seq->length;
    @exons = sort { $a->length <=> $b->length } @{$qtran->get_all_Exons};
  }
  my $min_exon_length = $exons[0]->length;
  if ( $protein ){
    if ( $qtran->translateable_seq=~/TAA$|TGA$|TAG$/i ){
      $qtranscript_length = ($qtranscript_length - 3)/3;
    }
    else{
      $qtranscript_length /= 3;
    }
    $min_exon_length /= 3;
  }  


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
   for(my $i=0; $i<$#clusters; $i++){
     my $gap = $cluster_hstarts[$i+1] - $cluster_hends[$i] - 1;
     print STDERR "gap: $gap, min_exon_length = $min_exon_length\n";
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
  my $seqname = $f->seqname;
  my ($chr_start,$chr_end);
  if ( $seqname =~ /(\S+)\.(\d+)-(\d+)/ ){
    $chr_start = $2 + $f->start - 1;
    $chr_end   = $2 + $f->end   - 1;
  }
  else{
    $chr_start = $f->start;
    $chr_end   =  $f->end;
  }
  
  my ($chr_hstart, $chr_hend);
  my $hseqname = $f->hseqname;
  if ( $hseqname =~ /(\S+)\.(\d+)-(\d+)/ ){
    $chr_hstart = $2 + $f->hstart - 1;
    $chr_hend   = $2 + $f->hend   - 1;
  }
  else{
    $chr_hstart =  $f->hstart;
    $chr_hend   =  $f->hend;
  }
  my $string =
    $f->seqname."\t".
      $chr_start."-".$chr_end."\t".
	$f->start."-".$f->end."\t".
	  ($f->end - $f->start + 1)."\t".
	    $f->strand."\t".
	      $f->hseqname."\t".
		$chr_hstart."-".$chr_hend."\t".
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

sub gap_penalty{
    my ($self,$value) = @_;
    if ( $value ){
	$self->{_gap_penalty} = $value;
    }
    return $self->{_gap_penalty};
}

############################################################


sub blast_genomic_isoforms{
  my ( $self,$tran1, $tran2, $coding_exons , $gene_id1, $gene_id2) = @_;
  
  my $verbose = 1;

  #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($tran1);
  #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($tran2);

  my $padding = 2000;
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
  print STDERR "1: start $chr_start, end: $chr_end\n";
  my $seqname;
  foreach my $exon ( @exons ){
    $seqname = $exon->seqname;
    last if $seqname;
  }
  $seqname =~/(\S+)\.\d*-\d*/;
  my $chr_name  = $1;
  my $slice_adaptor = $tran1->adaptor->db->get_SliceAdaptor;
  my ($seq1,$slice1);
  if ( $exons[0]->strand == 1 ){
      $slice1 = $slice_adaptor
	  ->fetch_by_chr_start_end($chr_name,$chr_start-$padding,$chr_end+$padding);
      #$seq1 = $slice1;
      $seq1 = $slice1->get_repeatmasked_seq(['RepeatMask'],1);
  }
  else{
      my $seq = $slice_adaptor
	  ->fetch_by_chr_start_end($chr_name,$chr_start-$padding,$chr_end+$padding);
      
      my $id  = "reverse_".$seq->display_id;
      $slice1   = $seq->invert;
      #$seq1     = $slice1;
      $seq1     = $slice1->get_repeatmasked_seq(['RepeatMask'],1);
      $seq1->display_id($id);
      $seq1->id($id);
      #$seq1->name($id);
      $seq1->desc('');
  }
  
  my $chr_start2 = $exons2[0]->start;
  my $chr_end2   = $exons2[-1]->end;
  print STDERR "2: start $chr_start2, end: $chr_end2\n";

  my $seqname2;
  foreach my $exon ( @exons2 ){
    $seqname2 = $exon->seqname;
    last if $seqname2;
  }
  $seqname2 =~/(\S+)\.\d*-\d*/;
  my $chr_name2  = $1;
  my $slice_adaptor2 = $tran2->adaptor->db->get_SliceAdaptor;
  my ($seq2,$slice2);
  if ( $exons2[0]->strand == 1 ){
      $slice2 = $slice_adaptor2->fetch_by_chr_start_end($chr_name2,$chr_start2,$chr_end2);
      #$seq2 = $slice2;
      $seq2 = $slice2->get_repeatmasked_seq(['RepeatMask'],1); 
  }
  else{
    my $seq = $slice_adaptor2->fetch_by_chr_start_end($chr_name2,$chr_start2,$chr_end2);
    my $id  = "reverse_".$seq->display_id;
    $slice2 = $seq->invert;
    $seq2   = $slice2->get_repeatmasked_seq(['RepeatMask'],1);
    #$seq2   = $slice2;
    
    $seq2->display_id($id);
    $seq2->id($id);
    #$seq2->name($id);
    $seq2->desc('');
  }
  my $length1 = $seq1->length;
  my $length2 = $seq2->length;
  print STDERR "length1 = $length1\n";
  print STDERR "length2 = $length2\n";
  
  #print STDERR "seq1: ".$seq1->subseq(1,1000);
  #print STDERR "seq2: ".$seq2->subseq(1,1000);
  
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
  
  
  my $blastz =  Bio::EnsEMBL::Pipeline::Runnable::Blastz
    ->new ('-query'     => $seq1,
	   '-database'  => $database,
	   '-options'   => 'B=0 C=2 K=1000 T=1',
	  );
  

  $blastz->run();
  
  my @featurepairs = $blastz->output();

  if ( $verbose ){
    foreach my $fp (@featurepairs) {
      print STDERR $self->print_Feature($fp);
    }
  }
  #print STDERR "##############################\n";
  unlink ( $database );
  unlink( $database."csq" );
  unlink( $database."nhd" );
  unlink( $database."ntb" );
  
  ############################################################
  # map the exons into the feature coordinates:
  
  my $copy1 = Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_clone_Transcript($tran1);
  my $copy2 = Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_clone_Transcript($tran2);
  my $transcript1 = $self->map_to_slice( $copy1,$slice1 );
  my $transcript2 = $self->map_to_slice( $copy2,$slice2 );
  
  my @exon_pairs;
  my ($exon_object_map, $exon_map) = $self->get_exon_pairs( $transcript1, $transcript2, \@featurepairs ,$coding_exons, $gene_id1, $gene_id2);
  
  #print STDERR "got exon pairs\n";

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

############################################################
# this method get the exons from a transcript
# all of them or only coding ones.
# It correct teh exon coordinates to
# bridge over small frameshifts

sub get_Exons{
  my ( $self, $trans , $coding) = @_;
  my @exons;
  my @newexons;
  my $strand = $trans->start_Exon->strand;

  my $verbose = 1;

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

sub get_exon_pairs{
    my ($self, $tran1, $tran2, $features, $coding , $gene_id1, $gene_id2) = @_;
    
    my @features = @$features;
    print STDERR "number of features: ".scalar(@features)."\n";

    my @exons1 = @{$self->get_Exons( $tran1,$coding)};
    my @exons2 = @{$self->get_Exons( $tran2,$coding)};
    
    my $verbose = 1;
    print STDERR "get_exon_pairs()\n" if $verbose;
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
	    if ( $i == $#exons1 ){
		#print STDERR "changing phase of exon1($i) from ".$exon1->end_phase." to "
		#    . (($exon1->phase + $exon1->length) %3 )."\n";
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
		    if ( $j == $#exons2 ){
			#print STDERR "changing phase of exon2($j) from ".$exon2->end_phase." to "
			#    . (($exon2->phase + $exon2->length) %3 )."\n";
			$exon2->end_phase( ($exon2->phase + $exon2->length) %3 );
		      }
		  }
		  if ( $exon2->start > $feat->hend || $exon2->end < $feat->hstart ){
		    next EXON2;
		  }
		  if ( $exon2->start <= $feat->hend && $exon2->end >= $feat->hstart ){
		    
		    ############################################################
		    # check whether these exons has been matched already
		    # for the time being, we avoid 2-to-1/1-to-2 cases
		    if ( keys( %{ $exon_map{$i} } ) ){
		      print STDERR "exon1($i) has been paired-up already - skipping\n" if $verbose;
		      next EXON1;
		    }
		    if ( $self->mouse_exon_has_human_match( \%exon_map, $j ) ){
		      print STDERR "exon2($j) has been already paired-up - skipping\n" if $verbose;
		      $start = $j+1;
		      next EXON2;
		    }



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
				      
				      $start = $j + 1;
				      # if there are previous exon_maps we reject this one:
				  
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
					  
					  print STDERR "exon_map($i)($j): ".$exon_map{$i}{$j}."\n" if $verbose;
					  $exon_pointer_map{$exon1}{$exon2} = $exon_map{$i}{$j};
				      }
				      
				      next EXON1;
				  }
			      }
			      ############################################################
			      # if one has the start outside the match state
			      # or both are but we are at the first feature:
			      elsif( $end1 && $end2 && ( $start1 || $start2 || $k==0 ) ){ 
				my $mismatch1 = ($s1 - $exon1->start);
				my $mismatch2 = ($s2 - $exon2->start);
				
				if ( $block eq $blocks[0] ){
				    my $left = ($mismatch1 - $mismatch2 );
				    if ($left){
					my $mismatch = ( $self->max($mismatch1,$mismatch2) -
							 $self->min($mismatch1,$mismatch2) );
					$exon_map{$i}{$j} .=      $mismatch."I" if $left>0;
					$exon_map{$i}{$j} .= abs($mismatch)."D" if $left<0;
				    }
				}
				# potentially alignable bases
				# not aligned by blastz
				my $extra_bases = max( 0, $self->min($mismatch1,$mismatch2) );
				
				$extra_bases = 0 unless ( $block eq $blocks[0] );

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
				  print STDERR "exon_map($i)($j): ".$exon_map{$i}{$j}."\n" if $verbose;
				}
				$start = $j + 1;
				$in_exon1 = 0;
				$in_exon2 = 0;
				next EXON1;
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
				  $start = $j + 1;
				  $in_exon1 = 0;
				  $in_exon2 = 0;
				  
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
				    print STDERR "exon_map($i)($j): ".$exon_map{$i}{$j}."\n" if $verbose;
				    $exon_pointer_map{$exon1}{$exon2} = $exon_map{$i}{$j};
				  }



				  next EXON1;
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
				  
				  if ( $block eq $blocks[-1] ){
				      my $right = ($mismatch1 - $mismatch2 );
				      if ($right){
					  my $mismatch = ( $self->max($mismatch1,$mismatch2) -
							   $self->min($mismatch1,$mismatch2) );
					  $exon_map{$i}{$j} .=      $mismatch."I" if $right>0;
					  $exon_map{$i}{$j} .= abs($mismatch)."D" if $right<0;
				      }
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
			      # and we are in the only exons
			      elsif( !($start1||$end1||$start2||$end2) 
				     && 
				     $k==0 && $k==$#features 
				     &&
				     ( $i==0 && $i== $#exons1 )
				     &&
				     ( $j==0 && $j== $#exons2 )
				     ){
				
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
				next EXON1;
			      }
			      
			  }
			  
			  ############################################################
			  # insert state
			  elsif ( $block =~ /I$/ ){
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
			  elsif ( $block =~ /D$/ ){
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
		      
		      print STDERR "finished looking at i=$i j=$j\n" if $verbose;
		      if ( $exon_map{$i}{$j} ){
			my $s = $exon_map{$i}{$j};
			print STDERR "checking exon_map($i)($j) for matches: ".$s."\n" if $verbose;
			my $matches    = ( $s =~ /M/g );
			print STDERR "$matches matches found\n" if $verbose;
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
    
    
    if ($verbose){
	#print STDERR "pairs found:\n";
	foreach my $i ( keys %exon_map ){
	    #print STDERR "matches for exon $i\n";
	    foreach my $j ( keys %{$exon_map{$i}} ){
		
		my $match    = ( $exon_map{$i}{$j} =~ /(\d*M)/g );
		my $mismatch = ( $exon_map{$i}{$j} =~ /(\d*[DI])/g );
		my $flag = '';
		if ( $match && !$mismatch ){
		    $flag = "match";
		}
		elsif( $mismatch ){
		    $flag = "mismatch";
		}
		print STDERR "exon($i) ".$self->exon_string($exons1[$i]).
	    "\t<---->\t".
		"exon($j) " .$self->exon_string($exons2[$j])."\t".$exon_map{$i}{$j}."\t[$flag]\n";
	    }
	}
    }

  ############################################################
  # get alignment of exons:
  my ($human_list,$mouse_list) = $self->get_exon_pair_alignment(\@exons1,\@exons2,\%exon_map,\%exon_pointer_map);
    
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
      my $mouse_phase = $mouse_list->[$i]->phase;
      my $mouse_end_phase = ( $mouse_length + $mouse_phase ) %3;
      my $human_phase = $human_list->[$i]->phase;
      my $human_end_phase = ( $human_length + $human_phase ) %3;
      if ( $mouse_phase == $human_phase
	   &&
	   $mouse_end_phase == $human_end_phase ){
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


sub get_exon_pair_alignment{
  my ($self,$human_list, $mouse_list, $exon_map, $exon_pointer_map ) = @_;
  
  my $verbose = 1;
  print STDERR "get_exon_pair_alignment()\n" if $verbose;
  my %exon_map = %$exon_map;
  my @human_list = @$human_list;
  my @mouse_list = @$mouse_list;

  print STDERR "human_list: ".scalar(@human_list)."\n" if $verbose;
  print STDERR "mouse_list: ".scalar(@mouse_list)."\n" if $verbose;
  
  
  my $human_pos = 0;
  my $mouse_pos = 0;
  my @human_aligned;
  my @mouse_aligned;
  
  while ( $human_pos <= $#human_list && $mouse_pos <= $#mouse_list ){
      
      print STDERR "human_pos: $human_pos\tmouse_pos: $mouse_pos\n" if $verbose;

      if ( @human_aligned && @mouse_aligned ){
	  print STDERR "human_list: ".$self->list_string(\@human_aligned)."\n" if $verbose;
	  print STDERR "mouse_list: ".$self->list_string(\@mouse_aligned)."\n" if $verbose;
      }
      
      if ( defined( $exon_map{$human_pos}{$mouse_pos} ) 
	   &&
	   $exon_map{$human_pos}{$mouse_pos} ){
	  
	  push (@human_aligned, $human_list[$human_pos] );
	  push (@mouse_aligned, $mouse_list[$mouse_pos] );
	  $human_pos++;
	  $mouse_pos++;
	  next;
	  
	  #if ( scalar( keys( %{ $exon_map{$human_pos} } ) ) == 1 
	  #	   &&
	  #	   !( $self->mouse_exon_has_another_human_match( \%exon_map, $human_pos, $mouse_pos ) )
	  #	 ){
	  #		
	  #}
	  
      }
      elsif ( !defined( $exon_map{$human_pos}{$mouse_pos} ) 
	      &&
	      $self->mouse_exon_has_another_human_match( \%exon_map, $human_pos, $mouse_pos ) 
	      ){
	  push( @human_aligned, $human_list[$human_pos] );
	  push( @mouse_aligned, "gap" );
	  $human_pos++;
	  next;
      }
      elsif ( !defined( $exon_map{$human_pos}{$mouse_pos} )
	      &&
	      keys( %{ $exon_map{$human_pos} } )
	      ){
	  push( @human_aligned, "gap" );
	  push( @mouse_aligned, $mouse_list[$mouse_pos] );
	  $mouse_pos++;
	  next;
      }
      else{
	  push( @human_aligned, "gap" );
	  push( @mouse_aligned, $mouse_list[$mouse_pos] );
	  push( @human_aligned, $human_list[$human_pos] );
	  push( @mouse_aligned, "gap" );
	  $human_pos++;
	  $mouse_pos++;
	  next;
      }
  }
  if ( $human_pos < $#human_list ){
      for ( my $i = $human_pos; $i<= $#human_list; $i++ ){
	  push( @human_aligned, $human_list[$i] );
	  push( @mouse_aligned, "gap" );
      }
      $human_pos = $#human_list;
  }

  if ( $mouse_pos < $#mouse_list ){
    for ( my $i = $mouse_pos; $i<= $#mouse_list; $i++ ){
      push( @mouse_aligned, $mouse_list[$i] );
      push( @human_aligned, "gap" );
    }
    $mouse_pos = $#mouse_list;
  }
  if ( @human_aligned && @mouse_aligned ){
      print STDERR "returning human_list: ".$self->list_string(\@human_aligned)."\n" if $verbose;
      print STDERR "returning mouse_list: ".$self->list_string(\@mouse_aligned)."\n" if $verbose;
    }
  return( \@human_aligned, \@mouse_aligned );
} 

############################################################

# not used anymore, it is only valid for 1-to-1 correspondances

sub old_get_exon_pair_alignment{
  my ($self,$human_list, $mouse_list, $exon_map, $exon_pointer_map ) = @_;

  my $verbose = 0;
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
	 !( $self->mouse_exon_has_human_match( \%exon_map, $mouse_length-1 ) )
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

############################################################

sub mouse_exon_has_human_match{
  my ($self,$map,$mouse_pos) = @_;
  my %exon_map = %$map;
  foreach my $i ( keys %exon_map ){
    foreach my $j ( keys %{$exon_map{$i}} ){
	#print STDERR "mouse_pos=$mouse_pos, j=$j\n"; 
	if ( $j == $mouse_pos && defined( $exon_map{$i}{$mouse_pos} ) && $exon_map{$i}{$mouse_pos} ){
	    #print STDERR "exon_map($i)($mouse_pos) = ".$exon_map{$i}{$mouse_pos}."\n";
	  return 1;
	}
    }
}
  return 0;
}

############################################################


sub mouse_exon_has_another_human_match{
  my ($self,$map,$human_pos, $mouse_pos) = @_;
  my %exon_map = %$map;
  foreach my $i ( keys %exon_map ){
    next unless $i>$human_pos;
    foreach my $j ( keys %{$exon_map{$i}} ){
      next unless $j == $mouse_pos;
      if ( defined $exon_map{$i}{$j} && $exon_map{$i}{$j} ){
	  #print STDERR "exon_map($i)($mouse_pos) = ".$exon_map{$i}{$mouse_pos}."\n";
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

############################################################

sub min{
    my ($self, $min, @values ) = @_;
    for (my $i=0; $i<@values; $i++ ){
	$min = $values[$i] if $values[$i]<$min;
    }
    return $min;
}

############################################################




############################################################
#
# Other methods
#
############################################################

=head2 put_Transcripts()

  function to include one or more transcripts in the cluster.
  Useful when creating a cluster. It takes as argument an array of transcripts, it returns nothing.

=cut

sub put_Transcripts {
    my ($self, @new_transcripts)= @_;
    
    if ( !$new_transcripts[0]->isa('Bio::EnsEMBL::Transcript') ){
	$self->throw( "Can't accept a [ $new_transcripts[0] ] instead of a Bio::EnsEMBL::Transcript");
    }
    
    push ( @{ $self->{'_transcript_array'} }, @new_transcripts );
    
}

#########################################################################

=head2 get_Transcripts()

  it returns the array of transcripts in the GeneCluster object

=cut
    
sub get_Transcripts {
    my $self = shift @_;
    return $self->{'_transcript_array'};
}

############################################################

sub score{
    my ($self,$score) = @_;
    if ( $score ){
	$self->{_score} = $score;
    }
    return $self->{_score};
}

#########################################################################

=head2 to_String()

  it returns a string containing the information about the transcripts in the TranscriptCluster object

=cut

sub to_String {
  my $self = shift @_;
  my $data='';
  foreach my $tran ( @{ $self->{'_transcript_array'} } ){
    my @exons = @{$tran->get_all_Exons};
    my $id;
    if ( $tran->stable_id ){
      $id = $tran->stable_id;
    }
    else{
      $id = $tran->dbID;
    }
 
    $data .= sprintf "Id: %-16s"             , $id;
    $data .= sprintf "Contig: %-21s"         , $exons[0]->contig->id;
    $data .= sprintf "Exons: %-3d"           , scalar(@exons);
    my ($start, $end) = $self->_get_start_end($tran);
    $data .= sprintf "Start: %-9d"           , $start;
    $data .= sprintf "End: %-9d"             , $end;
    $data .= sprintf "Strand: %-3d"          , $exons[0]->strand;
    $data .= sprintf "Exon-density: %3.2f\n", $self->exon_Density($tran);
  }
  return $data;
}

#########################################################################

sub exon_Density{
  my ($self, $transcript) = @_;  
  my $density;
  my $exon_span;
  my @exons = @{$transcript->get_all_Exons};
  @exons = sort { $a->start <=> $b->start } @exons;
  my $transcript_length = $exons[$#exons]->end - $exons[0]->start;
  foreach my $exon ( @exons ){
    $exon_span += $exon->length;
  }
  $density = $exon_span/$transcript_length;
  return $density;
}

=head2 _get_start_end()

 function to get the start and end positions - written as one method
 for efficiency

=cut

sub _get_start_end {
  my ($self, $transcript) = @_;
  my $start;
  my $end;
 
  my $start_Exon = $transcript->start_Exon;
  my $end_Exon = $transcript->end_Exon;
 
  if ($start_Exon->strand == 1) {
    $start = $start_Exon->start;
    $end   = $end_Exon->end;
  } else {
    $end   = $start_Exon->end;
    $start = $end_Exon->start;
  }
  return ($start, $end);
}    


############################################################

sub blast_CDSs{
    my ( $self,$tran1,$tran2) = @_;
    
    my $verbose = 1;

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
    eval{
      $seq1    = $tran1->translate;
      unless ( $seq1->display_id ){
	$seq1->display_id($id1);
      }
      $seq2    = $tran2->translate;
      unless ( $seq2->display_id ){
	$seq2->display_id($id2);
      }
    };
    unless ( $seq1 && $seq2 ){
      $self->warn("could not retrieve all the proteins");
      return (0,undef);
    }
    
    my $length1 = $seq1->length;
    my $length2 = $seq2->length;
    print STDERR "\t$id1 length = $length1\n";
    print STDERR "\t$id2 length = $length2\n";
        
    
    ############################################################
    # create database
    my $file = 'seq_'.$$.'.fa';
    my $database = "/tmp/".$file;
    open( DB_SEQ,">$database") || die("Could not open $database $!");
    
    my $seqout = Bio::SeqIO->new('-format' => 'Fasta',
				 '-fh'     => \*DB_SEQ);
    
    $seqout->write_seq($seq2);
    close( DB_SEQ );
    
    system("setdb $database > /dev/null 2>&1");
    
    my $options = "W=5 -warnings";
    my $blast =  
      Bio::EnsEMBL::Pipeline::Runnable::Blast->new ('-query'          => $seq1,
						    '-program'        => 'wublastp',
						    '-database'       => $database,
						    '-threshold_type' => "PVALUE",
						    #'-threshold'      => 1e-10,
						    '-options'        => $options,
						   );
    
    $blast->add_regex($file,'(\S+)');
    $blast->run();
    
    my @featurepairs = $blast->output();
    #foreach my $fp (sort {$a->hstart <=> $b->hstart} @featurepairs) {
    #  print $fp->gffstring . "\n";
    #}
    

    unlink ( $database );

        
    ############################################################
    # separate by strands

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
    
    my $best_score = 0;
    my $best_features;
    for (my $i=0; $i<4; $i++ ){
      unless ( $features[$i] && @{$features[$i]} ){
	next;
      }
      
      ############################################################
      # compute the score
      my $score = 0;
      foreach my $fp (sort {$a->hstart <=> $b->hstart} @{$features[$i]}) {
        $score += $fp->score;
        #print $fp->gffstring . "\n";
      }
      if ( $score > $best_score ){
        $best_score    = $score;
        $best_features = $features[$i];
      }
     }
    
    my $coverage = 0;
    my $spliced = 0;
    my $perc_id=0;
    my ($query_coverage, $query_spliced);
    my ($target_coverage, $target_spliced);
     if ( $best_features ){

      ############################################################
      # calculate coverage
      # we use query/target as in feature pairs the target=seqname and query=hseqname
       my $coding_exons = 1;
       ($query_coverage,  $query_spliced)  = 
	 $self->process_query( $best_features, $tran2 , $coding_exons,1);
       ($target_coverage, $target_spliced) = 
	 $self->process_target( $best_features, $tran1, $coding_exons,1);
       $coverage = ( $query_coverage + $target_coverage )/2;
       
       if ($verbose){
	 foreach my $f ( sort { $a->start <=> $b->start } @{$best_features} ){
	   $self->print_Feature($f,1);
	 }
       }
       if ( $query_spliced || $target_spliced ){
	 $spliced = 1;
       }
       
      foreach my $f ( @$best_features ){
	$perc_id += $f->percent_id;
      }
      $perc_id = sprintf "%.2f", ( $perc_id/scalar(@$best_features) );
      
      print STDERR "\tquery:$id1 coverage:$query_coverage spliced:$query_spliced\n";
      print STDERR "\ttarget:$id2 coverage:$target_coverage spliced:$target_spliced\n";
      print STDERR "\taveraged percent id: $perc_id\n";
      
      $best_score = 0 if ( $target_coverage + $query_coverage < 100 );
    }
    
    $best_score = 0 if $spliced;
    
    return ( $best_score, $best_features, $target_coverage, $query_coverage, $perc_id );
  }

############################################################

############################################################
# this method aligns two transcript sequences
# which are not attached to any ensembl transcript object

sub blast_unmapped_transcripts{
  my ( $self,$seq1,$seq2) = @_;
  
  my $verbose = 1;
  
  # query
  my $id1 = $seq1->display_id;
  my $id2 = $seq2->display_id;
  my $length1 = $seq1->length;
  my $length2 = $seq2->length;
  #print STDERR "\tquery : $id1 length = $length1\n";
  #print STDERR "\ttarget: $id2 length = $length2\n";
  
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
  #my $options = "-nogap W=5";
  
  # Ian's parameters:
  #my $options = "W=5 M=1 N=-1 Q=3 R=3";
  
  #my $options = "W=5";
  my $options = "";
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
  
  my @featurepairs = $blast->output();
  unlink ( $database );
        
  print STDERR "\t$id1 length = $length1\n";
  print STDERR "\t$id2 length = $length2\n";
  
  ############################################################
  # separate by strands
  
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
  
  my $best_score = 0;
  my $best_features;
  for (my $i=0; $i<4; $i++ ){
    unless ( $features[$i] && @{$features[$i]} ){
      next;
    }
    
    ############################################################
    # compute the score
    my $score = 0;
    foreach my $fp (sort {$a->hstart <=> $b->hstart} @{$features[$i]}) {
      $score += $fp->score;
      #print $fp->gffstring . "\n";
    }
    if ( $score > $best_score ){
      $best_score    = $score;
      $best_features = $features[$i];
    }
  }   
  my $coverage = 0;
  my $spliced = 0;
  my $perc_id=0;
  my ($query_coverage, $query_spliced, $max_query_gap);
  my ($target_coverage, $target_spliced, $max_target_gap);
  

  ############################################################
  # must check whether there are big splits in the alignment
  if ( $best_features ){
    
    ############################################################
    # calculate coverage
    # we use query/target as in feature pairs the target=seqname and query=hseqname
    ($query_coverage,  $query_spliced, $max_query_gap)  = 
      $self->process_unmapped_query( $best_features, $seq2);
    ($target_coverage, $target_spliced, $max_target_gap) = 
      $self->process_unmapped_target( $best_features, $seq1);
    $coverage = ( $query_coverage + $target_coverage )/2;
    
    if ($verbose){
      foreach my $f ( sort { $a->start <=> $b->start } @{$best_features} ){
	$self->print_Feature($f);
      }
    }
    if ( $query_spliced || $target_spliced ){
      $spliced = 1;
    }
    
    ############################################################
    # calculate the perc id
    foreach my $f ( @$best_features ){
      $perc_id += $f->percent_id;
    }
    $perc_id = sprintf "%.2f", ( $perc_id/scalar(@$best_features) );
    
    print STDERR "\tquery:$id2 coverage:$query_coverage spliced:$query_spliced\n";
    print STDERR "\ttarget:$id1 coverage:$target_coverage spliced:$target_spliced\n";
    print STDERR "\taveraged percent id: $perc_id\n";
    $best_score = 0 if ( $target_coverage + $query_coverage < 100 );
  }
  
  $best_score = 0 if $spliced;
  return ( $best_score, $best_features, $target_coverage, $max_target_gap, $query_coverage, $max_query_gap, $perc_id );
}


############################################################
# the unmapped target
############################################################

sub process_unmapped_target{
    my ($self,$feat, $seq, $coding_exons, $protein) = @_;
    print STDERR "target is: ".$seq->display_id."\n";
    my $length = $seq->length;
    my $min_exon_length = 50;
    if ( $protein ){
      $min_exon_length = 16;
    }
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
      #print STDERR "feature: ".$f->start."-".$f->end."\n";
      #print STDERR "cluster: ".$cluster_starts[$count]."-".$cluster_ends[$count]."\n";
      
      if (!($f->end < $cluster_starts[$count] || $f->start > $cluster_ends[$count])) {      
	push(@$cluster,$f);
	
	# re-adjust size of cluster
	if ($f->start < $cluster_starts[$count]) {
	  $cluster_starts[$count] = $f->start;
	}
	if ($f->end  > $cluster_ends[$count]) {
	  $cluster_ends[$count]   = $f->end;
	}
	#print STDERR "added - reset cluster: ".$cluster_starts[$count]."-".$cluster_ends[$count]."\n";
      }
      else{
	# else, start create a new cluster with this feature
	$count++;
	$cluster = [];
	push (@$cluster, $f);
	$cluster_starts[$count] = $f->start;
	$cluster_ends[  $count] = $f->end;
	#print STDERR "create new cluster:".$cluster_starts[$count]."-".$cluster_ends[$count]."\n";
	
	# store it in the list of clusters
	push(@clusters,$cluster);
      }
    }
    
    my $max_gap = 0;
    ############################################################
    # check whether the transcript has one or more exons unaligned
    print STDERR "In process_unmapped_target(): number of clusters: ".scalar( @clusters )."\n";
    if ( scalar( @clusters ) == 1 ){
      print STDERR "only one cluster\n";
      $is_spliced = 0;
    }
    else{
      # compute the size of the 'gaps'
      my @gaps;
      $is_spliced = 0;
      for (my $i=0; $i<$#clusters; $i++){
	#print STDERR "cluster[$i]  : ".$cluster_starts[$i]."-".$cluster_ends[$i]."\n";
	#print STDERR "cluster[$i+1]: ".$cluster_starts[$i+1]."-".$cluster_ends[$i+1]."\n";
	my $gap = $cluster_starts[$i+1] - $cluster_ends[$i] - 1;
	print STDERR "gap: $gap, min_exon_length = $min_exon_length\n";
	if ( $gap >= $min_exon_length ){
	  $is_spliced = 1;
	  print STDERR "is spliced\n";
	}
	if ($gap > $max_gap ){
	  $max_gap = $gap;
	}
      }
    }
    
    ############################################################
    # calculate the coverage of the transcript
    my $feature_length = 0;
    for(my $i=0; $i<=$#clusters; $i++){
      #print STDERR "target cluster $i: $cluster_starts[$i] - $cluster_ends[$i] (".( $cluster_ends[$i] - $cluster_starts[$i] + 1).")\n";
      $feature_length += $cluster_ends[$i] - $cluster_starts[$i] + 1;
      #print STDERR "feature length: $feature_length\n";
    }
    print STDERR "feature length   : $feature_length\n";
    print STDERR "transcript length: $length\n";
    my $coverage = sprintf "%.2f", 100*$feature_length/$length;
    #$self->print_exons_in_transcript($tran);
    return ($coverage,$is_spliced,$max_gap);
}

############################################################
# the unmapped_query 
############################################################
sub process_unmapped_query{
  my ($self,$feat, $seq, $coding_exons, $protein) = @_;
  print STDERR "query is: ".$seq->display_id."\n";
  my $transcript_length = $seq->length;
  my $min_exon_length = 50;
  if ( $protein ){
    $min_exon_length = 16;
  }  
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
  
  my $max_gap = 0;
  ############################################################
  # check whether the transcript has one or more exons unaligned
  if ( scalar( @clusters ) == 1 ){
    #print STDERR "only one cluster\n";
    $is_spliced = 0;
  }
  else{
    # compute the size of the 'gaps'
    $is_spliced = 0;
    for(my $i=0; $i<$#clusters; $i++){
      my $gap = $cluster_hstarts[$i+1] - $cluster_hends[$i] - 1;
      print STDERR "gap: $gap, min_exon_length = $min_exon_length\n";
      if ( $gap >= $min_exon_length ){
	$is_spliced = 1;
	#print STDERR "is spliced\n";
      }
      if ( $gap > $max_gap ){
	$max_gap = $gap;
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
  my $coverage = sprintf "%.2f", 100*$feature_length/$transcript_length;
  #print STDERR "coverage = $feature_length / $qtranscript_length = $coverage\n";
  
  #$self->print_exons_in_transcript($qtran);
  
  return ($coverage,$is_spliced, $max_gap);
}

############################################################


############################################################

sub blast_unmapped_proteins{
    my ( $self,$seq1,$seq2) = @_;
    
    my $verbose = 1;
    my $id1 = $seq1->display_id;
    my $id2 = $seq2->display_id;

    unless ( $seq1 && $seq2 ){
      $self->warn("could not retrieve all the proteins");
      return (0,undef);
    }
    
    my $length1 = $seq1->length;
    my $length2 = $seq2->length;
    print STDERR "\t$id1 length = $length1\n";
    print STDERR "\t$id2 length = $length2\n";
    
    
    ############################################################
    # create database
    my $file = 'seq_'.$$.'.fa';
    my $database = "/tmp/".$file;
    open( DB_SEQ,">$database") || die("Could not open $database $!");
    
    my $seqout = Bio::SeqIO->new('-format' => 'Fasta',
				 '-fh'     => \*DB_SEQ);
    
    $seqout->write_seq($seq2);
    close( DB_SEQ );
    
    system("setdb $database > /dev/null 2>&1");
    
    my $options = " -warnings";
    my $blast =  
      Bio::EnsEMBL::Pipeline::Runnable::Blast->new ('-query'          => $seq1,
						    '-program'        => 'wublastp',
						    '-database'       => $database,
						    '-threshold_type' => "PVALUE",
						    #'-threshold'      => 1e-10,
						    '-options'        => $options,
						   );
    
    $blast->add_regex($file,'(\S+)');
    $blast->run();
    
    my @featurepairs = $blast->output();
    #foreach my $fp (sort {$a->hstart <=> $b->hstart} @featurepairs) {
    #  print $fp->gffstring . "\n";
    #}
    

    unlink ( $database );

        
    ############################################################
    # separate by strands

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
    
    my $best_score = 0;
    my $best_features;
    for (my $i=0; $i<4; $i++ ){
      unless ( $features[$i] && @{$features[$i]} ){
	next;
      }
      
      ############################################################
      # compute the score
      my $score = 0;
      foreach my $fp (sort {$a->hstart <=> $b->hstart} @{$features[$i]}) {
        $score += $fp->score;
        #print $fp->gffstring . "\n";
      }
      if ( $score > $best_score ){
        $best_score    = $score;
        $best_features = $features[$i];
      }
     }
    
    my $coverage = 0;
    my $spliced = 0;
    my $perc_id=0;
    my ($query_coverage, $query_spliced, $max_query_gap);
    my ($target_coverage, $target_spliced, $max_target_gap);
    if ( $best_features ){
      
      ############################################################
      # calculate coverage
      # we use query/target as in feature pairs the target=seqname and query=hseqname
      my $coding_exons = 1;
      ($query_coverage,  $query_spliced, $max_query_gap)  = 
	$self->process_unmapped_query( $best_features, $seq2 , $coding_exons,1);
      ($target_coverage, $target_spliced, $max_target_gap) = 
	$self->process_unmapped_target( $best_features, $seq1, $coding_exons,1);
      $coverage = ( $query_coverage + $target_coverage )/2;
      
      if ($verbose){
	foreach my $f ( sort { $a->start <=> $b->start } @{$best_features} ){
	  $self->print_Feature($f,1);
	}
      }
      if ( $query_spliced || $target_spliced ){
	$spliced = 1;
      }
      
      foreach my $f ( @$best_features ){
	$perc_id += $f->percent_id;
      }
      $perc_id = sprintf "%.2f", ( $perc_id/scalar(@$best_features) );
      
      print STDERR "\ttarget:$id1 coverage:$target_coverage spliced:$query_spliced\n";
      print STDERR "\tquery :$id2 coverage:$query_coverage spliced:$target_spliced\n";
      print STDERR "\taveraged percent id: $perc_id\n";
      
      $best_score = 0 if ( $target_coverage + $query_coverage < 100 );
    }
    
    $best_score = 0 if $spliced;
    
    return ( $best_score, $best_features, $target_coverage, $max_target_gap, $query_coverage, $max_query_gap, $perc_id );
  }

############################################################

1;
