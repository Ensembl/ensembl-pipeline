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

    #print STDERR "\t$id1 length = $length1\n";
    #print STDERR "\t$id2 length = $length2\n";
    
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
      
      #foreach my $f ( @{$features[$i]} ){
      #	$self->print_Feature($f);
      #}

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
      #print STDERR "\tquery:$id1 coverage:$query_coverage spliced:$query_spliced\n";
      #print STDERR "\ttarget:$id2 coverage:$target_coverage spliced:$target_spliced\n";
    }
    #return ($max_score,$pair);
    $ave_score = 0 if $spliced;
    return ( $ave_score, $pair );
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

sub gap_penalty{
    my ($self,$value) = @_;
    if ( $value ){
	$self->{_gap_penalty} = $value;
    }
    return $self->{_gap_penalty};
}



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


1;
