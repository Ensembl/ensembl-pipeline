#!/usr/local/ensembl/bin/perl

use strict;  
use Bio::EnsEMBL::DBSQL::DBAdaptor;
#use Bio::EnsEMBL::Pipeline::GeneComparison::GeneComparison;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::GeneComparison::ObjectMap;
use Bio::EnsEMBL::Pipeline::GeneComparison::GenePair;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Bio::EnsEMBL::Pipeline::Runnable::Blast;
use Getopt::Long;


# human_db
my $human_dbname = 'homo_sapiens_core_13_31';
my $human_dbhost = 'ecs2f';
my $human = 'Homo sapiens';

# human_dnadb
my $human_dnadbname = 'homo_sapiens_core_13_31';
my $human_dnadbhost = 'ecs2f';

# mouse_db
my $mouse_dbname = 'mus_musculus_core_13_30';
my $mouse_dbhost = 'ecs2f';
my $mouse = 'Mus musculus';

# mouse_dnadb
my $mouse_dnadbname = 'mus_musculus_core_13_30';
my $mouse_dnadbhost = 'ecs2f';

# compara_db
my $compara_dbname = 'ensembl_compara_13_1';
my $compara_dbhost = 'ecs2f';


my $from_file;
my $gap_penalty;

# options
&GetOptions( 
	    'from_file'     => \$from_file,
	    'gap_penalty:n' => \$gap_penalty,
	   );

unless( $gap_penalty ){
  $gap_penalty = -100;
}

my $human_dnadb;
my $mouse_dnadb;
my $mouse_db;
my $human_db;
my $compara_db;

my %gene_pair;

############################################################
# connect to the databases 
$human_dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $human_dnadbhost,
						  -user  => 'ensro',
						  -dbname=> $human_dnadbname,
						 );

$human_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $human_dbhost,
					       -user  => 'ensro',
					       -dbname=> $human_dbname,
					       -dnadb => $human_dnadb,
					      );


$mouse_dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $mouse_dnadbhost,
						  -user  => 'ensro',
						  -dbname=> $mouse_dnadbname,
						 );

$mouse_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $mouse_dbhost,
					       -user  => 'ensro',
					       -dbname=> $mouse_dbname,
					       -dnadb => $mouse_dnadb,
					      );


my $human_adaptor = $human_db->get_GeneAdaptor;
my $mouse_adaptor = $mouse_db->get_GeneAdaptor;

my $gene_pair = Bio::EnsEMBL::Pipeline::GeneComparison::GenePair->new();


unless( $from_file ){

  print STDERR "Using compara_db\n";
  #my $compara_config =  '/nfs/acari/eae/ensembl/ensembl-compara/modules/Bio/EnsEMBL/Compara/Compara.conf';
  $compara_db = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(
							     -user      => 'ensro',
							     -dbname    => $compara_dbname,
							     -host      => $compara_dbhost,
							     #-conf_file => $compara_config,
							    );
  
  $compara_db->add_db_adaptor( $human_db);
  $compara_db->add_db_adaptor( $mouse_db);
  
  #my $mouse_db = $compara_db->get_db_adaptor($target_species,'NCBIM30');
  #my $human_db  = $compara_db->get_db_adaptor($focus_species ,'NCBI31');
  
  
  my $homol_adaptor = $compara_db->get_HomologyAdaptor;
  my @human_ids     = $homol_adaptor->list_stable_ids_from_species('Homo_sapiens');
  
  foreach my $human_id ( @human_ids ){
    
    my @mouse_homologs = 
      $homol_adaptor->fetch_homologues_of_gene_in_species('Homo_sapiens',$human_id,'Mus musculus');
    
    foreach my $homology ( @mouse_homologs ){
      
      print STDERR "comparing $human_id and ".$homology->stable_id."\n";
      my $human_gene = $human_adaptor->fetch_by_stable_id($human_id,1);
      my $mouse_gene = 
	  $mouse_adaptor->fetch_by_stable_id($homology->stable_id,1);
      $gene_pair->compare( $human_gene, $mouse_gene );
      push (@{$gene_pair{$human_id}},$homology->stable_id);
    }  
  }
}

if ( $from_file ){
  while(<>){
    chomp;
    my @entries = split;
    print STDERR "\n-------------------------------------------------------------\n\n";
    print STDERR "comparing $entries[0] and $entries[1]\n";

    my $human_gene = $human_adaptor->fetch_by_stable_id($entries[0],1);
    my $mouse_gene = $mouse_adaptor->fetch_by_stable_id($entries[1],1);
    $gene_pair->compare( $human_gene, $mouse_gene );
    push( @{$gene_pair{$entries[0]}}, $entries[1] );
  }   
}

############################################################

sub compare{
  my ($human_gene,$mouse_gene) = @_;
  
  my @human_transcripts = @{$human_gene->get_all_Transcripts};

  my @mouse_transcripts = @{$mouse_gene->get_all_Transcripts};
  
  
  ############################################################
  # we make a pair only if the transcripts align with gaps no longer than the smallest exon
  
  ############################################################
  # for the time being, take the best match for each human transcript
  # later on we can apply the stable marriage algorithm
  my $object_map = Bio::EnsEMBL::Pipeline::GeneComparison::ObjectMap->new();
  my @transcript_matches;
  foreach my $human_t ( @human_transcripts ){
    foreach my $mouse_t ( @mouse_transcripts ){
      print STDERR "blasting isoforms\n";
      my ($score,$pair) = &blast_isoforms( $human_t, $mouse_t );
      if ( $score && $pair ){
	$object_map->match($human_t, $mouse_t, $score );
      }
    }
  }
  my $best_pairs_object = $object_map->stable_marriage;
  
  ############################################################
  # pairs created:
  print STDERR "pairs created: ".scalar($best_pairs_object->list1)."\n";
  foreach my $element1 ( $best_pairs_object->list1 ){
    foreach my $partner ( $best_pairs_object->partners( $element1 ) ){
      # there should be only one
      print STDERR "Pair with score: ".$best_pairs_object->score( $element1, $partner )."\n";
      Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($element1);
      Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($partner);
    }
  }
  
  ############################################################
  # compare the exons for each pair
  print STDERR "comparing exons\n";
  foreach my $element1 ( $best_pairs_object->list1 ){
    foreach my $partner ( $best_pairs_object->partners( $element1 ) ){
      # there should be only one
      #print STDERR "Pair with score: ".$best_pairs_object->score( $element1, $partner )."\n";
      #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($element1);
      #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($partner);
      #my @score_matrix = 
      &compare_Exons( $element1, $partner, $gap_penalty);
      #&print_alignment( $element1, $partner, \@score_matrix);

    }
  }

}

############################################################


sub blast_isoforms{
    my ( $tran1,$tran2 ) = @_;
    
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
    
    print STDERR "comparing $id1 and $id2\n";
    
    my $seq1    = $tran1->seq;
    my $length1 = $seq1->length;
    unless ( $seq1->display_id ){
	$seq1->display_id($id1);
    }
    
    my $seq2    = $tran2->seq;
    my $length2 = $seq2->length;
    unless ( $seq2->display_id ){
	$seq2->display_id($id2);
    }
    
    
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
    my $blast =  Bio::EnsEMBL::Pipeline::Runnable::Blast->new ('-query'     => $seq1,
							       '-program'   => 'wublastn',
							       '-database'  => $database,
							       -threshold_type => "PVALUE",
							       '-threshold' => 1e-10,
							       #'-filter'    => $filter,
							       '-options'   => 'V=1000000'
							       );
    
    $blast->add_regex($file,'(\S+)');
    $blast->run();
    
    my @featurepairs = $blast->output();
    
    foreach my $fp (sort {$a->hstart <=> $b->hstart} @featurepairs) {
	&print_Feature($fp);
	#print $fp->gffstring . "\n";
    }
    print STDERR "$id1 length = $length1\n";
    print STDERR "$id2 length = $length2\n";
    
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
    for (my $i=0; $i<4; $i++ ){
      unless ( $features[$i] && @{$features[$i]} ){
	next;
      }
      # we use query/target as in feature pairs the target=seqname and query=hseqname
      my ($query_coverage,  $query_spliced)  = &process_query( $features[$i], $tran2 );
      my ($target_coverage, $target_spliced) = &process_target( $features[$i], $tran1 );
      my $score = ( $query_coverage + $target_coverage )/2;
      if ( $score > $max_score && $query_spliced == 0 && $target_spliced == 0 ){
	$max_score = $score;
	$pair = $features[$i];
      }
      print STDERR "query:$id1 coverage:$query_coverage spliced:$query_spliced\n";
      print STDERR "target:$id2 coverage:$target_coverage spliced:$target_spliced\n";
    }
    return ($max_score,$pair);
    
}

############################################################
# the target
############################################################

sub process_target{
    my ($feat, $tran) = @_;
    my $transcript_length = $tran->seq->length;
    my @exons = sort { $a->length <=> $b->length } @{$tran->get_all_Exons};
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
		print STDERR "is spliced\n";
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
    #&print_exons_in_transcript($tran);
    return ($coverage,$is_spliced);
   
}

############################################################
# the query 
############################################################
sub process_query{
 my ($feat, $qtran) = @_;
 my $qtranscript_length = $qtran->seq->length;
 my @exons = sort { $a->length <=> $b->length } @{$qtran->get_all_Exons};
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
	     print STDERR "is spliced\n";
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
 print STDERR "coverage = $feature_length / $qtranscript_length = $coverage\n";

 #&print_exons_in_transcript($qtran);
 
 return ($coverage,$is_spliced);
}

############################################################

sub print_Feature{
  my $f = shift;
  print STDERR
    $f->seqname."\t".
      $f->start."-".$f->end."\t".
	  ($f->end - $f->start + 1)."\t".
	    $f->strand."\t".
	      $f->hseqname."\t".
		$f->hstart."-".$f->hend."\t".
		    ($f->hend - $f->hstart + 1 )."\t".
		      $f->strand."\t".
			"score:".$f->score."\t".
			  "perc_id:".$f->percent_id."\n";
}
	    
############################################################

sub compare_Exons{
  my ($human_t, $mouse_t, $gap_penalty ) = @_;

  # get the exons 5' to 3'
  my @human_exons;
  if ( $human_t->start_Exon->strand == 1 ){
    @human_exons = sort {$a->start <=> $b->start} @{$human_t->get_all_Exons};
  }
  else{
    @human_exons = sort {$b->start <=> $a->start} @{$human_t->get_all_Exons};
  }
  my @mouse_exons;
  if ( $mouse_t->start_Exon->strand == 1 ){
    @mouse_exons = sort {$a->start <=> $b->start} @{$mouse_t->get_all_Exons};
  }
  else{
    @mouse_exons = sort {$b->start <=> $a->start} @{$mouse_t->get_all_Exons};
  }
  
  my @score_matrix;
  my %comparison_score;

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
      $comparison_score{$human_exons[$i-1]}{$mouse_exons[$j-1]} = &blast_Exons( $human_exons[$i-1], $mouse_exons[$j-1] );
      #print STDERR "comparison( ".$human_exons[$i-1]->stable_id."-".$mouse_exons[$j-1]->stable_id." ) = ".
      #$comparison_score{$human_exons[$i-1]}{$mouse_exons[$j-1]} ."\n";
      
      $score_matrix[$i][$j] = 
	&max( $score_matrix[$i-1][$j]   + $gap_penalty,
	      $score_matrix[$i][$j-1]   + $gap_penalty,
	      $score_matrix[$i-1][$j-1] + $comparison_score{$human_exons[$i-1]}{$mouse_exons[$j-1]} );
    }
  }
  #return $score_matrix[$human_length][$mouse_length];
  #return @score_matrix;
  
  my ($human_list, $mouse_list ) = 
    &get_alignment( \@human_exons, \@mouse_exons, \@score_matrix, \%comparison_score, $gap_penalty );

  for ( my $i=0; $i<scalar(@$human_list); $i++ ){
    my $human_string;
    my $mouse_string;
    if ( $human_list->[$i] eq 'gap'){
      $human_string = "             ####GAP####";
    }
    else{
      $human_string = &exon_string( $human_list->[$i] );
    }
    if ( $mouse_list->[$i] eq 'gap'){
      $mouse_string = "             ####GAP####";
    }
    else{
      $mouse_string = &exon_string( $mouse_list->[$i] );
    }
    my $score;
    if( !($human_string eq "gap" || $mouse_string eq "gap") ){
      $score = $comparison_score{$human_list->[$i]}{$mouse_list->[$i]};
    }
    else{
      $score = 0;
    }
    print STDERR $human_string."\t<---->\t".$mouse_string.
      "\t score= ".$score."\n";
  }  
}

############################################################

sub get_alignment{
  my ($human_list, $mouse_list, $matrix, $comparison, $gap_penalty) = @_;
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
       == $matrix[$human_length-1][$mouse_length-1] + $comparison{$human_last}{$mouse_last} ){
    pop @human_list;
    pop @mouse_list;
    my ( $human_list2, $mouse_list2) = 
      &get_alignment( \@human_list, \@mouse_list, $matrix, $comparison, $gap_penalty);
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
      &get_alignment( \@human_list, \@mouse_list, $matrix, $comparison, $gap_penalty);
    push ( @{$human_list2}, $human_last );
    push ( @{$mouse_list2}, "gap" );
    return ( $human_list2, $mouse_list2 );
  }
  ############################################################
  # last exons of the second list is paired up with a gap
  else{
    pop @mouse_list;
    my ( $human_list2, $mouse_list2) =
      &get_alignment( \@human_list, \@mouse_list, $matrix, $comparison, $gap_penalty);
    push ( @{$human_list2}, "gap" );
    push ( @{$mouse_list2}, $mouse_last );
    return ( $human_list2, $mouse_list2 );
  }
} 

############################################################

sub exon_string{
  my ($exon) = @_;
  my $string = $exon->seqname.":".$exon->start."-".$exon->end.
    " (".($exon->end - $exon->start + 1 ).")".
      " strand:".$exon->strand.
	" phase:".$exon->phase.
	  " endphase:".$exon->end_phase;
  
}    

############################################################

sub max{
  my ($max, @others ) = @_;
  foreach my $other (@others){
    $max = $other if $other > $max;
  }
  return $max;
}

############################################################

sub blast_Exons{
  my ($exon1, $exon2) =@_;
  
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

  #print STDERR "comparing $id1 and $id2\n";

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
  my $blast =  Bio::EnsEMBL::Pipeline::Runnable::Blast->new ('-query'     => $seq1,
							     '-program'   => 'wutblastx',
							     '-database'  => $database,
							     -threshold_type => "PVALUE",
							     '-threshold' => 1e-10,
							     #'-filter'    => $filter,
							     '-options'       => 'V=200 B=200 W=5 E=0.01 E2=0.01',
							     #'-options'   => 'V=1000000'
							    );
  
  
  $blast->add_regex($file,'(\S+)');
  $blast->run();
  
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
    
    my @feat_by_score = sort { $b->score <=> $a->score } @featurepairs;
    return $feat_by_score[0]->score;
  }
  else{
    return 0;
  }
}
  
############################################################

sub print_exons_in_transcript{
    my $tran = shift;
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


































