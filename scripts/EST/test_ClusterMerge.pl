#!/usr/local/ensembl/bin/perl -w


use strict;
use diagnostics;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Runnable::ClusterMerge;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Exon;
use Getopt::Long;

my $dbhost;
my $dbname;
my $dbuser   = 'ensro';
my $dnadbhost;
my $dnadbname;


my $input_id;
my $chr;
my $chrstart;
my $chrend;
my $genes    = 1;
my $genetype;
my $geneIDs;
my $transcriptIDs;

my $filter = 'NM_';
my $verbose = 0;
my $info;

$| = 1;


&GetOptions( 'dbhost:s'  => \$dbhost,
	     'dbname:s'  => \$dbname,
	     'dnadbhost:s'=> \$dnadbhost,
	     'dnadbname:s'=> \$dnadbname,
	     'input_id:s'=> \$input_id,
	     'genetype:s'=> \$genetype,
	   );

unless ($input_id){
  print STDERR "USAGE: $0 -dbhost -dbname  -input_id [[-dnadhost -dnadbname ... ]\n";
  exit(0);
}

if ( defined($geneIDs) && defined($transcriptIDs) ){
  print STDERR "You can only define one of them: -geneIDs or -transcriptIDs\n";
  exit(0);
}

# default: print transcript IDs
if ( !defined($geneIDs) && !defined($transcriptIDs) ){
  $transcriptIDs = 1;
}

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host             => $dbhost,
					    -user             => $dbuser,
					    -dbname           => $dbname,
					    );


if ( $dnadbname && $dnadbhost ){
  my $dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host             => $dnadbhost,
						 -user             => $dbuser,
						 -dbname           => $dnadbname,
						);
  $db->dnadb($dnadb);
}

my $sa = $db->get_SliceAdaptor();
$input_id =~/(\S+)\.(\d+)-(\d+)/;
$chr = $1;
$chrstart = $2;
$chrend   = $3;
print STDERR "fetching slice $chr . $chrstart - $chrend\n";

my $slice = $sa->fetch_by_chr_start_end( $chr,
					 $chrstart,
					 $chrend,
				       );


my @these_genes    = @{$slice->get_all_Genes}   if ($genes);
print STDERR scalar(@these_genes)." found\n";
my @genes;
my %evidence;
if ( $filter ){
 GENE:
  foreach my $gene ( @these_genes ){
    my $is_in = 0;
    foreach my $tran ( @{$gene->get_all_Transcripts} ){
      my @evidence = &get_evidence($tran);
      foreach my $evi ( @evidence ){
	if ( $evi =~/$filter/ ){
	  $evidence{ $tran } = $evi;
	  $is_in =1;
	  push ( @genes, $gene );
	  next GENE;
	}
      }
    }
  }
  print STDERR scalar(@genes)." genes after filtering\n";
}
else{
  print STDERR "no filter applied\n";
  @genes = @these_genes;
}

my @transcripts;

############################################################
# put the transcripts into a file
  
my $input_file = "input_transcripts.$input_id.gff";
open( OUT, ">$input_file") or die("cannot open file $input_file");  

my $exon_count = 0;
my $trans_count = 0;
foreach my $gene ( @genes ){
  my @trans = @{$gene->get_all_Transcripts};
  push ( @transcripts, @trans );
  foreach my $tran ( @trans ){
    $trans_count++;
    my $trans_id = $evidence{ $tran }."_$trans_count";
    foreach my $exon ( @{$tran->get_all_Exons} ){
      $exon_count++;
      my $strand = "+";
      if ($exon->strand == -1) {
	$strand = "-";
      }
      my $phase = ".";
      my $score = 100;
      my $g_type = 'input';
      print OUT $exon_count."\t".
	$g_type."\t".
	  "exon\t". 
	    ($exon->start) . "\t" . 
	      ($exon->end) . "\t" .
		$score. "\t" . 
		  $strand . "\t" . 
		    $phase . "\t" . 
		      $trans_id. "\n";
      
    }
  }
}

close (OUT);

############################################################
# chop the transcripts into pieces
my @trans_bits;
my $rand1 = rand;

GENE: 
foreach my $gene (@genes) {
  if ( $genetype ){
    next unless ( $gene->type eq $genetype );
  }  
  
  my $g_type = $gene->type;
  
  ############################################################
  # split the transcripts into small overlapping transcripts
  ############################################################
 TRAN:
  foreach my $tran (@{$gene->get_all_Transcripts}){

    # size of the overlap between the transcript bits is pseudo-random
    # for each transcript:
    my $rand2 = rand;
    my $rand = ($rand1 + $rand2)/2;
    
    my $overlap = int(20 + $rand*30);
    my $est_length = int(200 + $rand*200);
    print STDERR "overlap = $overlap - est_length = $est_length\n" if $verbose;

    my @this_trans_bits;
    my $length = $tran->length;
    my $strand = $tran->start_Exon->strand;
    print STDERR "transcript: (strand = $strand)\n" if $verbose;
    foreach my $exon ( @{$tran->get_all_Exons} ){
      print STDERR $exon->start."-".$exon->end."(".($exon->length).") " if $verbose;
    }
    print STDERR "\n" if $verbose;
    print STDERR "transcript length: $length\n" if $verbose;

    ############################################################
    # divide into ($est_length)bp pieces if the length is larger than that
    my $pieces;
    if ( $length <= $est_length ){
      $pieces = 1;
    }
    else{
      $pieces = int( $length/$est_length + 1 );
    }
    print STDERR "number of pieces: $pieces\n" if $verbose;

    ############################################################
    # calculate the cuts:
    my @cuts;
    my $start = 1;
    my $end   = $est_length + $overlap;
    if ( $pieces > 1 ){
      for(my $i=1; $<=$pieces; $i++){
	print STDERR "piece: ( $start, $end )\n" if $verbose;
	push( @cuts, [$start,$end] );
	$start = $start + $est_length;
	$end   = $start + $est_length + $overlap - 1;
	if ( $start > $length ){
	  last;
	}
	if ( $end > $length ){
	  $end = $length;
	}
      }
    }
    else{
      push( @cuts, [ 1, $length] );
    }
    
    my $parent_object_id = $tran->dbID;
    my $exon_count = 0;
    my @exons = sort { $a->start <=> $b->start } @{$tran->get_all_Exons};
  
    ############################################################
    # go over each cut and take the exons and
    # exon-pieces included in this cut
  CUT:
    foreach my $cut ( @cuts ){
      print STDERR "cut: [".$cut->[0].",".$cut->[1]."]\n" if $verbose;
      my @included_exons;
      my $exon_start = 1;
      my $exon_end;
      foreach my $exon ( @exons ){
	$exon_end = $exon_start + $exon->length - 1;
	
	print STDERR "checking exon ".$exon_start."-".$exon_end." (".($exon_end-$exon_start+1).")\n" if $verbose;
	# completely included
	if ( $exon_start >= $cut->[0] 
	     &&
	     $exon_end <= $cut->[1] 
	   ){
	  my $new_exon = Bio::EnsEMBL::Exon->new();
	  $new_exon->start( $exon->start );
	  $new_exon->end( $exon->end );
	  $new_exon->strand( $exon->strand );
	  push ( @included_exons, $new_exon );
	  print STDERR "completely included --> created: ".$new_exon->start."-".$new_exon->end."\n" if $verbose;
	}
	# prefix is included
	elsif( $exon_start >= $cut->[0] 
	       && 
	       $exon_start <= $cut->[1]
	       &&
	       $exon_end > $cut->[1]
	     ){
	  my $new_exon = Bio::EnsEMBL::Exon->new();
	  $new_exon->start( $exon->start );
	  $new_exon->end( $exon->end - ( $exon_end - $cut->[1] ) );
	  $new_exon->strand( $exon->strand );
	  push ( @included_exons, $new_exon );
	  print STDERR "prefix included --> created: ".$new_exon->start."-".$new_exon->end."\n" if $verbose;
	}
	# suffix is included
	elsif( $exon_start < $cut->[0] 
	       &&
	       $exon_end >= $cut->[0]
	       &&
	       $exon_end <= $cut->[1]
	     ){
	  my $new_exon = Bio::EnsEMBL::Exon->new();
	  $new_exon->start( $exon->start + ( $cut->[0] - $exon_start ));
	  $new_exon->end( $exon->end );
	  $new_exon->strand( $exon->strand );
	  push ( @included_exons, $new_exon );
	  print STDERR "suffix included --> created: ".$new_exon->start."-".$new_exon->end."\n" if $verbose;
	}
	# internal block is included
	elsif( $exon_start < $cut->[0] 
	       &&
	       $exon_end > $cut->[1]
	     ){
	  my $new_exon = Bio::EnsEMBL::Exon->new();
	  $new_exon->start( $exon->start + ( $cut->[0] - $exon_start ));
	  $new_exon->end( $exon->end - ( $exon_end - $cut->[1] ) );
	  $new_exon->strand( $exon->strand );
	  push ( @included_exons, $new_exon );
	  print STDERR "block included --> created: ".$new_exon->start."-".$new_exon->end."\n" if $verbose;
	}
	$exon_start = $exon_end + 1;
      }

      my $trans_bit = Bio::EnsEMBL::Transcript->new();
      foreach my $included_exon ( @included_exons ){
	$trans_bit->add_Exon( $included_exon );
      }
      push (@this_trans_bits, $trans_bit);
      
    }
    
    ############################################################
    # check:
    print STDERR "transcript bits created\n" if $verbose;
    foreach my $tran ( @this_trans_bits ){
      foreach my $exon ( @{$tran->get_all_Exons} ){
	print STDERR $exon->start."-".$exon->end."  " if $verbose;
      }
      print STDERR "\n" if $verbose;
    }
    push ( @trans_bits, @this_trans_bits );

  }
}

############################################################
# put the pieces into one file
  
my $bits_file = "transcript_pieces.$input_id.gff";
open( OUT, ">$bits_file") or die("cannot open file $bits_file");  

TRANS_BITS:
foreach my $tran ( @trans_bits ){
  $trans_count++;
  $tran->dbID($trans_count);
  foreach my $exon ( @{$tran->get_all_Exons} ){
    $exon_count++;
    $exon->dbID( $exon_count);
    my $strand = "+";
    if ($exon->strand == -1) {
      $strand = "-";
    }
    my $phase = ".";
    my $score = 100;
    my $g_type = 'bits';
    print OUT $exon_count."\t".
      $g_type."\t".
	"exon\t". 
	  ($exon->start) . "\t" . 
	    ($exon->end) . "\t" .
	      $score. "\t" . 
		$strand . "\t" . 
		  $phase . "\t" . 
		    $tran->dbID. "\n";
    
  }
}

close (OUT);

############################################################
# run ClusterMerge

print STDERR "\nRunning ClusterMerge algorithm\n\n";
my $merge_object 
  = Bio::EnsEMBL::Pipeline::Runnable::ClusterMerge->new(
							-transcripts      => \@trans_bits,
							-comparison_level => 3,
							-splice_mismatch  => 0,
							-intron_mismatch  => 0,
							-exon_match       => 0,
							-minimum_order    => 1,
						       );

$merge_object->run;
my @merged_transcripts = $merge_object->output;


###############################
# put the result in a new file
  
my $result_file = "merged_pieces.$input_id.gff";
open( OUT, ">$result_file") or die("cannot open file $result_file");  

TRANS_BITS:
foreach my $tran ( @merged_transcripts ){
  $trans_count++;
  $tran->dbID($trans_count);
  foreach my $exon ( @{$tran->get_all_Exons} ){
    $exon_count++;
    my $strand = "+";
    if ($exon->strand == -1) {
      $strand = "-";
    }
    my $phase = ".";
    my $score = 100;
    my $g_type = 'merged';
    print OUT $exon_count."\t".
      $g_type."\t".
	"exon\t". 
	  ($exon->start) . "\t" . 
	    ($exon->end) . "\t" .
	      $score. "\t" . 
		$strand . "\t" . 
		  $phase . "\t" . 
		    $trans_count. "\n";
    
  }
}

close (OUT);




############################################################

my $cluster_into_genes = 1;

if ( $cluster_into_genes ){
  print STDERR "original transcripts = ".scalar( @transcripts )."\n";
  my $clusters = &cluster_Transcripts( \@transcripts );
  print STDERR scalar(@$clusters)." transcript-clusters found\n";
  my $multi_cluster = 0;
  foreach my $c ( @$clusters ){
    if ( scalar( @{$c->get_Transcripts} ) > 1 ){
      $multi_cluster++;
    }
  }
  print STDERR "$multi_cluster clusters with more than one transcript\n";
  
  print STDERR "merged transcripts = ".scalar( @merged_transcripts )."\n";
  my $clusters2 = &cluster_Transcripts( \@merged_transcripts );
  print STDERR scalar(@$clusters2)." merged_transcript_clusters found\n";
  my $multi_cluster2 = 0;
  foreach my $c ( @$clusters2 ){
    if ( scalar( @{$c->get_Transcripts} ) > 1 ){
      $multi_cluster2++;
    }
  }
  print STDERR "$multi_cluster2 clusters with more than one transcript\n";
  my $count = 0;
  foreach my $c ( @$clusters2 ){
    $count++;
    print STDERR "cluster $count:\n";
    foreach my $t ( @{$c->get_Transcripts} ){
      Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript($t);
      #print STDERR $t->dbID." ";
    }
    print STDERR "\n";
    
  }
}

############################################################

sub get_transcript_start_end_strand {
  my ($t) = @_;
  my @exons = sort { $a->start <=> $b->start } @{$t->get_all_Exons}; 
  my $start = $exons[0]->start;
  my $end   = $exons[-1]->end;
  return ($start, $end, $exons[0]->strand);
}

sub by_transcript_high {
  my $alow;
  my $blow;
  my $ahigh;
  my $bhigh;
  
  if ($a->start_Exon->strand == 1) {
    $alow  = $a->start_Exon->start;
    $ahigh = $a->end_Exon->end;
  } 
  else {
    $alow  = $a->end_Exon->start;
    $ahigh = $a->start_Exon->end;
  }

  if ($b->start_Exon->strand == 1) {
    $blow  = $b->start_Exon->start;
    $bhigh = $b->end_Exon->end;
  } 
  else {
    $blow  = $b->end_Exon->start;
    $bhigh = $b->start_Exon->end;
  }

  if ($ahigh != $bhigh) {
    return $ahigh <=> $bhigh;
  } 
  else {
    return $alow <=> $blow;
  }
}

sub by_transcript_low {
  my @aexons = sort { $a->start <=> $b->start } @{$a->get_all_Exons};
  my @bexons = sort { $a->start <=> $b->start } @{$b->get_all_Exons};
  my $alow = $aexons[0]->start;
  my $blow = $bexons[0]->start;
  return $alow <=> $blow;
}

############################################################

sub cluster_Transcripts {
  my $trans = shift;
  my $forward_trans;
  my $reverse_trans;
  foreach my $t ( @$trans ){
    if ( $t->start_Exon->strand == 1 ){
      push ( @{$forward_trans}, $t );
    }
    else{
      push ( @{$reverse_trans}, $t );
    }
  }
  my @clusters;
  if ( $forward_trans ){
    my $f_clusters = &cluster_Transcripts_by_strand($forward_trans);
    push( @clusters, @{$f_clusters} );
    print STDERR scalar( @{$f_clusters} ). " clusters on forward strand\n";
  }
  if ( $reverse_trans ){
    my $r_clusters = &cluster_Transcripts_by_strand($reverse_trans);
    push( @clusters, @{$r_clusters} );
    print STDERR scalar( @{$r_clusters} ). " clusters on reverse strand\n";
  }
  return \@clusters;
}

sub cluster_Transcripts_by_strand {

  my $transcripts_unsorted = shift;
  my @transcripts_unsorted = @$transcripts_unsorted;

  my @transcripts = sort by_transcript_low @transcripts_unsorted;
  
  my @clusters = ();

  my $verbose = 0;

  # clusters transcripts by whether or not any exon overlaps with an exon in
  foreach my $tran (@transcripts) {
    my @matching_clusters;
    my ($trans_start, $trans_end, $trans_strand) = get_transcript_start_end_strand($tran);
    
    print STDERR "transcript limits: $trans_start $trans_end \n" if $verbose;

  CLUSTER: 
    foreach my $cluster (@clusters) {
      
      print STDERR "Testing against cluster with limits " . $cluster->start ." to " . $cluster->end . "\n" if $verbose;
      
      if (!($trans_start > $cluster->end || $trans_end < $cluster->start) &&
	  $trans_strand == $cluster->strand) {
	
	print STDERR "In range\n" if $verbose;
        foreach my $cluster_transcript (@{$cluster->get_Transcripts()}) {
          foreach my $exon1 (@{$tran->get_all_Exons}) {
	    
            foreach my $cluster_exon (@{$cluster_transcript->get_all_Exons}) {
              if ($exon1->overlaps($cluster_exon) &&
                  $exon1->strand == $cluster_exon->strand) {
                push (@matching_clusters, $cluster);
                next CLUSTER;
              }
            }
          }
        }
      }
    }
    
    if (scalar(@matching_clusters) == 0) {
      print STDERR "Found new cluster for " . $tran->dbID . "\n" if $verbose;
      my $newcluster = new Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster;
      $newcluster->put_Transcripts($tran);
      push(@clusters,$newcluster);
      
    } 
    elsif (scalar(@matching_clusters) == 1) {
      print STDERR "Adding to cluster for " . $tran->dbID . "\n" if $verbose;
      $matching_clusters[0]->put_Transcripts($tran);

    } 
    else {
      # Merge the matching clusters into a single cluster
      print STDERR "Merging clusters for " . $tran->dbID . "\n" if $verbose;
      my @new_clusters;
      my $merged_cluster = new Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster;
      foreach my $clust (@matching_clusters) {
        $merged_cluster->put_Transcripts(@{$clust->get_Transcripts});
      }
      $merged_cluster->put_Transcripts($tran);
      push @new_clusters,$merged_cluster;
      # Add back non matching clusters
      foreach my $clust (@clusters) {
        my $found = 0;
      MATCHING: 
	foreach my $m_clust (@matching_clusters) {
          if ($clust == $m_clust) {
            $found = 1;
            last MATCHING;
          }
        }
        if (!$found) {
          push @new_clusters,$clust;
        }
      }
      @clusters = @new_clusters;
    }
  }
  
  #Safety checks
  my $ntrans = 0;
  my %trans_check_hash;
  foreach my $cluster (@clusters) {
    $ntrans += scalar(@{$cluster->get_Transcripts});
    foreach my $trans (@{$cluster->get_Transcripts}) {
      if (defined($trans_check_hash{"$trans"})) {
        print STDERR ("Transcript " . $trans->dbID . " added twice to clusters\n") ;
      }
      $trans_check_hash{"$trans"} = 1;
    }
    if (!scalar(@{$cluster->get_Transcripts})) {
      print STDERR ("Empty cluster");
    }
  }
  
  if ($ntrans != scalar(@transcripts)) {
    print STDERR ("Not all transcripts have been added into clusters $ntrans and " . scalar(@transcripts). " \n");
  }
  #end safety checks

  #print "Started with " . scalar(@transcripts) . " transcripts and ended with " . scalar(@clusters) . " clusters\n";

  return \@clusters;
}

############################################################

sub get_evidence{
  my $t = shift;
  my %evidence;
  foreach my $exon ( @{$t->get_all_Exons} ){
    foreach my $evi ( @{$exon->get_all_supporting_features} ){
      $evidence{$evi->hseqname}++;
    }
  }
  return keys %evidence;
}
