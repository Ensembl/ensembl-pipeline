#!/usr/local/ensembl/bin/perl -w


use strict;
use diagnostics;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Runnable::ClusterMerge;
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


my @genes    = @{$slice->get_all_Genes}   if ($genes);
print STDERR scalar(@genes)." found\n";

############################################################
# put the transcripts into a file
  
my $input_file = "input_transcripts.$input_id.gff";
open( OUT, ">$input_file") or die("cannot open file $input_file");  

my $exon_count = 0;
my $trans_count = 0;
foreach my $gene ( @genes ){
  foreach my $tran ( @{$gene->get_all_Transcripts} ){
    $trans_count++;
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
		      $trans_count. "\n";
      
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
    my $est_length = int(200 + $rand*300);
    print STDERR "overlap = $overlap - est_length = $est_length\n";

    my @this_trans_bits;
    my $length = $tran->length;
    my $strand = $tran->start_Exon->strand;
    print STDERR "transcript: (strand = $strand)\n";
    foreach my $exon ( @{$tran->get_all_Exons} ){
      print STDERR $exon->start."-".$exon->end."(".($exon->length).") ";
    }
    print STDERR "\n";
    print STDERR "transcript length: $length\n";

    ############################################################
    # divide into ($est_length)bp pieces if the length is larger than that
    my $pieces;
    if ( $length <= $est_length ){
      $pieces = 1;
    }
    else{
      $pieces = int( $length/$est_length + 1 );
    }
    print STDERR "number of pieces: $pieces\n";

    ############################################################
    # calculate the cuts:
    my @cuts;
    my $start = 1;
    my $end   = $est_length + $overlap;
    if ( $pieces > 1 ){
      for(my $i=1; $<=$pieces; $i++){
	print STDERR "piece: ( $start, $end )\n";
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
      print STDERR "cut: [".$cut->[0].",".$cut->[1]."]\n";
      my @included_exons;
      my $exon_start = 1;
      my $exon_end;
      foreach my $exon ( @exons ){
	$exon_end = $exon_start + $exon->length - 1;
	
	print STDERR "checking exon ".$exon_start."-".$exon_end." (".($exon_end-$exon_start+1).")\n";
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
	  print STDERR "completely included --> created: ".$new_exon->start."-".$new_exon->end."\n";
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
	  print STDERR "prefix included --> created: ".$new_exon->start."-".$new_exon->end."\n";
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
	  print STDERR "suffix included --> created: ".$new_exon->start."-".$new_exon->end."\n";
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
	  print STDERR "block included --> created: ".$new_exon->start."-".$new_exon->end."\n";
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
    print STDERR "transcript bits created\n";
    foreach my $tran ( @this_trans_bits ){
      foreach my $exon ( @{$tran->get_all_Exons} ){
	print STDERR $exon->start."-".$exon->end."  ";
      }
      print STDERR "\n";
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

