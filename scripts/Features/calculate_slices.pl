#!/usr/local/ensembl/bin/perl -w

use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster;
use Bio::EnsEMBL::SeqFeature;
use Getopt::Long;


my $dbhost;
my $dbname;
my $chr_name;
my $outfile;


my $max_length = 200000;
my $max_feature_count = 200;


&GetOptions(
	    'chr_name:s'      => \$chr_name,
	    'dbhost:s'        => \$dbhost,
	    'dbname:s'        => \$dbname,
	    'outfile:s'       => \$outfile,
	    'max_length:s'    => \$max_length,
	    'max_feature_count' => \$max_feature_count,
	   );


unless ( $dbhost && $dbname && $chr_name && $outfile ){
  print STDERR  "script to calculate the cut points for the slices to run jobs\n";
  print STDERR  "based on the genes in the database\n";
  print STDERR  "Usage: $0 [-dbname -dbhost -outfile -chr_name]\n";
  exit(0);
}

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $dbhost,
					    '-user'   => 'ensro',
					    '-dbname' => $dbname,
					    );


my $slice = $db->get_SliceAdaptor->fetch_by_chr_name( $chr_name );
#my $slice = $db->get_SliceAdaptor->fetch_by_chr_start_end( $chr_name, 20500000, 21000000 );
my $chr_length = $slice->length;
print STDERR "got slice $chr_name.1-$chr_length\n";

my @all_genes = $slice->get_all_Genes;
unless ( @all_genes ){
  print STDERR "No genes found - exiting\n";
  exit(0);
}

my ($clusters,$feature_count) = &cluster_Genes(@all_genes);
print STDERR scalar(@$clusters)." clusters found\n";

my @clusters = sort { $a->start <=> $b->start } @$clusters;

open ( OUTPUT, ">$outfile" ) or die("cannot open $outfile");

my $start = 1;
my $count = 0;
my $end;
my $this_end;
my $order = 1;

CLUSTER:
for (my $i=0; $i<=$#clusters; $i++ ){
  
  my $extra_features = $feature_count->{$clusters[$i]};
  
  my $range_end = $chr_length;
  if ( $i<$#clusters ){
    $range_end = $clusters[$i]->end + int( ( $clusters[$i+1]->start - $clusters[$i]->end )/2 );
  }
  print OUTPUT "   this range:$chr_name.".($end+1)."-$range_end\tlength:".($range_end - $end)."\tfeature_count:$extra_features\n";
  
  if ( $i == 0 ){
    $end = $clusters[$i]->end + int( ( $clusters[$i+1]->start - $clusters[$i]->end )/2 );
  }
  
  if ( $i==$#clusters ){
    $this_end = $chr_length;    
  }
  else{
    $this_end = $clusters[$i]->end + int( ( $clusters[$i+1]->start - $clusters[$i]->end )/2 );
  }
  
  if ( ($this_end - $start + 1 ) > $max_length || ( $count + $extra_features ) > $max_feature_count ){
    if ( $order == 1 ){
      my $features = $count + $extra_features;
      print OUTPUT "$chr_name.$start-$this_end\tlength:".($this_end-$start+1)."\tfeature_count:$features (single range)\n";
      $start = $this_end + 1;
      $count = 0;
      $order = 1;
      next CLUSTER;
    }
    else{
      print OUTPUT "$chr_name.$start-$end\tlength:".($end-$start+1)."\tfeature_count:$count\n";
      $start = $end + 1;
      $count = $extra_features;
      $order = 1;
      next CLUSTER;
    }
  }
  else{
    $order++;
    $count += $extra_features;
    $end    = $this_end;
    next CLUSTER;
  }
}

if ( $end != $chr_length ){
  print OUTPUT "$chr_name.$start-$chr_length\tlength:".( $chr_length - $start + 1 )."\tfeature_count:$count\n";
}

close(OUTPUT);

############################################################

sub get_gene_start_end_strand {
  my ($g) = @_;
  my @exons = sort { $a->start <=> $b->start } @{$g->get_all_Exons}; 
  my $start = $exons[0]->start;
  my $end   = $exons[-1]->end;
  return ($start, $end, $exons[0]->strand);
}

############################################################

sub by_lower_coordinate {
    my @aexons = sort { $a->start <=> $b->start } @{$a->get_all_Exons};
    my @bexons = sort { $a->start <=> $b->start } @{$b->get_all_Exons};
    my $alow = $aexons[0]->start;
    my $blow = $bexons[0]->start;
    return $alow <=> $blow;
}

############################################################

sub get_feature_count{
    my $gene = shift;
    return scalar(@{$gene->get_all_Exons} );
}

#####################################################################

sub cluster_Genes{
  my ($genes) = @_;
  
  my %features;
  # first sort the genes
  my @genes = sort by_lower_coordinate @$genes;
  
  # create a new cluster 
  my $cluster = Bio::EnsEMBL::SeqFeature->new();
  
  my @clusters;
  
  my $first_gene = shift @genes;
  # put the first gene into these cluster
  my ($start,$end,$strand) = &get_gene_start_end_strand( $first_gene );
  
  $cluster->start( $start );
  $cluster->end($end);
  $features{$cluster} += &get_feature_count($genes[0]);
  # store the list of clusters
  push( @clusters, $cluster );
  
  # loop over the rest of the genes
 LOOP1:
  for (my $c=1; $c<=$#genes; $c++){
    
    my $next_gene = shift @genes;
    my ( $tstart,$tend,$tstrand ) = &get_gene_start_end_strand( $next_gene );
    if ( !( $tend   < $cluster->start 
	    ||
	    $tstart > $cluster->end 
	  ) 
       ){
      
      $features{$cluster} += &get_feature_count($next_gene);
      # re-adjust size of cluster
      if ( $tstart < $cluster->start ) {
	$cluster->start( $tstart );
      }
      if ( $tend > $cluster->end ) {
	$cluster->end($tend);
      }
    }
    else{
      # else, create a new cluster with this feature
      $cluster = Bio::EnsEMBL::SeqFeature->new();
      $cluster->start( $tstart );
      $cluster->end( $tend );
      $features{$cluster} += &get_feature_count($next_gene);
      push(@clusters,$cluster);
    }
  }
  return (\@clusters,\%features);
}


############################################################


