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


my $max_length = 100000;
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
#my $slice = $db->get_SliceAdaptor->fetch_by_chr_start_end( $chr_name, 84000000,191610523 );
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

my $current_start = 1;
my $current_end   = $clusters[0]->end + int( ( $clusters[1]->start - $clusters[0]->end )/2 );

my $extra_features = $feature_count->{$clusters[0]};
my $this_start  = $current_start;
my $this_end    = $current_end;
my $this_length = $this_end - $this_start + 1;
my $count = $extra_features;

print STDERR "   this range:$chr_name.$this_start-$this_end\tlength:".
    "$this_length\tfeature_count:$extra_features\n";

CLUSTER:
for (my $i=1; $i<=$#clusters; $i++ ){
    
    $extra_features = $feature_count->{$clusters[$i]};

    $this_end = $chr_length;
    if ( $i<$#clusters ){
	$this_end = $clusters[$i]->end + int( ( $clusters[$i+1]->start - $clusters[$i]->end )/2 );
    }
    $this_start  = $current_end + 1;
    $this_length = $this_end - $this_start + 1;
    print STDERR "   this range:$chr_name.$this_start-$this_end\tlength:".
	"$this_length\tfeature_count:$extra_features\n";
  
    my $new_length = $this_end - $current_start + 1;
    my $new_count  = $count + $extra_features;
    if ( $new_length > $max_length ||  $new_count > $max_feature_count ){
	print OUTPUT "$chr_name.$current_start-$current_end\tlength:".
	    ($current_end-$current_start +1)."\tfeature_count:$count\n";
	
	print STDERR "SLICE $chr_name.$current_start-$current_end\tlength:".
	    ($current_end-$current_start +1)."\tfeature_count:$count\n";
	
	$current_start = $current_end + 1;
	$current_end   = $this_end;
	$count = $extra_features;
    }
    else{
	$count = $new_count;
	$current_end = $this_end;
    }
    if ( $i==$#clusters ){
	print OUTPUT "$chr_name.$current_start-$chr_length\tlength:".( $chr_length - $current_start + 1 )."\tfeature_count:$count\n";
	print STDERR "SLICE $chr_name.$current_start-$chr_length\tlength:".( $chr_length - $current_start + 1 )."\tfeature_count:$count\n";
    }
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
  while( @genes ){
    
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


