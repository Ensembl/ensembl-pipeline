#!/usr/local/ensembl/bin/perl -w

use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster;
use Bio::EnsEMBL::SeqFeature;
use Getopt::Long;


my $dbhost;
my $dbname;
my $dbport;
my $region;
my $outfile;
my $v;

my $max_length = 100000;
my $max_feature_count = 200;



&GetOptions(
	    'region:s'          => \$region,
	    'dbhost:s'          => \$dbhost,
	    'dbname:s'          => \$dbname,
	    'dbport:s'          => \$dbport,
	    'max_length:s'      => \$max_length,
	    'max_feature_count' => \$max_feature_count,
	    'v!'                => \$v,
	   );


unless ( $dbhost && $dbname){
  print STDERR  "script to calculate the cut points for the slices to run jobs\n";
  print STDERR  "based on the genes in the database\n";
  print STDERR  "Usage: $0 [-dbname -dbhost -outfile -region]\n";
  exit(0);
}

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    -host   => $dbhost,
					    -user   => 'ensro',
					    -dbname => $dbname,
					    -port => $dbport,
					    );


my $slices; 

if (not defined $region){
   $slices = $db->get_SliceAdaptor->fetch_all('toplevel');
}
else{
    $slices = [$db->get_SliceAdaptor->fetch_by_region('toplevel',$region)];
}


if (!@$slices){
    print STDERR "cannot find slices -- bye.";
    exit(0); 
}

foreach my $slice (@$slices){

    my $length = $slice->length;

    my $slice_id = $slice->id;
    print STDERR "got slice $slice_id \n";

    my $all_genes = $slice->get_all_Genes;

    print STDERR scalar(@$all_genes)," genes found\n";
    
    next if ( ! (@$all_genes) );

    my ($clusters,$feature_count) = &cluster_Genes($all_genes); #cluster genes if they overlap over the genomic regions ('gene' refers to 'ensembl api gene')

    print STDERR scalar(@$clusters)." clusters found\n";

    my @clusters = sort { $a->start <=> $b->start } @$clusters;


    #group clusters if they satisfy the criteria ($max_length $max_feature_count)
    
    my $current_start = 1;
    my $current_end;

    my $this_cluster = shift (@clusters);

    my $count_feat = $feature_count->{$this_cluster};

    my ($this_start,$this_end);
 
    while (my $next_cluster = shift @clusters){
	 
	$this_start =  $this_cluster->start;
	$this_end   =  $this_cluster->end;	

	$current_end   = $this_cluster->end + int( ( $next_cluster->start - $this_cluster->end )/2 );
	
	print STDERR "\tslice:$slice_id st-ed:$this_start-$this_end feature_count:",$feature_count->{$this_cluster},"\n" if $v;
	 
	
	if (@clusters){
	    
	    my $after_next_cluster=$clusters[0];
	    
	    my $new_end = $next_cluster->end + int( ( $after_next_cluster->start - $next_cluster->end )/2 );
	    my $new_length = $new_end-$current_start+1;	
	    my $new_count = $count_feat + $feature_count->{$next_cluster};

	    if ( $new_length > $max_length ||  $new_count > $max_feature_count ){
		#print out group
		my $slice_to_print = $slice->sub_Slice($current_start,$current_end,1);

		print  $slice_to_print->name(),"\tlength:",($current_end-$current_start +1),"\tfeature_count:$count_feat\n";
		
		print STDERR "->GROUP:",$slice_to_print->name(),"\tlength:",($current_end-$current_start +1),"\tfeature_count:$count_feat\n" if $v; 
	    
		$current_start = $current_end + 1;       
		$count_feat=0;
	    }
	    
	}
	
	$this_cluster=$next_cluster;
	$count_feat += $feature_count->{$this_cluster};
	
    }
    print STDERR "\tslice:$slice_id st-ed:$this_start-$this_end feature_count:",$feature_count->{$this_cluster},"\n" if $v;
    $current_end = $length;
    my $slice_to_print = $slice->sub_Slice($current_start,$current_end,1);
    print  $slice_to_print->name(),"\tlength:",($current_end-$current_start +1),"\tfeature_count:$count_feat\n";
}


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

    my @exons = (@{$gene->get_all_Exons});
    return scalar(@exons);
    #return scalar(@{$gene->get_all_Exons} );
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

  $features{$cluster} += &get_feature_count($first_gene); 

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


