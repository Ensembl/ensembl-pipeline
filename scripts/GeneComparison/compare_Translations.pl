#!/usr/local/ensembl/bin/perl


=head1 NAME

comapre_Exons
 
=head1 DESCRIPTION

reads the config options from Bi::Ensembl::Pipeline::GeneComparison::GeneCompConf
and reads as input an input_id in the style of other Runnables, i.e. -input_id chr_name.chr_start-chr_end

=head1 OPTIONS

    -input_id  The input id: chrname.chrstart-chrend

=cut

use strict;  
use diagnostics;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::GeneComparison::GeneComparison;
use Getopt::Long;

## load all the parameters
use Bio::EnsEMBL::Pipeline::GeneComparison::GeneCompConf;


my $host1   = $DBHOST1;
my $dbname1 = $DBNAME1;
my $path1   = $PATH1;
my $type1   = $GENETYPES1;
my $user1   = $DBUSER1;

my $host2   = $DBHOST2;
my $dbname2 = $DBNAME2;
my $path2   = $PATH2;
my $type2   = $GENETYPES2;
my $user2   = $DBUSER2;

# get your refdb for the dna

my $refdb_host1 = $REF_DBHOST1;
my $refdb_name1 = $REF_DBNAME1;
my $refdb_path1 = $REF_PATH1;
my $refdb_user1 = $REF_DBUSER1;

my $refdb_host2 = $REF_DBHOST2;
my $refdb_name2 = $REF_DBNAME2;
my $refdb_path2 = $REF_PATH2;
my $refdb_user2 = $REF_DBUSER2;


my $runnable;
my $input_id;
my $write  = 0;
my $check  = 0;
my $params;
my $pepfile;

# can override db options on command line
&GetOptions( 
	     'input_id:s'  => \$input_id,

	     );
	
unless( $input_id){     
  print STDERR "Usage: run_GeneComparison.pl -input_id < chrname.chrstart-chrend >\n";
  exit(0);
}
    
# get genomic region 
my $chr      = $input_id;
$chr         =~ s/\.(.*)-(.*)//;
my $chrstart = $1;
my $chrend   = $2;

unless ( $chr && $chrstart && $chrend ){
       print STDERR "bad input_id option, try something like 20.1-5000000\n";
}

# connect to the databases 
my $dnadb1= new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $refdb_host1,
					       -user  => $refdb_user1,
					       -dbname=> $refdb_name1);

print STDERR "Connected to dna database $refdb_name1 : $refdb_host1 : $refdb_user1 \n";


my $db1= new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $host1,
					    -user  => $user1,
					    -dbname=> $dbname1,
					    -dnadb => $dnadb1,
					   );
    
print STDERR "Connected to database $dbname1 : $host1 : $user1 \n";

my $dnadb2= new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $refdb_host2,
					       -user  => $refdb_user2,
					       -dbname=> $refdb_name2);

print STDERR "Connected to dna database $refdb_name2 : $refdb_host2 : $refdb_user2 \n";


my $db2= new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $host2,
					    -user  => $user2,
					    -dbname=> $dbname2,
					    -dnadb => $dnadb2,
					   );
print STDERR "Connected to database $dbname2 : $host2 : $user2 \n";



# both databases should be in the same assembly, hence





# use different golden paths
$db1->static_golden_path_type($path1); 
$db2->static_golden_path_type($path2); 

my $sgp1 = $db1->get_StaticGoldenPathAdaptor;
my $sgp2 = $db2->get_StaticGoldenPathAdaptor;

# get a virtual contig with a piece-of chromosome #
my ($vcontig1,$vcontig2);

print STDERR "Fetching region $chr, $chrstart - $chrend\n";
$vcontig1 = $sgp1->fetch_VirtualContig_by_chr_start_end($chr,$chrstart,$chrend);
$vcontig2 = $sgp2->fetch_VirtualContig_by_chr_start_end($chr,$chrstart,$chrend);

# get the genes of type @type1 and @type2 from $vcontig1 and $vcontig2, respectively #
my (@genes1,@genes2);

foreach my $type ( @{ $type1 } ){
  print STDERR "Fetching genes of type $type\n";
  my @more_genes = $vcontig1->get_Genes_by_Type($type);
  push ( @genes1, @more_genes ); 
  print STDERR scalar(@more_genes)." genes found\n";
}

foreach my $type ( @{ $type2 } ){
  print STDERR "Fetching genes of type $type\n";
  my @more_genes = $vcontig2->get_Genes_by_Type($type);
  push ( @genes2, @more_genes ); 
  print STDERR scalar(@more_genes)." genes found\n";
}

# get a GeneComparison object 



my $gene_comparison = 
  Bio::EnsEMBL::Pipeline::GeneComparison::GeneComparison->new(
      
		      '-annotation_genes' => \@genes1,
							      
							      
	             '-prediction_genes' => \@genes2,
							      
		      '-input_id'        => $input_id,

									       
									       
										     );
										     



#########################################################

## cluster the genes we have passed to $gene_comparison

my @gene_clusters    = $gene_comparison->cluster_Genes;
my @unclustered = $gene_comparison->unclustered_Genes;

## print out the results of the clustering:

my @unclustered1;
my @unclustered2;

UNCLUSTER:
foreach my $uncluster ( @unclustered ){
  my @gene = $uncluster->get_Genes;
  if ( scalar(@gene)>1 ){
    print STDERR "genes @gene are wrongly unclustered\n";
  }
  my $this_type = $gene[0]->type;
  foreach my $type ( @{ $type1 } ){
    if ( $this_type eq $type ){
      push( @unclustered1, $uncluster );
      next UNCLUSTER;
    }
  }
  foreach my $type ( @{ $type2 } ){
    if ( $this_type eq $type ){
      push( @unclustered2, $uncluster );
      next UNCLUSTER;
    }
  }
}

print STDERR scalar(@gene_clusters)." gene clusters formed\n";
print STDERR scalar(@unclustered1)." genes of type @$type1 left unclustered\n";
print STDERR scalar(@unclustered2)." genes of type @$type2 left unclustered\n";


$gene_comparison->compare_Translations(\@gene_clusters);


