#!/usr/local/ensembl/bin/perl

=head1 NAME

comapre_Exons
 
=head1 DESCRIPTION

reads the config options from Bio::Ensembl::Pipeline::GeneComparison::GeneCompConf
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

# annotation
my $host1   = $DBHOST1;
my $dbname1 = $DBNAME1;
my $path1   = $PATH1;
my $port1   = $PORT1 || 3306;
my $type1   = $GENETYPES1;
my $user1   = $DBUSER1;

# prediction
my $host2   = $DBHOST2;
my $dbname2 = $DBNAME2;
my $path2   = $PATH2;
my $port2   = $PORT2 || 3306;
my $type2   = $GENETYPES2;
my $user2   = $DBUSER2;

# reference db (where the dna is)
my $ref_host1   = $REF_DBHOST1;
my $ref_dbname1 = $REF_DBNAME1;
my $ref_path1   = $REF_PATH1;
my $ref_port1   = $REF_PORT1 || 3306;
my $ref_user1   = $REF_DBUSER1;

# reference db (where the dna is for the prediction)
my $ref_host2   = $REF_DBHOST2;
my $ref_dbname2 = $REF_DBNAME2;
my $ref_path2   = $REF_PATH2;
my $ref_port2   = $REF_PORT2 || 3306;
my $ref_user2   = $REF_DBUSER2;



my $input_id;
my $pepfile;
my $gff_file;
my $chr_name;
my $check;
my $coding_exons;

# options
&GetOptions( 
	    'chr_name:s'  => \$chr_name,
	    'input_id:s'  => \$input_id,
	    'gff_file:s'  => \$gff_file,
	    'check'       => \$check,
	    'coding_exons'=> \$coding_exons,
	   );

unless( $input_id || $chr_name ){     
  print STDERR "Usage: compare_Exons.pl -input_id < chrname.chrstart-chrend >\n";
  print STDERR "                        -chr_name <chr_name>  (optional instead of -input_id)\n";
  print STDERR "                        -gff_file <file_name> (optional)\n";
  print STDERR "                        -check\n";
  exit(0);
}

if ( $check ){
  exit(0);
}
    
# connect to the databases 
my $dna_db1= new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $ref_host1,
					       -user  => $ref_user1,
					       -port  => $ref_port1,
					       -dbname=> $ref_dbname1);
$dna_db1->assembly_type($ref_path1); 
print STDERR "Connected to dna database $ref_dbname1 : $ref_host1 : $ref_user1\n";



my $db1= new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $host1,
					    -user  => $user1,
					    -dbname=> $dbname1,
					    -port  => $port1,
					    -dnadb => $dna_db1,
					   );

print STDERR "Connected to database $dbname1 : $host1 : $user1 \n";


my $dna_db2= new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $ref_host2,
						-user  => $ref_user2,
					        -port  => $ref_port2,
						-dbname=> $ref_dbname2);
$dna_db2->assembly_type($ref_path2); 
print STDERR "Connected to dna database $ref_dbname2 : $ref_host2 : $ref_user2\n";




my $db2= new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $host2,
					    -user  => $user2,
					    -dbname=> $dbname2,
					    -port  => $port2,
					    -dnadb => $dna_db2,
					   );

print STDERR "Connected to database $dbname2 : $host2 : $user2 \n";



# use different golden paths

my $sa1 = $db1->get_SliceAdaptor;
my $sa2 = $db2->get_SliceAdaptor;

# get a slice with a piece-of chromosome #
my ($slice1,$slice2);

# get genomic region 
my ($chr,$chrstart,$chrend);
if ($input_id){
  $chr      = $input_id;
  $chr         =~ s/\.(.*)-(.*)//;
  $chrstart = $1;
  $chrend   = $2;

  unless ( $chr && $chrstart && $chrend ){
    print STDERR "bad input_id option, try something like 20.1-5000000\n";
    exit(0);
  }
  
  print STDERR "Fetching region $chr, $chrstart - $chrend\n";
  $slice1 = $sa1->fetch_by_chr_start_end($chr,$chrstart,$chrend);
  $slice2 = $sa2->fetch_by_chr_start_end($chr,$chrstart,$chrend);
}
elsif ($chr_name ){
  
  $slice1 = $sa1->fetch_by_chr_name($chr_name);
  $slice2 = $sa2->fetch_by_chr_name($chr_name);

  $input_id = $chr_name.".1-".$slice1->length."\n";
  print STDERR "fetched region $input_id\n";
}


# get the genes of type @type1 and @type2 from $slice1 and $slice2, respectively #
my (@genes1,@genes2);
my (@trascripts1,@transcripts2);

foreach my $type ( @{ $type1 } ){
  print STDERR "Fetching genes of type $type\n";
  my @more_genes = @{$slice1->get_all_Genes_by_type($type)};
  my @more_trans = ();
  foreach my $gene ( @more_genes ){
    push ( @more_trans, @{$gene->get_all_Transcripts} );
  }
  push ( @genes1, @more_genes ); 
  print STDERR scalar(@more_genes)." genes found\n";
  print STDERR "with ".scalar(@more_trans)." transcripts\n";
}

foreach my $type ( @{ $type2 } ){
  print STDERR "Fetching genes of type $type\n";
  my @more_genes = @{$slice2->get_all_Genes_by_type($type)};
  my @more_trans = ();
  foreach my $gene ( @more_genes ){
    push ( @more_trans, @{$gene->get_all_Transcripts} );
  }
  push ( @genes2, @more_genes ); 
  print STDERR scalar(@more_genes)." genes found\n";
  print STDERR "with ".scalar(@more_trans)." transcripts\n";
}

# get a GeneComparison object 
my $gene_comparison = 
  Bio::EnsEMBL::Pipeline::GeneComparison::GeneComparison->new(					     
							      '-annotation_genes' => \@genes1,
							      '-prediction_genes' => \@genes2,
							      '-input_id'         => $input_id,
							     );



#########################################################

## cluster the genes we have passed to $gene_comparison

my @gene_clusters    = $gene_comparison->cluster_Genes;

my ($unmatched1,$unmatched2) = $gene_comparison->get_unmatched_Genes;

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


### if you want to print to GFF, include a gff_file in command line options
if ($gff_file){
  $gene_comparison->annotation_db( $db1 );
  $gene_comparison->prediction_db( $db2 );
  $gene_comparison->gff_file( $gff_file );
}

# run the analysis

my $coding = $coding_exons ? 1:0;

print STDERR "coding: $coding\n";

$gene_comparison->compare_Exons(\@gene_clusters,$coding,'verbose');

# If you remove the string 'verbose', it'll print out less stuff
# If you want to look at coding exons only, change the 0 for any defined value



############################################################

  




























































