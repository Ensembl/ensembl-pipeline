#!/usr/local/ensembl/bin/perl -w

use strict;  
use Bio::EnsEMBL::DBSQL::DBAdaptor;
#use Bio::EnsEMBL::Pipeline::GeneComparison::GeneComparison;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::GeneComparison::ObjectMap;
use Bio::EnsEMBL::Pipeline::GeneComparison::GenePair;
use Bio::EnsEMBL::Pipeline::GeneComparison::ComparativeTools;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Bio::EnsEMBL::Pipeline::Runnable::Blast;
use Getopt::Long;

# human_db
my $human_dbname = 'human_estgenes_eae';
my $human_dbhost = 'ecs2b';
my $human        = 'Homo sapiens';

# human_dnadb
my $human_dnadbname = 'human_NCBI33_raw2';
my $human_dnadbhost = 'ecs2';

# mouse_db
my $mouse_dbname = 'mouse_estgenes_eae';
my $mouse_dbhost = 'ecs2b';
my $mouse        = 'Mus musculus';

# mouse_dnadb
my $mouse_dnadbname = 'mus_musculus_core_13_30_ro';
my $mouse_dnadbhost = 'ecs2';

# rat_db 
my $rat_dbname = 'rattus_norvegicus_core_13_2'; 
my $rat_dbhost = 'ecs2f'; 
my $rat        = 'Rattus norvegicus'; 
 
# rat_dnadb 
my $rat_dnadbname = 'rattus_norvegicus_core_13_2'; 
my $rat_dnadbhost = 'ecs2f'; 

# compara_db
my $compara_dbname = 'ensembl_compara_15_1';
my $compara_dbhost = 'ecs2d';
my $compara_config =  '/nfs/acari/eae/ensembl/ensembl-compara/modules/Bio/EnsEMBL/Compara/Compara.conf';

my $gene_id;
my $gap_penalty;
my $check;
my $coding_exons;

# options
&GetOptions( 
	    'gene_id:s'       => \$gene_id,
	    'gap_penalty:n'      => \$gap_penalty,
	    'check'              => \$check,
	    'coding_exons'       => \$coding_exons,
	   );


unless( $gap_penalty ){
  $gap_penalty = -100;
}

if ( $check ){
    exit(0);
}

unless( $gene_id ){
  print STDERR "$0 -gene_id [ -coding_exons -check -gap_penalty]\n";
  exit(0);
}


my $dnadb1;
my $dnadb2;
my $db1;
my $db2;

my %gene_pair;

############################################################
# connect to the databases 
$dnadb1 = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $human_dnadbhost,
					     -user  => 'ensro',
					     -dbname=> $human_dnadbname,
					     );

$db1 = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $human_dbhost,
					  -user  => 'ensro',
					  -dbname=> $human_dbname,
					  -dnadb => $dnadb1,
					  );


$dnadb2 = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $mouse_dnadbhost,
					     -user  => 'ensro',
					     -dbname=> $mouse_dnadbname,
					     );

$db2 = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $mouse_dbhost,
					  -user  => 'ensro',
					  -dbname=> $mouse_dbname,
					  -dnadb => $dnadb2,
					  );

############################################################
# get the genes

my $adaptor1 = $db1->get_GeneAdaptor;

my $gene1 = $adaptor1->fetch_by_stable_id($gene_id,1);
my $chr_name  = $gene1->chr_name;
my $chr_start = $gene1->start;
my $chr_end   = $gene1->end;


my $padding = 2000;
# get the slice of the gene
my $slice_adaptor1 = $db1->get_SliceAdaptor;
my $slice1 = $slice_adaptor1->fetch_by_chr_start_end( $chr_name, $chr_start - $padding, $chr_end + $padding );



my $compara_db = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(
							      -user      => 'ensro',
							      -dbname    => $compara_dbname,
							      -host      => $compara_dbhost,
							      -conf_file => $compara_config,
							     );
#my $target_db = $compara_db->get_db_adaptor('Mus musculus','NCBIM30');
#my $focus_db  = $compara_db->get_db_adaptor('Homo sapiens','NCBI33');

my $focus_db  = $db1;
my $target_db = $db2;
my $target_gene_adaptor = $target_db->get_GeneAdaptor;

my $slices2 = Bio::EnsEMBL::Pipeline::GeneComparison::ComparativeTools
  ->get_all_syntenic_slices($slice1, $focus_db, 'Homo sapiens', $compara_db, $target_db, 'Mus musculus' );

my $gene_id1 = $gene1->stable_id || $gene1->dbID;
unless ( @$slices2 ){
  print STDERR "NO_PAIR gene $gene_id1 has no syntenic region\n";
  exit(0);
}

if (@$slices2){
  
 SLICE:
  foreach my $slice2 ( @$slices2 ){
    
    my @genes2 = @{$slice2->get_all_Genes};
    
  GENE2:
    foreach my $gene2 ( @genes2 ){
      
      my $gene_id2 = $gene2->stable_id || $gene2->dbID;
      
      print STDERR "************************************************************\n";
      print STDERR "comparing $gene_id1 and $gene_id2\n";
      ############################################################
      # call the comparison method
      if ( $gene2->stable_id ){
	$gene2 = $target_gene_adaptor->fetch_by_stable_id($gene_id2,1);
      }
      else{
	$gene2 = $target_gene_adaptor->fetch_by_dbID($gene_id2,1);
      }

      my $gene_pair = Bio::EnsEMBL::Pipeline::GeneComparison::GenePair->new();
      
      print STDERR "comparing isoforms and aligning exons\n";
      my $result = $gene_pair->compare( $gene1, $gene2, $coding_exons);
      
      unless( $result ){
	print STDERR "DP_NO_PAIR gene $gene_id1 has no good match  with $gene_id2 in the syntenic region\n";
      }
      print STDERR "comparing orthologous $gene_id1 $gene_id2\n";
      print STDERR "finding exact matches between transcripts\n";
      my $result2 = $gene_pair->find_exact_matches( $gene1, $gene2, $coding_exons);
      
      unless( $result2 ){
	print STDERR "NO_PAIR gene $gene_id1 has no good match with $gene_id2 in the syntenic region\n";
      }
      ############################################################
    }
  }
}



























