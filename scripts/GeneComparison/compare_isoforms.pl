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


my $human_gene_id;
my $mouse_gene_id;
my $gap_penalty;
my $check;

# options
&GetOptions( 
	     'human_gene:s'       => \$human_gene_id,
	     'mouse_gene:s'       => \$mouse_gene_id,
	     'gap_penalty:n'      => \$gap_penalty,
	     'check'              => \$check,
	   );

unless( $gap_penalty ){
  $gap_penalty = -100;
}

if ( $check ){
    exit(0);
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

############################################################
# get the genes

my $human_adaptor = $human_db->get_GeneAdaptor;
my $mouse_adaptor = $mouse_db->get_GeneAdaptor;

my $human_gene = $human_adaptor->fetch_by_stable_id( $human_id,1);
my $mouse_gene = $mouse_adaptor->fetch_by_stable_id( $mouse_id,1);

############################################################
# call the comparison method

my $gene_pair = Bio::EnsEMBL::Pipeline::GeneComparison::GenePair->new();
$gene_pair->compare( $human_gene, $mouse_gene);

############################################################






























