#!/usr/local/ensembl/bin/perl -w

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
my $human_dbname = 'homo_sapiens_core_14_31';
my $human_dbhost = 'ecs2f';
my $human        = 'Homo sapiens';
#my $human_dbname = 'homo_sapiens_core_15_33';
#my $human_dbhost = 'ecs2f';


# human_dnadb
my $human_dnadbname = 'homo_sapiens_core_14_31';
my $human_dnadbhost = 'ecs2f';
#my $human_dnadbname = 'homo_sapiens_core_15_33';
#my $human_dnadbhost = 'ecs2f';

# mouse_db
my $mouse_dbname = 'mus_musculus_core_10_3a';
my $mouse_dbhost = 'ecs1d';
my $mouse        = 'Mus musculus';
#my $mouse_dbname = 'mus_musculus_core_15_30';
#my $mouse_dbhost = 'ecs2f';

# mouse_dnadb
my $mouse_dnadbname = 'mus_musculus_core_10_3a';
my $mouse_dnadbhost = 'ecs1d';
#my $mouse_dnadbname = 'mus_musculus_core_15_30';
#my $mouse_dnadbhost = 'ecs2f';

# rat_db
#my $rat_dbname = 'rattus_norvegicus_core_13_2';
#my $rat_dbhost = 'ecs2f';
my $rat        = 'Rattus norvegicus';
my $rat_dbname = 'rattus_norvegicus_core_15_2';
my $rat_dbhost = 'ecs2f';

# rat_dnadb
my $rat_dnadbname = 'rattus_norvegicus_core_15_2';
my $rat_dnadbhost = 'ecs2f';

# compara_db
my $compara_dbname = 'ensembl_compara_13_1';
my $compara_dbhost = 'ecs2f';

my $gene_id1;
my $gene_id2;
my $gap_penalty;
my $check;
my $coding_exons;

# options
&GetOptions( 
	    'gene_id1:s'       => \$gene_id1,
	    'gene_id2:s'       => \$gene_id2,
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

############################################################
# connect to the databases 
############################################################
# connect to the databases 
my $human_dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $human_dnadbhost,
						  -user  => 'ensro',
						  -dbname=> $human_dnadbname,
						 );

my $human_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $human_dbhost,
					       -user  => 'ensro',
					       -dbname=> $human_dbname,
					       -dnadb => $human_dnadb,
					      );

my $mouse_dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $mouse_dnadbhost,
						  -user  => 'ensro',
						  -dbname=> $mouse_dnadbname,
					     );

my $mouse_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $mouse_dbhost,
					       -user  => 'ensro',
					       -dbname=> $mouse_dbname,
					       -dnadb => $mouse_dnadb,
					      );


my $rat_dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $rat_dnadbhost,
						-user  => 'ensro',
						-dbname=> $rat_dnadbname,
					       );

my $rat_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $rat_dbhost,
					     -user  => 'ensro',
					     -dbname=> $rat_dbname,
					     -dnadb => $rat_dnadb,
					    );

############################################################

my $db1 = $human_db;
my $db2 = $mouse_db;



############################################################
# get the genes

my $adaptor1 = $db1->get_GeneAdaptor;
my $adaptor2 = $db2->get_GeneAdaptor;

my $gene1 = $adaptor1->fetch_by_stable_id( $gene_id1,1);
my $gene2 = $adaptor2->fetch_by_stable_id( $gene_id2,1);

############################################################
# call the comparison method

my $gene_pair = Bio::EnsEMBL::Pipeline::GeneComparison::GenePair->new();

print STDERR "comparing CDSs and aligning exons\n";
$gene_pair->compare_CDSs( $gene1, $gene2, $coding_exons );

#print STDERR "comparing isoforms and aligning exons\n";
#$gene_pair->compare( $gene1, $gene2, $coding_exons);

#print STDERR "comparing orthologous $gene_id1 $gene_id2\n";
print STDERR "finding exact matches between transcripts\n";
$gene_pair->find_exact_matches( $gene1, $gene2, $coding_exons);


############################################################






























