# Copyright GRL & EBI 2001
# Author: Val Curwen
# Creation: 02.05.2001

# configuration information for GeneBuild scripts
# give useful keynames to things

# I've left in sample entries for the various options to hopefully make this easier to use

BEGIN {
package main;

# general parameters for database connection

%db_conf = (
#	    'dbhost'      => 'ecs1b',
	    'dbhost'      => '',

#	    'dbname'      => 'simon_dec12',
	    'dbname'      => '',

#	    'dbuser'      => 'ensadmin',
	    'dbuser'      => '',

#	    'dbpass'      => 'ensembl',
	    'dbpass'      => '',	  

#           'golden_path' => 'UCSC',
	    'golden_path' => '',
);

# parameters for ensembl-pipeline/scripts/GeneBuild/*.pl

%scripts_conf = ( 
	    # general options
#	    'runner'      => '/work2/vac/ensembl-pipeline/scripts/test_RunnableDB',
	    'runner'      => '',

#	    'tmpdir'      => '/scratch3/ensembl/vac',
	    'tmpdir'      => '',

#	    'queue'       => 'acarilong',
	    'queue'       => '',

	    # prepare_proteome options
#	    'refseq'      => '/work2/vac/TGW/Dec_gp/human_proteome/refseq.fa',
	    'refseq'      => '',

#	    'sptr'        => '/work2/vac/TGW/Dec_gp/human_proteome/sptr_minus_P17013.fa',
	    'sptr'        => '',

#	    'pfasta'      => '/work2/vac/GeneBuild/script-test/human_proteome.fa',
	    'pfasta'      => '',

	    # pmatch from below

	    # pmatch_filter options
            # pfasta from above
#	    'pmatch'      => '/work2/vac/rd-utils/pmatch',
	    'pmatch'      => '',

#	    'pm_output'   => '/work2/vac/GeneBuild/script-test/',
	    'pm_output'   => '',

#	    'fpcdir'      => '/work2/vac/data/humangenome/',
	    'fpcdir'      => '',

	    # protein2cdna options
#	    'rf_gnp'      => '/work2/vac/TGW/Dec_gp/human_proteome/refseq.gnp',
	    'rf_gnp'      => '',

#	    'sp_swiss'    => '/work2/vac/TGW/Dec_gp/human_proteome/sptr.swiss',
	    'sp_swiss'    => '',

#	    'efetch'      => '/usr/local/ensembl/bin/efetch.new', 
	    'efetch'      => '', 

#	    'cdna_pairs'  => '/work2/vac/GeneBuild/script-test/cpp.out',
	    'cdna_pairs'  => '',

#	    'cdna_seqs'   => '/work2/vac/GeneBuild/script-test/cdna_seqs.fa',
	    'cdna_seqs'   => '',

	    # options specific to Targetted runnables
#	    'targetted_runnables'   => ['TargettedGeneWise', 'TargettedGeneE2G'],
	    'targetted_runnables'   => [''],

	    # options specific to length based runnables
#	    'length_runnables'      => ['CombinedGeneBuild'],
	    'length_runnables'      => ['CombinedGeneBuild'],

#           'size'        => '5000000',
            'size'        => '5000000',

	   );

# seqfetch parameters - location of getseqs indices used by the RunnableDBs
%seqfetch_conf = (

#	    location of seq_index for protein sequences; if not set, will use pfetch
#	    'protein_index' => '/data/blastdb/Ensembl/swall.all',
	    'protein_index' => '',

);	    

# TargettedGeneE2G parameters
%targetted_conf = (
#	    location of seq_index for cdna sequences; if not set, will use pfetch
#	    'cdna_index' => '/data/blastdb/Ensembl/cdna_seqs.fa',
	    'cdna_index' => '',

#	    location of seq_index for human_proteome sequences - swissprot plus refseq; 
#	    if not set, will use pfetch
#	    'protein_index' => '/data/blastdb/Ensembl/human_proteome.fa',
	    'protein_index' => '',
);

# FPC_BlastMiniGenewise parameters
%similarity_conf = (

);

# Riken_BlastMiniGenewise parameters
%riken_conf = (
#	    location of seq_index for riken protein sequences
#	    'riken_index' => '/data/blastdb/Ensembl/riken_prot',
	    'riken_index' => '',
);

# GeneBuild parameters
%genebuild_conf = (
	  
);
}

1;
