# Copyright GRL & EBI 2001
# Author: Val Curwen
# Creation: 02.05.2001

# configuration information for GeneBuild scripts
# give useful keynames to things

# I've left in sample entries for the various options to hopefully make this easier to use

BEGIN {
package main;

%GB_conf = ( 
	    # general options
#	    'runner'      => '/work2/vac/ensembl-pipeline/scripts/test_RunnableDB',
	    'runner'      => '',

#	    'tmpdir'      => '/scratch3/ensembl/vac',
	    'tmpdir'      => '',

#	    'dbhost'      => 'ecs1b',
	    'dbhost'      => '',

#	    'dbname'      => 'simon_dec12',
	    'dbname'      => '',

#	    'dbuser'      => 'ensadmin',
	    'dbuser'      => '',

#	    'dbuser'      => 'ensembl',
	    'dbpass'      => '',

#	    'golden_path'      => 'UCSC',
	    'golden_path'      => '',

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

#	    'pm_output'   => '/work2/vac/GeneBuild/script-test/pmf.out',
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

	    # options specific to Targetted runnables
#	    'targetted_runnables'   => ['TargettedGeneWise', 'TargettedGeneE2G'],
	    'targetted_runnables'   => [''],

	    # options specific to length based runnables
#	    'length_runnables'      => ['Riken_BlastMiniGenewise', 'FPC_BlastMiniGenewise', 'GeneBuilder'],
	    'length_runnables'      => [''],

#            'size'        => '5000000',
            'size'        => '',
	   );


}

1;
