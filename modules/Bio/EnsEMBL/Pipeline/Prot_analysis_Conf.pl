# Copyright GRL & EBI 2001
# Author: Emmanuel Mongin
# Creation: 03.10.2001

# configuration information for Protein Analysis scripts
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
);

# parameters for ensembl-pipeline/scripts/protein_pipeline/*.pl

%scripts_conf = ( 
	    # general options
	    #'runner'      => '/work2/vac/ensembl-pipeline/scripts/test_RunnableDB',
		  'runner'      => '',

#	    'tmpdir'      => '/scratch3/ensembl/vac',
		  'tmpdir'      => '',
		  
#	    'queue'       => 'acarilong',
		  'queue'       => 'acari',
		  
#Location of the different databases
		  'prints'      => '',
		  
		  'prosite'     => '',
		  
		  'profile'     => '',

		  'pfam'        => '',

		  'scanprosite'  => '',

#Chunk size for each analysis 1 (whole genome but 1 entry at a time) or 0 (give a full peptide dataset)
#If nothing is set, will consider that you don't want to run it...

		  'prints_chunk' => '',
		  
		  'prosite_chunk' => '',
		  
		  'profile_chunk' => '',

		  'pfam_chunk'    => '',
		  
		  'scanprosite_chunk' => '',

		  'tmhmm_chunk'       => '',
		  
		  'coils_chunk'       => '',

		  'signalp_chunk'     => '',
		  
		  'seg_chunk'         => '',

#paracel use 0 from no, 1 for yes. If nothing is set, suppose that the paracel is not used

		  'paracel'  => ''
 
}



1;
