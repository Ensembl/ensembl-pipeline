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

#	    'dbname'      => 'human_110',
	    'dbname'      => '',

#	    'dbuser'      => 'ensadmin',
	    'dbuser'      => '',

#	    'dbpass'      => 'ensembl',
	    'dbpass'      => '',	  
);

# parameters for ensembl-pipeline/scripts/protein_pipeline/*.pl

%scripts_conf = ( 
	    # general options
	    #'runner'      => '',
		  'runner'      => '',

#	    'tmpdir'      => '',
		  'tmpdir'      => '',
		  
#	    'queue'       => '',
		  'queue'       => '',
		  
#Location of the different databases
		  'prints'      => '',
		  
		  'prosite'     => '',
		  
		  'profile'     => '',

		  'pfam'        => '',

		  'scanprosite'  => '',

		  'pep_file'     => '',

#Put here all of the analysis which should be run, the analysis names are:
#prints, prosite, profile, pfam, scanprosite, tmhmm, coils, signalp, seg
		  #'2berun' => 'prints,prosite',

		  '2berun' => '',

#Chunk size for each analysis 1 (whole genome but 1 entry at a time), 0 (give a full peptide dataset) and 2 if want to use the chunks...promise I will change it

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

		  'paracel'  => '',

#Define the chunk size (how many protein are going to be in each chunks)
		  'chunk_size' => ''
		  )
}


1;
