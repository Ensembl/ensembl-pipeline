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
		  
#Location of the peptide file (file containing predicted peptides in fasta format)
		  #'pep_file'     => '/work4/mongin/dros/GeneBuild3/mapping/dros_pep.fa',

		  'pep_file'     => '',

#Put here all of the analysis which should be run, the analysis names are:
#Prints, Prosite, Profile, Pfam, Tmhmm, ncoils, Signalp, Seg
		  #'2berun' => 'prints,prosite',

		  '2berun' => '',

#Chunk size for each analysis 1 (whole genome but 1 entry at a time), 3 (give a full peptide dataset) and 2 if want to use the chunks...promise I will change it

		  #'Prints_chunk' => '2',
		  'Prints_chunk' => '',
		  
		  #'Prosite_chunk' => '3',
		  'Prosite_chunk' => '',

		  #'Profile_chunk' => '2',
		  'Profile_chunk' => '',

		  #'Pfam_chunk'    => '',
		  'Pfam_chunk'    => '',

		  #'Tmhmm_chunk'       => '2',
		  'Tmhmm_chunk'       => '',
		  
		  #'ncoils_chunk'       => '2',
		  'ncoils_chunk'       => '',

		  #'Signalp_chunk'     => '2',
		  'Signalp_chunk'     => '',

		  #'Seg_chunk'         => '3',
		  'Seg_chunk'         => '',

#paracel use 0 from no, 1 for yes. If nothing is set, suppose that the paracel is not used

		  'paracel'  => '',

#Define the chunk size (how many protein are going to be in each chunks)
		  #'chunk_size' => '100'
		  'chunk_size' => ''
		  )
}


1;
