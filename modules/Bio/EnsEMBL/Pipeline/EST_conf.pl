# Copyright GRL & EBI 2001
# Author: Val Curwen
# Creation: 19.07.2001

# configuration information for EST scripts
# give useful keynames to things

# I've left in sample entries for the various options to hopefully make this easier to use

BEGIN {
package main;

# options need by the various scripts
%scripts_conf = ( 
### general options
		 'runner'      => '/nfs/acari/eae/ensembl/ensembl-pipeline/scripts/run_EST_RunnableDB',
#		'runner'      => '',

		 'scriptdir'       => "/nfs/acari/eae/ensembl/ensembl-pipeline/scripts/EST/",
#		'scriptdir'       => "",

#		path to scratch directory where output subdirectories and files will be placed
#		 'tmpdir'      => '/scratch2/ensembl/eae/est/',
#		 'tmpdir'      => '/scratch4/ensembl/vac',
		 'tmpdir'      => 

#		'queue'          => "acarilong -m \"bcscab1 bcscab3 bcscab4 bcscab5 bcscab6 bcscab7 bcscab8 bcscab9\"",
#		'queue'          => "acari",
		'queue'          => "",

### for make_bsubs.pl - where to put the bsubs
#		'exonerate_bsubsfile'       => "/scratch4/ensembl/vac/exonerate_est.jobs",
#		'filter_bsubsfile'          => "/scratch4/ensembl/vac/filter_and_e2g.jobs",
#               'EST_GeneBuilder_bsubsfile' => "/work6a/eae.tmp/Mouse/jobs/EST_GeneBuilder.jobs",
		'exonerate_bsubsfile' => "",
		'filter_bsubsfile'    => "",
		'EST_GeneBuilder_bsubsfile' => "",

### for prepare_ests.pl

#		path to executable which will be used for splitting estfiles into chunks
		'filesplitter' => "/work2/gs2/gs2/bin/fastasplit",
#		'filesplitter' => "",

#		path to file with all the ESTs/cDNAs/whatever in it
#		'estfile'        => "/scratch3/ensembl/vac/ESTs/ests.fa",
		'estfile'        => "",

#	     path to dir where chunked ests are to be put
#            NB Sanger/EBI - this needs to be somewhere on acari!!!
#		'estfiledir'     => "/work2/vac/MGC/MGC_chunks",
#		'estfiledir'     => "/scratch3/ensembl/vac/ESTs",

#	     number of chunk files to prepare
#		'estchunknumber' => 350,
#		'estchunknumber' => 1,

#	     path to location of makeindex
		'makeindex'      => "/usr/local/ensembl/bin/makeseqindex",
#		'makeindex'      => "",


### for exonerate_ests.pl
#	     path to location of file with repeatmasked dusted genomic sequences
#	     *or* input id in form chrname.start-end
#            NB Sanger/EBI if a (large) file, this needs to be distributed across the farm or NFS will be an issue ...
#		'genomic'   => "/data/blastdb/Golden_Path_Archive/april_masked_golden_contigs.dust.fasta",
		'genomic'   => "",

### for filter_and_e2g.pl
#	     size of chunk to be processed in each job; best is 1Mb
		'filter_chunksize' => 1000000,

### for EST_GeneBuilder
#	     size of chunk to be processed in each job
		 'est_genebuilder_chunksize' => 1000000,
		);




%exonerate_conf = (
#	     path to exonerate executable
#		   'exonerate'           => "/work2/gs2/gs2/bin/exonerate-0.3d",
#		   'exonerate_args'      => "" ,
#		   'exonerate_runnable'  => "Bio::EnsEMBL::Pipeline::RunnableDB::ExonerateESTs",
#		   'exonerate_runnable'  => "",
		   
);

%est_genome_conf = (
#		    'est_genome_runnable' => "Bio::EnsEMBL::Pipeline::RunnableDB::FilterESTs_and_E2G",
#		    'est_genome_runnable' => "",

#		    location of est sequence file and index
#		    'est_index'           => "/scratch3/ensembl/vac/ESTs/ests.fa",
#		    source of ests - for putting into/retrieving from analysis
#		    'est_source'    => 'dbEST',

		    # specify the type of genes being read by EST_GeneBuilder
		    'genetype'            => "exonerate_e2g",

);

#######################################################
##
## We only allow exon-clusters with
## 
##        number_of_support_evidence > strict_lower_bound
## 
## we allow clusters with cDNA (or mRNA) evidence if they have 
## one or more exons, hence in this case 'strict_lower_bound' = 0, but we don't allow 
## single-exon clusters with only ESTs as supporting evidence, hence in this case 'strict_lower_boundt' = 1

%evidence_conf = ( 
		  		  
		  'strict_lower_bound' => 1,    # this is for ESTs only
		  #'strict_lower_bound' => 0,   # this is for cDNAs and mRNAs only

		  # specify the source_tag of the supporting evidence used in EST_GeneBuilder
		  'evidence_tag'       => "exonerate_e2g",
		 
);

#######################################################



%genomewise_conf = (
#		    'genomewise_runnable' => "Bio::EnsEMBL::Pipeline::RunnableDB::Genomewise",
		    'genomewise_runnable' => "",

		    # specify the type with which the genes are going to be written
		    'genetype'            => "genomewise",
);


# we use two databases in the EST build:
# ref_db = holds the static golden path, contig and dna information
# est_db = where we load up exonerate results into the feature table, build genes and write the exons out as features

%db_conf = (
#	    'refdbhost'      => "ecs1d",
	    'refdbhost'      => "",
	    
#	    'refdbname'      => "ens_UCSC_0801",        # where human dna lives
	    'refdbname'      => "",
	    
	    'refdbuser'      => "ensro",
#	    'refdbuser'      => "",
	    
	    'refdbpass'      => "",
	    
	    'estdbhost'      => "",
#	    'estdbhost'      => "ecs1f",

	    'estdbname'      =>
#	    'estdbname'      => "ens_UCSC_0801_est90",   # the latest est2genome results
	    
#	    'estdbuser'      => "ensadmin",
	    'estdbuser'      => "",
	    
#	    'estdbpass'      => "ensembl",
	    'estdbpass'      => "",
	    
#	    'golden_path'    => "UCSC",
	    'golden_path'    => "",
);

}

1;
