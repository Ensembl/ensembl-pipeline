# Copyright GRL & EBI 2001
# Author: Val Curwen
# Creation: 19.07.2001

# configuration information for EST scripts
# give useful keynames to things

# I've left in sample entries for the various options to hopefully make this easier to use

BEGIN {
package main;

%EST_conf = ( 
### general options
#	     'runner'      => '/work2/vac/dev/ensembl-pipeline/scripts/test_RunnableDB',
	     'runner'      => '',

#            'exonerate_runnable'  => "Bio::EnsEMBL::Pipeline::RunnableDB::ExonerateESTs",
	     'exonerate_runnable'  => "",

#            'est_genome_runnable' => "Bio::EnsEMBL::Pipeline::RunnableDB::FilterESTs_and_E2G",
	     'est_genome_runnable' => "",

#            'genomewise_runnable' => "Bio::EnsEMBL::Pipeline::RunnableDB::Genomewise",
	     'genomewise_runnable' => "",

#	     path to scratch directory where output subdirectories and files will be placed
	     'tmpdir'      => '/scratch4/ensembl/vac',

# we use two databases in the EST build:
# ref_db = holds the static golden path, contig and dna information
# est_db = where we load up exonerate results into the feature table, build genes and write the exons out as features

#	     'refdbhost'      => "ecs1d",
	     'refdbhost'      => "",

#	     'refdbname'      => "ens_apr01",
	     'refdbname'      => "",

#	     'refdbuser'      => "ensro",
	     'refdbuser'      => "",

#	     'estdbhost'      => "ecs1f",
	     'estdbhost'      => "",

#	     'estdbname'      => "exonerate_est",
	     'estdbname'      => "",

#	     'estdbuser'      => "ensadmin",
	     'estdbuser'      => "",

#	     'queue'          => "acarilong",
	     'queue'          => "",

### for make_bsubs.pl - where to put the bsubs
             'exonerate_bsubsfile' => "/scratch4/ensembl/vac/exonerate_est.jobs",
             'filter_bsubsfile'    => "/scratch4/ensembl/vac/filter_and_e2g.jobs",


### for prepare_ests.pl

#	     path to executable which will be used for splitting estfiles into chunks
#	     'filesplitter' => "/work2/gs2/gs2/bin/fastasplit",
	     'filesplitter' => "",

#	     path to file with all the ESTs/cDNAs/whatever in it
#	     'estfile'        => "/work2/vac/MGC/data/MGC_Hs.fa",
	     'estfile'        => "",

#	     path to dir where chunked ests are to be put
#            NB Sanger/EBI - this needs to be somewhere on acari!!!
#            'estfiledir'     => "/work2/vac/MGC/MGC_chunks",
             'estfiledir'     => "",

#	     number of chunk files to prepare
#	     'estchunknumber' => 350,
	     'estchunknumber' => 1,

#	     path to location of makeindex
#	     'makeindex'      => "/usr/local/ensembl/bin/makeindex",
	     'makeindex'      => "",


### for exonerate_ests.pl

#	     path to exonerate executable
#            'exonerate' => "/work2/gs2/gs2/bin/exonerate-0.3d",
             'exonerate' => "",

#	     path to location of file with repeatmasked dusted genomic sequences
#	     *or* input id in form chrname.start-end
#            NB Sanger/EBI if a (large) file, this needs to be distributed across the farm or NFS will be an issue ...
#	     'genomic'   => "/data/blastdb/Golden_Path_Archive/april_masked_golden_contigs.dust.fasta",
	     'genomic'   => "",


### for filter_and_e2g.pl
# where does the index to be used to fetch EST seqs live?
#	     'estindex'  => "/work2/vac/MGC/data/MGC_Hs.fa",
	     'estindex'  => "",

#	     size of chunk to be processed in each job; best is 1Mb
	     'filter_chunksize' => 1000000,
	   );


}

1;
