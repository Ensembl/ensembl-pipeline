package LowCoverageGeneBuildConf;

use strict;
use vars qw( %LowCoverageGeneBuildConf );

# Notes: this is mammal-specific

%LowCoverageGeneBuildConf = (
			     
	 	# User-specific configuration:
		  LC_USER		        => '', #ba1
		  LC_GeneBuilderID		=> '', #4
		  LC_DATE		        => '', #MMYYEnsembl
                  LC_ASSEMBLY_NAME              => '', #otoGar1, what broad calls the assembly
                  LC_ASSEMBLY_DATE              => '', #YYMM, when broad completed assembly
		  #perl5lib - this is the path to be added (temporarily) to the front of your perl5lib, make sure ensembl-config goes first
		  LC_PERL5LIB			=> '', #/ecs4/work1/sd3/rabbit1/ensembl-config/rabbit/rabbit_1.0:/ecs4/work1/sd3/rabbit1/ensembl/modules:/ecs4/work1/sd3/rabbit1/ensembl-pipeline/:/ecs4/work1/sd3/rabbit1/ensembl-pipeline/modules:/ecs4/work1/sd3/rabbit1/ensembl-analysis/:/ecs4/work1/sd3/rabbit1/ensembl-analysis/modules:/ecs4/work1/sd3/rabbit1/ensembl-compara/modules
		  LC_cvsDIR  			=> '',#/nfs/acari/ba1/PerlCode/
		  LC_workDIR 			=> '',#/ecs4/work3/ba1/test/
		  LC_scratchDIR			=> '', #/ecs2/scratch3/ba1/test/    Pipeline output directory
		  #LC_scratchDIR2		=> '',#/ecs4/scratch1/ba1/armadillo1/arma_1.0/
		  
		  
		  # Species-specific configuration:
		  LC_SPECIES 		=> '',#armadillo	
		  LC_BUILD_VERSION	=> '', #arma_1.0   in config dir eg. ensembl-config/armadillo/arma_1.0/	 
		  LC_GENOMICSEQS        => '',#/data/blastdb/Ensembl/Large/Dasypus_novemcinctus/BROAD/Dasypus_novemcinctus.BROAD.stdrep.fa
		  LC_CONFIG	        => '',#/nfs/acari/ba1/PerlCode/ensembl-config/
		  LC_PMATCH 		=> '', #/ecs4/work3/ba1/armadillo1/seqdata/arma_proteome.fa  location of fasta file to be written by new_prepare_proteome.pl and used in pmatch analysis
		  LC_KILL 	                => '',#/nfs/acari/ba1/PerlCode/ensembl-pipeline/scripts/GeneBuild/kill_list.txt	
		  LC_PROTEOME 			=> '',#/ecs4/work3/ba1/armadillo1/seqdata/arma_pep_all_fixed_header.fa
		  LC_PROT_INDEX 		=> '',#/ecs4/work3/ba1/armadillo1/seqdata/arma_proteome
		  LC_PEPFILE 			=> '', #/ecs4/work3/ba1/armadillo1/seqdata/ensembl_peptides/arma_pep.fa     path to dumped peptides
		  LC_CHUNKS_DIR 		=> '', #/ecs4/work3/ba1/armadillo1/seqdata/ensembl_peptides/chunks     path to directory chunks to be written
		  
		  
		  #data-specific cofiguration from BROAD and EBI
		  LC_AB_INITIO_LIB      => '', #/ecs4/work1/sd3/rabbit1/repeat_libraries/ab_initio.lib
		  LC_SUPP_LIB           => '', #/ecs4/work1/sd3/rabbit1/repeat_libraries/supplemental.lib
		  LC_ASSEMBLY_AGP	=> '', #/ecs4/work1/sd3/rabbit1/assembly/assembly.agp
		  LC_ASSEMBLY		=> '', #/ecs4/work1/sd3/rabbit1/assembly/assembly.bases #unzipped name
		  LC_SCAFFOLDS		=> '', #/ecs4/work1/sd3/rabbit1/assembly/scaffolds.fasta #unzipped name
		
		# Taxonomy
		  LC_DEFAULT	=> '', #RABBIT #use a short name - look at ensj-healthcheck/src/org/ensembl/healthcheck/Species.java for ideas
		  LC_NAME       => '', #Oryctolagus cuniculus #scientific name
                  TAXON_DBNAME  => 'ncbi_taxonomy',
		  TAXON_DBHOST  => 'ecs2',
		  TAXON_DBPORT  => '3365',

		 
	 	# database to put sequence and genes into
		  LC_DBNAME   => '',#ba1_test_code
		  LC_DBHOST   => '',#ia64g
		  LC_DBUSER   => '',#ensadmin
		  LC_DBPASS   => '',
		  LC_DBPORT   => '',#3306
		  LC_DBprefix => '',#ba1_test_
		  LC_DBro     => '',
		  LC_ANALYSIS => '',#/nfs/acari/ba1/Projects/Armadillo1/arma_1.0/analysis.conf
		  LC_RULE     => '',#/nfs/acari/ba1/Projects/Armadillo1/arma_1.0/rule.conf


		# if want the debug statements printed
		  LC_DEBUG   => 1,

		# Your choice of which RepeatMasks should be used for Pipeline
		# This is decided upon <after> all 3 repeatMask analyses are run, and <before> the 
		# RepeatMask-dependent analyses are run
		  LC_REPMASK_CHOICE		=>  [('RepeatMask','Supp_RepeatMask')], #[ ('RepeatMask','Supp_RepeatMask') ]

		# Groups of analyses
		  LC_REPEATMASKS		=> [
				 { 
				   logic_name 	=> 'RepeatMask',				  
				 },
				 { 
				   logic_name 	=> 'Supp_RepeatMask',				  
				 },
				 {
				   logic_name 	=> 'Ab_initio_RepeatMask',
				 }
				],
				 
		  LC_REPEATMASK_INDEPENDENT	=> [
		  		 {
				   logic_name 	=> 'CpG',
				 },
				 {
				   logic_name	=> 'Dust',
				 },
				 {
				   logic_name	=> 'Eponine',
				 },
				 {
				   logic_name	=> 'TRF',
				 },
				 {
				   logic_name	=> 'tRNAscan',
				 },
				 {
				 },
		  		],
		  LC_REPEATMASK_DEPENDENT	=> [
		  		 {
				   logic_name	=> 'Genscan',
				 },
				 {
				   logic_name	=> 'Unigene',
				 },
				 {
				   logic_name	=> 'Uniprot',
				 },
				 {
				   logic_name	=> 'Vertrna',
				 },
		  		],	  				 
			    );


sub import {
    my ($callpack) = caller(0); # Name of the calling package
    my $pack = shift; # Need to move package off @_

    # Get list of variables supplied, or else
    # all of GeneConf:
    my @vars = @_ ? @_ : keys( %LowCoverageGeneBuildConf );
    return unless @vars;

    # Predeclare global variables in calling package
    eval "package $callpack; use vars qw("
         . join(' ', map { '$'.$_ } @vars) . ")";
    die $@ if $@;

    foreach (@vars) {
	if ( defined $LowCoverageGeneBuildConf{ $_ } ) {
            no strict 'refs';
	    # Exporter does a similar job to the following
	    # statement, but for function names, not
	    # scalar variables:
	    *{"${callpack}::$_"} = \$LowCoverageGeneBuildConf{ $_ };
	} else {
	    die "Error: LowCoverageGeneBuildConf: $_ not known\n";
	}
    }
}

1;
