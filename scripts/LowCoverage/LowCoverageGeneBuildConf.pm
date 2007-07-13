package LowCoverageGeneBuildConf;

use strict;
use vars qw( %LowCoverageGeneBuildConf );

# Notes: this is mammal-specific

%LowCoverageGeneBuildConf = (
			     
	 	# User-specific configuration:
		  LC_USER		        => '', #ba1
		  LC_GeneBuilderID		=> '', #4
		  LC_DATE		        => '', #YYYY-MM-Ensembl
                  LC_ASSEMBLY_NAME              => '', #otoGar1, what broad calls the assembly
                  LC_ASSEMBLY_DATE              => '', #YYYY-MM, when broad completed assembly
		  #perl5lib - this is the path to be added (temporarily) to the front of your perl5lib, make sure ensembl-config goes first
		  LC_PERL5LIB			=> '', #/ecs4/work1/sd3/rabbit1/ensembl-config/rabbit/rabbit_1.0:/ecs4/work1/sd3/rabbit1/ensembl/modules:/ecs4/work1/sd3/rabbit1/ensembl-pipeline/modules:/ecs4/work1/sd3/rabbit1/ensembl-analysis/modules:/ecs4/work1/sd3/rabbit1/ensembl-analysis/scripts/buildchecks:/ecs4/work1/sd3/rabbit1/ensembl-compara/modules:/ecs4/work1/sd3/rabbit1/ensembl-pipeline/scripts
		  LC_cvsDIR  			=> '',#/nfs/acari/ba1/PerlCode/
		  LC_workDIR 			=> '',#/lustre/work1/ba1/test/
		  LC_scratchDIR			=> '', #/lustre/scratch1/ba1/test/    Pipeline output directory
		  
		  
		  # Species-specific configuration:
		  LC_SPECIES 		=> '',#armadillo	
		  LC_BUILD_VERSION	=> '', #arma_1.0   in config dir eg. ensembl-config/armadillo/arma_1.0/	 
		  
		  
		  #data-specific cofiguration from BROAD and EBI
		  LC_AB_INITIO_LIB      => '', #/ecs4/work1/sd3/rabbit1/repeat_libraries/ab_initio.lib
		  LC_SUPP_LIB           => '', #/ecs4/work1/sd3/rabbit1/repeat_libraries/supplemental.lib
		  LC_ASSEMBLY_AGP	=> '', #/ecs4/work1/sd3/rabbit1/assembly/assembly.agp
		  LC_ASSEMBLY		=> '', #/ecs4/work1/sd3/rabbit1/assembly/assembly.bases #unzipped name
		  LC_SCAFFOLDS		=> '', #/ecs4/work1/sd3/rabbit1/assembly/scaffolds.fasta #unzipped name
		
		# Taxonomy
		# Not used any more (assembly name used instead)LC_DEFAULT	=> '', #RABBIT #use a short name (used for assembly) - look at ensj-healthcheck/src/org/ensembl/healthcheck/Species.java for ideas
		  LC_NAME       => '', #Oryctolagus cuniculus #scientific name
                  TAXON_DBNAME  => 'ncbi_taxonomy',
		  TAXON_DBHOST  => 'ens-livemirror',
		  TAXON_DBPORT  => '3306',

		 
	 	# database to put sequence and genes into
		  LC_DBNAME   => '',#ba1_test_code
		  LC_DBHOST   => '',#genebuild1
		  LC_DBUSER   => '',#ensadmin
		  LC_DBPASS   => '',
		  LC_DBPORT   => '',#3306
		  LC_DBprefix => '',#ba1_test_
		  LC_DBro     => '',#ensro


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
