# Cared for by EnsEMBL (ensembl-dev@ebi.ac.uk)
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::EST_GeneBuilder_Conf - imports global variables used by EnsEMBL EST analysis

=head1 SYNOPSIS

    use Bio::EnsEMBL::Pipeline::EST_GeneBuilder_Conf;
 
=head1 DESCRIPTION

This class is a pure ripoff of humConf written by James Gilbert.

humConf is based upon ideas from the standard perl Env environment
module.

It imports and sets a number of standard global variables into the
calling package, which are used in many scripts in the human sequence
analysis system.  The variables are first declared using "use vars",
so that it can be used when "use strict" is in use in the calling
script.  Without arguments all the standard variables are set, and
with a list, only those variables whose names are provided are set.
The module will die if a variable which doesn\'t appear in its
C<%EST_GeneBuilder_Conf> hash is asked to be set.

The variables can also be references to arrays or hashes.

Edit C<%EST_GeneBuilder_Conf> to add or alter variables.

All the variables are in capitals, so that they resemble environment
variables.

=head1 CONTACT

=cut


package Bio::EnsEMBL::Pipeline::EST_GeneBuilder_Conf;

use strict;
use vars qw( %EST_GeneBuilder_Conf );

# Hash containing config info]
%EST_GeneBuilder_Conf = (
	    # general options for scripts

	    # path to run_ESTRunnableDB
	    	    
	    # for Rat
	    EST_INPUTID_REGEX => '(\S+)\.(\d+)-(\d+)',
	    
	    # path to run_EST_GeneBuilder.pl, script that launches EST_GeneBuilder.pm
	    EST_GENE_RUNNER   => '/nfs/acari/eae/ensembl/ensembl-pipeline/scripts/EST/run_EST_GeneBuilder.pl',
	    
	    	    
	    # where the result-directories are going to go	
    	    EST_TMPDIR        => '/ecs2/scratch3/ensembl/eae/NCBI_31/ests',
	    
	    # job-queue in the farm
	    EST_QUEUE         => 'acari',
	    
	    EST_GENEBUILDER_BSUBS => '/ecs2/scratch3/ensembl/eae/NCBI_31/ests/EST_GeneBuilder.jobs',
	    EST_EXPRESSION_BSUBS  => '/ecs2/scratch3/ensembl/eae/NCBI_31/ests/genes2ests.jobs', 
	    
	    ############################################################
	    # each runnable has an analysis
	    ############################################################
	    
	    EST_GENEBUILDER_RUNNABLE   => 'Bio::EnsEMBL::Pipeline::RunnableDB::EST_GeneBuilder',
	    EST_GENEBUILDER_ANALYSIS   => 'genomewise',

	    EST_EXPRESSION_RUNNABLE    => 'Bio::EnsEMBL::Pipeline::RunnableDB::MapGeneToExpression',
	    EST_EXPRESSION_ANALYSIS    => 'expression',

	    ############################################################
	    # EST_GeneBuilder
	    ############################################################
			 
	    EST_GENEBUILDER_CHUNKSIZE        => 1000000,      #  we use 1000000 (ie 1MB) chunks
		  
	    EST_GENEBUILDER_INPUT_GENETYPE => 'human_est',
	    EST_GENOMEWISE_GENETYPE        => 'genomewise',
	    EST_EVIDENCE_TAG               => 'human_est',
			 
			 
	                 ## you must choose one type of merging for cdnas/ests: 2 and 3 are the common ones
			 EST_GENEBUILDER_COMPARISON_LEVEL => 4,
			 
			 # for details see documentation 
			 # in Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptComparator
			 # 1 --> strict: exact exon matching (unrealistic). 
			 # 2 --> allow edge exon mismatches
			 # 3 --> allow internal mismatches
			 # 4---> allow intron mismatches
			 
			 # you can alow a mismatch in the splice sites
			 EST_GENEBUILDER_SPLICE_MISMATCH  => 10,
			 
			 # you can alow to bridge over small introns:
			 EST_GENEBUILDER_INTRON_MISMATCH  => 10,
			 
			 # you can choose whether you only want tw ests/cdnas to merge if
			 # they have the same number of exons
			 EST_GENEBUILDER_EXON_MATCH     => 0,
  	    	    
			 # how much discontinuity we allow in the supporting evidence
			 # this might be correlated with the 2-to-1 merges, so we
			 # usually put it =  EST_GENEBUILDER_INTRON_MISMATCH for ESTs
			 EST_MAX_EVIDENCE_DISCONTINUITY  => 10,
			 REJECT_SINGLE_EXON_TRANSCRIPTS  => 1,
			 GENOMEWISE_SMELL                => 1,
			 
			 # exons smaller than this will not be include in the merging algorithm
			 EST_MIN_EXON_SIZE               => 8,

			 # ests with intron bigger than this will not be incuded either
			 EST_MAX_INTRON_SIZE             => 200000,
			 
			 # this says to ClusterMerge what's the minimum
			 # number of ESTs/cDNAs that must be 'included' into a
			 # final transcript
			 CLUSTERMERGE_MIN_EVIDENCE_NUMBER => 1,
			 
	    # database config
	    # IMPORTANT: make sure that all databases involved in each analysis are
	    # not on the same mysql instance 
	    # database contention arises from having too many db conections open to the same database
            # if you have more than a couple of hundred jobs contacting the same database at the same 
            # time you will need multiple database but ifyou only have a few jobs you will probably be 
	    # able to get away with only 1 database.
	    
	    ############################################################
	    # ref_db - holds the static golden path, contig and dna information
	    ############################################################
	    
	    EST_REFDBNAME               => 'ens_NCBI_31',
	    EST_REFDBHOST               => 'ecs2b',
	    EST_REFDBUSER               => 'ensro',
	    EST_REFDBPASS               => '',

	    # est_e2g_db = where we write the genes we produce from the exonerate features
	    # this is in general the database where we read the mapped ests/cdnas
	    # from, in order to use them for the EST_GeneBuilder
			 
	    EST_E2G_DBNAME                  => 'ens_NCBI_31_est',
	    EST_E2G_DBHOST                  => 'ecs2f',
	    EST_E2G_DBUSER                  => 'ensro',
	    EST_E2G_DBPASS                  => '',

	    # est_gene_db = where we write the genes we produce from e2g transcripts
	    EST_GENE_DBNAME                  => 'ens_NCBI_31_estgene',
	    EST_GENE_DBHOST                  => 'ecs2c',
	    EST_GENE_DBUSER                  => 'ensadmin',
	    EST_GENE_DBPASS                  => 'ensembl',
	    
	    # if you want to use ests together with cdnas in EST_GeneBuilder
	    # and your cdnas are in a SEPARATE DATABASE, you can specify it here:
	    USE_cDNA_DB                  => 0,  # set it to a def/undef value if you do/don't want to use it
	    
	    cDNA_DBNAME                  => '',
	    cDNA_DBHOST                  => '',
	    cDNA_DBUSER                  => '',
	    cDNA_DBPASS                  => '',
	    cDNA_GENETYPE                => '',	  

	    
	    ############################################################
	    # parameters for the map of expression data. 
	    # Currently only available for human
	    ############################################################



	    # if you want to map expression data to a set of genes via ESTs, 
	    # ( see Bio::EnsEMBL::Pipeline::RunnableDB::MapGeneToExpression )
	    # set this to 1 if you want to map genes to ests after the genebuild
	    # this is only for human
	    MAP_GENES_TO_ESTS            => '1',

	    # you can specify here the database where those genes are
	    EST_TARGET_DBNAME            => 'ens_NCBI_31',
	    EST_TARGET_DBHOST            => 'ecs2b',
	    EST_TARGET_DBUSER            => 'ensro',
	    EST_TARGET_DBPASS            => '',
	    
	    # gene type to which we are going to map the ests
	    EST_TARGET_GENETYPE          => 'ensembl',
	   
	    	    
	    # you can specify here the database (non ensembl schema)
	    # where the expression database is.
	    # So far we have one adaptor only for SANBI's Stack database
	    EST_EXPRESSION_DBHOST        => 'ecs2a',
	    EST_EXPRESSION_DBNAME        => 'eae_human_expression2',
	    EST_EXPRESSION_DBUSER        => 'ensro',
	    EST_EXPRESSION_DBPASS        => '',

	    # the gene2est analysis is run with this script:
	    #(put your own path here)
	    EST_EXPRESSION_RUNNER        => '/nfs/acari/eae/ensembl/ensembl-pipeline/scripts/EST/map_ExpressionData.pl',
	    
	    # size of the chunks to run the gene to est map
	    # this should be the same one used in the GeneBuild as specified in GeneConf.pm
	    EST_EXPRESSION_CHUNKSIZE          => '5000000',		

	   );

sub import {
    my ($callpack) = caller(0); # Name of the calling package
    my $pack = shift; # Need to move package off @_

    # Get list of variables supplied, or else
    # all of EST_GeneBuilder_Conf:
    my @vars = @_ ? @_ : keys( %EST_GeneBuilder_Conf );
    return unless @vars;

    # Predeclare global variables in calling package
    eval "package $callpack; use vars qw("
         . join(' ', map { '$'.$_ } @vars) . ")";
    die $@ if $@;


    foreach (@vars) {
	if ( defined $EST_GeneBuilder_Conf{ $_ } ) {
            no strict 'refs';
	    # Exporter does a similar job to the following
	    # statement, but for function names, not
	    # scalar variables:
	    *{"${callpack}::$_"} = \$EST_GeneBuilder_Conf{ $_ };
	} else {
	    die "Error: EST_GeneBuilder_Conf: $_ not known\n";
	}
    }
}

1;
