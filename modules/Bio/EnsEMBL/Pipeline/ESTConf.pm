#
# BioPerl module for Bio::EnsEMBL::Analysis::ESTConf;
#
# Cared for by EnsEMBL (ensembl-dev@ebi.ac.uk)
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::ESTConf - imports global variables used by EnsEMBL EST analysis

=head1 SYNOPSIS

    use Bio::EnsEMBL::Pipeline::ESTConf;
    use Bio::EnsEMBL::Pipeline::ESTConf qw(  );

=head1 DESCRIPTION

ESTConf is a pure ripoff of humConf written by James Gilbert.

humConf is based upon ideas from the standard perl Env environment
module.

It imports and sets a number of standard global variables into the
calling package, which are used in many scripts in the human sequence
analysis system.  The variables are first declared using "use vars",
so that it can be used when "use strict" is in use in the calling
script.  Without arguments all the standard variables are set, and
with a list, only those variables whose names are provided are set.
The module will die if a variable which doesn\'t appear in its
C<%ESTConf> hash is asked to be set.

The variables can also be references to arrays or hashes.

Edit C<%ESTConf> to add or alter variables.

All the variables are in capitals, so that they resemble environment
variables.

=head1 CONTACT

=cut


package Bio::EnsEMBL::Pipeline::ESTConf;

use strict;
use vars qw( %ESTConf );

# Hash containing config info]
%ESTConf = (
	    # general options for scripts

	    # path to run_ESTRunnableDB
	    EST_INPUTID_REGEX => '(^\S+\.\S+)\.(\d+)-(\d+)',
	    EST_RUNNER                  => '/acari/work4a/lec/blade_code/main_trunk/ensembl-pipeline/scripts/run_EST_RunnableDB', 	    
	
	    # path to run_EST_GeneBuilder.pl, script that launches EST_GeneBuilder.pm
	    # we use a different on from EST_RUNNER to use different EST_DB's
	    EST_GENE_RUNNER			=> '/acari/work4a/lec/blade_code/main_trunk/ensembl-pipeline/scripts/EST/run_EST_GeneBuilder.pl',
	    
	    # path to ensembl-pipeline/scripts/EST
	    EST_SCRIPTDIR               => '/acari/work4a/lec/blade_code/main_trunk/ensembl-pipeline/scripts/EST/', 	    
	    # where the result-directories are going to go
	    #EST_TMPDIR                  => '/ecs2/scratch1/ensembl/eae/NCBI_30/cdnas',
	    EST_TMPDIR => '/acari/scratch4/ensembl/lec/new_out/',
	    # job-queue in the farm
	    EST_QUEUE                   => 'acari',
	    
	    # make_bsubs.pl options
	    EST_EXONERATE_BSUBS         => '/acari/work4a/lec/briggsae_sequence/exonerates.dat',
	    EST_FILTER_BSUBS            => '/acari/work4a/lec/briggsae_sequence/filter_and_e2g.jobs',
	    EST_GENEBUILDER_BSUBS       => '/acari/work4a/lec/briggsae_sequence/EST_GeneBuilder.jobs',	
	    
	    # for prepare_ests.pl
	    
	    # path to executable for chunking seqeuncefiles eg fastasplit
	    EST_FILESPLITTER            => '/acari/work2/gs2/gs2/bin/fastasplit', 	    
	    
	    # path to file containign ALL ESTs/cDNAs
	    EST_FILE                    => ' /usr/local/ensembl/data/blastdb/Worms/worm_mRNA.fa', 	   

	    # path to directory where EST chunks live
	    #EST_CHUNKDIR                => '/acari/work6a/eae.tmp/Mouse/Mouse_5.3.1/riken_mRNAs/second_rechunks',

	    #EST_CHUNKDIR                => '/work6a/eae.tmp/Human/NCBI_28/ESTs/chunks',
	    EST_CHUNKDIR                => '/acari/work4a/lec/briggsae_sequence/mRNA_chunks/',

 	    # how many chunks?
	    # for NCBI_28 we have 3690891 ests, at approx. 350 ests per chunk, we estimate
	    EST_CHUNKNUMBER             => 15, 	 
	       
	    # path to makeseqindex executable
	    EST_MAKESEQINDEX            => '/usr/local/ensembl/bin/makeseqindex', 	    
	    

	    # for exonerate_ests.pl
	    #EST_GENOMIC                 => '/data/blastdb/Ensembl/NCBI_30_dusted_masked_contigs.fa',
	    EST_GENOMIC                 => '/acari/work4a/lec/briggsae_sequence/briggsae_dusted_masked_12_08_02.fa',
	    #EST_GENOMIC                 => '/data/blastdb/Ensembl/NCBI_30_dusted_masked_all_contigs.fa',

	                                            # path to file with repeatmasked dusted genomic sequences
	                                            # NB this file is huge - distribute it across the farm or 
                                                    # be prepared to face the wrath of systems when the network 
                                                    # seizes up!
	    # for filter_and_e2g.pl
	    EST_FILTERCHUNKSIZE         => 1000000,      #  we use 1000000 (ie 1MB) chunks
	    
	    # ExonerateESTs options
	    EST_EXONERATE               => '/usr/local/ensembl/bin/exonerate',      # path to executable
	    EST_EXONERATE_ARGS          => '',
	    EST_EXONERATE_RUNNABLE     => 'Bio::EnsEMBL::Pipeline::RunnableDB::ExonerateESTs',      

	    # FilterESTs_and_E2G options
	    EST_FILTER_RUNNABLE         => 'Bio::EnsEMBL::Pipeline::RunnableDB::FilterESTs_and_E2G',

	    # new index, path where the directory of the index is
	    EST_INDEX                   => '/data/blastdb/Worms/worm_mRNAs',
	    
	                                            # name of EST file that has been indexed and is 
	                                            # accessible across the farm
	                                            # for analysis(process)
      	    EST_SOURCE                  => 'worm_mRNA',      
	    

	    # for EST_GeneBuilder
	    EST_GENEBUILDER_CHUNKSIZE   => 1000000,      #  we use 1000000 (ie 1MB) chunks
	                                                 # it should be equal to EST_FILTERCHUNKSIZE
	    EST_GENEBUILDER_INPUT_GENETYPE => 'exonerate_e2g',   
	    EST_GENOMEWISE_RUNNABLE     => 'Bio::EnsEMBL::Pipeline::RunnableDB::Genomewise',
	    EST_GENOMEWISE_GENETYPE     => 'genomewise',
	    EST_EVIDENCE_TAG            => 'exonerate_e2g',

	    EST_MAX_EVIDENCE_DISCONTINUITY => 10,
	    EST_MAX_INTRON_SIZE => 10000,
	    
	    # not used at this moment
	    EST_STRICT_LOWER_BOUND      => '', # 1 for ESTs only, 0 for cDNAs/mRNAs only
	    EST_MIN_EVIDENCE_SIMILARITY => 65,
	    EST_MIN_PERCENT_ID => 70,
	    EST_MIN_COVERAGE => 90,
	    # database config
	    # IMPORTANT: make sure that all databases involved in each analysis are
	    # not on the same mysql instance 
	    # database contention arises from having tomanny db conections open to the same databas if you have more than a couple of hundread jobs contacting the same database at the same time you will need multiple database but ifyou only have a few jobs you will probably beable to get away with only 1 database
	    # ref_db - holds the static golden path, contig and dna information
	    EST_REFDBNAME               => 'briggsae_newschema',
	    EST_REFDBHOST               => 'ecs1b',
	    EST_REFDBUSER               => 'ensro',
	    EST_REFDBPASS               => '',

	    EST_GOLDEN_PATH             => 'briggsae_170602',

	    # est_db = where we load up exonerate results into the feature table
	    EST_DBNAME                  => 'briggsae_est_newschema',
	    EST_DBHOST                  => 'ecs1c',
	    EST_DBUSER                  => 'ensadmin',
	    EST_DBPASS                  => 'ensembl',
	    
	    # est_e2g_db = where we write the genes we produce from the exonerate features
	    EST_E2G_DBNAME                  => 'briggsae_est_newschema', # this should be in a different place from the one above
	    EST_E2G_DBHOST                  => 'ecs1c',
	    EST_E2G_DBUSER                  => 'ensadmin',
	    EST_E2G_DBPASS                  => 'ensembl',
	    
	    # est_gene_db = where we write the genes we produce from e2g transcripts
	    EST_GENE_DBNAME                  => 'briggsae_est_newschema',
	    EST_GENE_DBHOST                  => 'ecs1c',
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

	    # if you want to map expression data to a set of genes via ESTs, 
	    # ( see Bio::EnsEMBL::Pipeline::RunnableDB::MapGeneToExpression )
	    # you can specify here the database where those genes are
	    EST_TARGET_DBNAME            => 'ens_NCBI_30',
	    EST_TARGET_DBHOST            => 'ecs1e',
	    EST_TARGET_DBUSER            => 'ensadmin',
	    EST_TARGET_DBPASS            => 'ensembl',
	    
	    # gene type to which we are going to map the ests
	    EST_TARGET_GENETYPE          => 'ensembl',	   
	   
	    # set this to 1 if you want to map genes to ests after the genebuild
	    # this is only for human
	    MAP_GENES_TO_ESTS            => '1',
	    
	    # you can specify here the database (non ensembl schema)
	    # where the expression database is.
	    # So far we have one adaptor only for SANBI's Stack database
	    EST_EXPRESSION_DBHOST        => 'ecs2a',
	    EST_EXPRESSION_DBNAME        => 'eae_human_expression2',
	    EST_EXPRESSION_DBUSER        => 'ensro',
	    EST_EXPRESSION_DBPASS        => '',

	    # the gene2est analysis is run with this script:
	    #(put your own path here)
	    EST_EXPRESSION_RUNNER        => '/nfs/acari/eae/ensembl-branch-121/ensembl-pipeline/scripts/EST/map_ExpressionData.pl',
	    
	    # size of the chunks to run the gene to est map
	    # this should be the same one used in the GeneBuild as specified in GeneConf.pm
	    EST_EXPRESSION_CHUNKSIZE          => '5000000',		

	   );

sub import {
    my ($callpack) = caller(0); # Name of the calling package
    my $pack = shift; # Need to move package off @_

    # Get list of variables supplied, or else
    # all of ESTConf:
    my @vars = @_ ? @_ : keys( %ESTConf );
    return unless @vars;

    # Predeclare global variables in calling package
    eval "package $callpack; use vars qw("
         . join(' ', map { '$'.$_ } @vars) . ")";
    die $@ if $@;


    foreach (@vars) {
	if ( defined $ESTConf{ $_ } ) {
            no strict 'refs';
	    # Exporter does a similar job to the following
	    # statement, but for function names, not
	    # scalar variables:
	    *{"${callpack}::$_"} = \$ESTConf{ $_ };
	} else {
	    die "Error: ESTConf: $_ not known\n";
	}
    }
}

1;
