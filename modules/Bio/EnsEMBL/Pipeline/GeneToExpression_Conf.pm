#
# BioPerl module for Bio::EnsEMBL::Pipeline::GeneToExpression_Conf;
#
# Cared for by EnsEMBL (ensembl-dev@ebi.ac.uk)
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::GeneToExpression_Conf

=head1 SYNOPSIS

    use Bio::EnsEMBL::Pipeline::GeneToExpression_Conf;
    use Bio::EnsEMBL::Pipeline::GeneToExpression_Conf qw(  );

=head1 DESCRIPTION

It is a pure ripoff of humConf written by James Gilbert.

humConf is based upon ideas from the standard perl Env environment
module.

It imports and sets a number of standard global variables into the
calling package, which are used in many scripts in the human sequence
analysis system.  The variables are first declared using "use vars",
so that it can be used when "use strict" is in use in the calling
script.  Without arguments all the standard variables are set, and
with a list, only those variables whose names are provided are set.
The module will die if a variable which doesn\'t appear in its
C<%GeneToExpression_Conf> hash is asked to be set.

The variables can also be references to arrays or hashes.

Edit C<%GeneToExpression_Conf> to add or alter variables.

All the variables are in capitals, so that they resemble environment
variables.

=head1 CONTACT

=cut


package Bio::EnsEMBL::Pipeline::GeneToExpression_Conf;;

use strict;
use vars qw( %GeneToExpression_Conf );

# Hash containing config info
%GeneToExpression_Conf = (
	    # general options for scripts


	    EST_INPUTID_REGEX => '(\S+)\.(\d+)-(\d+)',
	    EST_RUNNER        => '/nfs/acari/eae/ensembl/ensembl-pipeline/scripts/EST/map_ExpressionData.pl', 	   
	    
	    # where the result-directories are going to go	
    
	    EST_TMPDIR        => '/ecs2/scratch3/ensembl/eae/NCBI_31/ests/',
	    
	    # job-queue in the farm
	    EST_QUEUE         => 'acari',
	    
	    EST_EXPRESSION_BSUBS  => '/ecs2/scratch3/ensembl/eae/NCBI_31/ests/genes2ests.jobs', 
	    
	    
	    ############################################################
	    # each runnable has an analysis
	    ############################################################


	    EST_EXPRESSION_RUNNABLE    => 'Bio::EnsEMBL::Pipeline::RunnableDB::MapGeneToExpression',
	    EST_EXPRESSION_ANALYSIS    => 'expression',

	    EST_GENEBUILDER_INPUT_GENETYPE => 'human_est',

	    ############################################################
	    # ref_db - holds the static golden path, contig and dna information
	    ############################################################
	    
	    EST_REFDBNAME               => 'ens_NCBI_31',
	    EST_REFDBHOST               => 'ecs2b',
	    EST_REFDBUSER               => 'ensro',
	    EST_REFDBPASS               => '',

	    ############################################################
	    # est_db = where we have the ests mapped to the genome as transcrtips
	    ############################################################

	    EST_DBNAME                  => 'ens_NCBI_31_est',
	    EST_DBHOST                  => 'ecs2f',
	    EST_DBUSER                  => 'ensadmin',
	    EST_DBPASS                  => 'ensembl',

	    ############################################################
	    # parameters for the map of expression data. 
	    # Currently only available for human
	    ############################################################

	    # these two parameters will reduce the number of 
            # ensembl <-> ests maps:
			  
	    # this will reject all unspliced ests		  
	    SKIP_SINGLE_EXONS              => '0',
            
	    # this will take only the ests that are present in the expression db 
	    USE_ONLY_ESTS_IN_EXPRESSION_DB => '1',	  
			  
	    # this will put a limit to the number of ests linked to a transcript
            MAX_ESTS_PER_TRANSCRIPT        => 700,

	    ############################################################		  
	    # you can specify here the database where those genes are
            ############################################################
	    EST_TARGET_DBNAME            => 'ens_NCBI_31',
	    EST_TARGET_DBHOST            => 'ecs2b',
	    EST_TARGET_DBUSER            => 'ensadmin',
	    EST_TARGET_DBPASS            => 'ensembl',
	    
	    # gene type to which we are going to map the ests
	    EST_TARGET_GENETYPE          => 'ensembl',	   
	   
	    	    
	    # you can specify here the database (non ensembl schema)
	    # where the expression database is.
	    # So far we have one adaptor only for SANBI's Stack database
	    EST_EXPRESSION_DBHOST        => 'ecs2a',
	    EST_EXPRESSION_DBNAME        => 'eae_human_expression2',
	    EST_EXPRESSION_DBUSER        => 'ensadmin',
	    EST_EXPRESSION_DBPASS        => 'ensembl',

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
    # all of GeneToExpression_Conf
    my @vars = @_ ? @_ : keys( %GeneToExpression_Conf );
    return unless @vars;

    # Predeclare global variables in calling package
    eval "package $callpack; use vars qw("
         . join(' ', map { '$'.$_ } @vars) . ")";
    die $@ if $@;


    foreach (@vars) {
	if ( defined $GeneToExpression_Conf{ $_ } ) {
            no strict 'refs';
	    # Exporter does a similar job to the following
	    # statement, but for function names, not
	    # scalar variables:
	    *{"${callpack}::$_"} = \$GeneToExpression_Conf{ $_ };
	} else {
	    die "Error: GeneToExpression_Conf: $_ not known\n";
	}
    }
}

1;
