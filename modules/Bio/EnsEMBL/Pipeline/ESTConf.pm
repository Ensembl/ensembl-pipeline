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
analysis system.  The variables are first decalared using "use vars",
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

# Hash containing config info
%ESTConf = (
	    # general options for scripts

	    EST_RUNNER                  => '', 	    # path to run_ESTRunnableDB
	    EST_SCRIPTDIR               => '', 	    # path to ensembl-pipeline/scripts/EST
	    EST_TMPDIR                  => '',
	    EST_QUEUE                   => '',

	    # make_bsubs.pl options
	    EST_EXONERATE_BSUBS         => '',
	    EST_FILTER_BSUBS            => '',
	    EST_GENEBUILDER_BSUBS       => '',

	    # for prepare_ests.pl
	    EST_FILESPLITTER            => '', 	    # path to executable for chunking seqeuncefiles eg fastasplit
	    EST_FILE                    => '', 	    # path to file containign ALL ESTs/cDNAs
	    EST_CHUNKDIR                => '', 	    # path to directory where EST chunks live
	    EST_CHUNKNUMBER             => '', 	    # how many chunks?
	    EST_MAKESEQINDEX            => '', 	    # path to makeseqindex executable
	    
	    # for exonerate_ests.pl
	    EST_GENOMIC                 => '',      # path to file with repeatmasked dusted genomic sequences
	                                            # NB this file is huge - distribute it across the farm or 
                                                    # be prepared to face the wrath of systems when the network 
                                                    # seizes up!
	    # for filter_and_e2g.pl
	    EST_FILTERCHUNKSIZE         => '',      #  we use 1000000 (ie 1MB) chunks
	    
	    # ExonerateESTs options
	    EST_EXONERATE               => '',      # path to executable
	    EST_EXONERATE_ARGS          => '',
	    EST_EXONERATE_RUNNABLE     => 'Bio::EnsEMBL::Pipeline::RunnableDB::ExonerateESTs',      

	    # FilterESTs_and_E2G options
	    EST_FILTER_RUNNABLE         => 'Bio::EnsEMBL::Pipeline::RunnableDB::FilterESTs_and_E2G',
	    EST_INDEX                   => '',      # name of EST file that has been indexed and is 
	                                            # accessible across the farm
	    EST_SOURCE                  => '',      # for analysis(process)
	    

	    # for EST_GeneBuilder
	    EST_GENEBUILDER_CHUNKSIZE   => '',      #  we use 1000000 (ie 1MB) chunks
	    EST_GENEBUILDER_INPUT_GENETYPE => 'exonerate_e2g',   #
	    EST_GENOMEWISE_RUNNABLE     => 'Bio::EnsEMBL::Pipeline::RunnableDB::Genomewise',
	    EST_GENOMEWISE_GENETYPE     => 'genomewise',
	    EST_STRICT_LOWER_BOUND      => '', # 1 for ESTs only, 0 for cDNAs/mRNAs only
	    EST_EVIDENCE_TAG            => 'exonerate_e2g',
	    EST_MIN_EVIDENCE_SIMILARITY => '',
	    EST_MAX_EVIDENCE_DISCONTINUITY => '',

	    # database config
	    # ref_db - holds the static golden path, contig and dna information
	    # est_db = where we load up exonerate results into the feature table and build genes
	    EST_REFDBHOST               => '',
	    EST_REFDBNAME               => '',
	    EST_REFDBUSER               => '',
	    EST_REFDBPASS               => '',
	    EST_DBNAME                  => '',
	    EST_DBHOST                  => '',
	    EST_DBUSER                  => '',
	    EST_DBPASS                  => '',
	   );

sub import {
    my ($callpack) = caller(0); # Name of the calling package
    my $pack = shift; # Need to move package off @_

    # Get list of variables supplied, or else
    # all of GeneConf:
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
