#
# module for Bio::EnsEMBL::Analysis::GeneCombiner;
#
# Cared for by EnsEMBL (ensembl-dev@ebi.ac.uk)
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::Config::GeneBuild::GeneCombiner - imports global variables used by EnsEMBL gene building

=head1 SYNOPSIS

    use Bio::EnsEMBL::Pipeline::Config::GeneBuild::GeneCombiner;
    use Bio::EnsEMBL::Pipeline::Config::GeneBuild::GeneCombiner qw(  );

=head1 DESCRIPTION

It is a pure ripoff of humConf written by James Gilbert.

humConf is based upon ideas from the standard perl Env environment
module.

It imports and sets a number of standard global variables into the
calling package, which are used in many scripts in the human sequence
analysis system.  The variables are first decalared using "use vars",
so that it can be used when "use strict" is in use in the calling
script.  Without arguments all the standard variables are set, and
with a list, only those variables whose names are provided are set.
The module will die if a variable which doesn\'t appear in its
C<%GeneCombiner> hash is asked to be set.

The variables can also be references to arrays or hashes.

Edit C<%GeneConmbinerConf> to add or alter variables.

All the variables are in capitals, so that they resemble environment
variables.

=head1 CONTACT

=cut


package Bio::EnsEMBL::Pipeline::Config::GeneBuild::GeneCombiner;

use strict;
use vars qw( %GeneCombiner );

my $prefix='COB';

# Hash containing config info
%GeneCombiner = (
		     ############################################################
		     # EST GENES PROPERTIES
		     ############################################################
		     
		     # db from where we get the cdna-genes
		     ESTGENE_DBHOST              => '',
		     ESTGENE_DBUSER              => '', 
		     ESTGENE_DBNAME              => '',
		     ESTGENE_DBPASS              => '',
		     ESTGENE_DBPORT              => '',
		     # gene type for est-genes
		     ESTGENE_TYPE                => '',
		     
		     # in case you want to restrict the length of the introns
		     ESTGENE_MAX_INTRON_LENGTH    => 200001,
		     
		     # db from which we get ensembl genes
		     ENSEMBL_DBHOST              => '',
		     ENSEMBL_DBUSER              => '',
		     ENSEMBL_DBNAME              => '',
		     ENSEMBL_DBPASS              => '',
		     ENSEMBL_DBPORT              => '',
		     # gene type for the ensembl genes
		     ENSEMBL_TYPE                => 'ensembl',
		     		     
		     # refdb where the dna is, so that we do not need to have it everywhere
		     REF_DBHOST              => '',
		     REF_DBUSER              => '',
		     REF_DBNAME              => '',
		     REF_DBPASS              => '',
		     REF_DBPORT              => '',
		     
		     # different db for writing final genes to - to get round table locks
		     # this db needs to have clone & contig & static_golden_path tables populated
		     FINAL_DBHOST             => '',
		     FINAL_DBNAME             => '',
		     FINAL_DBUSER             => '',
		     FINAL_DBPASS             => '',
		     FINAL_DBPOT             => '',
		     
		     # final gene type
		     FINAL_TYPE               => 'ensembl',
		     
		     ############################################################
		     # general variables
		     ############################################################

		     # the input id should give chr_name . start - end
		     GENECOMBINER_INPUTID_REGEX => '(^\S+)\.(\d+)-(\d+)',



		     # path to run_GeneBuild_RunnableDB
		     RUNNER      => 'ensembl-pipeline/scripts/run_GeneCombiner.pl',
		     GENECOMBINER_RUNNABLES    =>  [
						    {
						     runnable => 'GeneCombiner',
						     analysis => 'ensembl',
						    },
						   ],
		     
		     # directory where the output files will go
		     OUTPUT_DIR  => '',
		     
		     # directory where the jobs files will go
		     JOBS_DIR    => '',
		     
		     # LSF queue plus any options you want to use
		     QUEUE       => 'acari',
		     
		     # size of slices to use in length based build
		     SLICE_SIZE                  => '5000000',
		     
		     # GeneBuilder parameters
		     GC_VCONTIG              => 1,
		     GC_SKIP_BMG             => 0,


		     # Post gene build integrity checking script parameters
		     GB_MINSHORTINTRONLEN    => 7, 
		     GB_MAXSHORTINTRONLEN    => 50, 
		     GB_MINLONGINTRONLEN     => 100000, 
		     GB_MINSHORTEXONLEN      => 3, 
		     GB_MAXSHORTEXONLEN      => 10, 
		     GB_MINLONGEXONLEN       => 50000, 
		     GB_MINTRANSLATIONLEN    => 10, 
		     GB_MAX_EXONSTRANSCRIPT  => 150, 
		     GB_MAXTRANSCRIPTS       => 10, 
		     GB_IGNOREWARNINGS       => 1, 
		     
		    );

sub import {
    my ($callpack) = caller(0); # Name of the calling package
    my $pack = shift; # Need to move package off @_

    # Get list of variables supplied, or else
    # all of GeneConf:
    my @vars = @_ ? @_ : keys( %GeneCombiner );
    return unless @vars;

    # Predeclare global variables in calling package
    eval "package $callpack; use vars qw("
         . join(' ', map { '$'.$_ } @vars) . ")";
    die $@ if $@;


    foreach (@vars) {
	if ( defined $GeneCombiner{ $_ } ) {
            no strict 'refs';
	    # Exporter does a similar job to the following
	    # statement, but for function names, not
	    # scalar variables:
	    *{"${callpack}::$_"} = \$GeneCombiner{ $_ };
	} else {
	    die "Error: GeneCombiner: $_ not known\n";
	}
    }
}

1;



