#
# module for Bio::EnsEMBL::Analysis::GeneCombinerConf;
#
# Cared for by EnsEMBL (ensembl-dev@ebi.ac.uk)
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::GeneCombinerConf - imports global variables used by EnsEMBL gene building

=head1 SYNOPSIS

    use Bio::EnsEMBL::Pipeline::GeneCombinerConf;
    use Bio::EnsEMBL::Pipeline::GeneCombinerConf qw(  );

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
C<%GeneCombinerConf> hash is asked to be set.

The variables can also be references to arrays or hashes.

Edit C<%GeneConmbinerConf> to add or alter variables.

All the variables are in capitals, so that they resemble environment
variables.

=head1 CONTACT

=cut


package Bio::EnsEMBL::Pipeline::GeneCombinerConf;

use strict;
use vars qw( %GeneCombinerConf );

my $prefix='COB';

# Hash containing config info
%GeneCombinerConf = (
		     # db from where we get the est-genes
		     ESTGENE_DBHOST              => 'ecs1f',
		     ESTGENE_DBUSER              => 'ensro', 
		     ESTGENE_DBNAME              => 'ens_NCBI_30_cdna_genes',
		     ESTGENE_DBPASS              => 'ensembl',
		     
		     # gene type for est-genes
		     ESTGENE_TYPE                => 'genomewise',
		     
		     # db from which we get ensembl genes
		     ENSEMBL_DBHOST              => 'ecs1a',
		     ENSEMBL_DBUSER              => 'ensro',
		     ENSEMBL_DBNAME              => 'ens_NCBI_30_final_build',
		     ENSEMBL_DBPASS              => 'ensembl',
		     
		     # gene type for the ensembl genes
		     ENSEMBL_TYPE                => 'ensembl',
		     
		     # refdb where the dna is, so that we do not need to have it everywhere
		     REF_DBHOST              => 'ecs1e',
		     REF_DBUSER              => 'ensro',
		     REF_DBNAME              => 'ens_NCBI_30',
		     REF_DBPASS              => 'ensembl',
		     
		     # different db for writing final genes to - to get round table locks
		     # this db needs to have clone & contig & static_golden_path tables populated
		     FINAL_DBHOST             => 'ecs1c',
		     FINAL_DBNAME             => 'genecombiner_test_NCBI_30',
		     FINAL_DBUSER             => 'ensadmin',
		     FINAL_DBPASS             => 'ensembl',
		     
		     # final gene type
		     FINAL_TYPE               => 'final_ensembl',
		     
		     # general variables
		     # path to run_GeneBuild_RunnableDB
		     RUNNER      => '/nfs/acari/eae/ensembl-branch-121/ensembl-pipeline/scripts/run_GeneCombiner.pl',
		     RUNNABLE    => 'GeneCombiner',
		     
		     # directory where the output files will go
		     OUTPUT_DIR  => '/ecs2/scratch1/ensembl/eae/NCBI_30/test_genecombiner',
		     
		     # directory where the jobs files will go
		     JOBS_DIR    => '/ecs2/scratch1/ensembl/eae/NCBI_30/test_genecombiner',
		     
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
		     
		     # old id generating variables
		     EXON_ID_SUBSCRIPT       => $prefix.'E',
		     EXON_ID_DIGITS          => 11,
		     TRANSCRIPT_ID_SUBSCRIPT => $prefix.'T',
		     TRANSCRIPT_ID_DIGITS    => 11,
		     GENE_ID_SUBSCRIPT       => $prefix.'G',
		     GENE_ID_DIGITS          => 11,
		     PROTEIN_ID_SUBSCRIPT    => $prefix.'P',
		     PROTEIN_ID_DIGITS       => 11 ,
		     
		    );

sub import {
    my ($callpack) = caller(0); # Name of the calling package
    my $pack = shift; # Need to move package off @_

    # Get list of variables supplied, or else
    # all of GeneConf:
    my @vars = @_ ? @_ : keys( %GeneCombinerConf );
    return unless @vars;

    # Predeclare global variables in calling package
    eval "package $callpack; use vars qw("
         . join(' ', map { '$'.$_ } @vars) . ")";
    die $@ if $@;


    foreach (@vars) {
	if ( defined $GeneCombinerConf{ $_ } ) {
            no strict 'refs';
	    # Exporter does a similar job to the following
	    # statement, but for function names, not
	    # scalar variables:
	    *{"${callpack}::$_"} = \$GeneCombinerConf{ $_ };
	} else {
	    die "Error: GeneCombinerConf: $_ not known\n";
	}
    }
}

1;
