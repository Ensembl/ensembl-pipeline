#
# BioPerl module for Bio::EnsEMBL::Analysis::GeneConf;
#
# Cared for by EnsEMBL (ensembl-dev@ebi.ac.uk)
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::GeneConf - imports global variables used by EnsEMBL gene building

=head1 SYNOPSIS

    use Bio::EnsEMBL::Pipeline::GeneConf;
    use Bio::EnsEMBL::Pipeline::GeneConf qw(  );

=head1 DESCRIPTION

GeneConf is a pure ripoff of humConf written by James Gilbert.

humConf is based upon ideas from the standard perl Env environment
module.

It imports and sets a number of standard global variables into the
calling package, which are used in many scripts in the human sequence
analysis system.  The variables are first decalared using "use vars",
so that it can be used when "use strict" is in use in the calling
script.  Without arguments all the standard variables are set, and
with a list, only those variables whose names are provided are set.
The module will die if a variable which doesn\'t appear in its
C<%GeneConf> hash is asked to be set.

The variables can also be references to arrays or hashes.

Edit C<%GeneConf> to add or alter variables.

All the variables are in capitals, so that they resemble environment
variables.

=head1 CONTACT

=cut


package Bio::EnsEMBL::Pipeline::GeneConf;

use strict;
use vars qw( %GeneConf );

my $prefix='COB';

# Hash containing config info
%GeneConf = (
	     # database specific variables
	     GB_DBHOST                  => '',
#	     GB_DBNAME                  => '',
	     GB_DBNAME                  => '',
	     GB_DBUSER                  => '',
	     GB_DBPASS                  => '',

	     #db for writing final genewise genes to - to get round table locks
	     # this db needs to have clone & contig & static_golden_path tables populated
	     GB_FINALDBHOST             => '',
	     GB_FINALDBNAME             => '',

	     GB_GOLDEN_PATH             => '',

	     # general variables
	     # path to run_GeneBuild_RunnableDB
	     GB_RUNNER      => '',
	     GB_TMPDIR      => '',
	     # LSF queue plus any options you want to use
	     GB_QUEUE       => '',
	     GB_TBLASTN     => '',
	     
	     # pmatch related variables - for Targetted build
	     # path to refseq fasta file
	     GB_REFSEQ      => '',
	     # path to swissprot fasta file
	     GB_SPTR        => '',
	     # path to file where we'll write cleaned up  proteome data
	     GB_PFASTA      => '',
	     # path pmatch executable
	     GB_PMATCH      => '',
	     # path to directory where fpc/chromosoaml sequences are 
	     GB_FPCDIR      => '',
	     # directory to write pmatch results
	     GB_PM_OUTPUT   => '',

	     # eg TargettedGeneE2G
	     GB_TARGETTED_RUNNABLES   => [''],
	     # eg FPC_TargettedGeneE2G
	     GB_LENGTH_RUNNABLES      => [''],
	     # size of chunk to use in length based build
	     GB_SIZE                  => '5000000',

	     # location of sequence indices
	     GB_PROTEIN_INDEX           => '',
	     # species specific protein index
	     GB_TARGETTED_PROTEIN_INDEX => '',
	     GB_TARGETTED_CDNA_INDEX    => '',

	     # similairity genewise specific parameters
	     GB_SIMILARITY_TYPE      => 'swall',
	     GB_SIMILARITY_THRESHOLD => 200,

	     # GeneBuilder parameters
	     GB_VCONTIG              => 1,
	     GB_SKIP_BMG             => 0,

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
    my @vars = @_ ? @_ : keys( %GeneConf );
    return unless @vars;

    # Predeclare global variables in calling package
    eval "package $callpack; use vars qw("
         . join(' ', map { '$'.$_ } @vars) . ")";
    die $@ if $@;


    foreach (@vars) {
	if ( defined $GeneConf{ $_ } ) {
            no strict 'refs';
	    # Exporter does a similar job to the following
	    # statement, but for function names, not
	    # scalar variables:
	    *{"${callpack}::$_"} = \$GeneConf{ $_ };
	} else {
	    die "Error: GeneConf: $_ not known\n";
	}
    }
}

1;
