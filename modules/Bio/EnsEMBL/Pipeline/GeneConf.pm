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
	     GB_DBHOST                  => 'ecs1c',
	     GB_DBNAME                  => 'briggsae_intermediate_newschema',
	     GB_DBUSER                  => 'ensadmin',
	     GB_DBPASS                  => 'ensembl',

	     #db for writing final genewise genes to - to get round table locks
	     # this db needs to have clone & contig & static_golden_path tables populated
	     GB_FINALDBHOST             => 'ecs1c',
	     GB_FINALDBNAME             => 'briggsae_genebuild_newschema',
	     GB_INPUTID_REGEX => '',
	     # general variables
	     # path to run_GeneBuild_RunnableDB
	     GB_RUNNER      => '/acari/work4a/lec/blade_code/main_trunk/ensembl-pipeline/scripts/run_GeneBuild_RunnableDB',
	     GB_TMPDIR      => '/acari/scratch4/ensembl/lec/new_out',
	     # LSF queue plus any options you want to use
	     GB_QUEUE       => 'acari',
	     GB_TBLASTN     => '',
	     
	     # pmatch related variables - for Targetted build
	     # path to refseq fasta file
	     GB_REFSEQ      => '/usr/local/ensembl/data/blastdb/Worms/wormpep82.pep',
	     # path to swissprot fasta file
	     GB_SPTR        => '',
	     # path to swissprot "evidence kill list" file
	     GB_KILL_LIST   => '',
	     # path to file where we'll write cleaned up  proteome data
	     #GB_PFASTA      => '/acari/work4a/lec/cDNAs/prepare_cdnas.fa',
	     GB_PFASTA      => '/acari/work4a/lec/briggsae_sequence/proteome82.fa',
	     # path pmatch executable
	     GB_PMATCH      => '/usr/local/ensembl/bin/pmatch',
	     # path to directory where fpc/chromosoaml sequences are 
	     GB_FPCDIR      => '/acari/scratch4/ensembl/lec/new_contigs/',
	     # directory to write pmatch results
	     GB_PM_OUTPUT   => '/acari/scratch4/ensembl/lec/pm_out/',

	     # eg TargettedGeneE2G
	     GB_TARGETTED_RUNNABLES   => [''],
	     # eg FPC_TargettedGeneE2G
	     GB_LENGTH_RUNNABLES      => ['FPC_TargettedGeneWise'],
					  #, 'FPC_BlastMiniGenewise','Combine_Genewises_and_E2Gs', 'Gene_Builder'],
	     # size of chunk to use in length based build
	     GB_SIZE                  => '1000000',

	     # targetted genewise/geneE2G specific parameters
	     # species specific protein index
	     GB_TARGETTED_PROTEIN_INDEX => '/data/blastdb/Worms/proteome82.fa',
	     #GB_TARGETTED_PROTEIN_INDEX => '/acari/work4a/lec/briggsae_sequence/proteome82.fa',
	     GB_TARGETTED_CDNA_INDEX    => '',
	     # minimum required coverage for multiexon predictions
	     GB_TARGETTED_MULTI_EXON_COVERAGE      => '25',
	     # minimum required coverage for single predictions
	     GB_TARGETTED_SINGLE_EXON_COVERAGE     => '80',
	     # maximum allowed size of intron in Targetted gene
	     GB_TARGETTED_MAX_INTRON               => '20000',
	     # minimum coverage required to prevent splitting on long introns - keep it high!
	     GB_TARGETTED_MIN_SPLIT_COVERAGE       => '95',
	     # genetype for Targetted_GeneWise
	     #GB_TARGETTED_GW_GENETYPE              => 'targetted_cdna',
	     GB_TARGETTED_GW_GENETYPE              => 'TGE_gw',
	     #for est2genome runnabledb
	     GB_EST_DATABASES => [
                                  # fill in one complete hash for each database from which blast 
                                  # features are to be retrieved
                                  {                                  
                                   'type'       => '',
                                   'threshold'  => '',
                                   'index'      => ''
                                  },
                
                                  #                                         {
                                  #                                          'type'       => 'dbEST', #logic name of analysis type
                                  #                                          'threshold'  => '100',
                                  #                                          'index'      => '/full/path/to/dbEST'
                                  #                                         },
                                        ],
	     GB_EST_GENETYPE => 'est2genome',
	     # similairity genewise specific parameters
	     GB_SIMILARITY_DATABASES => [
					 # fill in one complete hash for each database from which blast 
					 # features are to be retrieved
					 {				  
					  'type'       => 'Swall',
					  'threshold'  => '200',
					  'index'      => '/usr/local/ensembl/data/blastdb/Ensembl/swall_270602'
					 },
# example:
					 {
					  'type'       => 'Wormpep',
					  'threshold'  => '200',
					  'index'      => '/usr/local/ensembl/data/blastdb/Worms/wormpep'
					 },
					],
	     
	     # minimum required parent protein coverage
	     GB_SIMILARITY_COVERAGE           => 70,
	     # maximum allowed size of intron 
	     GB_SIMILARITY_MAX_INTRON         => 10000,
	     # minimum coverage required to prevent splitting on long introns - keep it high!
	     GB_SIMILARITY_MIN_SPLIT_COVERAGE => 90,
	     # low complexity threshold - transcripts whose translations have low
	     # complexity % higher than GB_MAX_LOW_COMPLEXITY will be discarded
	     GB_SIMILARITY_MAX_LOW_COMPLEXITY => 60,
	     # gene type for FPC_BlastMiniGenewise
	     GB_SIMILARITY_GENETYPE           => 'similarity_genewise',

	     # Combine Genewises_and_E2Gs specific parameters
	     # gene type for Combine_Genewises_and_E2Gs
	     GB_COMBINED_GENETYPE           => 'combined_gw_e2g',
	     GB_COMBINED_MAX_INTRON         => 100000,
	     #GB_COMBINED_GENETYPE           => 'combined_cdnaorf_e2g',
	     # GeneBuilder parameters
	     GB_VCONTIG              => 1,
	     GB_SKIP_BMG             => 0,
	     GB_MIN_GENSCAN_EXONS    => 4,
	     GB_GENSCAN_MAX_INTRON   => 15000,
	     GB_FINAL_GENETYPE       => 'ensembl',

	     # Post gene build integrity checking script parameters
	     GB_MINSHORTINTRONLEN    => 7, 
	     GB_MAXSHORTINTRONLEN    => 50, 
	     GB_MINLONGINTRONLEN     => 10000, 
	     GB_MINSHORTEXONLEN      => 3, 
	     GB_MAXSHORTEXONLEN      => 10, 
	     GB_MINLONGEXONLEN       => 5000, 
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
