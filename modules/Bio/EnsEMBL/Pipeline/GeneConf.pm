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
	     # DATABASE SPECIFIC VARIABLES:

	     # database containing the blast and pmatch features
	     GB_DBHOST                  => 'ecs1e',
	     GB_DBNAME                  => 'ens_NCBI_30',
	     GB_DBUSER                  => 'ensadmin',
	     GB_DBPASS                  => 'ensembl',

	     # database containing the genewise genes (TGE_gw,similarity_genewise)
	     GB_GW_DBHOST                  => 'ecs1b',
	     GB_GW_DBNAME                  => 'ens_NCBI_30_genewises',
	     GB_GW_DBUSER                  => 'ensadmin',
	     GB_GW_DBPASS                  => 'ensembl',

	     # database where the combined_gw_e2g genes will be stored
	     # IMPORTANT: we should have copied the genewise genes to this db before hand:
	     GB_COMB_DBHOST                  => 'ecs1f',
	     GB_COMB_DBNAME                  => 'ens_NCBI_30_combined_genes',
	     GB_COMB_DBUSER                  => 'ensro',
	     GB_COMB_DBPASS                  => '',
	     
	     # database containing the cdnas mapped, to be combined with the genewises
	     # by putting this info here, we free up ESTConf.pm so that two analysis can
	     # be run at the same time
	     GB_cDNA_DBHOST                  => 'ecs1a',
	     GB_cDNA_DBNAME                  => 'ens_NCBI_30_cdna',
	     GB_cDNA_DBUSER                  => 'ensro',
	     GB_cDNA_DBPASS                  => '',
	     
	     #db for writing final genewise genes to - to get round table locks
	     # this db needs to have clone & contig & static_golden_path tables populated
	     GB_FINALDBHOST             => 'ecs1a',
	     GB_FINALDBNAME             => 'ens_NCBI_30_final_build',
	     
	     # general variables
	     # path to run_GeneBuild_RunnableDB
	     GB_RUNNER      => '/nfs/acari/eae/ensembl-branch-121/ensembl-pipeline/scripts/run_GeneBuild_RunnableDB',
	     GB_TMPDIR      => '/ecs2/scratch1/ensembl/eae/NCBI_30/',
	     # LSF queue plus any options you want to use
	     GB_QUEUE       => 'acari',
	     GB_TBLASTN     => '',
	     
	     # pmatch related variables - for Targetted build
	     # path to refseq fasta file
	     GB_REFSEQ      => '/acari/work6a/eae.tmp/Human/NCBI_29/proteome/refseq.fa',
	     # path to swissprot fasta file
	     GB_SPTR        => '/acari/work6a/eae.tmp/Human/NCBI_29/proteome/sptr.fa',
	     # path to file where we'll write cleaned up  proteome data
	     GB_PFASTA      => '/acari/work6a/eae.tmp/Human/NCBI_29/proteome/proteome.fa',
	     # path pmatch executable
	     GB_PMATCH      => '/nfs/acari/eae/rd-utils/pmatch',
	     # path to directory where fpc/chromosomal sequences are 
	     GB_FPCDIR      => '/acari/scratch5/ensembl/eae/NCBI_29/genome/',
	     # directory to write pmatch results
	     GB_PM_OUTPUT   => '/acari/scratch5/ensembl/eae/NCBI_29/pmatch/',

	     # eg TargettedGeneE2G
	     GB_TARGETTED_RUNNABLES   => [''],
	     # eg FPC_TargettedGeneE2G
	     GB_LENGTH_RUNNABLES      => ['FPC_TargettedGeneWise', 'FPC_BlastMiniGenewise','Combine_Genewises_and_E2Gs', 'Gene_Builder'],
	     # size of chunk to use in length based build
	     GB_SIZE                  => '5000000',

	     # location of sequence indices
	     GB_PROTEIN_INDEX           => '',
	     # species specific protein index
	     GB_TARGETTED_PROTEIN_INDEX => '/data/blastdb/Ensembl/NCBI_30_proteome.fa.jidx',
	     
	     # one of the available modules: 
	     # Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch
	     # Bio::EnsEMBL::Pipeline::SeqFetcher::Getseqs
	     # Bio::EnsEMBL::Pipeline::SeqFetcher::OBDAIndexSeqFetcher
	     #
	     GB_TARGETTED_PROTEIN_SEQFETCHER => '',

	     #GB_TARGETTED_PROTEIN_INDEX => '/acari/work6a/eae.tmp/Human/NCBI_29/proteome/NCBI_29_proteome.fa',
	     GB_TARGETTED_CDNA_INDEX    => '',

	     # targetted genewise/geneE2G specific parameters
	     # minimum required coverage for multiexon predictions
	     GB_TARGETTED_MULTI_EXON_COVERAGE      => '25',
	     # minimum required coverage for single predictions
	     GB_TARGETTED_SINGLE_EXON_COVERAGE     => '80',
	     # maximum allowed size of intron in Targetted gene
	     GB_TARGETTED_MAX_INTRON               => '250000',
	     # minimum coverage required to prevent splitting on long introns - keep it high!
	     GB_TARGETTED_MIN_SPLIT_COVERAGE       => '95',
	     # genetype for Targetted_GeneWise
	     GB_TARGETTED_GW_GENETYPE              => 'TGE_gw',

	     # similairity genewise specific parameters
	     GB_SIMILARITY_DATABASES => [
					 # fill in one complete hash for each database from which blast 
					 # features are to be retrieved
					 {				  
					  'type'       => 'swall',
					  'threshold'  => '200',
					  'index'      => '/data/blastdb/Ensembl/swall_indicate_index/',
					  'seqfetcher' => 'Bio::EnsEMBL::Pipeline::SeqFetcher::OBDAIndexSeqFetcher'
 					 },
					],
	     
	     # minimum required parent protein coverage
	     GB_SIMILARITY_COVERAGE           => 70,
	     # maximum allowed size of intron 
	     GB_SIMILARITY_MAX_INTRON         => 150000,
	     # minimum coverage required to prevent splitting on long introns - keep it high!
	     GB_SIMILARITY_MIN_SPLIT_COVERAGE => 90,
	     # gene type for FPC_BlastMiniGenewise
	     GB_SIMILARITY_GENETYPE           => 'similarity_genewise',

	     # low complexity threshold - transcripts whose translations have low
	     # complexity % higher than GB_MAX_LOW_COMPLEXITY will be discarded
	     GB_SIMILARITY_MAX_LOW_COMPLEXITY => 100,
	     
	     # Combine Genewises_and_E2Gs specific parameters
	     # gene type for Combine_Genewises_and_E2Gs
	     GB_COMBINED_GENETYPE           => 'combined_gw_e2g',
	     GB_COMBINED_MAX_INTRON         => 100000,

	     # GeneBuilder parameters
	     GB_VCONTIG              => 1,
	     GB_SKIP_BMG             => 0,
	     GB_MIN_GENSCAN_EXONS    => 4,
	     GB_GENSCAN_MAX_INTRON   => 150000,
	     GB_FINAL_GENETYPE       => 'ensembl',

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
