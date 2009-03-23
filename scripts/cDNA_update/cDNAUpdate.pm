package cDNAUpdate;

=pod

=head1 NAME

cDNAupdate.pm

=head1 DESCRIPTION

This file contains the necessary configuration for running cDNA_update
for both human and mouse cDNAs. You should read this POD along with
the POD for cDNA_update for an overview of the cDNA_update process.

As in all configurations there is a default section containg all
variables that are required, followed by the individual settings for
human and mouse. Which species to use is set when cDNA_update.pl is
run.

Note that the import subroutine is slightly different to the method
usually used for importing the configurations as the cDNA_update
procedure itself doesn't follow a regular genebuild pipeline.

Some variables are pre-defined in the default section, these should
not require minimal modicifations once they have been set. The user
section will change every time the procedure is run.

=cut
use strict;
use vars qw(%Config);

%Config = (
    DEFAULT => {
        # Admin rights are required
        WB_DBUSER     => 'ensadmin',
        WB_DBPASS     => 'ensembl',
        WB_REF_DBPORT => 3306,

        # Path to gss file
        GSS_PATH => '/ensembl-personal/genebuilders/cDNA_update/gss_acc.txt',

        # Various scripts required by the process.
        FASTA_SPLIT             => '/nfs/acari/searle/progs/production_code/ensembl-trunk_1106/ensc-core/src/Programs/fastasplit',
        POLYA_CLIPPING_PATH     => '/ensembl-pipeline/scripts/EST/new_polyA_clipping.pl',
        FIND_N_PATH             => '/ensembl-pipeline/scripts/cDNA_update/find_N.pl',
        STORE_UNMAPPED_PATH     => '/ensembl-pipeline/scripts/cDNA_update/store_unmapped_cdnas.pl',
        UNMAPPED_REASONS_PATH   => '/ensembl/misc-scripts/unmapped_reason/unmapped_reason.txt',
        LOAD_TAX_PATH           => '/ensembl-pipeline/scripts/load_taxonomy.pl',

        # Exonerate specifications
        PROGRAM_NAME    => "exonerate",
        PROGRAM_VERSION => "0.9.0",
        PROGRAM_FILE    => "/usr/local/ensembl/bin/exonerate-0.9.0",
        MODULE_NAME     => "Exonerate2Genes",

        # Source data files details
        SOURCE_HOST => 'cbi4',
        SOURCE_DIR  => '/data/blastdb/',

        # Taxonomy db for loading meta_table - should not need to change
        TAXONDBNAME => 'ncbi_taxonomy',
        TAXONDBHOST => 'ens-livemirror',
        TAXONDBPORT => 3306,

        # For the comparison only
        OLD_FEATURE_NAME => 'cDNA_update',

        # These variables need to be set, please do so in the
        # corresponding hash tables for the species you want to
        # analyse.

        # User details
        USER         => '',
        HOST         => '',
        GENEBUILD_ID => undef,

        # Reference db (current build)
        WB_REF_DBNAME => '',
        WB_REF_DBHOST => '',

        # New source db (PIPELINE)
        WB_PIPE_DBNAME => '',
        WB_PIPE_DBHOST => '',
        WB_PIPE_DBPORT => 3306,

        # New target db (ESTGENE)
        WB_TARGET_DBNAME => '',
        WB_TARGET_DBHOST => '',
        WB_TARGET_DBPORT => 3306,

        # Older cDNA db (needed for comparison only) -
        # check schema is up to date!!!!!!
        WB_LAST_DBNAME => '',
        WB_LAST_DBHOST => '',
        WB_LAST_DBPORT => undef,

        # Reference db (last build, needed for comparison only)
        WB_LAST_DNADBNAME => '',
        WB_LAST_DNADBHOST => '',
        WB_LAST_DNADBPORT => undef,

        # Path to your cvs directory
        CVS_DIR => '',

        # Where the output files should go
        DATA_DIR => '',

        # Path to the genomic sequence
        GENOMICSEQS => '',

        # Chunk size recommendations: 5500 for human
        # otherwise get AWOL jobs in first run
        CHUNK => undef,

        # Sequence files
        VERTRNA        => '',
        VERTRNA_UPDATE => '',
        REFSEQ         => '',

        # Species information
        COMMON_SPECIES_NAME => '',
        SPECIES             => '',
        TAX_ID              => undef,
    },

    human => {
        USER         => 'amonida',
        HOST         => 'bc-9-1-01',
        GENEBUILD_ID => 25,

        WB_REF_DBNAME => 'amonida_human_core_55',
        WB_REF_DBHOST => 'genebuild7',

        WB_PIPE_DBNAME => 'amonida_homo_cdna0509_ref',
        WB_PIPE_DBHOST => 'genebuild1',
        WB_PIPE_DBPORT => 3306,

        WB_TARGET_DBNAME => 'amonida_homo_cdna0509_update',
        WB_TARGET_DBHOST => 'genebuild1',
        WB_TARGET_DBPORT => 3306,

        WB_LAST_DBNAME => 'homo_sapiens_cdna_54_36p',
        WB_LAST_DBHOST => 'ens-livemirror',
        WB_LAST_DBPORT => 3306,

        WB_LAST_DNADBNAME => 'homo_sapiens_core_54_36p',
        WB_LAST_DNADBHOST => 'ens-livemirror',
        WB_LAST_DNADBPORT => 3306,

        CVS_DIR => "$ENV{CVSDIR}",

        DATA_DIR => "$ENV{WORK}",

        # You shouldn't need to change the settings below but do check #
        # that they are correct.                                       #
        ################################################################

        # Path to the genomic sequence
        GENOMICSEQS => '/data/blastdb/Ensembl/Human/GRCh37/genome/softmasked/softmasked_dusted.fa',

        # Chunk size recommendations: 5500 for human
        # otherwise get AWOL jobs in first run
        CHUNK => 5500,

        # Sequence files
        VERTRNA        => 'embl_vertrna-1',
        VERTRNA_UPDATE => 'emnew_vertrna-1',

        # Using human sequence
        REFSEQ => 'hs.fna',

        # Species information
        COMMON_SPECIES_NAME => 'human',
        SPECIES             => 'Homo sapiens',
        TAX_ID              => 9606,
    },

    mouse => {
        USER         => 'amonida',
        HOST         => 'bc-9-1-03',
        GENEBUILD_ID => 25,

        WB_REF_DBNAME => 'amonida_mouse_core_53',
        WB_REF_DBHOST => 'genebuild4',

        WB_PIPE_DBNAME => 'amonida_mus_test_ref',
        WB_PIPE_DBHOST => 'genebuild4',
        WB_PIPE_DBPORT => 3306,

        WB_TARGET_DBNAME => 'amonida_mus_test_update',
        WB_TARGET_DBHOST => 'genebuild4',
        WB_TARGET_DBPORT => 3306,

        WB_LAST_DBNAME => 'mus_musculus_cdna_53_37f',
        WB_LAST_DBHOST => 'ensdb-archive',
        WB_LAST_DBPORT => 5304,

        WB_LAST_DNADBNAME => 'mus_musculus_core_53_37f',
        WB_LAST_DNADBHOST => 'ens-livemirror',
        WB_LAST_DNADBPORT => 3306,

        CVS_DIR => '/nfs/acari/amonida/projects/cdna_update/mouse/',

        DATA_DIR => "$ENV{SCRATCH}",

        # You shouldn't need to change the settings below but do check #
        # that they are correct.                                       #
        ################################################################

        # Path to the genomic sequence
        GENOMICSEQS => '/data/blastdb/Ensembl/Mouse/NCBIM37/genome/softmasked_dusted/toplevel_sequence.fa',

        # Chunk size recommendations: 1500 for mouse
        # otherwise get AWOL jobs in first run
        CHUNK => 1500,

        # Sequence files
        VERTRNA        => 'embl_vertrna-1',
        VERTRNA_UPDATE => 'emnew_vertrna-1',
        REFSEQ         => 'mouse.fna',

        # Species information
        COMMON_SPECIES_NAME => 'mouse',
        SPECIES             => 'Mus musculus',
        TAX_ID              => 10090,
    }
);

sub import {
  my ($callpack) = caller(0); # Name of the calling package
  my $pack = shift; # Need to move package off @_

  # Get list of variables supplied, or else
  # all of General:
  my @vars = @_ ? @_ : keys( %Config );
  return unless @vars;

  # Predeclare global variables in calling package
  eval "package $callpack; use vars qw("
    . join(' ', map { '$'.$_ } @vars) . ")";
    die $@ if $@;

    foreach my $var (@vars) {
        if ( defined $Config{ $var } ) {
            no strict 'refs';
            foreach my $key (keys %{$Config{$var}}) {
            # Exporter does a similar job to the following
            # statement, but for function names, not
            # scalar variables:
                *{"${callpack}::$key"} = \$Config{$var}{$key};
            }
        } else {
            die "Error: Config: $var not known\n";
        }
    }
}

1;

