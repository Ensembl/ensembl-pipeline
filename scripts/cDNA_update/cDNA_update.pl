#!/usr/local/ensembl/bin/perl

#$Id: cDNA_update.pl,v 1.50 2010-06-21 12:50:40 at6 Exp $

# Original version cDNA_update.pl for human cDNAs
# Adapted for use with mouse cDNAs - Sarah Dyer 13/10/05

# Uses new polyAclipping, stop list for gene trap cdnas (gss) and
# cdnas which hit too many times (kill list).
#
# Will make three different logic names, eg mouse_cDNA_update,
# mouse_cDNA_update_2 and mouse_cDNA_update_3 mouse_cDNA_update_2 and
# mouse_cDNA_update_3 are the same.
#
# Mouse jobs can run in BatchQueue normal, human needs BatchQueue long
#
# TODO:
# 1/ signal capturing would be useful when a user interrupts a step.
# 2/ use Mole for getting the sequences instead of downloading the
# fasta files directly.
#
# The following errors are reported, but can be ignored:
#
# "OSTYPE: Undefined variable"
#
# "You should also add a rule that has SubmitcDNAChunk as its goal,
# or this rule will never have its conditions fulfilled."
#
#
# "-------------------- WARNING ----------------------
# MSG: Some of your analyses don t have entries in the input_id_type_analysis table
# FILE: Pipeline/Utils/PipelineSanityChecks.pm LINE: 99
# CALLED BY: Pipeline/Utils/PipelineSanityChecks.pm  LINE: 73
# ---------------------------------------------------"
#
# "MSG: Could not find fasta file for 'MT in '/data/blastdb/Ensembl/Human/NCBI35/softmasked_dusted'"

=head1 NAME

cDNA_setup

Procedure for getting the latest human or mouse cDNAs aligned to an
existing build in an automated fashion.

=head1 SYNOPSIS

cDNA_update.pl [ prepare | test-run | run | compare | clean ]

=head1 OPTIONS

The options should be run in the above order.

prepare     Read the configuration file cDNAUpdate.pm.
            Pipeline and target databases are created and configured. Sequence
            files are also downloaded and kill-list objects removed from
            the set. Analyses are added and input_ids made.

test-run    (optional but recommended) Run a test_Runnable.

run         Start the analysis.

compare     Compare the number of alignments per chromosome with the
            previous cDNA_update.

clean       Remove temporary files, output directories, sql dumps etc.

=head1 DESCRIPTION

Please read the configuration POD along with this POD.

This is a set-up script for generating new human or mouse cDNA
alignments as an isolated step of the build process with the Ensembl
pipeline. We are using the current build as a basis to add the latest
cDNA information as an additional track. This is done using exonerate
in the pipeline, fasta files and repeat-masked genome files. The
results are dna_align_features in an ESTGENE-type database.

Check out the latest code to match the database to be updated, for
example:

   cvs co -r branch-ensembl-53 ensembl
   cvs co -r branch-ensembl-53 ensembl-pipeline
   cvs co -r branch-ensembl-53 ensembl-analysis
   cvs co -r branch-ensembl-53 ensembl-compara
   cvs co -r branch-ensembl-53 ensembl-killlist

*MAKE SURE THAT YOU HAVE THE LATEST ensembl/sql/table.sql FILE*

=head2 The steps the script performs

  1. config_setup: check config variables & files.
  2. DB_setup: partly copy current db (PIPELINE database),
     insert analysis etc., create OUTPUT database, synchronise.
  3. fastafiles: get & read input files, kill-list objects are removed
     at this stage.
  4. chop off the polyA tails and chunk the fasta files into smaller
     pieces.
  5. (optional) test-run: run a test_Runnable.
  5. run_analysis: run exonerate using the pipeline.
  6. rerun cDNAs which did not align using different Exonerate
     parameters this step is usually repeated twice.
  7. find_many_hits: identify those cDNAs which aligned to many places
     in the genome.
  8. why_cdnas_missed: compile a list of reasons why hits were rejected.
  9. comparison: health-check by comparing the results to previous
     alignments.
     quicker version: get number of alignments for each chromosome
     from previous and new db version.
     extended version: look for new hits, track missing hits.
 10. cleanup: post-process result DB, restore config files, remove tmp
     files and dbs.

Note that run_analysis is called 3 times:
  . First it runs jobs with normal exonerate parameters and the
    specified chunk size
  . Second it runs exonerat with more exhaustive parameters
  . Last it rechunks into smaller chunk s and runs exonerate with
    the same parameters as second run

=head2 Check points

The check points (progress_status) are introduced at various steps for
easy tracking of which steps have been done. The progress_status can
be reset to a previous step if the process was interrupted and rerun
from that point onwards.

The check points have been set after the following steps:

  1. Databases have been setup
  2. Fasta files have been downloaded and kill-list objects
     removed.
  3. Test-run
  4. The first run
  5. Finding missing cDNAs
  6. Resetting the configuration and remaking the fasta files
  7. AWOL jobs
  8. Finding cDNAs with many hits to genome
  9. Storing cDNAs as unmapped objects
  10. Updating the meta coordinate table
  11. Updating the meta table
  12. Comparison with previous analysis

You would need to remove the progress_status entry from the meta table
before handind the database over.

=head2 What YOU will need to do

  1. Fill in the configuration file (cDNAUpdate.pm)
     for the species that you will be using along with
     /Bio/EnsEMBL/Analysis/Config/GeneBuild/KillListFilter.pm.
  2. Ask systems to push the genome files across the farm
     (after the prepare step) if necessary.
  3. Run it; check the set-up and re-run if there are errors.
  4. Check the results directly by querying the target db.
  5. Hand over target-database (patch to new version if necessary).

=head2 Potential error handling

If there is an error and the script dies, the original config files
are restored without removing the data files and databases, allowing
the re-run of the script. The tracking system in place as well which
keeps record of which steps have been run can be used for rerunning
any steps that did not finish properly.

=head2 Consumed time

The setup of scripts and databases runs for ~ 15 min, the exonerate
pipeline needs around 24 h for the first round of human cDNAs,
depending on farm sage.

Change specifications manually in BatchQueue.pm and re-run the
pipeline command if jobs fail or take too long:
  resource => 'select[mem>2500] rusage[mem=2500]',
  queue    => 'bugmem'

=head2 Run the healthchecks

Run the healthchecks as below:

  run-healthcheck.sh -d <user>_cDNA_update -output problem -species homo_sapiens -type cdna post_genebuild

Hand-over target db.

=head1 CONTACT

ensembl-dev@ebi.ac.uk

=cut


use strict;
use warnings;
use Carp;
use File::Copy;
use Data::Dumper;
use DBI;
use DBD::mysql;

#use Bio::EnsEMBL::KillList::KillList;
use Bio::EnsEMBL::KillList::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw ( get_input_arg );
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBConnection;
use cDNAUpdate; 

# We need the Net::SSH module from somewhere:
use lib '/software/perl-5.8.8/lib/site_perl/5.8.8/';
use Net::SSH qw(sshopen2);

my $species;

print "Which species are you running, human or mouse? ";
$species = <STDIN>;
chomp $species;

if ( ($species ne 'human') and ($species ne 'mouse') ){
    print("Please specify either human or mouse!\n");
    exit;
} 

eval {
       import cDNAUpdate ('DEFAULT', $species )
     };



$| = 1;

# Variables from configuration

my $POLYA_CLIPPING      = $CVS_DIR . $POLYA_CLIPPING_PATH;
my $FIND_N              = $CVS_DIR . $FIND_N_PATH;
my $STORE_UNMAPPED      = $CVS_DIR . $STORE_UNMAPPED_PATH;
my $UNMAPPED_REASONS    = $CVS_DIR . $UNMAPPED_REASONS_PATH;
my $LOAD_TAX            = $CVS_DIR . $LOAD_TAX_PATH;
my $GSS;

if ( defined($GSS_PREFIX) ) {
    # In case gss file is not under the same directory as the code
    # checkouts.
    $GSS = $GSS_PREFIX . $GSS_PATH;
} else {
    $GSS = $CVS_DIR . $GSS_PATH;
}

# When comparing to a previously updated cdna db
# $oldFeatureName = $newFeatureName
my $newFeatureName  = "cDNA_update"; # the analysis name!

my %saved_files;
my $cmd;
my $status;

# Temp. dirs & files:
my $config_file     = $CVS_DIR."/ensembl-pipeline/scripts/cDNA_update/config_files.txt";
my $newfile         = "cdna_update";
my $configDIR       = $DATA_DIR."/configbackup";
my $chunkDIR        = $DATA_DIR."/chunks";
my $outDIR          = $DATA_DIR."/output";
my @masked_genome   = $GENOMICSEQS;
my $submitName      = "SubmitcDNAChunk";

my $genomelist      = join "\',\'", @masked_genome;

my %configvars = (
                 "MIN_LENGTH"       => $MIN_LENGTH ,         # from cDNAUpdate
                 "CVS_DIR"          => $CVS_DIR,             # from cDNAUpdate
                 "DATA_DIR"         => $DATA_DIR,            # from cDNAUpdate
                 "VERTRNA"          => $VERTRNA,             # from cDNAUpdate
                 "VERTRNA_UPDATE"   => $VERTRNA_UPDATE,      # from cDNAUpdate
                 "REFSEQ"           => $REFSEQ,              # from cDNAUpdate
                 "FASTASPLIT"       => $FASTA_SPLIT,         # from cDNAUpdate
                 "POLYA_CLIPPING"   => $POLYA_CLIPPING,      # from cDNAUpdate
                 "DBUSER"        => $DBUSER,           # from cDNAUpdate
                 "DBPASS"        => $DBPASS,           # from cDNAUpdate
                 "REF_DBNAME"    => $REF_DBNAME,       # from cDNAUpdate
                 "REF_DBHOST"    => $REF_DBHOST,       # from cDNAUpdate
                 "REF_DBPORT"    => $REF_DBPORT,       # from cDNAUpdate
                 "PIPE_DBNAME"   => $PIPE_DBNAME,      # from cDNAUpdate
                 "PIPE_DBHOST"   => $PIPE_DBHOST,      # from cDNAUpdate
                 "PIPE_DBPORT"   => $PIPE_DBPORT,      # from cDNAUpdate
                 "OUTPUT_DBNAME" => $OUTPUT_DBNAME,    # from cDNAUpdate
                 "OUTPUT_DBHOST" => $OUTPUT_DBHOST,    # from cDNAUpdate
                 "OUTPUT_DBPORT" => $OUTPUT_DBPORT,    # from cDNAUpdate
                 "taxonomy_id"      => $TAX_ID,              # from cDNAUpdate
                 "PROGRAM_NAME"     => $PROGRAM_NAME,        # from cDNAUpdate
                 "PROGRAM_VERSION"  => $PROGRAM_VERSION,     # from cDNAUpdate
                 "PROGRAM_FILE"     => $PROGRAM_FILE,        # from cDNAUpdate
                 "MODULE_NAME"      => $MODULE_NAME,         # from cDNAUpdate
                 "SOURCE_DIR"       => $SOURCE_DIR,          # from cDNAUpdate
                 "chunkDIR"         => $chunkDIR,
                 "outDIR"           => $outDIR,
                 "configDIR"        => $configDIR,
                 "newfile"          => $newfile,
                 "config_file"      => $config_file,
                 "masked_genome"    => $genomelist,
                 "newFeatureName"   => $newFeatureName );

# Fasta chunk specifications:
my $maxseqlength      = 17000;
my $num_missing_cdnas = 0;
my $rerun_flag        = 0;

my $option = $ARGV[0]; 

if ( !$option or
            ( $option ne "prepare"
          and $option ne "prepare"
          and $option ne "test-run"
          and $option ne "run"
          and $option ne "clean"
          and $option ne "compare" ) ) { 
   print "\n\nYou have to supply an option to the cDNA_update.pl command - valid options are : \n\n\t";
   print  "perl cDNA_update.pl [ prepare | test-run | run | compare | clean ]\n\n" ; 
   sleep(1) ; 
   exec('perldoc', $0);
}


my $pipe_db = connect_db( $PIPE_DBHOST, $PIPE_DBPORT,
                     $PIPE_DBNAME, $DBUSER,
                     $DBPASS );

my $progress_status = undef;

if ( $option eq "prepare" ) {

    if( exists($ENV{STY}) ) {
        print "Good, you're running this in screen\n";
    } else {
        print(   "The program will exit now. "
               . "Restart it again in a screen session.\n" );
        unclean_exit();
    }

    print("\nStarting cDNA-update procedure.\n");

    if ( !defined($progress_status) ) {

        my $ref_db = connect_db( $REF_DBHOST, $REF_DBPORT,
                                 $REF_DBNAME, "ensro" );

        # Check existence of source databases
        if ( !$ref_db ) {
            croak("could not find $REF_DBNAME.");
        }

        # Disconnect from the database
        $ref_db->dbc->disconnect_if_idle();

        print("Shall we set the configuration? (y/n) ");
        if ( get_input_arg() ) {
            config_setup();
        }

        print "\nSet-up the databases?(y/n) ";

        if ( get_input_arg() ) {
            if ( !DB_setup() ) {
                unclean_exit();
            }

            my $progress = get_status($pipe_db->dbc());
            if ( !defined($progress) ) {
                print "No progress_status, setting it to 1!\n";
                $progress_status = 1;

                # Set progress_status
                my $sql = "INSERT INTO meta (meta_key, meta_value) "
                        . "VALUES ('progress_status', $progress_status)";
                my $sth = $pipe_db->dbc->prepare($sql);
                $sth->execute;
            } else {
                croak(   "Something is wrong, progress_status shouldn't be "
                       . "set yet, rerun the prepare step again. "
                       . "progress_status = $progress\n" );
            }
        } else {

            my $pipe_db = connect_db( $PIPE_DBHOST, $PIPE_DBPORT,
                                      $PIPE_DBNAME, $DBUSER,
                                      $DBPASS );

            my $target_db = connect_db( $OUTPUT_DBHOST, $OUTPUT_DBPORT,
                                        $OUTPUT_DBNAME, "ensro" );

            if ( $pipe_db->dbc->dbname() && $target_db->dbc->dbname() ) {
                print("\nDatabases exist, good! so we continue with the process.\n");

                $target_db->dbc->disconnect_if_idle();

                my $progress = get_status($pipe_db->dbc());
                if ( defined($progress) ) {
                    $progress_status = $progress;
                    print("progress_status is already set ($progress_status), "
                        . "jumping to the next step.\n");
                } else {
                    croak("The progress_status should have already been set. "
                        . "Perhaps you didn't set the databases correctly?\n"
                        . "\nPlease rerun the prepare step now.\n");
                }
            } else {
                carp(  "The databases have not been set up, "
                     . "please do this now.\n" );
                unclean_exit();
            }
        } ## end else [ if ( get_input_arg() )
    } ## end if ( !defined($progress_status...

    $progress_status = get_status($pipe_db->dbc());
    print("\n\nAfter DB setup progress_status is: $progress_status\n\n");

    if ( $progress_status == 1 ) {
        print(   "\nGet fasta files? "
               . "kill-list objects will be removed as well. (y/n) " );
        if ( get_input_arg() ) {
            if ( !fastafiles() ) {
                unclean_exit();
            }
        } else {
            print "\nYou said don't fetch files, checking if files exist...\n\n";
            if (    -e $DATA_DIR . "/" . $REFSEQ
                 && -e $DATA_DIR . "/" . $VERTRNA ) {
                print(   "Files are present in the directory, "
                       . "removing sequences on the kill-list.\n\n" );

                my $trim_file = remove_kill_list_object();
                polya_clipping($trim_file);
            } else {
                carp("\nThe files don't exist, exiting the program.\n\n");
                unclean_exit();
            }
        }

        # Create analysis, set the rule and make the input_ids
        setup_analysis();

        # Synchronise reference and target databases
        print( "\nHave now added analyses and rules to the pipeline database "
               . "so should synchronise the target database now...\n" );

        sync_databases();

        $progress_status = 2;
        update_progress_status($progress_status);

        print(   "\n\nFinished setting up the analysis.\n\n"
               . "NEXT: If you want to run a test use the 'test-run' "
               . "option. Alternatively, start the analysis itself "
               . "using the 'run' option.\n\n" );

    } ## end if ( $progress_status ...
} ## end if ( $option eq "prepare")

elsif ( $option eq "test-run" ) {

    $progress_status = get_status($pipe_db->dbc());
    print("\tprogress_status is: $progress_status\n\n");

    if ( $progress_status == 2 ) {

        print("\nDo we need to set/re-set the configs? (y/n) ");
        if ( get_input_arg() ) {
            config_setup();
        }

        test_run();

        $progress_status = 3;
        update_progress_status($progress_status);
    }
    print("NEXT: start the analysis using the 'run' option.\n\n");
} ## end elsif ( $option eq "test-run")

elsif ( $option eq "run" ) {

    print("\nDo we need to set/re-set the configs? (y/n) ");
    if ( get_input_arg() ) {

        # Is this needed? The options are set by the script, so user doesn't
        # have to change the configuration.
        print(   "Are you sure you want to reset the configs?\n"
               . "\tIf you're in the second stage you need to\n"
               . "\tedit BatchQueue.pm as well as Exonerate2Genes.pm!\n"
               . "Shall we set the configs? (y/n) " );

        if ( get_input_arg() ) {
            config_setup();
        } else {
            print "Not resetting config.\n\n" ;
        }
    }

    $progress_status = get_status($pipe_db->dbc());
    if ( $progress_status >= 2 && $progress_status < 6 ) {

        if ( $progress_status == 2 || $progress_status == 3 ) {
            if ( !run_analysis() ) {

                # Set status to 2 (fasta files & kill list removal are done).
                $progress_status = 2;
                update_progress_status($progress_status);
                unclean_exit();
            }

            # Have run analysis and it was successful so update the
            # status.
            $progress_status = 4;
            update_progress_status($progress_status);
        }

        $progress_status = get_status($pipe_db->dbc());

        if ( $progress_status == 4 ) {

            $num_missing_cdnas = find_missing_cdnas();

            $progress_status = 5;
            update_progress_status($progress_status);

            print(   "\n$num_missing_cdnas have not aligned to the genome.\n"
                . "Would you like to rerun these cdnas with "
                . "adjusted Exonerate parameters:\n"
                . "\tmax_intron     = 400,000 and\n"
                . "\tsoftmasktarget = FALSE? (y/n) " );
        }

        # Rerunning the analysis
        if ( $progress_status >= 5 && get_input_arg() ) {
            $rerun_flag = 1;

            # Change the logic_name and directories

            # To show different params
            $newFeatureName               = $newFeatureName . "_2";

            $configvars{"newFeatureName"} = $newFeatureName;
            $chunkDIR                     = $DATA_DIR . "/chunks2";
            $configvars{"chunkDIR"}       = $chunkDIR;
            $outDIR                       = $DATA_DIR . "/output2";
            $configvars{"outDIR"}         = $outDIR;

            $progress_status = get_status($pipe_db->dbc());

            # Remaking fasta files
            if ( $progress_status == 5 ) {

                config_setup();
                remake_fasta_files();

                $progress_status = 6;
                update_progress_status($progress_status);

                print("\nSet databases for next run? (y/n) ");
                if ( get_input_arg() ) {
                    $rerun_flag = 1;
                    if ( !DB_setup() ) {
                        unclean_exit();
                    }
                    print(   "\n\nFinished setting up the analysis "
                           . "for next round.\n\n" );
                }
                print("Should we start the 2nd run ? (y/n) ");
                if ( get_input_arg() ) {
                    if ( !run_analysis() ) {

                        # If run_analysis is not successful, return to
                        # the step before the config setup, need to
                        # recalulate the missing cdnas.
                        $progress_status = 4;
                        update_progress_status($progress_status);
                        carp("Analysis didn't go well, "
                            . "resetting progress_status to step before.\n"
                            . "progress_status is now set to: $progress_status. "
                            . "You'll need to rerun the analysis.\n\n");
                        unclean_exit();
                    }
                }
            # maybe need to increment the status here?
            }
        }
    }

    elsif ( $progress_status >= 6 ) {

        if ( $progress_status == 6 ) {
            print("Do you want to check for AWOL jobs ? (y/n) ");
            if ( get_input_arg() ) {
                print "checking for AWOL jobs...\n";
                chase_jobs();

                $progress_status = 7;
                update_progress_status($progress_status);
            }
        }

        # Checking cDNAs with many hits to the genome
        if ( $progress_status == 7 ) {

            print(   "Would you like to check for cDNAs which "
                   . "have hit many places in the genome? (y/n) " );

            if ( get_input_arg() ) {

                find_many_hits();

                $progress_status = 8;
                update_progress_status($progress_status);
            }
        }

        # Storing cDNAs as unmapped_objects
        if ( $progress_status == 8 ) {

            print(   "Would you like to store the cDNAs which "
                   . "have not aligned as unmapped_objects? (y/n) " );

            if ( get_input_arg() ) {

                why_cdnas_missed();

                $progress_status = 9;
                update_progress_status($progress_status);
            }
        }

        # Updating the meta_coord table
        if ( $progress_status == 9 ) {
            print("updating meta_coord table...\n");

            update_metacoord();

            $progress_status = 10;
            update_progress_status($progress_status);
        }

        # Sorting out the meta table
        if ( $progress_status == 10 ) {

            print("sorting out the meta table...\n");

            fix_metatable();

            $progress_status = 11;
            update_progress_status($progress_status);
        }

        print(   "\n\nNOTE!!! You should now change the analysis_ids "
                . "of the cdnas in your database\n"
                . "so that they all have the same logic name, "
                . "otherwise the comparison script won't work;\n"
                . "you will need to change both the gene and "
                . "dna_align_feature tables.\n\n" );
    } ## end elsif ( $progress >= 6 )
}

elsif ( $option eq "compare" ) {
    print("\nHave you changed the analysis_id for the runs? (y/n) ");

    if ( get_input_arg() ) {
        $progress_status = get_status($pipe_db->dbc());

        if ( $progress_status == 11 ) {
            print(   "\nRunning checks after cDNA-update.\n"
                . "checking through alignments & genes.\n" );

            check_vars();
            compare();

            $progress_status = 12;
            update_progress_status($progress_status);
        }
    } else {
        croak("\nChange the analysis_ids first before running "
            ." the comparison step.\n\n");
    }
}

elsif ( $option eq "clean" ) {

    $progress_status = get_status($pipe_db->dbc());
    if ( $progress_status == 12 ) {
        print("\nCleaning up after cDNA-update.\n");
        clean_up(0);
    }
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Write required config files for analysis.
# Stores content and location of original files for later restoration.
# The required config values are written into placeholders in the config skeleton file (config_files.txt)
#  using the variables defined above.

sub config_setup {
    $status = 0;
    my $filecount = 0;
    my ( $header, $filename, $path );

    # Set env var to avoid warnings
    if ( !defined $ENV{"OSTYPE"} ) {
        $ENV{"OSTYPE"} = "";
    }
    # import function, to be included in all config files
    my $import_sub = '
        sub import {
            my ($callpack) = caller(0);
            my $pack = shift;
            my @vars = @_ ? @_ : keys(%Config);
            return unless @vars;
            eval "package $callpack; use vars qw(".
            join(\' \', map { \'$\'.$_ } @vars) . ")";
            die $@ if $@;
            foreach (@vars) {
            if (defined $Config{ $_ }) {
                no strict \'refs\';
                *{"${callpack}::$_"} = \$Config{ $_ };
            }else {
                die "Error: Config: $_ not known\n";
            }
            }
        }
        1;
        ';
    check_vars();

    # Check existence of source databases
    if ( !connect_db( $REF_DBHOST, $REF_DBPORT,
                      $REF_DBNAME, "ensro" ) ) {
        croak("Could not find $REF_DBNAME.");
    }

    # Go through config info to create defined files,
    # back-up the original, write version with filled-in variables
    open( RP, "<", "$config_file" )
      or croak("Can't open config file definitions $config_file.");

    local $/ = '>>';
    <RP>;

    while ( my $content = <RP> ) {

        # Get specific config-file name
        $content =~ s/([\w\/_\-\.]+)\n//;
        $header = $1;

        $header =~ m/(.+\/)([\w\._\-]*)$/g;
        $path     = $1;
        $filename = $2;
        $content =~ s/>>//;

        # Replace variables in config file
        foreach my $configvarref ( keys %configvars ) {
            $content =~ s/\<$configvarref\>/$configvars{$configvarref}/g;
        }

        # If rerunning seqs change Exonerate options:
        if (   ( $filename =~ /Exonerate2Genes/ )
            && ( $rerun_flag == 1 ) ) {
            # Because only interested in those which don't align after new parameters

            # $content =~ s/-verbosity => 0/-verbosity => 1/;  # old regex, didn't work

            # Due to config hash formatting, there are usually many trailing
            # whitespaces between the hash key and values. Must include the
            # spaces in the regex to carry out the substitution as intended:
            $content =~ s/\-verbosity\s+=>\s0/-verbosity => 1/;

            my $substitute = '--maxintron 400000 --bestn 10 --softmasktarget FALSE';
            $content =~ s/--softmasktarget TRUE/$substitute/;
        }

        # Backup file if exists
        if ( -e $CVS_DIR . "/" . $header ) {
            print $CVS_DIR. "/" . $header . " exists\n";
            eval {
                move( $CVS_DIR      . "/" . $header,
                      $configDIR    . "/" . $filecount );
            };
            if ($@) {
                croak( "Failed to move "
                     . $CVS_DIR     . "/" . $header
                     . " to "
                     . $configDIR   . "/" . $filecount
                     . ". $@" );
            }
        }

        # Store file location
        $saved_files{$filecount} = $header;

        # Write modified file
        $filename = $CVS_DIR . "/" . $path . $filename;

        open( WP, ">", $filename )
          or croak("Can't create new config file $filename.\n");
        print WP $content . "\n" . $import_sub;
        close WP;
        $filecount++;
    } ## end while ( my $content = <RP>)

    close(RP);

    # Save data dump with config paths
    $Data::Dumper::Purity = 1;
    open( WP, "> config_paths.perldata" )
      or croak("\nCan't create file for data dumping!\n");
    print WP Data::Dumper->Dump( [ \%saved_files ], ['*saved_files'] );
    close(WP);
    $/ = "\n";
    print(   "\nCreated backup of current config files, "
           . "new config files written.\n" );
} ## end sub config_setup


# Check files & directories, create if necessary
sub check_vars {
    foreach my $configvarref ( keys %configvars ) {
        if ( !$configvarref ) {
            croak(   "Please define all configuration variables! "
                   . "[$configvarref]\n" );
        }

        if ( $configvarref =~ m/.+DIR.*/ ) {
            if ( !-e $configvars{$configvarref} ) {
                if ( system("mkdir $configvars{$configvarref}") ) {
                    croak(   "Could not create directory! "
                           . "[$configvars{$configvarref}]. $!\n" );
                }
            }
            if ( !-r $configvars{$configvarref} ) {
                croak(   "Directory not accessible! "
                       . "[$configvars{$configvarref}]. $!\n" );
            }
        }
    }
}

# Delete old content if any for a given directory
# (recursively, as some tmp-dir get pretty crowded...)
sub checkdir {
    my $dirname = shift;
    my $option  = shift;

    # Go through dirs recursively
    unless ( opendir( DIR, $dirname ) ) {
        closedir(DIR);
        print "\nCan't open $dirname.\n";
        return 0;
    }

    my @files = grep ( !/^\.\.?$/, readdir(DIR) );
    foreach my $file (@files) {
        if ( -f "$dirname/$file" ) {
            system("rm $dirname/$file");
        }
        if ( -d "$dirname/$file" ) {
            checkdir("$dirname/$file");
        }
    }
    closedir(DIR);
    if ( !$option ) {
        eval { system("rmdir $dirname"); };
    }
    return 1;
} ## end sub checkdir


# Fetch fasta files, combine them, chop them up
sub fastafiles {
    $status = 0;

    my $vertrna_ver     = 1;
    my $vertrna_upd_ver = 1;
    my $refseq_ver      = 1;

    # Set to 1 if want to redo clipping and chunks regardless of
    # changes in vertna etc
    my $update = 0;

    my @filestamp;

    eval {
        # Check file versions, copy only if changed
        $cmd = "cd $SOURCE_DIR; ls -n "
                . $VERTRNA . " "
                . $VERTRNA_UPDATE . " "
                . $REFSEQ;

        sshopen2( "$USER\@$HOST", *READER, *WRITER, "$cmd" ) || croak "ssh: $!";
        while (<READER>) {
            @filestamp = split( " ", $_ );
            my $stampA = join( "-", @filestamp[ 5 .. 7 ] );
            $cmd = "cd " . $DATA_DIR . "; " . "ls -n " . $filestamp[8];
            @filestamp = split( " ", `$cmd` );
            my $stampB = join( "-", @filestamp[ 5 .. 7 ] );
            if ( $stampA eq $stampB ) {
                # No changes...
                if ( $filestamp[8] eq $VERTRNA ) {
                    $vertrna_ver = 0;
                } elsif ( $filestamp[8] eq $VERTRNA_UPDATE ) {
                    $vertrna_upd_ver = 0;
                } elsif ( $filestamp[8] eq $REFSEQ ) {
                    $refseq_ver = 0;
                }
            }
        } ## end while (<READER>)
        close(READER);
        close(WRITER);

        # Copy files
        if ($vertrna_ver) {
            $cmd = "scp -p " . $SOURCE_HOST . ":"
                . $SOURCE_DIR . "/" . $VERTRNA . " "
                . $DATA_DIR   . "/" . $VERTRNA;

            print $cmd, "\n";

            $status += system($cmd);
        }
        if ($vertrna_upd_ver) {
            $cmd = "scp -p " . $SOURCE_HOST . ":"
                . $SOURCE_DIR . "/" . $VERTRNA_UPDATE . " "
                . $DATA_DIR   . "/" . $VERTRNA_UPDATE;

            print $cmd, "\n";
            $status += system($cmd);
        }
        if ($refseq_ver) {
            $cmd = "scp -p " . $SOURCE_HOST . ":"
                . $SOURCE_DIR . "/" . $REFSEQ . " "
                . $DATA_DIR   . "/" . $REFSEQ;

            print $cmd, "\n";
            $status += system($cmd);
        }
        if ($status) { croak("Error while copying files.\n"); }
        print "Copied necessary files.\n";
        my $file = $DATA_DIR . $newfile;
        if ( $vertrna_upd_ver or $vertrna_ver or $refseq_ver ) {
            $update = 1;
            write_to_file();
        }

        my $newfile2 = remove_kill_list_object();
        if ($update) {
            polya_clipping($newfile2);
        }
    };
    if ($@) {
        print STDERR "\nERROR: $@";
        return 0;
    }
    return 1;
} ## end sub fastafiles

sub write_to_file {
    my $update;

#    $status = 0;

    my %EMBL_ids;
    my $header;
    my @filestamp;

    # Get entries for species of interest, combine base file and
    # update file.
    # Read update file.
    local $/ = "\n>";
    open( RP, "<", $DATA_DIR . "/" . $VERTRNA_UPDATE )
      or croak("can't read $VERTRNA_UPDATE\n");
    open( WP, ">", $DATA_DIR . "/" . $newfile )
      or croak("can't create $newfile\n");
    while ( my $entry = <RP> ) {
        # Need this to include the first record when using $/='\n>'
        $entry =~ s/^>//;
        if ( $entry =~ m/$SPECIES/ ) {
            # Extract & save id
            $entry =~ s/^([\w\.\d]+)\s.*\n{1}?/$1\n/;
            if ( !$1 ) {
                croak(   "\n$VERTRNA_UPDATE: "
                       . "unmatched id pattern:\n$entry\n" );
            }
            $EMBL_ids{$1} = 1;
            # Re-write fasta entry
            $entry =~ s/\>//g;
            print WP '>' . $entry;
        }
    }
    close(RP);
    print("\nRead update $VERTRNA_UPDATE EMBL file.\n");

    my @embl_vertrna_files;
    opendir( DIR, "$DATA_DIR" || croak "Cannot opendir $DATA_DIR $!" );
    while ( my $filename = readdir(DIR) ) {
        if ( $filename =~ /$VERTRNA/ ) {
            push @embl_vertrna_files, $filename;
        }
    }

    # now loop through all embl_vertna files
    #read base file
    foreach my $embl_vert (@embl_vertrna_files) {
        open( RP, "<", $DATA_DIR . "/" . $embl_vert )
          or croak("Can't read $embl_vert\n");
        #<RP>;
        while ( my $entry = <RP> ) {
            # Need this to include the first record when using $/='\n>'
            $entry =~ s/^>//;
            if ( $entry =~ m/$SPECIES/ ) {
                # Extract & save id
                $entry =~ s/^([\w\.\d]+)\s.*\n{1}?/$1\n/;
                if ( !$1 ) {
                    croak("\n$embl_vert: unmatched id pattern:\n$entry\n");
                }
                if ( !defined( $EMBL_ids{$1} ) ) {
                    # Add fasta entry for unchanged id
                    $entry =~ s/\>//g;
                    print WP '>' . $entry;
                }
            }
        }
        close(RP);
        print "Read base $embl_vert EMBL file.\n";
    }

    # Read RefSeq file
    open( RP, "<", $DATA_DIR . "/" . $REFSEQ )
      or croak("Can not read $REFSEQ.\n");
    while ( my $entry = <RP> ) {

        # Need this to include the first record when using $/='\n>'
        # we're not using 'predicted' XM entries for now
        $entry =~ s/^>//;
        if ( $entry =~ m/^gi.+ref\|(NM_.+)\| $SPECIES.*/ ) {
            $header = $1;
        } elsif ( $entry =~ m/^gi.+ref\|(NR_.+)\| $SPECIES.*/ ) {
            $header = $1;
        } else {
            next;
        }
        $entry =~ s/\>//g;
        if ($header) {
            # Reduce header to accession number
            $entry =~ s/^gi.+\n{1}?/$header\n/g;
            print WP '>' . $entry;
        }
    }
    print "read RefSeq file.\n";
    close(RP);
    close(WP);
    local $/ = "\n";

} ## end sub write_to_file

# Now get the kill_list
# Config found at /Bio/EnsEMBL/Pipeline/Config/GeneBuild/KillListFilter.pm
sub remove_kill_list_object {
    require Bio::EnsEMBL::KillList::KillList;
    my $kill_list_object =
      Bio::EnsEMBL::KillList::KillList->new( -TYPE => 'cDNA_update' );
    my %kill_list = %{ $kill_list_object->get_kill_list() };

    open( LIST, "<", $GSS ) or croak("can't open gss list $GSS");
    my %gss;
    while (<LIST>) {
        my @tmp = split /\s+/, $_;
        $gss{ $tmp[1] } = 1;
    }
    close LIST;

    # Go through file removing any seqs which appear on the kill list
    local $/ = "\n>";
    if ( !( -e $DATA_DIR . "/" . $newfile ) ) {
        print("\tNewfile not here so need to create it first.\n");
        write_to_file();
    }
    my $newfile2 = $newfile . ".seqs";
    open( SEQS, "<", $DATA_DIR . "/" . $newfile )
      or croak("Can't open seq file $newfile");
    open( OUT, ">", $DATA_DIR . "/" . $newfile2 )
      or croak("Can't open seq file $newfile2");
    while (<SEQS>) {
        s/>//g;

        my @tmp = split /\n/, $_;

        # Store the accession number
        my $acc;
        if ( $tmp[0] =~ /(\w+)\./ ) {
            $acc = $1;
        }
        if ( ( !exists $kill_list{$acc} ) && ( !exists $gss{$acc} ) ) {
            print OUT ">$_";
        }
    }

    local $/ = "\n";
    close OUT;
    close SEQS;

    return $newfile2;
} ## end sub remove_kill_list_object


# Have already made the sequence file from the previously clipped
# seqs, just need to rechunk it.
sub remake_fasta_files {

    my $file = $DATA_DIR . "/missing_cdnas.fasta";

    # How many files do we want? automatically adjust chunk_num,
    # don't want >20 seqs/file because softmasktarget = false
    my $chunk_num = int( $num_missing_cdnas/20 );

    # Split fasta files, store into new CHUNKDIR
    print("Splitting new fasta file in to chuncks.\n");

    $cmd = "$FASTA_SPLIT $file $chunk_num $chunkDIR";
    if ( system($cmd) ) {
       croak("Couldn't split file.$@\n");
    }

    # Isolate biggest sequences
    check_chunksizes();

    print "\nChopped up file.\n";
}


# Find the really big sequences & put them into separate chunks
sub check_chunksizes {
    local $/ = '>';
    my $allseqs;
    my $file;
    my $toolongs;
    my $seqname;
    my $newfile;

    unless ( opendir( DIR, $chunkDIR ) ) {
        croak("Can't read $chunkDIR");
    }

    foreach ( readdir(DIR) ) {
        if ( ( $_ =~ /^\.+$/ ) || ( $_ =~ /^newchunk.+$/ ) ) { next; }
        $file     = $chunkDIR . "/" . $_;
        $toolongs = 0;
        $allseqs  = "";

        open( CHUNKFILE, "<$file" ) or croak("can t open file $file.");

        # Skipping the first one as it just contains ">"
        <CHUNKFILE>;

        while ( my $seq = <CHUNKFILE> ) {
            $seq =~ s/\>//;
            $seq =~ m/(.+)\n/;
            $seqname = $1;
            if ( length($seq) > $maxseqlength ) {
                print "\nTOO LONG: $seqname";
                if ( !$toolongs ) {
                    $toolongs = 1;
                }
                $newfile = $chunkDIR . "/newchunk_" . $seqname;
                open( NEWFILE, ">$newfile" )
                  or croak("Can't create new fasta file $newfile!");
                print NEWFILE ">" . $seq;
                close(NEWFILE);

            } else {
                $allseqs .= ">" . $seq;
            }
        }

        close(CHUNKFILE);

        if ($toolongs) {
            open( CHUNKFILE, ">$file" ) or croak("Can't open file $file.");
            print CHUNKFILE $allseqs;
            close(CHUNKFILE);
        }
    } ## end foreach ( readdir(DIR) )

    closedir(DIR);
    local $/ = "\n";
} ## end sub check_chunksizes


# Prepare required databases; pipeline db, target db and fill required
# tables with data.
sub DB_setup {
    $status = 0;

    eval {
        if ( $rerun_flag == 0 ) {

            # Create the pipeline database
            # 1 = pipeline db
            $status = create_db( $PIPE_DBNAME, 1 );

            my $tables_to_delete =  'analysis '           . 'assembly '
                                  . 'assembly_exception ' . 'attrib_type '
                                  . 'coord_system '       . 'meta '
                                  . 'meta_coord '         . 'seq_region_attrib ';

            # Delete unnecessary tables
            # 1 = pipeline db
            delete_unwanted_tables( $tables_to_delete, 1 );

            # Copy defined db tables from current build
            my $tables_to_dump =  'assembly '       . 'assembly_exception '
                                . 'attrib_type '    . 'coord_system '
                                . 'meta '           . 'seq_region '
                                . 'seq_region_attrib ';

            # Dump the tables to import_tables1.sql
            # 1 = file count
            $status += dump_tables( $tables_to_dump, 1);

            # Import tables from import_tables1
            # 1 = file count
            # 1 = pipeline db
            $status += import_tables(1, 1);

            # Copy dna table from current build and import it to pipeline
            my $dna_table = 'dna';

            # Dump the table to import_tables2
            # 2 = file count
            $status += dump_tables( $dna_table, 2 );

            # 2 = file count
            # 1 = pipeline
            $status += import_tables(2, 1);

            # Create the target database
            $status += create_db($OUTPUT_DBNAME);

            # Remove unwanted tables from above plus the seq_region table from
            # target db.
            $tables_to_delete .= 'seq_region ';

#            print( "\nTables to delete in targe db:\n" . $tables_to_delete,
#                   "\n\n" );
            delete_unwanted_tables( $tables_to_delete );

            # Load the external_db and unmapped_reason tables from
            # /ensembl/misc-scripts.
            load_misc_script_files();

            if ($status) { croak("Could not create databases!\n"); }

            print("\nCreated databases.\n\n");

        } else {
            # If rerunning without rebuilding databases - clear out jobs tables first:

            my $pipe_db = connect_db( $PIPE_DBHOST, $PIPE_DBPORT,
                                      $PIPE_DBNAME, $DBUSER,
                                      $DBPASS );

            my $tables_to_remove =  'job '              .   'job_status '
                                  . 'rule_goal '        .   'rule_conditions '
                                  . 'input_id_analysis '.   'input_id_type_analysis ';

            my @tables = split(/ /, $tables_to_remove);
            foreach my $table ( @tables ) {
                print "\tDeleting table $table\n";

                my $sql = "DELETE FROM $table";
                my $sth = $pipe_db->dbc->do($sql);
            }

            $pipe_db->dbc->do( "DELETE from analysis where logic_name = '$submitName'");
            $pipe_db->dbc->do("DELETE from analysis where logic_name = '$newFeatureName'");

            # Create analysis, set the rule and make the input_ids
            setup_analysis();

        } ## end else [ if ( $rerun_flag == 0 )

        # Synchronise pipeline and target databases
        print("\nSynchronising the pipeline and target databases...\n");
        sync_databases();

    };    ## end eval

    if ($@) {
        print STDERR "\nERROR: $@";
        return 0;
    }
    return 1;
} ## end sub DB_setup


# Running a test first
sub test_run {
    print("\nRunning the test-RunnablDB.\n");

    # Get one input id for testing
    my $db = connect_db( $PIPE_DBHOST, $PIPE_DBPORT,
                         $PIPE_DBNAME, $DBUSER,
                         $DBPASS );

    my $sql = "SELECT input_id "
            . "FROM input_id_analysis i, analysis a "
            . "WHERE i.analysis_id = a.analysis_id "
            . "AND a.logic_name = '" . $submitName
            . "' LIMIT 1;";

    my $sth = $db->dbc->prepare($sql)
        or croak("Sql error getting an input-id!\n$!");

    $sth->execute();
    my ($input_id) = $sth->fetchrow_array;

    print "\n$input_id\n";

    if ( !$input_id ) {
        croak(  "\nCould not get an input id from database!\n"
              . "Query used: $sql\n\n" );
    }

    $cmd = "perl " . $CVS_DIR . "/ensembl-analysis/scripts/test_RunnableDB"
         . " -dbhost "     . $PIPE_DBHOST
         . " -dbport "     . $PIPE_DBPORT
         . " -dbuser "     . $DBUSER
         . " -dbpass "     . $DBPASS
         . " -dbname "     . $PIPE_DBNAME
         . " -input_id "   . $input_id
         . " -logic_name " . $newFeatureName
         . " -verbose"
         . " -nowrite";

    print $cmd . "\n";
    system($cmd);

    print( "\nDid the test_RunnableDB cmd:\n\n"
         . $cmd
         . "\n\nproduce any results?\n\n"
         . "If you answer no to this question the script will exit (y/n) " );

    if ( !get_input_arg() ) {
        print( $cmd
             . "\n\nDid not give any results, this probably means "
             . "there are issues with your set-up. "
             . "This script will exit you need to investigate.\n\n" );
        exit;
    }
} ## end sub test_run


# Call rulemanager to start the exonerate runs, leaving the set-up
# script.
sub run_analysis {
    eval {

        print("\n\nShall we start the actual analysis? (y/n) ");

        if ( get_input_arg() ) {
            $cmd = "perl " . $CVS_DIR . "/ensembl-pipeline/scripts/rulemanager.pl"
                 . " -dbhost "   . $PIPE_DBHOST
                 . " -dbport "   . $PIPE_DBPORT
                 . " -dbuser "   . $DBUSER
                 . " -dbpass "   . $DBPASS
                 . " -dbname "   . $PIPE_DBNAME;

            print(  "\nSTARTING THE PIPELINE.\n"
                  . "Using command:\n\n"
                  . $cmd
                  . "\n\nPlease monitor the output with:\n\n"
                  . "perl " . $CVS_DIR . "/ensembl-pipeline/scripts/monitor"
                  . " -dbhost "   . $PIPE_DBHOST
                  . " -dbport "   . $PIPE_DBPORT
                  . " -dbuser "   . $DBUSER
                  . " -dbpass "   . $DBPASS
                  . " -dbname "   . $PIPE_DBNAME
                  . " -current_summary"
                  . " -finishedpercent\n\n" );

            print "\n\nIf you stop the rule-manager due to inactivity / hanging you can restart the\n";
            print "cDNA update procedure using :\n\n\t\t\tperl cDNA_update.pl run\n\n" ; 

            system($cmd);
        }
        else {
            print("\nAre you sure you want to interrupt the analysis? (y/n) ");
            if (get_input_arg() ) {
                print("\nProcess interrupted. Not running pipeline.\n\n");
            }
            else {
                if ( !run_analysis() ) {
                    unclean_exit();
                }
            }
        }
    };  ## end eval

    if ($@) {
        print STDERR "\nERROR: $@";
        return 0;
    }
    return 1;
} ## end sub run_analysis

# Identify cdnas which did not align to the genome.
sub find_missing_cdnas {
    # Find all the cdnas which have hits in the database.
    my $db = connect_db( $OUTPUT_DBHOST, $OUTPUT_DBPORT,
                         $OUTPUT_DBNAME, $DBUSER,
                         $DBPASS );

    my $sql = ("SELECT distinct hit_name FROM dna_align_feature");

    my $q1 = $db->dbc->prepare($sql) or croak("Sql error.$!");
    $q1->execute();

    # Make list of cdnas with hits in the database.
    my (%cdna_hits);
    while ( my $cdna = $q1->fetchrow_array ) {
        $cdna_hits{$cdna} = 1;
    }

    # Now go through clipped sequence file and extract those sequences
    # which do not have any hits in the database.
    open( OUT, ">" . $DATA_DIR . "/missing_cdnas.fasta" )
      or croak("Can't open file missing_cdnas.fasta");

    local $/ = "\n>";
    my $cdna_file = $DATA_DIR . "/" . $newfile . ".seqs.clipped";

    open( IN, "<$cdna_file" ) or croak("Can't open file $cdna_file!\n$!");
    while (<IN>) {
        my $seq = $_;

        if ( $seq =~ /(\w+\.\d+)\n/ ) {
            if ( !exists $cdna_hits{$1} ) {
                $seq =~ s/>//g;
                print OUT ">$seq\n";
            }
        }
    }
    close IN;
    close OUT;

    $num_missing_cdnas = `grep -c ">" $DATA_DIR/missing_cdnas.fasta`;
    chomp $num_missing_cdnas;
    return $num_missing_cdnas;

} ## end sub find_missing_cdnas

# Run a check to see if there are any unfinished jobs in the database.
sub chase_jobs {

    # In case have skipped previous sections, reset variables:
    $rerun_flag = 1;
    $chunkDIR   = $DATA_DIR . "/chunks2";

    my $db = connect_db( $PIPE_DBHOST, $PIPE_DBPORT,
                         $PIPE_DBNAME, $DBUSER,
                         $DBPASS );

    # Want to find the list of input files which did not finish
    # running in the db.
    my $sql = (  "SELECT input_id FROM job as j, job_status as s "
               . "WHERE j.job_id = s.job_id && s.is_current = 'y'" );

    my $q1 = $db->dbc->prepare($sql) or croak("Sql error$!\n");
    $q1->execute();

    my %chunks;
    while ( my $file = $q1->fetchrow_array ) {
        $chunks{$file} = 1;
    }

    my $chunk_numbers = keys %chunks;
    if ($chunk_numbers) {
        print "$chunk_numbers chunks did not finish.\n";

        # Store the chunks into a single file:
        open( OUT, ">" . $DATA_DIR . "/single_file.out" )
          or croak("Can't open file single_file.out");

        # Store the list incase need to rerun
        open( LIST, ">" . $chunkDIR . "/went_awol.txt" )
          or croak("Can't open file went_awol.txt");

        my $seq_count = 0;
        foreach my $file ( keys %chunks ) {
            print LIST "$file\n";
            open(IN, $chunkDIR . "/" . $file)
              or croak("Can't open " . $chunkDIR . "/" . $file ."!$!\n");

            while (<IN>) {
                if ( $_ =~ />/ ) { $seq_count++; }
                print OUT "$_";
            }
            close IN;
        }
        close OUT;
        close LIST;

        print( "\nThere were $seq_count cdnas in the files which didn't run.\n"
             . "Would you like to try with smaller chunk files? (y/n) " );

        my $ans = "";
        if ( get_input_arg() ) {
            print(  "Please specify number of chunk files to make "
                  . "(maximum = 1 cdna per file): " );

                my $ans = "";
                chomp( $ans = <STDIN> );

                if ( $ans > $seq_count ) {
                carp("This would give less than 1 sequence per file.\n");
                exit;
            } elsif ( $ans <= $chunk_numbers ) {
                carp("This is the same number of chunk files as last time - "
                    . "it would be better to increase the number of files\n");
                exit;
            } else {
                if ( $newFeatureName =~ /_2/ ) {
                    # To show different run
                    $newFeatureName =~ s/_2/_3/;
                } else {
                    # If have restarted from point after the second run
                    $newFeatureName = $newFeatureName . "_3";
                }
                $configvars{"newFeatureName"} = $newFeatureName;
                $chunkDIR                     = $DATA_DIR . "/chunks3";
                $configvars{"chunkDIR"}       = $chunkDIR;
                $outDIR                       = $DATA_DIR . "/output3";
                $configvars{"outDIR"}         = $outDIR;

                # Check that the new exonerate parameters are set
                config_setup();

                $CHUNK    = $ans;
                $chunkDIR = $DATA_DIR . "/chunks3";
                print("\nSplitting into $CHUNK chunks.\n");
                
                $cmd = "$FASTA_SPLIT $DATA_DIR/single_file.out "
                      . $CHUNK . " " . $chunkDIR; 

                print $cmd . "\n" ; 

                if ( system($cmd) ) {
                    croak("Couldn't split file.$@\n");
                }

                # Isolate biggest sequences
                check_chunksizes();

                print("\nChopped up file.\n\n"
                    . "Set databases for next run? (y/n) ");

                if ( get_input_arg() ) {
                    $rerun_flag = 1;
                    if ( !DB_setup() ) { unclean_exit(); }
                    print("\n\nFinished setting up the analysis.\n");
                }

                if ( !run_analysis() ) {
                    unclean_exit();
                }

                print(  "\nYou should check for any AWOL jobs now, "
                      . "hopefully there won't be any.\n\n" );
            } ## end else
        } ## end if ( get_input_arg() )
    } ## end if ($chunk_numbers)
} ## end sub chase_jobs

# Check the database for those cDNAS which hit many times - might be
# worth adding these to the kill list depending on what they are eg
# LINEs.
sub find_many_hits {
    # Mysql queries involving temporary tables
    my $db = connect_db( $OUTPUT_DBHOST, $OUTPUT_DBPORT,
                         $OUTPUT_DBNAME, $DBUSER,
                         $DBPASS );

    # Make a table containing each transcript matching a cDNA
    my $sql1 =
      (   "CREATE temporary table tmp1 "
        . "SELECT hit_name, exon_transcript.transcript_id "
        . "FROM dna_align_feature, supporting_feature, exon_transcript "
        . "WHERE dna_align_feature_id = feature_id "
        . "AND supporting_feature.exon_id = exon_transcript.exon_id "
        . "GROUP by hit_name, exon_transcript.transcript_id" );

    my $q1 = $db->dbc->prepare($sql1) or croak("Sql error 1.\n$!");
    $q1->execute();

    # Group these to find the number of hits per cDNA
    my $sql2 = (  "CREATE temporary table tmp2 SELECT hit_name, "
                . "count(*) as hits FROM tmp1 GROUP by hit_name" );

    my $q2 = $db->dbc->prepare($sql2) or croak("Sql error 2.\n$!");
    $q2->execute();

    # Examine those which hit more than 20 places
    my $sql3 = ("SELECT * FROM tmp2 WHERE hits > 20 ORDER by hits desc");

    my $q3 = $db->dbc->prepare($sql3) or croak("Sql error 3\n$!");
    $q3->execute();

    my $many_hits_flag = 0;
    while ( my ( $cdna, $hits ) = $q3->fetchrow_array ) {
        print "$cdna\t$hits\n";
        $many_hits_flag = 1;
    }

    if ($many_hits_flag) {
        print(  "\nIt might be worth investigating these sequences to "
              . "see whether these are likely to be genuine hits.\n"
              . "If we don't want them in the database you "
              . "should add them to the kill list\n\n" );
    }
} ## end sub find_many_hits

# Run the script to parse output from ExonerateTranscriptFilter to
# identify reasons for failures.
sub why_cdnas_missed {

    # Make a file containing all of the relevant lines from
    # ExonerateTranscriptFilter.pm outputs.
    my @output = ( "output2", "output3" );
    my $file = $DATA_DIR . "/failed_hits.out";

    open( OUT, ">$file" ) or croak("Can not open file $file");
    close OUT;

    for my $output (@output) {
        `find $DATA_DIR/$output/. | xargs -l1 grep "rpp" >> $file`;
        `find $DATA_DIR/$output/. | xargs -l1 grep "only" >> $file`;
        `find $DATA_DIR/$output/. | xargs -l1 grep "reject" >> $file`;
        `find $DATA_DIR/$output/. | xargs -l1 grep "max_coverage" >> $file`;
    }

    # Need to pass all the variables to the script:
    $cmd = "perl "              . $STORE_UNMAPPED
        . " -gss "              . $GSS
        . " -seq_file "         . $DATA_DIR . "/missing_cdnas.fasta"
        . " -user "             . $DBUSER
        . " -pass "             . $DBPASS
        . " -host "             . $OUTPUT_DBHOST
        . " -port "             . $OUTPUT_DBPORT
        . " -dbname "           . $OUTPUT_DBNAME
        . " -species \""        . $SPECIES  . "\""
        . " -vertrna "          . $DATA_DIR . "/" . $VERTRNA
        . " -refseq "           . $DATA_DIR . "/" . $REFSEQ
        . " -vertrna_update "   . $DATA_DIR . "/" . $VERTRNA_UPDATE
        . " -infile "           . $file
        . " -findN_prog "       . $FIND_N
        . " -reasons_file "     . $UNMAPPED_REASONS;

    if ( system($cmd) ) {
        carp("Erros when running $STORE_UNMAPPED!\n$cmd\n");
    }
    else {
        print("\nUnmapped objects stored.\n");
    }
} ## end sub why_cdnas_missed

# Update the meta-coord table
sub update_metacoord {
    my @table_names = qw(gene
                         exon
                         dna_align_feature
                         transcript);

    my $db = connect_db( $OUTPUT_DBHOST, $OUTPUT_DBPORT,
                         $OUTPUT_DBNAME, $DBUSER,
                         $DBPASS );

    my $sql = "TRUNCATE meta_coord";
    my $sth = $db->dbc->prepare($sql);
    $sth->execute;
    $sth->finish;

    foreach my $table_name (@table_names) {
        my $sql = "INSERT into meta_coord "
                . "SELECT '$table_name', s.coord_system_id, "
                . "max(t.seq_region_end-t.seq_region_start+1) "
                . "FROM $table_name t, seq_region s "
                . "WHERE t.seq_region_id = s.seq_region_id "
                . "GROUP by s.coord_system_id";
        my $sth = $db->dbc->prepare($sql);
        $sth->execute;
        $sth->finish;
    }
    print STDERR "Finished updating meta_coord table\n";
} ## end sub update_metacoord


# Setting various meta table entries.
sub fix_metatable {

    my $db = connect_db( $OUTPUT_DBHOST, $OUTPUT_DBPORT,
                         $OUTPUT_DBNAME, $DBUSER,
                         $DBPASS );

    # Remove previous entries.
    my $sql = "DELETE FROM meta WHERE meta_key like 'genebuild%" . "id'";
    my $sth = $db->dbc->prepare($sql);
    $sth->execute;
    $sth->finish;

    # Set the genebuild id.
    $sql = "INSERT into meta (meta_key, meta_value) "
         . "VALUES ('genebuild_id', '$GENEBUILD_ID' )";
    $sth = $db->dbc->prepare($sql);
    $sth->execute;
    $sth->finish;

    # Set all gene and transcript statuses to putative.
    $sql = "UPDATE gene set status = 'PUTATIVE'";
    $sth = $db->dbc->prepare($sql);
    $sth->execute;
    $sth->finish;

    $sql = "UPDATE transcript set status = 'PUTATIVE'";
    $sth = $db->dbc->prepare($sql);
    $sth->execute;
    $sth->finish;

    # Reload the taxonomy to make sure it's up to date.
    my $cmd = "perl "       . $LOAD_TAX
      . " -name \""         . $SPECIES . "\""
      . " -taxondbhost "    . $TAXONDBHOST
      . " -taxondbport "    . $TAXONDBPORT
      . " -taxondbname "    . $TAXONDBNAME
      . " -lcdbhost "       . $OUTPUT_DBHOST
      . " -lcdbport "       . $OUTPUT_DBPORT
      . " -lcdbname "       . $OUTPUT_DBNAME
      . " -lcdbuser "       . $DBUSER
      . " -lcdbpass "       . $DBPASS;

    if ( system($cmd) ) {
        carp("Errors occurred when running $LOAD_TAX!\n$cmd\n");
    } else {
        print("Meta table fixed.\n");
    }
} ## end sub fix_metatable


# Remove files and database leftovers after analysis,
# restore original config files.
sub clean_up {
    my $option  = shift;
    $status     = 0;

    # Read data dump
    open( RP, "< config_paths.perldata" ) or $status = 1;
    if ( !$status ) {
        undef $/;
        eval <RP>;
        if ($@) {
            $/ = "\n";
            croak("\nCan't recreate data dump.\n$@\n");
        }
        close(RP);
        $/ = "\n";
    } else {
        carp("\nCan't open data dumping file! Already cleaned?\n");
        $status = 0;
    }

    if ( !$option ) {
        # Remove files (fasta, chunks, sql)
        if ( -e $DATA_DIR . "/" . $VERTRNA ) {
            $cmd     = "rm " . $DATA_DIR . "/" . $VERTRNA;
            $status += system($cmd);
        }
        if ( -e $DATA_DIR . "/" . $VERTRNA_UPDATE ) {
            $cmd     = "rm " . $DATA_DIR . "/" . $VERTRNA_UPDATE;
            $status += system($cmd);
        }
        if ( -e $DATA_DIR . "/" . $REFSEQ ) {
            $cmd     = "rm " . $DATA_DIR . "/" . $REFSEQ;
            $status += system($cmd);
        }
        print("\n\nShould we remove the clipped fasta files? (y/n) ");
        if ( get_input_arg() ) {
            if ( -e $DATA_DIR . "/" . $newfile ) {
                $cmd     = "rm " . $DATA_DIR . "/" . $newfile;
                $status += system($cmd);
            }
            if ( -e $DATA_DIR . "/$newfile" . ".seqs" ) {
                $cmd     = "rm " . $DATA_DIR . "/$newfile" . ".seqs";
                $status += system($cmd);
            }
            if ( -e $DATA_DIR . "/$newfile" . ".seqs.clipped" ) {
                $cmd     = "rm " . $DATA_DIR . "/$newfile" . ".seqs.clipped";
                $status += system($cmd);
            }
        }
        if ( -e $DATA_DIR . "/import_tables1.sql" ) {
            $cmd     = "rm " . $DATA_DIR . "/import_tables1.sql";
            $status += system($cmd);
        }
        if ( -e $DATA_DIR . "/import_tables2.sql" ) {
            $cmd     = "rm " . $DATA_DIR . "/import_tables2.sql";
            $status += system($cmd);
        }
        if ( -e $DATA_DIR . "/import_tables3.sql" ) {
            $cmd     = "rm " . $DATA_DIR . "/import_tables3.sql";
            $status += system($cmd);
        }
        # Clean output directories
        if ( !checkdir( $chunkDIR, 0 ) ) {
            carp("Could not prepare directory! [" . $chunkDIR . "]");
        }
        if ( !checkdir( $chunkDIR . "2", 0 ) ) {
            carp("Could not prepare directory! [" . $chunkDIR . "2]");
        }
        if ( !checkdir( $chunkDIR . "3", 0 ) ) {
            carp("Could not prepare directory! [" . $chunkDIR . "3]");
        }
        if ( !checkdir( $outDIR, 0 ) ) {
            carp("Could not prepare directory! [" . $outDIR . "]");
        }
        if ( !checkdir( $outDIR . "2", 0 ) ) {
            carp("Could not prepare directory! [" . $outDIR . "2]");
        }
        if ( !checkdir( $outDIR . "3", 0 ) ) {
            carp("Could not prepare directory! [" . $outDIR . "3]");
        }

        #remove the other temporary output files:

        if ( -e $DATA_DIR . "/missing_cdnas.fasta" ) {
            $cmd     = "rm " . $DATA_DIR . "/missing_cdnas.fasta";
            $status += system($cmd);
        }
        if ( -e $DATA_DIR . "/single_file.out" ) {
            $cmd     = "rm " . $DATA_DIR . "/single_file.out";
            $status += system($cmd);
        }
        if ( -e $DATA_DIR . "/failed_hits.out" ) {
            $cmd     = "rm " . $DATA_DIR . "/failed_hits.out";
            $status += system($cmd);
        }

        if ($status) {
            carp("Error deleting files.\n");
            $status = 0;
        }

        # Remove dbs
        print("\n\nShould we remove the pipeline database? (y/n) ");
        if ( get_input_arg() ) {
            drop_create_database( $PIPE_DBNAME, $PIPE_DBHOST,
                                  $PIPE_DBPORT );
        }

        print "Cleaned out databases, removed temporary files.\n";
    } ## end if ( !$option )

    if (%saved_files) {

        # Restore original config files
        foreach my $config_file ( keys %saved_files ) {
            if ( $saved_files{$config_file} ) {

                print $saved_files{$config_file}, "\n";

                $cmd = "mv "
                  . $DATA_DIR
                  . "/configbackup/"
                  . $config_file . " "
                  . $CVS_DIR . "/"
                  . $saved_files{$config_file};
            } else {
                $cmd = "rm " . $configDIR . "/" . $config_file;
            }
            $status += system($cmd);
        }
    }
    if ($status) { carp("Error restoring config files.\n") }
    print "restored original config files.\n\n";
    if (     ( -e "config_paths.perldata" )
         and ( system("rm config_paths.perldata") ) )
    {
        carp("\ncould not remove perldata file.\n");
    }
} ## end sub clean_up


# Partly clean-up after error
sub unclean_exit {
    clean_up(1);
    print(   "\nRestored original config files.\n"
           . "Check for errors and restart script.\n\n" );
    exit 1;
}


# Compare results to previous data as a health-check
# can also bsub a further function call for every chromosome
sub compare {
    my ( %chromosomes_1,    %chromosomes_2 );
    my ( $sql,              $sql2, $sth1, $sth2 );
    my ( %hits_per_chrom_1, %hits_per_chrom_2 );
    my $hitcount1 = 0;
    my $hitcount2 = 0;

    # Should we exclude all the NT_x-regions?
    my $exclude_NT = 1;

    # Pld alignments
    my $db1 = connect_db( $LAST_DBHOST, $LAST_DBPORT,
                          $LAST_DBNAME, "ensro" );
    # New alignments
    my $db2 = connect_db( $OUTPUT_DBHOST, $OUTPUT_DBPORT,
                          $OUTPUT_DBNAME, "ensro" );

    # Get chromosome names / ids
    $sql = "SELECT coord_system_id "
         . "FROM coord_system "
         . "WHERE name='chromosome' && attrib='default_version'";
    $sth1 = $db1->dbc->prepare($sql) or croak("Sql error!$!\n");
    $sth1->execute();
    my ($coord_system_id1) = $sth1->fetchrow_array;
    $sth2 = $db2->dbc->prepare($sql) or croak("Sql error!$!\n");
    $sth2->execute();
    my ($coord_system_id2) = $sth2->fetchrow_array;

    $sql = "SELECT seq_region_id, name "
         . "FROM seq_region "
         . "WHERE coord_system_id = " . $coord_system_id1;
    if ($exclude_NT) {
        $sql .= " and name not like '%NT%'";
    }
    $sth1 = $db1->dbc->prepare($sql) or croak("Sql error!$!\n");
    $sth1->execute();
    while ( my ( $seq_region_id, $name ) = $sth1->fetchrow_array ) {
        $chromosomes_1{$name} = $seq_region_id;
    }

    $sql = "SELECT seq_region_id, name "
         . "FROM seq_region "
         . "WHERE coord_system_id = " . $coord_system_id2;
    if ($exclude_NT) {
        $sql .= " and name not like '%NT%'";
    }
    $sth2 = $db2->dbc->prepare($sql) or croak("Sql error!$!\n");
    $sth2->execute();
    while ( my ( $seq_region_id, $name ) = $sth2->fetchrow_array ) {
        $chromosomes_2{$name} = $seq_region_id;
    }

    print "Do you want to start the detailed analysis? (y/n) ";
    if ( get_input_arg() ) {

        # Create LSF jobs for in-depth analysis
        print "\nSubmitting jobs for detailed analysis.\n\n";
        foreach my $chromosome ( keys %chromosomes_1 ) {

            $cmd = "bsub -q normal "
                 . "-o "                . $DATA_DIR . "/" . $chromosome . ".out "
                 . "perl "              . $CVS_DIR . "/ensembl-pipeline/scripts/cDNA_update/comparison.pl "
                 . " -chrom "           . $chromosome
                 . " -oldname "         . $OLD_FEATURE_NAME
                 . " -newname "         . $newFeatureName
                 . " -dir "             . $DATA_DIR
                 . " -olddbhost "       . $LAST_DBHOST
                 . " -olddbport "       . $LAST_DBPORT
                 . " -olddbname "       . $LAST_DBNAME
                 . " -newdbhost "       . $OUTPUT_DBHOST
                 . " -newdbport "       . $OUTPUT_DBPORT
                 . " -newdbname "       . $OUTPUT_DBNAME
                 . " -olddnadbhost "    . $LAST_DNADBHOST
                 . " -olddnadbport "    . $LAST_DNADBPORT
                 . " -olddnadbname "    . $LAST_DNADBNAME
                 . " -newdnadbhost "    . $PIPE_DBHOST
                 . " -newdnadbport "    . $PIPE_DBPORT
                 . " -newdnadbname "    . $PIPE_DBNAME;

            if ( system($cmd) ) {
                carp("Error occurred when submitting job!\n$cmd\n\n");
            }
        } ## end foreach my $chromosome ( keys...
    } ## end if ( get_input_arg() )

    print "\nGetting hits per chromosome\n" . "\told\tnew\tdiff\n";

    # Check hits per chromosome
    $sql = "SELECT COUNT(distinct hit_name) "
         . "FROM dna_align_feature daf, analysis a "
         . "WHERE a.logic_name = '" . $OLD_FEATURE_NAME
         . "' and a.analysis_id = daf.analysis_id "
         . "and daf.seq_region_id = ?";
    $sth1 = $db1->dbc->prepare($sql) or croak("Sql error!$!\n");
    $sql = "SELECT COUNT(distinct hit_name) "
         . "FROM dna_align_feature daf, analysis a "
         . "WHERE a.logic_name = '" . $newFeatureName
         . "' and a.analysis_id = daf.analysis_id "
         . "and daf.seq_region_id = ?";

    $sth2 = $db2->dbc->prepare($sql) or croak("Sql error!$!\n");

    my @sorted_chromosomes = sort bychrnum keys %chromosomes_1;
    foreach my $chromosome (@sorted_chromosomes) {
        $sth1->execute( $chromosomes_1{$chromosome} );
        $sth2->execute( $chromosomes_2{$chromosome} );
        $hits_per_chrom_1{$chromosome} = $sth1->fetchrow_array;
        $hits_per_chrom_2{$chromosome} = $sth2->fetchrow_array;
        my $diff =
          $hits_per_chrom_2{$chromosome} - $hits_per_chrom_1{$chromosome};
        print "\n$chromosome:" . "\t"
          . $hits_per_chrom_1{$chromosome} . "\t"
          . $hits_per_chrom_2{$chromosome} . "\t"
          . $diff;
        $hitcount1 += $hits_per_chrom_1{$chromosome};
        $hitcount2 += $hits_per_chrom_2{$chromosome};
    }

    print "\n\nsum:" . "\t" . $hitcount1 . "\t" . $hitcount2 . "\n\n";
} ## end sub compare


# Sort chroms by name
sub bychrnum {
    my @awords = split /_/, $a;
    my @bwords = split /_/, $b;

    my $anum = $awords[0];
    my $bnum = $bwords[0];

    $anum =~ s/chr//;
    $bnum =~ s/chr//;

    if ( $anum !~ /^[0-9]*$/ ) {
        if ( $bnum !~ /^[0-9]*$/ ) {
            return $anum cmp $bnum;
        } else {
            return 1;
        }
    }
    if ( $bnum !~ /^[0-9]*$/ ) {
        return -1;
    }

    if ( $anum <=> $bnum ) {
        return $anum <=> $bnum;
    } else {
        if ( $#awords == 0 ) {
            return -1;
        } elsif ( $#bwords == 0 ) {
            return 1;
        } else {
            return $awords[1] cmp $bwords[1];
        }
    }
} ## end sub bychrnum

# Update the progress status in meta table
sub update_progress_status {
    my $progress_status = shift;

    my $sql = "UPDATE meta set meta_value = '"
            . $progress_status
            . "' WHERE meta_key = 'progress_status'";
    my $sth = $pipe_db->dbc->prepare($sql);
    $sth->execute;
    $sth->finish;
    print(   "\n\nprogress_status has been updated to "
           . "$progress_status.\n\n" );
}

# Create a generic databse
sub create_db {
    my ( $db_name, $is_pipeline ) = ( $_[0], $_[1] );
    my $status = 0;
    my $cmd;

    my $mysql_cmd;
    my ($host, $port);
    if ($is_pipeline) {
        print("\nCreating a pipeline db...\n");
        $host = $PIPE_DBHOST;
        $port = $PIPE_DBPORT;

        $mysql_cmd =  "mysql "
            . "-h $host "
            . "-P $port "
            . "-u $DBUSER "
            . "-p$DBPASS ";

    } else {
        print("\nCreating a target db...\n");
        $host = $OUTPUT_DBHOST;
        $port = $OUTPUT_DBPORT;

        $mysql_cmd =  "mysql "
            . "-h $host "
            . "-P $port "
            . "-u $DBUSER "
            . "-p$DBPASS ";
    }

    # Drop the database if exists and then create it
    # 1 = drop AND create
    drop_create_database($db_name, $host, $port, 1);

    $cmd     = $mysql_cmd . " " . $db_name . " < "
             . $CVS_DIR . "/ensembl/sql/table.sql";
    $status += system($cmd);

    # If it's a pipeline db, add the pipeline tables.
    if ($is_pipeline) {
        $cmd     = $mysql_cmd . "$db_name < "
                 . $CVS_DIR . "/ensembl-pipeline/sql/table.sql";
        $status += system($cmd);
    }

    if($status) { croak("Couldn't create databases!\n"); }

    return $status;
}

# Delete unwanted tables from a database
sub delete_unwanted_tables {
    my ($tables_to_delete, $is_pipeline) = ($_[0], $_[1]);

    # Get an array separator
    my @tables = split(/ /, $tables_to_delete);
    $" = "\n";
    print "\nDeleting the following tables:\n"
        . "@tables\n";

    my $db;
    if ($is_pipeline) {
        $db = connect_db( $PIPE_DBHOST, $PIPE_DBPORT,
                          $PIPE_DBNAME, $DBUSER,
                          $DBPASS );
    } else {
        $db = connect_db( $OUTPUT_DBHOST, $OUTPUT_DBPORT,
                          $OUTPUT_DBNAME, $DBUSER,
                          $DBPASS );
    }

    foreach my $table ( @tables ) {
        $db->dbc->do("DELETE FROM $table");
    }
}

# Dump tables from reference or pipeline db
sub dump_tables {
    my ( $tables_to_dump, $count, $is_pipeline ) = ( $_[0], $_[1], $_[2] );

    my $mysql_dump;
    if ($is_pipeline) {
        $mysql_dump = "mysqldump "
                    . "-h $PIPE_DBHOST "
                    . "-P $PIPE_DBPORT "
                    . "-u $DBUSER "
                    . "-p$DBPASS "
                    . "-t $PIPE_DBNAME ";
        print("\nDumping from pipeline db.\n");
    }
    else {
        # Dumping from the reference db
        $mysql_dump = "mysqldump "
                    . "-h $REF_DBHOST "
                    . "-P $REF_DBPORT "
                    . "-u $DBUSER "
                    . "-p$DBPASS "
                    . "-t $REF_DBNAME ";
        print("\nDumping from ref db.\n");
    }

    # Copy defined db tables from current build
    # (removed analysis table)
    my $file_name = $DATA_DIR . "/import_tables" . $count . ".sql";

    my $cmd;
    if ( -e $file_name ) {
        print(   "\n$file_name exists, "
               . "do you still want to dump the tables? (y/n) " );

        if ( get_input_arg() ) {
            $cmd = $mysql_dump . $tables_to_dump . " > " . $file_name;
            print $cmd, "\n";
            $status += system($cmd);
        }
        else {
            print("\nNot dumping the tables again.\n\n");
        }
    }
    else {
        $cmd = $mysql_dump . $tables_to_dump . " > " . $file_name;
        print $cmd, "\n";
        $status += system($cmd);
    }

    return $status;
}

# Import tables to pipeline or target db
sub import_tables {
    my ( $count, $is_pipeline ) = ( $_[0], $_[1] );

    $status = 0;

    my $mysql_query;
    if ($is_pipeline) {
        $mysql_query = "mysql "
            . "-h $PIPE_DBHOST "
            . "-P $PIPE_DBPORT "
            . "-u $DBUSER "
            . "-p$DBPASS "
            . "-D $PIPE_DBNAME";
        print("\nImporting tables to pipeline db.\n");
    } else {
        $mysql_query = "mysql "
            . "-h $OUTPUT_DBHOST "
            . "-P $OUTPUT_DBPORT "
            . "-u $DBUSER "
            . "-p$DBPASS "
            . "-D $OUTPUT_DBNAME";
        print("\nImporting tables to target db.\n");
    }

    print "\nLoading the tables from file...\n";
    my $cmd = $mysql_query . " < "
            . $DATA_DIR . "/import_tables" . $count . ".sql";

    print $cmd, "\n";
    $status += system($cmd);

    return $status;
}

# Load the external_dbs and unmapped_reason tables from
# /ensembl/misc-scripts.
sub load_misc_script_files {

    my $dir = $CVS_DIR . "/ensembl/misc-scripts/";
    my %table_files = (
            'external_db'     => 'external_db/external_dbs.txt',
            'unmapped_reason' => 'unmapped_reason/unmapped_reason.txt'
    );

    my $target_db = connect_db( $OUTPUT_DBHOST, $OUTPUT_DBPORT,
                                $OUTPUT_DBNAME, $DBUSER,
                                $DBPASS );

    foreach my $table ( keys %table_files ) {
        my $prepare =  "LOAD DATA LOCAL INFILE \'"
                     . $dir . $table_files{$table}
                     . "\' INTO TABLE " . $table;

        my $target_sth = $target_db->dbc->prepare($prepare);
        $target_sth->execute();
        $target_sth->finish();
    }
}

# Add analysis, rule and create the input ids
sub setup_analysis {
    print("Adding analyses and making input_ids...\n\n");

    $status = 0;
    my $sql_pipe = " -dbname " . $PIPE_DBNAME
                 . " -dbhost " . $PIPE_DBHOST
                 . " -dbuser " . $DBUSER
                 . " -dbpass " . $DBPASS
                 . " -dbport " . $PIPE_DBPORT;

    my $cmd = "perl " . $CVS_DIR . "/ensembl-pipeline/scripts/add_Analysis "
            . $sql_pipe
            . " -logic_name "      . $newFeatureName
            . " -program "         . $PROGRAM_NAME
            . " -program_version " . $PROGRAM_VERSION
            . " -program_file "    . $PROGRAM_FILE
            . " -module "          . $MODULE_NAME
            . " -module_version "  . "1 "
            . " -gff_source "      . "Exonerate "
            . " -gff_feature "     . "similarity "
            . " -input_id_type "   . "FILENAME";

    $status += system($cmd);
    print $cmd, "\n";

    $cmd = "perl " . $CVS_DIR . "/ensembl-pipeline/scripts/add_Analysis "
         . $sql_pipe
         . " -logic_name "      . $submitName
         . " -module "          . "dummy"
         . " -input_id_type "   . "FILENAME";

    $status += system($cmd);
    print $cmd, "\n";

    $cmd =  "perl " . $CVS_DIR . "/ensembl-pipeline/scripts/RuleHandler.pl "
          . $sql_pipe
          . " -goal "       . $newFeatureName
          . " -condition "  . $submitName
          . " -insert";

    $status += system($cmd);
    print $cmd, "\n";

    $cmd =  "perl ".$CVS_DIR."/ensembl-pipeline/scripts/make_input_ids "
          . $sql_pipe
          . " -file "
          . " -dir "        . $chunkDIR
          . " -logic_name " . $submitName;

    print $cmd, "\n";
    $status += system($cmd);

    if ($status) { croak("Error while setting up the database.\n"); }
    print "\nDatabase is set up now.\n\n";
} ## end sub setup_analysis

# Drop and create a database
sub drop_create_database {
    my ($db_name, $host, $port, $create) = @_;

    my $serverh = DBI->connect( "DBI:mysql:mysql"
                              . ";host="      . $host
                              . ":port="      . $port
                              . ";user="      . $DBUSER
                              . ";password="  . $DBPASS )
                            or croak($DBI::errstr);

    my $drop_sql = "DROP DATABASE IF EXISTS " . $db_name;
    $serverh->do($drop_sql);

    if ($create) {
        my $create_sql = "CREATE DATABASE " . $db_name;
        $serverh->do($create_sql);
    }

    $serverh->disconnect();
}

# Get the progress_status
sub get_status {
    my ($dbc) = @_;

    my $status_sql = "SELECT meta_value FROM meta "
        . "WHERE meta_key = 'progress_status'";
    my $sth = $dbc->prepare($status_sql);
    $sth->execute;
    my $status = $sth->fetchrow_array;

    return $status;
}

# Clipping the polyA tail
sub polya_clipping {
    my ($trim_file) = @_;

    # Clip ployA tails
    print("\nPerforming polyA clipping...\n");
    my $newfile3 = $DATA_DIR . "/" . $trim_file. ".clipped";
    $cmd = "perl " . $POLYA_CLIPPING . " " ; 
    if ( $MIN_LENGTH ) { 
       $cmd.="-min_length $MIN_LENGTH "; 
    } 
      $cmd .=  $DATA_DIR . "/" . $trim_file . " " . $newfile3;

    print $cmd, "\n";

    if ( system($cmd) ) {
        croak("Couldn't clip file.$@\n");
    }

    # Split fasta files, store into CHUNKDIR
    print("Splitting fasta file.\n");
    $cmd = "$FASTA_SPLIT $newfile3 $CHUNK $chunkDIR";
    if ( system($cmd) ) {
        croak("Couldn't split file.$@\n");
    }

    # Isolate biggest sequences
    check_chunksizes();

    print "\nChopped up the file.\n";
}

# Synchronise reference and target databases
sub sync_databases {

    $status = 0;

    # Copy analysis entries (and others, just to make sure)
    my $tables_to_copy =  'analysis '           . 'assembly '
                        . 'assembly_exception ' . 'attrib_type '
                        . 'coord_system '       . 'meta '
                        . 'meta_coord '         . 'seq_region '
                        . 'seq_region_attrib';

    # Dump the tables from ref db
    # 3 = file count
    # 1 = pipeline db
    $status += dump_tables( $tables_to_copy, 3, 1 );

    # Delete the unwanted tables in target db
    delete_unwanted_tables( $tables_to_copy );

    # Import the dumped tables to target db
    # 3 = count of the output file
    $status += import_tables(3);

    if ($status) { croak("Error while synchronising databases.\n") }

    print "\n\nDatabases in sync.\n";
}

# Connect to a given database, optional with attached DNA db
sub connect_db {
    my $host   = shift;
    my $port   = shift;
    my $dbname = shift;
    my $user   = shift;
    my $pass   = shift;
    my $dnadb  = shift;
    my $dbObj;

    if ($dnadb) {
        $dbObj =
          new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $host
                                                -port => $port,
                                              -dbname => $dbname,
                                              -user   => $user,
                                              -pass   => $pass,
                                              -dnadb  => $dnadb );
    } else {
        $dbObj =
          new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $host,
                                              -port   => $port,
                                              -dbname => $dbname,
                                              -user   => $user,
                                              -pass   => $pass );
    }
    if ( !$dbObj ) {
        return 0;
    }
    return $dbObj;
} ## end sub connect_db


1;


__END__
