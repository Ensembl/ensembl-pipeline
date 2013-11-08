#!/usr/bin/env perl

=pod

=head1 NAME

test_single_analysis.pl.  This is a script which will test a single
analysis inside the pipeline system.  This is not an atomic test but
instead tests the whole process by using the job_submission script.

=head1 SYNOPSIS

The test_single_analysis script invokes first a TestDB object which sets
up a test database then a RunTest object which knows how to fill the
database with data and how to invoke the job_submission script to run
the specified analysis.  At the end of a run the RunTest can compare the
results to a reference set of results if specifed.

Note currently the script must be run in this directory
(ensembl-pipeline/test_system).  When the script starts
running the test, PERL5LIB is automatically altered to include
ensembl-pipeline/test_system/config as the first place to look for Perl
modules.  If -blastdb is used on the commandline to specify a particular
directory for BLASTDB location, BLASTDB will be altered too for the
duration of the test.  Once the test has finished running, the test
system will automatically reset the PERL5LIB and BLASTDB settings.  All
these changes in the environment variables (PERL5LIB and BLASTDB) were
done by ./modules/Environment.pm while using the "run_single_analysis"
method in RunTest.pm.

=head1 OPTIONS

The compulsory options are "-logic_name" and "-feature_table".  The rest
are optional.


  -species
        The species you are running the test on and this tells the
        TestDB object which zip file in reference_data to use. By
        default, "homo_sapiens" is used.

  -verbose
        Toggle to indicate whether to be verbose.  See also notes on the
        "-run_comparison" flag below.

  -logic_name
        The logic name of the analysis you wish to run.

  -feature_table
        The name of the table which should be filled by this analysis.
        The specified table should always have a column named
        "analysis_id" as the system relies on the analysis_id to work
        out which rows in the table were actually generated during the
        test and then reports the count of number of results from the
        test.  This helps the system to distinguish test-generated
        results from those which had been loaded from text files
        prior to the test.  Most of the appropriate tables for this
        option will have names like xxxx_feature, e.g. repeat_feature,
        simple_feature, marker_feature, etc.

        For RepeatMask, use "repeat_feature" (not "repeat_consensus").
        For Genscan, use "prediction_transcript", (not
        "prediction_exon").

  -table_to_load
        The names of the tables which should be filled before
        this analysis can run.  Note that for every analysis,
        the core and pipeline tables are loaded by default (see
        TestDB::load_core_tables method).  Also, there is a method
        "table_groups" which contains a list of tables to fill for
        specific analyses based on logic_name.  If your analysis
        appears in "table_groups", you don't need to give any tables
        in 'table_to_load' either.  Use this option ONLY IF any other
        tables need to be filled (i.e. not core, not pipeline, not in
        "table_groups" method).  This option can appear multiple times
        on the command line, e.g., -table_to_load table1 -table_to_load
        table2

  -output_dir
        The directory the job output will be written to.  Otherwise,
        the DEFAULT_OUTPUT_DIR from BatchQueue.pm will be used.  Note
        that if 'output_dir' is specified on the command line, output
        will NOT be sent to analysis-specific output directories
        specified in BatchQueue.pm.  However, DEFAULT_OUTPUT_DIR and
        analysis-specific directories in BatchQueue.pm will still be
        created by the system if they don't exist.  See documentation in
        config/Bio/EnsEMBL/Pipeline/Config/BatchQueue.pm for more info.

  -queue_manager
        The BatchSubmission module to use.  If not specified, the
        QUEUE_MANAGER from BatchQueue.pm is used.

  -run_comparison
        Toggle to indicate to run comparison with reference data
        set.  Please note in some cases, where the analysis generates
        thousands of features in the results, the comparsion will take
        quite a while as every single query feature has to be checked
        against the reference.  Also, if the -verbose flag is also used,
        all the query-target matching pairs (a large number of them)
        will be printed on screen.  Unless the details of each matching
        pair are required, or else it's recommended to turn verbose off.

  -conf_file
        The name of the conf file to use when setting up the test DB to
        run the analysis.  By default TestDB.conf is used.

  -comparison_conf
        The name of the conf file to use when setting up the database
        containing reference data for comparison.  By default,
        RefDB.conf is used.

  -dont_cleanup
        A toggle to indicate not to delete the database, the unpacked
        reference data directory or test output directories.

  -blastdb
        A string to change the BLASTDB to.

  -job_submission_verbose
        Whether to make the job submission script verbose.

=head1 EXAMPLES

  ./test_single_analysis.pl -logic_name CpG -run_comparison

  Will get the correct tables from the table_groups method:
  ./test_single_analysis.pl -logic_name Genscan

=head1 SEE ALSO

  test_whole_pipeline.pl
  docs/running_tests.txt
  ensembl-pipeline/scripts/job_submission.pl

  any questions pass to dev@ensembl.org

=cut

use strict;
use warnings;

use lib './modules';
use lib './config';

use Bio::EnsEMBL::Utils::Exception qw(throw warning verbose);
use Bio::EnsEMBL::Pipeline::Config::BatchQueue;
use TestDB;
use Environment;
use RunTest;
use Getopt::Long;

my $blastdb         = $ENV{BLASTDB} || '/data/blastdb/Supported';
my $comparison_conf = 'RefDB.conf';
my $species         = 'homo_sapiens';

my $conf_file;
my $dont_cleanup;
my $feature_table;
my $job_submission_verbose;
my $logic_name;
my $output_dir;
my $queue_manager;
my $run_comparison;
my @tables_to_load;
our $verbose;

my $help  = 0;
my $local = 0;

if ( !GetOptions( 'species=s'               => \$species,
                  'verbose!'                => \$verbose,
                  'logic_name=s'            => \$logic_name,
                  'feature_table=s'         => \$feature_table,
                  'table_to_load=s@'        => \@tables_to_load,
                  'output_dir=s'            => \$output_dir,
                  'run_comparison!'         => \$run_comparison,
                  'comparison_conf'         => \$comparison_conf,
                  'conf_file=s'             => \$conf_file,
                  'dont_cleanup!'           => \$dont_cleanup,
                  'local'                   => \$local,
                  'blastdb=s'               => \$blastdb,
                  'job_submission_verbose!' => \$job_submission_verbose,
                  'help!'                   => \$help,
     ) ||
     $help )
{
  perldoc();
}

if ( !defined($logic_name) || !defined($feature_table) ) {
  throw( "Must specific which analysis " .
         "you wish to run and what table its " .
         "results are written to with " .
         "-logic_name and -feature_table" );
}

if ( @tables_to_load == 0 ) {
  @tables_to_load = @{ table_groups($logic_name) };
}

print("Tables to load are: @tables_to_load\n");

my $testdb = TestDB->new( -SPECIES   => $species,
                          -VERBOSE   => $verbose,
                          -CONF_FILE => $conf_file,
                          -LOCAL     => $local, );

my $environment = Environment->new( $testdb, $verbose );

# To be used to modify PERL5LIB in $environment
my $extra_perl = sprintf( "%s/config", $testdb->curr_dir() );

my $test_runner = RunTest->new( -TESTDB        => $testdb,
                                -ENVIRONMENT   => $environment,
                                -OUTPUT_DIR    => $output_dir,
                                -EXTRA_PERL    => $extra_perl,
                                -BLASTDB       => $blastdb,
                                -TABLES        => \@tables_to_load,
                                -QUEUE_MANAGER => $QUEUE_MANAGER,
                                -DONT_CLEANUP  => $dont_cleanup,
                                -VERBOSE       => $verbose, );

if ($run_comparison) {
  $test_runner->comparison_conf($comparison_conf);
}

$test_runner->run_single_analysis( $logic_name, $feature_table,
                                   $job_submission_verbose );

sub table_groups {
  my ($logic_name) = @_;

  $logic_name = lc($logic_name);

  my %tables = (
         '__default' => [ 'core', 'pipeline' ],
         'genscan' =>
           [ 'core', 'pipeline', 'repeat_feature', 'repeat_consensus' ],
         'uniprot' => [ 'core',            'pipeline',
                        'repeat_feature',  'repeat_consensus',
                        'prediction_exon', 'prediction_transcript'
         ],
         'vertrna' => [ 'core',            'pipeline',
                        'repeat_feature',  'repeat_consensus',
                        'prediction_exon', 'prediction_transcript'
         ],
         'unigene' => [ 'core',            'pipeline',
                        'repeat_feature',  'repeat_consensus',
                        'prediction_exon', 'prediction_transcript'
         ],
         'firstef' =>
           [ 'core', 'pipeline', 'repeat_feature', 'repeat_consensus' ],
         'marker' => [ 'core',           'pipeline',
                       'marker',         'map',
                       'marker_synonym', 'marker_map_location'
         ],
         'bestpmatch' =>
           [ 'core', 'pipeline', 'protein_align_feature.Pmatch' ], );

  if ( exists( $tables{$logic_name} ) ) {
    return $tables{$logic_name};
  }
  else {
    return $tables{'__default'};
  }
} ## end sub table_groups

sub perldoc {
  exec( 'perldoc', $0 );
  exit(0);
}
