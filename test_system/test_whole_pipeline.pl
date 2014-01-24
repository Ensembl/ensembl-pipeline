#!/usr/bin/env perl

=pod

=head1 NAME

This is a script which will test running a whole pipeline start to end.
It does this by running the rulemanager script.

=head1 SYNOPSIS

The test_whole_pipeline script invokes first a TestDB object which sets
up a test database, and then a RunTest object which knows how to fill
the database with data and how to invoke the rulemanager script to run
the pipeline.

Note currently the script must be run in this directory
(ensembl-pipeline/test_system).  When the script starts
running the test, PERL5LIB is automatically altered to include
ensembl-pipeline/test_system/config as the first place to look for
Perl modules.  If -blastdb is used on the commandline to specify a
particular dir for BLASTDB location, BLASTDB will be altered too for
the duration of the test.  Once the test has finished running, the test
system will automatically reset the PERL5LIB and BLASTDB settings.  All
these changes in the environment variables (PERL5LIB and BLASTDB) were
done by modules/Environment.pm while using the "run_pipeline" method in
RunTest.pm.

=head1 OPTIONS

  -species
      This is the species you are running the test on and this tells the
      TestDB object which zip file in reference_data to use.

  -verbose
      Toggle to indicate whether to be verbose.

  -output_dir
      The directory the job output will be written to.  Otherwise
      the DEFAULT_OUTPUT_DIR from BatchQueue.pm will be used.  Note
      if output_dir is specified on the command line, output will
      NOT be sent to analysis-specific output directories specified
      in BatchQueue.pm either.  However, DEFAULT_OUTPUT_DIR and
      analysis-specific directories in BatchQueue.pm will still be
      created by the system if they don't exist.  See documentation in
      config/Bio/EnsEMBL/Pipeline/Config/BatchQueue.pm for more info.

  -queue_manager
      The BatchSubmission module to use otherwise the QUEUE_MANAGER from
      BatchQueue.pm.

  -conf_file
      The name of the conf file to use when running the analysis.  By
      default TestDB.conf is used.

  -dont_cleanup
      A toggle to indicate not to delete the database, the unzipped
      reference data directory and test output directories.

  -blastdb
      A string to change the BLASTDB to.

  -rulemanager_verbose
      Whether to make the job submission script verbose.

=head1 SEE ALSO

  test_single_analysis.pl
  docs/running_tests.txt
  ensembl-pipeline/scripts/job_submission.pl

  any questions pass to http://lists.ensembl.org/mailman/listinfo/dev

=cut

use strict;
use warnings;

use lib ( './config', './modules' );

use Bio::EnsEMBL::Utils::Exception qw(throw warning verbose);
use Bio::EnsEMBL::Pipeline::Config::BatchQueue;
use TestDB;
use Environment;
use RunTest;
use Getopt::Long;

my $species;
our $verbose;
my $feature_table;
my @tables_to_load = ( 'core',           'pipeline',
                       'marker',         'map',
                       'marker_synonym', 'marker_map_location' );
my $output_dir;
my $queue_manager;
my $conf_file;
my $dont_cleanup;
my $blastdb;
my $job_submission_verbose;
my $run_comparison;
my $local = 0;
my $comparison_conf;
my $help;

&GetOptions( 'species:s'            => \$species,
             'verbose!'             => \$verbose,
             'output_dir:s'         => \$output_dir,
             'queue_manager:s'      => \$queue_manager,
             'conf_file:s'          => \$conf_file,
             'dont_cleanup!'        => \$dont_cleanup,
             'blastdb:s'            => \$blastdb,
             'rulemanager_verbose!' => \$job_submission_verbose,
             'comparison_conf'      => \$comparison_conf,
             'run_comparison!'      => \$run_comparison,
             'local'                => \$local,
             'help!'                => \$help,
  ) or
  perldoc();

perldoc() if ($help);

$output_dir      = $DEFAULT_OUTPUT_DIR     if ( !$output_dir );
$queue_manager   = $QUEUE_MANAGER          if ( !$queue_manager );
$species         = 'homo_sapiens'          if ( !$species );

if ( !defined($blastdb) ) {
  $blastdb = $ENV{'BLASTDB'} || '/data/blastdb/Supported';
}

$comparison_conf = 'RefDB.conf'            if ( !$comparison_conf );

my $testdb = TestDB->new( -SPECIES   => $species,
                          -VERBOSE   => $verbose,
                          -CONF_FILE => $conf_file,
                          -LOCAL     => $local, );

my $environment = Environment->new( $testdb, $verbose );

my $extra_perl = $testdb->curr_dir . "/config";

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

$test_runner->run_pipeline($job_submission_verbose);

sub perldoc {
  exec( 'perldoc', $0 );
  exit(0);
}
