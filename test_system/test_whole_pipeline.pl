#!/usr/local/ensembl/bin/perl -w

use lib './modules';
use lib './config';

use strict;
use Bio::EnsEMBL::Utils::Exception qw(throw warning verbose);
use Bio::EnsEMBL::Pipeline::Config::BatchQueue;
use TestDB;
use Environment;
use RunTest;
use Getopt::Long;

my $species;
our $verbose;
my $feature_table;
my @tables_to_load = ('core', 'pipeline');
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

&GetOptions(
            'species:s' => \$species,
            'verbose!' => \$verbose,
            'output_dir:s' => \$output_dir,
            'queue_manager:s' => \$queue_manager,
            'conf_file:s' => \$conf_file,
            'dont_cleanup!' => \$dont_cleanup,
            'blastdb:s' => \$blastdb,
            'rulemanager_verbose!' => \$job_submission_verbose,
            'comparison_conf' => \$comparison_conf,
            'run_comparison!' => \$run_comparison,
            'local' => \$local,
            'help!' => \$help,
           ) or perldoc();




perldoc() if($help);

$output_dir = $DEFAULT_OUTPUT_DIR if(!$output_dir);
$queue_manager = $QUEUE_MANAGER if(!$queue_manager);
$species = 'homo_sapiens' if(!$species);
$blastdb = '/data/blastdb/Ensembl' if(!$blastdb);
$comparison_conf = 'RefDB.conf' if(!$comparison_conf);

my $testdb = TestDB->new(
                         -SPECIES => $species, 
                         -VERBOSE => $verbose,
                         -CONF_FILE => $conf_file,
                         -LOCAL => $local,
                        );

my $environment = Environment->new($testdb, $verbose);

my $extra_perl = $testdb->curr_dir."/config";

my $test_runner = RunTest->new(
                               -TESTDB => $testdb,
                               -ENVIRONMENT => $environment,
                               -OUTPUT_DIR => $DEFAULT_OUTPUT_DIR,
                               -EXTRA_PERL => $extra_perl,
                               -BLASTDB => $blastdb,
                               -TABLES => \@tables_to_load,
                               -QUEUE_MANAGER => $QUEUE_MANAGER,
                               -DONT_CLEANUP => $dont_cleanup,
                               );

if($run_comparison){
  $test_runner->comparison_conf($comparison_conf);
}

$test_runner->run_pipeline($job_submission_verbose);

sub perldoc{
	exec('perldoc', $0);
	exit(0);
}

=pod

=head1 NAME

test_whole_pipeline.pl, This is a script which will test running a 
whole pipeline start to end. It does this by running the rulemanager
script

=head1 SYNOPSIS

The test_whole_pipeline script invokes first a TestDB object which sets up
a test database then a RunTest object which knows how to fill the
database with data and how to invoke the rulemanager script to run
the pipeline. 

Note currently the script must be run in this directory. This script
also has to alter your perl5lib and blastdb to run the test properly
they should always be returned to their previous state after the script
is run though

=head1 OPTIONS

  -species, this is the species you are running the test on and this tells
   the TestDB object which zip file in reference_data to us
  -verbose toggle to indicate whether to be verbose
  -output_dir the directory the job output will be written to
   otherwise the DEFAULT_OUTPUT_DIR from BatchQueue.pm is used
  -queue_manager the BatchSubmission module to use otherwise the 
   QUEUE_MANAGER from BatchQueue.pm
  -conf_file the name of the conf file to use when running the analysis
   by default TestDB.conf is used
  -dont_cleanup a toggle to indicate not to delete the database and output
   directories
  -blastdb what string to change the blastdb to
  -job_submission_verbose whether to make the job submission script verbose


=head1 SEE ALSO

test_single_analysis.pl
docs/running_tests.txt
ensembl-pipeline/scripts/job_submission.pl

any questions pass to ensembl-dev@ebi.ac.uk

=cut
