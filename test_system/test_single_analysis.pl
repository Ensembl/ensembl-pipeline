=pod

=head1 NAME

test_single_analysis.pl. This is a script which will test a single analysis
inside the pipeline system. This is not an atomic test but instead
tests the whole process by using the job_submission script.

=head1 SYNOPSIS

The test_single_analysis script invokes first a TestDB object which sets up
a test database then a RunTest object which knows how to fill the
database with data and how to invoke the job_submission script to run
the specified analysis. At the end of a run the RunTest can compare
the results to a reference set of results if specifed.

Note currently the script must be run in this directory. To run the test properly,
you will have to alter your perl5lib and blastdb settings. They should always be 
returned to their previous state after the script is run though.

=head1 OPTIONS

  -species              the species you are running the test on and this tells
   the TestDB object which zip file in reference_data to use

  -verbose              toggle to indicate whether to be verbose

  -logic_name           logic_name of the analysis you wish to run

  -feature_table        name of the table which should be filled by this analysis                     #at6 (8 Dec): what happens if >1 table filled!?

  -table_to_load        names of the tables which should be filled before this
        analysis can run. As standard the core and pipeline tables
        are filled (see TestDB::load_core_tables method). If other tables need
        to be filled they need to be specified using this option which can
        appear multiple times on the commandline. Note there is a method
        table_groups which contains a list of tables to fill for specific
        analyses based on logic_name. If your analysis appears in this method
        you don't need to give any tables in 'table_to_load'

  -output_dir           the directory the job output will be written to. 
   Otherwise the DEFAULT_OUTPUT_DIR from BatchQueue.pm is used

  -queue_manager        the BatchSubmission module to use otherwise the
   QUEUE_MANAGER from BatchQueue.pm is used

  -run_comparison       toggle to indicate to run comparison with reference data
   set

  -conf_file            the name of the conf file to use when setting up the test DB to
   run the analysis.  By default TestDB.conf is used

  -comparison_conf      the name of the conf file to use when setting up the
   database containing reference data for comparison. RefDB.conf is used by default

  -dont_cleanup         a toggle to indicate not to delete the database and output
   directories

  -blastdb              a string to change the blastdb to

  -job_submission_verbose whether to make the job submission script verbose


=head1 EXAMPLES

./test_single_analysis -logic_name CpG -run_comparison

./test_single_analysis -logic_name Genscan (will get the correct tables
from the table_groups method)


=head1 SEE ALSO

test_whole_pipeline.pl
docs/running_tests.txt
ensembl-pipeline/scripts/job_submission.pl

any questions pass to ensembl-dev@ebi.ac.uk

=cut

#!/usr/local/ensembl/bin/perl -w

use lib './modules';  
use lib './config';  

use strict;
use Bio::EnsEMBL::Utils::Exception qw(throw warning verbose);
use Bio::EnsEMBL::Pipeline::Config::BatchQueue;  # BatchQueue.pm from "./config" dir
use TestDB;
use Environment; 
use RunTest;
use Getopt::Long;

my $species;
our $verbose;
my $logic_name;
my $feature_table;
my @tables_to_load;
my $output_dir;
my $queue_manager;
my $conf_file;
my $dont_cleanup;
my $blastdb;
my $job_submission_verbose;
my $run_comparison;
my $comparison_conf;
my $help;
my $local = 0;
&GetOptions(
            'species:s' => \$species,
            'verbose!' => \$verbose,
            'logic_name:s' => \$logic_name,
            'feature_table:s' => \$feature_table,
            'table_to_load:s@' => \@tables_to_load,
            'output_dir:s' => \$output_dir,
            'run_comparison!' => \$run_comparison,
            'comparison_conf' => \$comparison_conf,
            'conf_file:s' => \$conf_file,
            'dont_cleanup!' => \$dont_cleanup,
            'local' => \$local,
            'blastdb:s' => \$blastdb,
            'job_submission_verbose!' => \$job_submission_verbose,
            'help!' => \$help,
           ) or perldoc();



if(@tables_to_load == 0){
  @tables_to_load = @{table_groups($logic_name)};
}
print "Tables to load are: @tables_to_load\n";
perldoc() if($help);
if(!$logic_name || !$feature_table){
  throw("Must specific which analysis you wish to run and what table its ".
        "results are written to with -logic_name and -feature_table");
}

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

my $extra_perl = $testdb->curr_dir."/config"; #to be used to modify PERL5LIB in $environment

my $test_runner = RunTest->new(
                               -TESTDB => $testdb,
                               -ENVIRONMENT => $environment,
                               -OUTPUT_DIR => $output_dir,
                               -EXTRA_PERL => $extra_perl,
                               -BLASTDB => $blastdb,
                               -TABLES => \@tables_to_load,
                               -QUEUE_MANAGER => $QUEUE_MANAGER,
                               -DONT_CLEANUP => $dont_cleanup,
                               );


if($run_comparison){
  $test_runner->comparison_conf($comparison_conf);
}
$test_runner->run_single_analysis($logic_name, $feature_table, 
                                  $job_submission_verbose);


sub table_groups{
  my ($logic_name) = @_;
  $logic_name = lc($logic_name);
  my %tables;
  $tables{'genscan'} = ['core', 'pipeline', 'repeat_feature',
                        'repeat_consensus'];
  $tables{'uniprot'} = ['core', 'pipeline', 'repeat_feature',
                        'repeat_consensus', 'prediction_exon',
                        'prediction_transcript'];
  $tables{'vertrna'} = ['core', 'pipeline', 'repeat_feature',
                      'repeat_consensus', 'prediction_exon',
                      'prediction_transcript'];
  $tables{'unigene'} = ['core', 'pipeline', 'repeat_feature',
                        'repeat_consensus', 'prediction_exon',
                        'prediction_transcript'];
  $tables{'firstef'} = ['core', 'pipeline', 'repeat_feature',
                         'repeat_consensus'];
  $tables{'marker'} = ['core', 'pipeline','marker', 
                        'marker_synonym', 'marker_map_location'];
  $tables{'bestpmatch'} = ['core', 'pipeline', 
                            'protein_align_feature.Pmatch'];
  if($tables{$logic_name}){
    return $tables{$logic_name};
  }else{
    return ['core', 'pipeline'];
  }
}


sub perldoc{
	exec('perldoc', $0);
	exit(0);
}


