#!/usr/local/ensembl/bin/perl -w

use strict;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use Bio::EnsEMBL::Pipeline::Utils::PipelineSanityChecks;
use Bio::EnsEMBL::Pipeline::RuleManager;


$| = 1;

my $dbhost    = $ENV{'ENS_DBHOST'};
my $dbname    = $ENV{'ENS_DBNAME'};
my $dbuser    = $ENV{'ENS_DBUSER'};
my $dbpass    = $ENV{'ENS_DBPASS'};
my $dbport    = $ENV{'ENS_DBPORT'} || 3306;
my $help;
my $verbose;
my $queue_manager;
my $max_job_time;
my $min_job_time;
my $sleep_per_job;
my $runner;
my $output_dir;
my $job_limit;
my $mark_awol = 1;
my $rename_on_retry = 1;
my $once;
my @starts_from;
my @analyses_to_run;
my @analyses_to_skip;
my @types_to_run;
my @types_to_skip;
my $ids_to_run;
my $ids_to_skip;
my $config_sanity = 1;
my $accumulator_sanity = 1;
my $rules_sanity = 1;
my $db_sanity = 1;
my $kill_jobs;
my $killed_time;
my $kill_file;
my $rerun_sleep = 300;

GetOptions(
           'dbhost=s'      => \$dbhost,
           'dbname=s'      => \$dbname,
           'dbuser=s'      => \$dbuser,
           'dbpass=s'      => \$dbpass,
           'dbport=s'      => \$dbport,
           'help!' => \$help,
           'verbose!' => \$verbose,
           'queue_manager=s' => \$queue_manager,
           'max_job_time=s' => \$max_job_time,
           'min_job_job=s' => \$min_job_time,
           'sleep_per_job=s' => \$sleep_per_job,
           'runner=s' => \$runner,
           'output_dir=s' => \$output_dir,
           'job_limit=s' => \$job_limit,
           'mark_awol!' => \$mark_awol,
           'rename_on_retry' => \$rename_on_retry,
           'starts_from=s@' => \@starts_from,
           'analysis=s@' => \@analyses_to_run,
           'skip_analysis=s@' => \@analyses_to_skip,
           'input_id_type=s@' => \@types_to_run,
           'skip_input_id_type=s@' => \@types_to_skip,
           'input_id_file=s' => \$ids_to_run,
           'skip_input_id_file=s' => \$ids_to_skip,
           'config_sanity!' => \$config_sanity,
           'accumulator_sanity!' => \$accumulator_sanity,
           'db_sanity!' => \$db_sanity,
           'rules_sanity!' => \$rules_sanity,
           'kill_jobs!' => \$kill_jobs,
           'killed_time=s' => \$killed_time,
           'kill_file=s' => \$kill_file,
           'rerun_sleep=s' => \$rerun_sleep,
           ) or throw("Can't get rulemanager commandline arguments");


my $db = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new(
                                                       -host   => $dbhost,
                                                       -dbname => $dbname,
                                                       -user   => $dbuser,
                                                       -pass   => $dbpass,
                                                       -port   => $dbport,
                                                      );

my $sanity = Bio::EnsEMBL::Pipeline::Utils::PipelineSanityChecks->new
  (
   -DB => $db,
  );
my $rulemanager = Bio::EnsEMBL::Pipeline::RuleManager->new
  (
   -DB => $db,
   -QUEUE_MANAGER => $queue_manager,
   -MARK_AWOL => $mark_awol,
   -RENAME_ON_RETRY => $rename_on_retry,
   -MAX_JOB_SLEEP => $max_job_time,
   -MIN_JOB_SLEEP => $min_job_time,
   -SLEEP_PER_JOB => $sleep_per_job,
   -VERBOSE => $verbose,
   -RUNNER => $runner,
   -OUTPUT_DIR => $output_dir,
   -JOB_LIMIT => $job_limit,
   );
   
