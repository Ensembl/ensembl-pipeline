package RunTest;

use lib './modules';
use lib './config';
use strict;
use warnings;
use File::Path;
use Bio::EnsEMBL::Utils::Exception qw(throw warning verbose);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::Config::BatchQueue;
use TestDB;
use FeatureComparison;
use Bio::EnsEMBL::Pipeline::Monitor;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Root);



=head2 new

  Arg [1]   : RunTest
  Arg [2]   : TestDB object
  Arg [3]   : Environment object
  Arg [4]   : string, path to output directory
  Arg [5]   : string, extra info for perl5lib
  Arg [6]   : arrayref, list of tables to import
  Arg [7]   : string, name of Bio::EnsEMBL::Pipeline::BatchSubmission
  module to require
  Arg [8]   : int, toggle whether to delete the databases and directories
  created
  Function  : create a RunTest object 
  Returntype: RunTest
  Exceptions: throws if not passed in a TestDB object or it isnt a TestDB
  object, if not passed in an Environment object or it isnt an Environemt
  object, throws if it cant require the queue manager, or if not passed any
  tables to import or if the variable passed isnt an array ref
  Example   : my $runtest = RunTest->new
  (
   -TESTDB => $testdb,
   -ENVIRONMENT => $environment,
   -OUTPUT_DIR => $DEFAULT_OUTPUT_DIR,
   -EXTRA_PERL => $extra_perl,
   -BLASTDB => $blastdb,
   -TABLES => \@tables_to_load,
   -QUEUE_MANAGER => $QUEUE_MANAGER,
   -DONT_CLEANUP => $dont_cleanup,
  );

=cut




sub new{
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  &verbose('WARNING');
  my ($testdb, $environment, $output_dir,
      $extra_perl, $tables, $queue_manager,
      $cleanup, $verbose, $conf, $blastdb) =
        rearrange(['TESTDB', 'ENVIRONMENT', 'OUTPUT_DIR',
                   'EXTRA_PERL', 'TABLES', 'QUEUE_MANAGER',
                   'DONT_CLEANUP', 'VERBOSE', 'COMPARISON_CONF',
                   'BLASTDB'], @args);
  if(!$testdb || !$testdb->isa('TestDB')){
    throw("Can't run without a TestDB object or a with a ".$testdb);
  }
  if(!$environment || !$environment->isa('Environment')){
    throw("Can't run without an Environment object or with a ".
          $environment);
  }
  if(!$output_dir){
    $output_dir = $DEFAULT_OUTPUT_DIR;
  }
  if(!$queue_manager){
    $queue_manager = $QUEUE_MANAGER; #found in BatchQueue.pm
  }
  my $batch_q_module = 
    "Bio::EnsEMBL::Pipeline::BatchSubmission::$queue_manager";
  my $file = "$batch_q_module.pm";
  $file =~ s{::}{/}g;
  eval {
    require "$file";
  };
  if($@){
    $self->exception("Can't find $file [$@]");
  }
  if(!$batch_q_module->can('job_stats')){
    $self->exception($batch_q_module." doesn't have the job_stats method ".
                     "this won't work");
  }
  if(!$tables || ref($tables) ne 'ARRAY'){
    $self->exception("Must define some tables or table groups to load and".
                     " this must passed in as an array ref $tables");
  }
  $self->testdb($testdb);
  $self->environment($environment);
  $self->output_dir($output_dir);
  $self->extra_perl($extra_perl);
  $self->tables($tables);
  $self->queue_manager($batch_q_module);
  $self->dont_cleanup_tests($cleanup);
  $self->verbosity($verbose);
  $self->comparison_conf($conf);
  $self->blastdb($blastdb);
  return $self;
}



=head2 containers

  Arg [1]   : RunTest
  Arg [2]   : variable, normally a string
  Function  : Containers for the variables passed into the constructor
  Returntype: variable stored in container
  Exceptions: none
  Example   : my $testdb = $self->testdb;

=cut



sub testdb{
  my $self = shift;
  $self->{'testdb'} = shift if(@_);
  return $self->{'testdb'};
}

sub ref_testdb{
  my $self = shift;
  $self->{'ref_testdb'} = shift if(@_);
  return $self->{'ref_testdb'};
}

sub blastdb{
  my $self = shift;
  $self->{'blastdb'} = shift if(@_);
  return $self->{'blastdb'};
}

sub environment{
  my $self = shift;
  $self->{'environment'} = shift if(@_);
  return $self->{'environment'};
}
sub output_dir{
  my $self = shift;
  $self->{'output_dir'} = shift if(@_);
  return $self->{'output_dir'};
}
sub extra_perl{
  my $self = shift;
  $self->{'extra_perl'} = shift if(@_);
  return $self->{'extra_perl'};
}
sub tables{
  my $self = shift;
  $self->{'tables'} = shift if(@_);
  return $self->{'tables'};
}
sub queue_manager{
  my $self = shift;
  $self->{'queue_manager'} = shift if(@_);
  return $self->{'queue_manager'};
}
sub dont_cleanup_tests{
  my $self = shift;
  $self->{'cleanup_tests'} = shift if(@_);
  return $self->{'cleanup_tests'};
}

sub cleanup_dir{
  my $self = shift;
  $self->{'cleanup_dir'} = shift if(@_);
  return $self->{'cleanup_dir'};
}

sub verbosity{
  my $self = shift;
  $self->{'verbose'} = shift if(@_);
  return $self->{'verbose'};
}
sub comparison_conf{
  my $self = shift;
  $self->{'comparison_conf'} = shift if(@_);
  return $self->{'comparison_conf'};
}


=head2 job_submission_command

  Arg [1]   : RunTest
  Arg [2]   : string, logic_name of analysis to be run
  Arg [3]   : int, toggle for script verbosity
  Function  : constructs a job_submission_command with the standard 
  options for database, logic_name and -force so it ignores the rules
  system
  Returntype: string, the job_submission command
  Exceptions: throws if the job_submission command isnt executable 
  Example   : 

=cut



sub job_submission_command{
  my ($self, $logic_name, $verbose) = @_;

  my $db_conf = $self->testdb->conf_hash;
  my $dbport = $db_conf->{'port'};
  my $dbhost = $db_conf->{'host'};
  my $dbpass = $db_conf->{'pass'};
  my $dbuser = $db_conf->{'user'};
  my $dbname = $db_conf->{'dbname'};

  my $job_submission = "../scripts/job_submission.pl";
  if(! -e $job_submission){
    $self->execption("Can't run $job_submission if it doesn't exist");
  }
  my $db_args = $self->database_args($self->testdb);
  my $cmd = "perl ".$job_submission." ";
  $cmd .= $db_args." ";
  $cmd .= "-logic_name $logic_name ";
  $cmd .= "-force";
  $cmd .= "-verbose" if($verbose);
  return $cmd;
}

=head2 rulemanager_command

  Arg [1]   : RunTest
  Arg [2]   : int, toggle for script verbostiy
  Function  : constructs a rulemanager command with standard database
  options
  Returntype: string, the rulemanager command
  Exceptions: throws if the rulemanager command isnt executable 
  Example   : 

=cut

sub rulemanager_command{
  my ($self, $verbose) = @_;

  my $db_conf = $self->testdb->conf_hash;
  my $dbport = $db_conf->{'port'};
  my $dbhost = $db_conf->{'host'};
  my $dbpass = $db_conf->{'pass'};
  my $dbuser = $db_conf->{'user'};
  my $dbname = $db_conf->{'dbname'};

  my $job_submission = "../scripts/rulemanager.pl";
  if(! -e $job_submission){
    $self->exception("Can't run $job_submission if it doesn't exist");
  }
  my $db_args = $self->database_args($self->testdb);
  my $cmd = "perl ".$job_submission." ";
  $cmd .= $db_args." ";
  $cmd .= "-verbose" if($verbose);
  return $cmd;
}

=head2 cleanup

  Arg [1]   : RunTest
  Arg [2]   : int, toggle as whether to delete the output directory or
  not
  Function  : calling the cleanup methods on the other objects to return
  the status quo and delete the output directoy if required
  Returntype: none
  Exceptions: none
  Example   : 

=cut


sub cleanup{
  my ($self, $testdb) = @_;
  if($self->cleanup_dir){
    print "Deleting directory tree ".$self->output_dir."\n" 
      if($self->verbosity);
    rmtree($self->output_dir);
  }
  if(!$testdb){
    $testdb = $self->testdb;
  }
  $testdb->cleanup if($testdb);
  $self->environment->return_environment if($self->environment);
  return;
}



=head2 analysis_stats

  Arg [1]   : RunTest
  Arg [2]   : string, logic_name
  Arg [3]   : string, table_name
  Function  : prints out statistics about the run, how much should of been
  run, what has run sucessfully, what is left in the job table and how
  many results there are in total
  Returntype: none
  Exceptions: none
  Example   : 

=cut


sub analysis_stats{
  my ($self, $logic_name, $table) = @_;
  my $db = $self->testdb->db;
  my $aa = $db->get_AnalysisAdaptor;
  my $ra = $db->get_RuleAdaptor;
  my $sic = $db->get_StateInfoContainer;
  my $ja = $db->get_JobAdaptor;
  my $analysis = $aa->fetch_by_logic_name($logic_name);
  my $rule = $ra->fetch_by_goal($analysis);
  my $conditions = $rule->list_conditions;
  my $input_id_type;
 COND:foreach my $condition(@$conditions){
    my $condition_analysis =  $aa->fetch_by_logic_name($condition);
    $input_id_type = $condition_analysis->input_id_type;
    if($input_id_type ne 'ACCUMULATOR'){
      last COND;
    }
  }
  my $total_input_ids = $sic->list_input_ids_by_type($input_id_type);
  my $analysis_input_ids = $sic->list_input_ids_by_analysis
    ($analysis->dbID);
  my $results_count = $self->count_rows_by_analysis($db, $table, 
                                                    $analysis->dbID);
  print "\n\nThere were ".@$total_input_ids." input ids to analyse\n";
  print @$analysis_input_ids." input_ids where analysed sucessfully\n";
  print "This produced ".$results_count." results\n\n";
  $self->job_details($ja);
}



=head2 whole_pipeline_stats

  Arg [1]   : RunTest
  Function  : produce stats about the pipeline run, This is done using
  Bio::EnsEMBL::Pipeline::Monitor
  Returntype: none
  Exceptions: none
  Example   : 

=cut


sub whole_pipeline_stats{
  my ($self) = @_;
  my $db = $self->testdb->db;
  my $sic = $db->get_StateInfoContainer;
  my $input_ids = $sic->get_all_input_id_analysis_sets;
  print "\nTotal input_ids by type\n";
  foreach my $type(keys(%$input_ids)){
    my @input_ids = keys(%{$input_ids->{$type}});
    print $type." ".@input_ids."\n";
  }
  print "\n";
  my $monitor = new Bio::EnsEMBL::Pipeline::Monitor(-dbobj => $db);
  $monitor->show_current_status_summary;
  $monitor->show_finished_summary(1);
}

=head2 job_details

  Arg [1]   : RunTest
  Arg [2]   : Bio::EnsEMBL::Pipeline::DBSQL::JobAdaptor
  Function  : prints out information about the jobs in the job table
  Returntype:
  Exceptions: 
  Example   : 

=cut


sub job_details{
  my ($self, $job_adaptor) = @_;
  my @jobs = $job_adaptor->fetch_all;
  my %status_count;
  foreach my $j(@jobs){
    if(!$status_count{$j->current_status->status}){
      $status_count{$j->current_status->status} = 1;
    }else{
      $status_count{$j->current_status->status}++;
    }
  }
  if(keys(%status_count)){
    print "Unsucessful job details\n";
    foreach my $status(keys(%status_count)){
      print "Job status ".$status." ".$status_count{$status}." jobs\n";
    }
    print "\n";
  }
}


=head2 count_rows_by_analysis

  Arg [1]   : RunTest
  Arg [2]   : DBI
  Arg [3]   : string, table_name
  Arg [4]   : int, analysis id
  Function  : provides a count of results in the specified table
  with the provided analysis id
  Returntype: int, count
  Exceptions: none
  Example   : 

=cut



sub count_rows_by_analysis {
  my $self = shift;
  my $db = shift;
  my $tablename = shift;
  my $analysis_id = shift;

  my $sth = $db->prepare( "select count(*) from $tablename where ".
                          "analysis_id = $analysis_id" );
  $sth->execute();
  my ( $count ) = $sth->fetchrow_array();
  return $count;
}


=head2 setup_databases

  Arg [1]   : RunTest
  Arg [2]   : TestDB, if not provided takes the testdb the object holds
  Arg [3]   : arrayref, list of tables or table groups, if not provided
  tables the are taken from the object container
  Function  : To insert to the listed tables into the given database
  the mehtod first checks to see if the name if a table_group as TestDB
  has a method for inserting all the tables in the group. If TestDB doesnt
  have a method which fits the load_tablename_tables the tablename is 
  assumed to be an actual table an put on a list to be passed to the
  TestDB::load_tables method
  Returntype: none
  Exceptions: none
  Example   : 

=cut



sub setup_database{
  my ($self, $testdb, $tables) = @_;
  if(!$testdb){
    $testdb = $self->testdb;
  }
  if(!$tables){
    $tables = $self->tables;
  }
  my @unloaded;
  foreach my $table(@{$tables}){
    my $method = "load_".$table."_tables";
    if($testdb->can($method)){
      $testdb->$method;
    }else{
      push(@unloaded, $table);
    }
  }
  if(@unloaded >= 1){
    my $data_dir = $self->testdb->curr_dir."/".$self->testdb->species;
    $testdb->load_tables(\@unloaded, $data_dir);
  }
}



=head2 check_output_dir

  Arg [1]   : RunTest
  Function  : checks if directory specifed in output_dir exists. If it
  does it returns 0 so the dir wont be deleted otherwise the directory
  is created
  Returntype: int
  Exceptions: throws if directory creation fail
  Example   : 

=cut


sub check_output_dir{
  my ($self) = @_;
  if(-d $self->output_dir){
    return 0;
  }else{
    eval{
      mkdir($self->output_dir);
    };
    if($@){
      $self->exception("Failed to create ".$self->output_dir." $@");
    }
    $self->cleanup_dir(1);
    return 1;
  }
}


=head2 run_single_analysis

  Arg [1]   : RunTest
  Arg [2]   : string, logic_name
  Arg [3]   : string, table the analysis results will be stored in
  Function  : sets up the environment, runs job_submission, checks
  if jobs are running. Once the system contains no more jobs
  a set of stats are produced about the analysis run using the 
  analysis_stats method. If a comparison database conf if passed in the
  results would also be compared to the data in the reference database
  Returntype: none
  Exceptions: throws if fails run the job_submission command
  Example   : $runtest->run_single_analysis('cpg', 'simple_feature');

=cut



sub run_single_analysis{
  my ($self, $logic_name, $table_to_fill, $verbose) = @_;

  $self->check_output_dir;
  $self->environment->add_to_perl5lib($self->extra_perl);
  $self->environment->setup_paths($logic_name, $self->testdb->species);
  $self->environment->change_blastdb($self->blastdb);
  $self->setup_database;

  my $cmd = $self->job_submission_command($logic_name, $verbose);

  print $cmd."\n" if($self->verbosity || $verbose);

  system($cmd) == 0 or $self->exception("Failed to run ".$cmd);

  my $run = 1;

 RUNNING: while ($run == 1) {
    my $jobs = $self->queue_manager->job_stats;
    if (keys(%$jobs) == 0) {
      $self->analysis_stats($logic_name, $table_to_fill);
      $run = 0;
    } else {
      sleep($self->testdb->conf_hash->{'job_stats_sleep'});
    }
  }

 if($self->comparison_conf){
    my $ref_testdb = TestDB->new(
                                 -SPECIES => $self->testdb->species, 
                                 -VERBOSE => $self->verbosity,
                                 -CONF_FILE => $self->comparison_conf,
                                );
    my $tables_to_fill = $self->tables;
    push(@$tables_to_fill, $table_to_fill);
    $self->setup_database($ref_testdb, $tables_to_fill);
    $self->ref_testdb($ref_testdb);
    my $method = "compare_".$table_to_fill;
    if($self->can($method)){
      $self->$method($ref_testdb, $logic_name);
    }else{
      print "No comparison can be made as ".$method." doesn't exist\n";
    }
    $ref_testdb->cleanup unless($self->dont_cleanup_tests);
    if($self->dont_cleanup_tests){
      $self->cleanup_command($ref_testdb);
    }
  }
  $self->cleanup() unless($self->dont_cleanup_tests);
  if($self->dont_cleanup_tests){
    $self->cleanup_command;
    $self->environment->return_environment;
  }
}


=head2 run_pipeline

  Arg [1]   : RunTest
  Arg [2]   : int, toggle for script verbosity
  Function  : sets up environment, sets up databases creates rulemanager
  commandline then runs rulemanager until sucessful completion
  Returntype: none
  Exceptions: none
  Example   : 

=cut


sub run_pipeline{
  my ($self, $verbose) = @_;
  my $cleanup_dir = $self->check_output_dir;
  $self->environment->add_to_perl5lib($self->extra_perl);
  $self->environment->change_blastdb($self->blastdb);
  $self->environment->setup_all_paths($self->testdb->species);
  $self->setup_database;
  my $cmd = $self->rulemanager_command();
  print $cmd."\n" if($self->verbosity || $verbose);
  system($cmd) == 0 or $self->exception("Failed to run ".$cmd);
  $self->whole_pipeline_stats();
  if($self->comparison_conf){
    my @comparison_tables = ('repeat_feature',  'prediction_transcript', 
                              'marker_feature', 'dna_align_feature', 
                              'protein_align_feature', 'simple_feature');
    my $ref_testdb = TestDB->new(
                                 -SPECIES => $self->testdb->species, 
                                 -VERBOSE => $self->verbosity,
                                 -CONF_FILE => $self->comparison_conf,
                                );
    my $tables_to_fill = $self->tables;
    push(@$tables_to_fill, @comparison_tables);
    $self->setup_database($ref_testdb, $tables_to_fill);
    $self->ref_testdb($ref_testdb);
    $self->feature_comparison->pipeline_compare
      (\@comparison_tables, $self->testdb->db, $ref_testdb->db);
    $ref_testdb->cleanup unless($self->dont_cleanup_tests);
    if($self->dont_cleanup_tests){
      $self->cleanup_command($ref_testdb);
    }
  }
  $self->cleanup($cleanup_dir, $self->testdb) 
    unless($self->dont_cleanup_tests);
  if($self->dont_cleanup_tests){
    $self->cleanup_command;
    $self->environment->return_environment;
  }
}


=head2 cleanup_command

  Arg [1]   : RunTest
  Arg [2]   : TestDB
  Function  : creates a command and prints it to screen to cleanup after
  a test run. This is to allow easy removal of test databases and output
  directories if you wish to keep the output for a little while for 
  investigation
  Returntype: none
  Exceptions: 
  Example   : 

=cut


sub cleanup_command{
  my ($self, $testdb) = @_;
  $testdb = $self->testdb unless($testdb);
  my $db_args = $self->database_args($testdb);
  my $data_dir = $testdb->curr_dir."/".$testdb->species;
  my $cleanup_command = "cleanup_output.pl ";
  $cleanup_command .= $db_args." ";
  $cleanup_command .= " -output_dir ".$self->output_dir;
    $cleanup_command .= " -sql_data_dir ".$data_dir;
  print "You have specifed -dont_cleanup when running your test \n".
    "If you want to delete your output you can run this script ".
      "ensembl-pipeline/test_system/cleanup_output.pl\n".
        "this is the command you should use \n".$cleanup_command."\n".
            "If you don't want any of the data sets deleted remove either".
              " -dbname, -sql_data_dir or -output_dir options from the ".
                "commandline\n";
}





=head2 database_args

  Arg [1]   : RunTest
  Arg [2]   : TestDB
  Function  : create a string of the standard database args for
  ensembl-pipeline scripts for the script RunTest runs
  Returntype: string
  Exceptions: none
  Example   : 

=cut


sub database_args{
  my ($self, $testdb) = @_;
  $testdb = $self->testdb if(!$testdb);
  my $db_conf = $testdb->conf_hash;
  my $dbport = $db_conf->{'port'};
  my $dbhost = $db_conf->{'host'};
  my $dbpass = $db_conf->{'pass'};
  my $dbuser = $db_conf->{'user'};
  my $dbname = $db_conf->{'dbname'};
  my $db_args = " -dbhost ".$dbhost." -dbuser ".$dbuser;
  $db_args .= " -dbpass ".$dbpass if($dbpass);
  $db_args .= " -dbport ".$dbport if($dbport);
  $db_args .= " -dbname ".$dbname." ";
  return $db_args;
}



=head2 compare_tables

  Arg [1]   : RunTest
  Arg [2]   : TestDB for reference database
  Arg [3]   : logic_name for analysis to compare
  Function  : These are a series of methods used for comparing results
  between the test database and a reference database using the 
  FeatureComparison module
  Returntype: none
  Exceptions: throws if it has no input_ids on which to base its fetching
  Example   : 

=cut




sub compare_simple_feature{
  my ($self, $ref_db, $logic_name) = @_;
  my ($analysis, $ref_ana, $input_ids) = $self->
    get_analyses_and_input_ids($ref_db, $logic_name);
  my $test_id = $input_ids->[0];
  my $method;
  my $query_fs;
  my $target_fs;
  if(!$test_id){
    $self->exeception("Something is wrong analysis ".$logic_name.
                      " of type ".$analysis->input_id_type.
                      " has produced no input_ids");
  }else{
    my @array = split(/:/,$test_id);
    if(scalar(@array) < 3 || scalar(@array) > 6) {
      $query_fs = $self
        ->fetch_features_by_dbID($self->testdb->db,
                                 "get_SimpleFeatureAdaptor",
                                 $logic_name);
      $target_fs = $self
        ->fetch_features_by_dbID($ref_db->db,
                                 "get_SimpleFeatureAdaptor",
                                 $logic_name);
    }else{
      $query_fs = $self
        ->fetch_features_by_slice_name($input_ids, $self->testdb->db,
                                      "get_SimpleFeatureAdaptor",
                                      $logic_name);
      $target_fs = $self
        ->fetch_features_by_slice_name($input_ids, $ref_db->db,
                                       "get_SimpleFeatureAdaptor",
                                       $logic_name);
    }
    my $feature_comparison = $self->feature_comparison;
    $feature_comparison->query($query_fs);
    $feature_comparison->target($target_fs);
    $feature_comparison->compare;
  }
}


sub compare_pmatch_feature{
  my ($self, $ref_db, $logic_name) = @_;
  my $data_dir = $self->testdb->curr_dir."/".$self->testdb->species;
  $ref_db->load_tables(['protein'], $data_dir);
  my $test_pfa = $self->testdb->db->get_PmatchFeatureAdaptor;
  my @query_fs = $test_pfa->fetch_by_logic_name($logic_name);
  my $ref_pfa = $ref_db->db->get_PmatchFeatureAdaptor;
  my @target_fs = $ref_pfa->fetch_by_logic_name($logic_name);
  my $feature_comparison = $self->feature_comparison;
  $feature_comparison->query(\@query_fs);
  $feature_comparison->target(\@target_fs);
  $feature_comparison->pmatch_compare;
}

sub compare_prediction_transcript{
  my ($self, $ref_db, $logic_name) = @_;
  my ($analysis, $ref_ana, $input_ids) = $self->
    get_analyses_and_input_ids($ref_db, $logic_name);
  my $data_dir = $self->testdb->curr_dir."/".$self->testdb->species;
  $ref_db->load_tables(['prediction_exon'], $data_dir);
  my $test_id = $input_ids->[0];
  my $method;
  my $query_fs;
  my $target_fs;
  if(!$test_id){
    $self->exception("Something is wrong analysis ".$logic_name.
                     " of type ".$analysis->input_id_type.
                     " has produced no input_ids");
  }else{
    my @array = split(/:/,$test_id);
    if(scalar(@array) < 3 || scalar(@array) > 6) {
      $query_fs = $self
        ->fetch_features_by_dbID($self->testdb->db,
                                 "get_PredictionTranscriptAdaptor",
                                 $logic_name);
      $target_fs = $self
        ->fetch_features_by_dbID($ref_db->db,
                                 "get_PredictionTranscriptAdaptor",
                                 $logic_name);
    }else{
      $query_fs = $self
        ->fetch_features_by_slice_name($input_ids, $self->testdb->db,
                                      "get_PredictionTranscriptAdaptor",
                                      $logic_name);
      $target_fs = $self
        ->fetch_features_by_slice_name($input_ids, $ref_db->db,
                                       "get_PredictionTranscriptAdaptor",
                                       $logic_name);
    }
    my $feature_comparison = $self->feature_comparison;
    $feature_comparison->query($query_fs);
    $feature_comparison->target($target_fs);
    $feature_comparison->compare;
  }
}


sub compare_repeat_feature{
  my ($self, $ref_db, $logic_name) = @_;
  my ($analysis, $ref_ana, $input_ids) = $self->
    get_analyses_and_input_ids($ref_db, $logic_name);
  my $data_dir = $self->testdb->curr_dir."/".$self->testdb->species;
  $ref_db->load_tables(['repeat_consensus'], $data_dir);
  my $test_id = $input_ids->[0];
  my $method;
  my $query_fs;
  my $target_fs;
  if(!$test_id){
    $self->exception("Something is wrong analysis ".$logic_name.
                     " of type ".$analysis->input_id_type.
                     " has produced no input_ids");
  }else{
    my @array = split(/:/,$test_id);
    if(scalar(@array) < 3 || scalar(@array) > 6) {
      $query_fs = $self
        ->fetch_features_by_dbID($self->testdb->db,
                                 "get_RepeatFeatureAdaptor",
                                 $logic_name);
      $target_fs = $self
        ->fetch_features_by_dbID($ref_db->db,
                                 "get_RepeatFeatureAdaptor",
                                 $logic_name);
    }else{
      $query_fs = $self
        ->fetch_features_by_slice_name($input_ids, $self->testdb->db,
                                      "get_RepeatFeatureAdaptor",
                                      $logic_name);
      $target_fs = $self
        ->fetch_features_by_slice_name($input_ids, $ref_db->db,
                                       "get_RepeatFeatureAdaptor",
                                       $logic_name);
    }
    my $feature_comparison = $self->feature_comparison;
    $feature_comparison->query($query_fs);
    $feature_comparison->target($target_fs);
    $feature_comparison->compare;
  }
}


sub compare_dna_align_feature{
  my ($self, $ref_db, $logic_name) = @_;
  my ($analysis, $ref_ana, $input_ids) = $self->
    get_analyses_and_input_ids($ref_db, $logic_name);
  my $data_dir = $self->testdb->curr_dir."/".$self->testdb->species;
  $ref_db->load_tables(['prediction_exon'], $data_dir);
  my $test_id = $input_ids->[0];
  my $method;
  my $query_fs;
  my $target_fs;
  if(!$test_id){
    $self->exception("Something is wrong analysis ".$logic_name.
                     " of type ".$analysis->input_id_type.
                     " has produced no input_ids");
  }else{
    my @array = split(/:/,$test_id);
    if(scalar(@array) < 3 || scalar(@array) > 6) {
      $query_fs = $self
        ->fetch_features_by_dbID($self->testdb->db,
                                 "get_DnaAlignFeatureAdaptor",
                                 $logic_name);
      $target_fs = $self
        ->fetch_features_by_dbID($ref_db->db,
                                 "get_DnaAlignFeatureAdaptor",
                                 $logic_name);
    }else{
      $query_fs = $self
        ->fetch_features_by_slice_name($input_ids, $self->testdb->db,
                                      "get_DnaAlignFeatureAdaptor",
                                      $logic_name);
      $target_fs = $self
        ->fetch_features_by_slice_name($input_ids, $ref_db->db,
                                       "get_DnaAlignFeatureAdaptor",
                                       $logic_name);
    }
    my $feature_comparison = $self->feature_comparison;
    $feature_comparison->query($query_fs);
    $feature_comparison->target($target_fs);
    $feature_comparison->compare;
  }
}

sub compare_protein_align_feature{
  my ($self, $ref_db, $logic_name) = @_;
  my ($analysis, $ref_ana, $input_ids) = $self->
    get_analyses_and_input_ids($ref_db, $logic_name);
  my $data_dir = $self->testdb->curr_dir."/".$self->testdb->species;
  $ref_db->load_tables(['prediction_exon'], $data_dir);
  my $test_id = $input_ids->[0];
  my $method;
  my $query_fs;
  my $target_fs;
  if(!$test_id){
    $self->exception("Something is wrong analysis ".$logic_name.
                     " of type ".$analysis->input_id_type.
                     " has produced no input_ids");
  }else{
    my @array = split(/:/,$test_id);
    if(scalar(@array) < 3 || scalar(@array) > 6) {
      $query_fs = $self
        ->fetch_features_by_dbID($self->testdb->db,
                                 "get_ProteinAlignFeatureAdaptor",
                                 $logic_name);
      $target_fs = $self
        ->fetch_features_by_dbID($ref_db->db,
                                 "get_ProteinAlignFeatureAdaptor",
                                 $logic_name);
    }else{
      $query_fs = $self
        ->fetch_features_by_slice_name($input_ids, $self->testdb->db,
                                      "get_ProteinAlignFeatureAdaptor",
                                      $logic_name);
      $target_fs = $self
        ->fetch_features_by_slice_name($input_ids, $ref_db->db,
                                       "get_ProteinAlignFeatureAdaptor",
                                       $logic_name);
    }
    my $feature_comparison = $self->feature_comparison;
    $feature_comparison->query($query_fs);
    $feature_comparison->target($target_fs);
    $feature_comparison->compare;
  }
}

sub compare_marker_feature{
  my ($self, $ref_db, $logic_name) = @_;
  my ($analysis, $ref_ana, $input_ids) = $self->
    get_analyses_and_input_ids($ref_db, $logic_name);
  my $data_dir = $self->testdb->curr_dir."/".$self->testdb->species;
  my $test_id = $input_ids->[0];
  my $method;
  my $query_fs;
  my $target_fs;
  if(!$test_id){
    $self->exception("Something is wrong analysis ".$logic_name.
                     " of type ".$analysis->input_id_type.
                     " has produced no input_ids");
  }else{
    my @array = split(/:/,$test_id);
    if(scalar(@array) < 3 || scalar(@array) > 6) {
      $query_fs = $self
        ->fetch_features_by_dbID($self->testdb->db,
                                 "get_MarkerFeatureAdaptor",
                                 $logic_name);
      $target_fs = $self
        ->fetch_features_by_dbID($ref_db->db,
                                 "get_MarkerFeatureAdaptor",
                                 $logic_name);
    }else{
      $query_fs = $self
        ->fetch_features_by_slice_name($input_ids, $self->testdb->db,
                                      "get_MarkerFeatureAdaptor",
                                      $logic_name);
      $target_fs = $self
        ->fetch_features_by_slice_name($input_ids, $ref_db->db,
                                       "get_MarkerFeatureAdaptor",
                                       $logic_name);
    }
    my $feature_comparison = $self->feature_comparison;
    $feature_comparison->query($query_fs);
    $feature_comparison->target($target_fs);
    $feature_comparison->compare;
  }
}



=head2 feature_comparison

  Arg [1]   : RunTest
  Arg [2]   : FeatureComparision
  Function  : a container for the FeatureComparison object it will
  create an empty FeatureComparison object if one isn't passed in
  Returntype: FeatureComparison
  Exceptions: throws if passed a second argument which isn't a 
  FeatureComparison
  Example   : 

=cut



sub feature_comparison{
  my ($self, $feature_comparison) = @_;
  if($feature_comparison){
    $self->exception("Must pass in a FeatureComparison not a "
                     .$feature_comparison) 
      unless($feature_comparison->isa('FeatureComparison'));
    $self->{'feature_comparison'} = $feature_comparison;
  }
  if(!$self->{'feature_comparison'}){
    $self->{'feature_comparison'} = FeatureComparison->new
      (
       -VERBOSE=>$self->verbosity,
      );
  }
  return $self->{'feature_comparison'};
}


=head2 get_analyses_and_input_ids

  Arg [1]   : RunTest
  Arg [2]   : TestDB for reference database
  Arg [3]   : string, logic name for analysis to fetch
  Function  : fetches the analysis object from both the reference database
  and the test database and fetched the input_ids from the test database
  Returntype: Bio::EnsEMBL::Pipeline::Analysis, 
  Bio::EnsEMBL::Pipeline::Analysis, hashref
  Exceptions: none
  Example   : 

=cut



sub get_analyses_and_input_ids{
  my ($self, $ref_db, $logic_name) = @_;
  my $analysis = $self->testdb->db->get_AnalysisAdaptor
    ->fetch_by_logic_name($logic_name);
  my $ref_ana = $ref_db->db->get_AnalysisAdaptor
    ->fetch_by_logic_name($logic_name);
  my $input_ids = $self->testdb->db->get_StateInfoContainer
    ->list_input_ids_by_type($analysis->input_id_type);
  return($analysis, $ref_ana, $input_ids);
}



=head2 fetch_features_by_slice

  Arg [1]   : RunTest
  Arg [2]   : arrayref of slice names
  Arg [3]   : dbadaptor
  Arg [4]   : string method to fetch the appropriate adaptor
  Arg [5]   : logic_name of the analysis results desired
  Function  : fetches features on a slice by slice basis
  using the adaptors fetch_all_by_Slice method
  Returntype: arrayref of Bio::EnsEMBL::Features
  Exceptions: 
  Example   : 

=cut



sub fetch_features_by_slice_name{
  my ($self, $names, $db, $fetch_adaptor_method, $logic_name) = @_;
  my @features;
  my $sa = $db->get_SliceAdaptor;
  my $fa = $db->$fetch_adaptor_method;
  foreach my $name(@$names){
    my $slice = $sa->fetch_by_name($name);
    my $features = $fa->fetch_all_by_Slice($slice, $logic_name);
    push(@features, @$features);
  }
  return \@features;
}



=head2 fetch_features_by_dbID

  Arg [1]   : RunTest
  Arg [2]   : dbadaptor
  Arg [3]   : string method to fetch the appropriate adaptor
  Arg [4]   : logic_name of the analysis results desired
  Function  : when the input_ids for a particular input_id_type
  doesnt match the standard slice name then all the features are
  fetched out of the database on the basis of dbID and then the features
  of appropriate analysis are filtered out
  Returntype: arrayref of Bio::EnsEMBL::Features
  Exceptions: 
  Example   : 

=cut



sub fetch_features_by_dbID{
  my ($self, $db, $fetch_adaptor_method, $logic_name) = @_;
  my @features;
  my $sa = $db->get_SliceAdaptor;
  my $fa = $db->$fetch_adaptor_method;
  my $dbIDs = $fa->list_dbIDs;
  if(@$dbIDs <= 1000){
    @features = @{$fa->fetch_by_dbID_list($dbIDs)};
  }else{
    my $num_in_list = 1000;
    my @list_of_lists;
    while(@$dbIDs){
      my @list = splice(@$dbIDs, 0, $num_in_list);
      push(@list_of_lists, \@list);
    }
    foreach my $list(@list_of_lists){
      my @chunk = @{$fa->fetch_by_dbID_list($list)};
      push(@features, @chunk);
    }
  }
  my @analysis_features;
  foreach my $f(@features){
    if($f->analysis->logic_name eq $logic_name){
      push(@analysis_features, $f);
    }
  }
  return(@analysis_features);
}



sub before_throw{
  my ($self) = @_;
  if(!$self->dont_cleanup_tests){
    $self->cleanup;
    if($self->ref_testdb){
      $self->cleanup($self->ref_testdb);
    }
  }else{
    $self->cleanup_command;
    if($self->ref_testdb){
      $self->cleanup_command($self->ref_testdb);
    }
  }
}



sub exception{
  my ($self, $msg) = @_;
  $self->before_throw;
  throw($msg);
}

1;
