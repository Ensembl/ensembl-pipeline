#!/usr/local/ensembl/bin/perl -w

use strict;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use Bio::EnsEMBL::Pipeline::Utils::PipelineSanityChecks;
use Bio::EnsEMBL::Pipeline::RuleManager;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;

$| = 1;

my $term_sig =  0;
my $rst_sig  =  0;

$SIG{USR1} = \&sighandler;
$SIG{TERM} = \&termhandler;
$SIG{INT} = \&termhandler;

#database arguments
my $dbhost    = $ENV{'ENS_DBHOST'};
my $dbname    = $ENV{'ENS_DBNAME'};
my $dbuser    = $ENV{'ENS_DBUSER'};
my $dbpass    = $ENV{'ENS_DBPASS'};
my $dbport    = $ENV{'ENS_DBPORT'} || 3306;

#command line args
my $help; #get docs about script
my $verbose; #print statements about the running script
my $queue_manager; #what Bio::EnsEMBL::BatchSubmission module to use
#The rulemanger sleeps when a predefined number of jobs has been reached
my $max_job_time; # the longest time the rulemanger should sleep for
my $min_job_time; # the shortest time the rulemanger sleeps for
my $sleep_per_job; #the length of time to sleep per job when limit is 
#reached
my $runner; #the runner script to use when running jobs (this will be
#over ridden by anything in config
my $output_dir; #the output_dir to use when running jobs (this won't be
#over ridden when running jobs
my $job_limit; #the maximum number of jobs of a perdefined status that are
#in the system. The predefined statuses are set in the BatchQueue.pm config
#module and currently refer to LSF
my $mark_awol = 1; #Flag as whether to mark jobs which have gone missing 
#from the system
my $rename_on_retry = 1; #Whether to rename jobs stdout/err which are 
#being retried
my $once; #Only run the loop once
my @starts_from; #array of logic_names which are used to get the starting
#input_id hash
my @analyses_to_run; #array of logic_names of analyses to run
my @analyses_to_skip; #array of logic_names of analyses to skip
my @types_to_run; #array of input_id_types to consider for the run
my @types_to_skip; #array of input_id_types to not consider for the run
my $ids_to_run; #filepath to file of input_ids to run
my $ids_to_skip; #filepath to file of input_ids to skip
my $config_sanity = 1; #flag as to whether to check configuration sanity
my $accumulator_sanity = 1; #as above for accumulators
my $rules_sanity = 1;#same for rules
my $db_sanity = 1;#same for the database tables
my $kill_jobs; #whether to kill jobs which have been running for to long
my $killed_time;#how long to let jobs have before being killed
my $kill_file; #a path to a file to record what jobs have been killed
my $rerun_sleep = 300; #how long to sleep for at the end of a loop if no 
#jobs get submitted
my $utils_verbosity = 'WARNING'; #how verbose do you want the 
#Bio::EnsEMBL::Utils::Exceptions module to be by default it is set to
#WARNING as this gives warning and throws but not deprecates or infos
my $shuffle;
my $accumulators = 1;
my $reread_input_ids = 0; #toggle whether to reread input_id each time the
#script loops
my $reread_rules = 0; #toggle whether to reread rules each time the script
#loops

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
           'utils_verbosity=s' => \$utils_verbosity,
           'shuffle!' => \$shuffle,
           'accumulators!' => \$accumulators,
           'reread_input_ids!' => \$reread_input_ids,
           'reread_rules!' => \$reread_rules,
           'once!' => \$once,
           ) or throw("Can't get rulemanager commandline arguments");

verbose($utils_verbosity);
unless ($dbhost && $dbname && $dbuser) {
    print STDERR "Must specify database with -dbhost, -dbname, -dbuser and -dbpass\n";
    print STDERR "Currently have -dbhost $dbhost -dbname $dbname ".
      "-dbuser $dbuser -dbpass $dbpass -dbport $dbport\n";
    exit 1;
}
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
   

my $done;
my $reset;

if($config_sanity){
  $sanity->config_sanity_check;
}
if($db_sanity){
  $sanity->db_sanity_check;
}

my %analyses_to_run = %{$rulemanager->logic_name2dbID(\@analyses_to_run)};
my %analyses_to_skip = %{$rulemanager->
                           logic_name2dbID(\@analyses_to_skip)};



if ($ids_to_run && ! -e $ids_to_run) {
  throw("Must be able to read $ids_to_run");
}
if($ids_to_skip && ! -e $ids_to_skip){
  throw("Must be able to read $ids_to_skip");
}

my $all_rules = $rulemanager->rules;

if($accumulators && $accumulator_sanity){
  $accumulators = $sanity->accumulator_sanity_check($all_rules, 
                                                    $accumulators);
}
if($rules_sanity){
  $sanity->rule_type_sanity($all_rules, $verbose);
}

$rulemanager->add_created_jobs_back; # a method to ensure unsubmitted jobs
#are readded to the system

my %always_incomplete_accumulators;
my %accumulator_analyses;

foreach my $rule (@$all_rules) {
  if ($rule->goalAnalysis->input_id_type eq 'ACCUMULATOR') {
    $accumulator_analyses{$rule->goalAnalysis->logic_name}
      = $rule->goalAnalysis;
  }
}

setup_pipeline(\%analyses_to_run, \%analyses_to_skip, $all_rules, 
               \%accumulator_analyses, \%always_incomplete_accumulators, 
               $ids_to_run, $ids_to_skip, \@types_to_run, \@types_to_skip,
               \@starts_from);

my %completed_accumulator_analyses;
while(1){
  print "Reading IDs \n" if $verbose;
  my $submitted = 0;
  %completed_accumulator_analyses = %{$rulemanager->fetch_complete_accumulators};
  my %incomplete_accumulator_analyses = %always_incomplete_accumulators;
  if($reread_input_ids){
    $rulemanager->input_id_setup($ids_to_run, $ids_to_skip, 
                                 \@types_to_run, \@types_to_skip,
                                 \@starts_from); 
  }
  if($reread_rules){
    $rulemanager->rules_setup(\%analyses_to_run, \%analyses_to_skip, 
                              $all_rules, \%accumulator_analyses, 
                              \%always_incomplete_accumulators);
  }
  my $id_hash = $rulemanager->input_ids;
 INPUT_ID_TYPE:foreach my $type(keys(%$id_hash)){
    next INPUT_ID_TYPE if ($type eq 'ACCUMULATOR');
    my @id_list = keys(%{$id_hash->{$type}});
    @id_list = shuffle(@id_list) if $shuffle;
   INPUT_ID:foreach my $input_id(@id_list){
       if ($term_sig) {
         print "Got term signal\n" if($verbose);
         $done = 1;
         last INPUT_ID_TYPE;
       }
       
       if ($rst_sig) {
         $done = 0;
         $reset = 1;
         print "Got reset signal\n" if($verbose);
         setup_pipeline(\@analyses_to_run, \@analyses_to_skip, $all_rules, 
                        \%accumulator_analyses, 
                        \%always_incomplete_accumulators, 
                        $ids_to_run, $ids_to_skip, \@types_to_run, \@types_to_skip,
                        \@starts_from);
         last INPUT_ID_TYPE;
       }
       my %analHash;
       my @anals = @{$rulemanager->stateinfocontainer->
                       fetch_analysis_by_input_id($input_id)};
       foreach my $rule(@{$rulemanager->rules}){
         my $anal = $rule->check_for_analysis
           (\@anals, $type, \%completed_accumulator_analyses);
         if(UNIVERSAL::isa($anal,'Bio::EnsEMBL::Pipeline::Analysis')){
           $analHash{$anal->dbID} = $anal;
         }else{
           if ($rule->goalAnalysis->input_id_type eq 'ACCUMULATOR' &&
               $rule->has_condition_of_input_id_type($type) ) {
             $incomplete_accumulator_analyses
               {$rule->goalAnalysis->logic_name} = 1;
           }
         }
       }
     ANAL: for my $anal (values %analHash) {
         if($rulemanager->can_job_run($input_id, $anal)){
           $submitted++;
         }
       }
     }
   }
  if(!$done && !$reset){
    if($accumulators){
      %completed_accumulator_analyses = $rulemanager->
          fetch_complete_accumulators;
      foreach my $logic_name (keys %accumulator_analyses) {
        if (!exists($incomplete_accumulator_analyses{$logic_name}) &&
            !exists($completed_accumulator_analyses{$logic_name})) {
          if($rulemanager->can_job_run('ACCUMULATOR', 
                                       $accumulator_analyses{$logic_name}
                                      )){
            $submitted++;
          } elsif (exists($incomplete_accumulator_analyses
                         {$logic_name})) {
            print "Accumulator type analysis $logic_name ".
              "conditions unsatisfied\n" if $verbose;
          } else {
            print "Accumulator type analysis $logic_name ".
              "already run\n" if $verbose;
          }
        }
      }
    }
  }
  $rulemanager->job_stats($job_limit);
  if(!$done){
    if(!$rulemanager->check_if_done){
      $done = 1;
    }else{
      $done = 0;
    }
  }
  $rulemanager->cleanup_waiting_jobs();
  if($done || $once){
    $rulemanager->db->pipeline_unlock;
    exit 0;
  }
  sleep($rerun_sleep) if $submitted == 0;
  print "Waking up and run again!\n" if $verbose;
}


sub setup_pipeline{
  my ($analyses_to_run, $analyses_to_skip, $all_rules, 
      $accumulator_analyses, $always_incomplete_accumulators, 
      $ids_to_run, $ids_to_skip, $types_to_run, $types_to_skip,
      $starts_from) = @_;

  $rulemanager->rules_setup($analyses_to_run, $analyses_to_skip, 
                            $all_rules, $accumulator_analyses, 
                            $always_incomplete_accumulators);
  $rulemanager->input_id_setup($ids_to_run, $ids_to_skip, 
                               $types_to_run, $types_to_skip, 
                               $starts_from);
}

sub shuffle {
    my (@in) = @_;
    my @out;

    srand;

    push @out, splice(@in, rand @in, 1) while @in;

    return @out;
}


sub termhandler {
    $term_sig = 1;
}


# handler for SIGUSR1
sub sighandler {

    $rst_sig = 1;
    $SIG{SIG1} = \&sighandler;
};
