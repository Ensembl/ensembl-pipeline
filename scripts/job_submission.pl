#!/usr/local/ensembl/bin/perl 

use strict;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use Bio::EnsEMBL::Pipeline::Utils::PipelineSanityChecks;
use Bio::EnsEMBL::Pipeline::RuleManager;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Utils::InputIDFactory;
$| = 1;

my $term_sig =  0;

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
my $runner; #the runner script to use when running jobs (this will be
#over ridden by anything in config
my $output_dir; #the output_dir to use when running jobs (this won't be
#over ridden when running jobs
my $mark_awol = 1; #Flag as whether to mark jobs which have gone missing 
#from the system
my $rename_on_retry = 1; #Whether to rename jobs stdout/err which are 
#being retried
my $config_sanity = 1; #flag as to whether to check configuration sanity
my $utils_verbosity = 'WARNING'; #how verbose do you want the 
#Bio::EnsEMBL::Utils::Exceptions module to be by default it is set to
#WARNING as this gives warning and throws but not deprecates or infos
my $logic_name; #the logic_name of the analysis to be run
my $force; #force the analysis to run regardless of the rules
my $ids_to_run; #filepath to file of input_ids to run
my $perldoc;
my @command_args = @ARGV;
my $make_input_ids;
my $slice_size;
my $slice_overlap;
my $coord_system;
my $coord_system_version;
my $slice;
my $input_id_type;
my $file;
my $dir;
my $regex;
my $single;
my $translation_id;
my $name = 'genome';
my $seq_level;
my $top_level;
my $insert_analysis;
my $store_input_ids;
GetOptions(
           'dbhost=s'      => \$dbhost,
           'dbname=s'      => \$dbname,
           'dbuser=s'      => \$dbuser,
           'dbpass=s'      => \$dbpass,
           'dbport=s'      => \$dbport,
           'help!' => \$help,
           'verbose!' => \$verbose,
           'queue_manager=s' => \$queue_manager,
           'runner=s' => \$runner,
           'output_dir=s' => \$output_dir,
           'mark_awol!' => \$mark_awol,
           'rename_on_retry' => \$rename_on_retry,
           'config_sanity!' => \$config_sanity,
           'utils_verbosity=s' => \$utils_verbosity,
           'perldoc!' => \$perldoc,
           'analysis=s' => \$logic_name,
           'force!' => \$force,
           'input_id_file=s' => \$ids_to_run,
           'make_input_ids!' => \$make_input_ids,
           'coord_system:s'       => \$coord_system,
           'coord_system_version:s' => \$coord_system_version,
           'slice'        => \$slice,
           'slice_size:s' => \$slice_size,
           'slice_overlap:s' => \$slice_overlap,
           'logic_name:s' => \$logic_name,
           'input_id_type:s' => \$input_id_type,
           'file'         => \$file,
           'dir:s'        => \$dir,
           'file_regex:s' => \$regex,
           'single'       => \$single,
           'single_name:s'=> \$name,
           'translation_ids' => \$translation_id,
           'seq_level' => \$seq_level,
           'top_level' => \$top_level,
           'insert_analysis!' => \$insert_analysis,
           'store_input_ids!' => \$store_input_ids,
           ) or useage(\@command_args);





perldoc() if $perldoc;
verbose($utils_verbosity);
unless ($dbhost && $dbname && $dbuser) {
    print STDERR "Must specify database with -dbhost, -dbname, -dbuser and -dbpass\n";
    print STDERR "Currently have -dbhost $dbhost -dbname $dbname ".
      "-dbuser $dbuser -dbpass $dbpass -dbport $dbport\n";
    $help = 1;
}
if(!$logic_name){
  print "Can't run without an analysis to run on \n";
  print "specific analysis with -analysis logic_name";
  $help = 1;
}
if($help){
  useage(\@command_args);
}

if(!defined($mark_awol)){
  $mark_awol = 0;
}
if($make_input_ids && !$force){
  print STDERR "Setting force to on as if you are making the input ids ".
    "the odds are any rule won't work/n";
  $force = 1;
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
   -VERBOSE => $verbose,
   -RUNNER => $runner,
   -OUTPUT_DIR => $output_dir,
   );
   
if($config_sanity){
  $sanity->config_sanity_check;
}

my $analysis = $rulemanager->analysis_adaptor->
  fetch_by_logic_name($logic_name);

if(!$analysis || !$analysis->input_id_type){
  throw("Must have an analysis object $logic_name $analysis and analysis ".
        "must have an input_id_type\n");
}
if($analysis->input_id_type eq 'ACCUMULATOR'){
  throw("Can't use this script to run accumulators");
}
if ($ids_to_run && ! -e $ids_to_run) {
  throw("Must be able to read $ids_to_run");
}
my $input_ids = setup_input_ids($analysis, $rulemanager, $ids_to_run, 
                                $make_input_ids, $slice, $file, 
                                $translation_id, $single, $slice_size, 
                                $slice_overlap, $dir, $regex, $name, 
                                $seq_level, $top_level, $verbose, 
                                $logic_name, $input_id_type, 
                                $insert_analysis, $store_input_ids);

my @input_ids = keys(%{$input_ids->{$analysis->input_id_type}});

my $rule = $rulemanager->rule_adaptor->fetch_by_goal($analysis);
my %completed_accumulator_analyses = %{$rulemanager->fetch_complete_accumulators};
if(@input_ids == 0){
  throw("Can't do anything I have no input ids");
}
if($force){
  warning("You are forcing this job to be run without checking if ".
          "the rules allow it or if it has already run/is running ".
          "are you sure you want to do this\n");
}
print STDERR "Trying to submit jobs for ".$analysis->logic_name."\n" 
  if($verbose);
INPUT_ID:foreach my $input_id(@input_ids){
  print $input_id."\n" if($verbose);
  if ($term_sig) {
    print "Got term signal\n" if($verbose);
    last INPUT_ID;
  }
  if($force){
   my $job = $rulemanager->create_and_store_job($input_id, $analysis);
   $job->batch_runRemote;
  }else{
    my @anals = @{$rulemanager->stateinfocontainer->
                    fetch_analysis_by_input_id($input_id)};
    my $anal = $rule->check_for_analysis
      (\@anals, $analysis->input_id_type, 
       \%completed_accumulator_analyses, $verbose);
    if(UNIVERSAL::isa($anal,'Bio::EnsEMBL::Pipeline::Analysis')){
      $rulemanager->can_job_run($input_id, $analysis);
    }
  }
}
$rulemanager->cleanup_waiting_jobs();
$rulemanager->db->pipeline_unlock;
sub setup_input_ids{
  my ($analysis, $rulemanager, $ids_to_run, $make_input_ids,
      $slice, $file, $translation_id, $single, $slice_size, 
      $slice_overlap, $dir, $regex, $name, $seq_level, $top_level, 
      $verbose, $logic_name, $input_id_type, $insert_analysis,
      $store_input_ids) = @_;
  
  my $id_hash;
  if($ids_to_run){
    $id_hash = $rulemanager->read_id_file($ids_to_run);
    my @types = keys(%$id_hash);
    if(scalar(@types) != 1){
      throw("You have passed in a file with ".@types." input ".
            "id types something funny is going on");
    }
    if($types[0] ne $analysis->input_id_type){
      throw("If your input_ids aren't the same tupe are your analysis ".
            $types[0]." compared to ".$analysis->input_id_type." this ".
            "won't work");
    }
    return $id_hash;
  }elsif($make_input_ids){
    print STDERR "Making input ids\n";
    $id_hash = make_input_ids($slice, $file, $translation_id, $single, 
                          $slice_size, $slice_overlap, $dir, 
                          $regex, $name, $seq_level, $top_level, 
                          $verbose, $logic_name, $input_id_type, 
                          $insert_analysis, $store_input_ids);
  }else{
    $id_hash->{$analysis->input_id_type} = {};
    my @ids = @{$rulemanager->stateinfocontainer
                  ->list_input_ids_by_type($analysis->input_id_type)};
    foreach my $id(@ids){
      $id_hash->{$analysis->input_id_type}->{$id} = 1;
    }
    $rulemanager->input_ids($id_hash);
    return $id_hash;
  }
}


sub termhandler {
    $term_sig = 1;
}


sub make_input_ids{
  my ($slice, $file, $translation_id, $single, $slice_size, 
      $slice_overlap, $dir, $regex, $name, $seq_level, $top_level, 
      $verbose, $logic_name, $input_id_type, $insert_analysis,
      $store_input_ids) = @_;
  my $inputIDFactory = new Bio::EnsEMBL::Pipeline::Utils::InputIDFactory
    (
     -db => $db,
     -slice => $slice,
     -single => $single,
     -file => $file,
     -translation_id => $translation_id,
     -seq_level => $seq_level,
     -top_level => $top_level,
     -dir => $dir,
     -regex => $regex,
     -single_name => $name,
     -verbose => $verbose,
     -logic_name => $logic_name,
     -input_id_type => $input_id_type,
     -insert_analysis => $insert_analysis,
     -coord_system => $coord_system,
     -coord_system_version => $coord_system_version,
     -slice_size => $slice_size,
     -slice_overlaps => $slice_overlap,
    );



  my $ids = $inputIDFactory->generate_input_ids;
  if($store_input_ids){
    $inputIDFactory->store_input_ids;
  }
  return $inputIDFactory->get_id_hash;
}
