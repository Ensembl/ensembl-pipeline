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
           'logic_name=s' => \$logic_name,
           'force!' => \$force,
           'input_id_file=s' => \$ids_to_run,
           'make_input_ids!' => \$make_input_ids,
           'coord_system:s'       => \$coord_system,
           'coord_system_version:s' => \$coord_system_version,
           'slice'        => \$slice,
           'slice_size:s' => \$slice_size,
           'slice_overlap:s' => \$slice_overlap,
           'logic_name:s' => \$logic_name,
           'file'         => \$file,
           'dir:s'        => \$dir,
           'file_regex:s' => \$regex,
           'single'       => \$single,
           'single_name:s'=> \$name,
           'translation_ids' => \$translation_id,
           'seq_level' => \$seq_level,
           'top_level' => \$top_level,
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
                                $logic_name);

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
      $verbose, $logic_name) = @_;
  
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
                          $verbose, $logic_name, $input_id_type);
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
      $verbose, $logic_name, $store_input_ids) = @_;
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
     -coord_system => $coord_system,
     -coord_system_version => $coord_system_version,
     -slice_size => $slice_size,
     -slice_overlaps => $slice_overlap,
    );



  my $ids = $inputIDFactory->generate_input_ids;

  return $inputIDFactory->get_id_hash;
}

=pod 

=head1 NAME

job_submission.pl

=head1 SYNOPSIS

job_submission.pl a script for submitting a single analysis' jobs
to the pipeline.

=head1 DESCRIPTION

this script will run a single analyses jobs through the pipeline
it will check rules and retry failed jobs but it doesn't have too.
It can either read the input_ids from a file or create the
input_ids using the InputIDFactory

=head OPTIONS
 
DB Connection Details


   -dbhost     The host where the pipeline database is.
   -dbport     The port.
   -dbuser     The user to connect as.
   -dbpass     The password to use.
   -dbname     The database name.

Analysis Details

  -logic_name the logic_name of the analysis you want to run. You 
  must have this analysis already in the analysis table and it
  must has an input_id_type and module specified. This logic
  name is also used to fetch the rule. This analysis should be
  the goal of the rule you want executed
  -force this forces the script to ignore the rule and just create
  and submit the jobs with no regard to whether they should be run
  or are already running. If this option is specifed jobs already 
  in the pipeline running this analysis are ignored and failed
  jobs of this type are not retried
  

RuleManager details

  -utils_verbosity, this affects the amount of chatter you recieve from
  the core module Bio::EnsEMBL::Utils::Exception. By default this is set
  to WARNING which means you see the prints from warnings and throws but
  not from deprecate and info calls. See the modules itself for more 
  information
  -verbose, toggles whether some print statements are printed in this
  script and in the RuleManager
  -config_sanity this is a test to check certain values are defined in
  your General.pm and BatchQueue.pm config files

  Some of the follow options can intially be set in either the General.pm
  of BatchQueue.pm config files see docs of  rulemanager.pl for more 
  details
  
  -queue_manager this specifies which 
  Bio::EnsEMBL::Pipeline::BatchSubmission module is used by 
  Bio::EnsEMBL::Pipeline::RuleManager
  -output_dir the path to an output directory when the jobs stderr and
  stdout will be redirected. This alwaysoverrides the values specified in
  BatchQueue
  -mark_awol toggle to specify whether to mark jobs as awol if lost from
  the submission system this can apply strain to the job table. It is on
  by default and can be switched off using the -nomark_awol flag
  -runner path to a default runner script. This will override what is
  set in General.pm but will be overidden by any analyses specific
  settings found in BatchQueue.pm
  -rename_on_retry a toggle to specify whether to rename stderr/stdout 
  file when a job is retried as otherwise the submission system just
  cats them all together


Input id details

  If you specify no options to do with input_ids it will just take 
  all the input_ids from the input_id_analysis table with an appropriate
  input_id_type as specified by the analysis object

  -input_id_file this is a text file in the format input_id input_id_type
  if used these are the only input_ids which are considered
  -make_input_ids this indicates you want to use the InputIDFactory
   to make the input_ids for you. If you specify this option it
   makes you use force too so the rules are ignored as if the
   analysis needs its input ids created it won't pass a rule check'
 
  These are options needed for the manufacture of input ids

  -slice  signals to insert slice type input ids using
  the format 
  coord_system:coord_system_version:seq_region_name:start:end:strand
  -coord_system the coordinate system you want slices in
  -coord_system_version the version of the coord system you want
  -slice_size the size to make the slice ids
  -slice_overlap the slice overlap (non-overlapping by default)
  -file     if the input_ids are to be a list of filenames from a directory
  -dir      the directory to read the filenames from
  -file_regex a regex to impose on the filenames before using them
  -single if you just want a single dummy input_id ie for genome wide 
  analyses
  -single_name , by default this is genome but you can specify something
  different here, for example for protein annotation jobs which use the 
  whole proteome this must be proteome
  -translation_ids if you want your input ids to be translation ids
  -verbose if you want more information about what the script is doing
  -input_id_type if you want to specific an input_id_type not already
  used by the analysis object
  -insert_analysis if you want to insert an analysis object if it doesn't
     already exist in the database'
  -seq_level if you want the ids for the seq_level seq_regions, can
  work with slice_size but the -slice options isn't required'
  -top_level this will fetch all the non_redundant pieces in
  the database this may produce ids which are a mixture of different
  coordinate systems, if -coord_system_version is specified it will
  be ignored
 


Misc options

 -help will print out the standard help
 -perldoc will print out these perl docs

=head1 CONTACT

Post general queries to <ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut
