#!/usr/local/ensembl/bin/perl -w

use strict;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use Bio::EnsEMBL::Pipeline::Utils::PipelineSanityChecks;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Utils::InputIDFactory;
use Bio::EnsEMBL::Pipeline::Config::General;
use Bio::EnsEMBL::Pipeline::Config::BatchQueue;
use Bio::EnsEMBL::Pipeline::Analysis;

$| = 1;



#database arguments

my $dbhost = $ENV{'ENS_DBHOST'};
my $dbname = $ENV{'ENS_DBNAME'};
my $dbuser = $ENV{'ENS_DBUSER'};
my $dbpass = $ENV{'ENS_DBPASS'};
my $dbport = $ENV{'ENS_DBPORT'} || 3306;

my $help;            #get docs about script
my $verbose;         #print statements about the running script
my $perldoc;
my @command_args = @ARGV;
my $utils_verbosity = 'WARNING'; #how verbose do you want the 
                                 #Bio::EnsEMBL::Utils::Exceptions module to be by default it is set to
                                 #WARNING as this gives warning and throws but not deprecates or infos
my $queue_manager;   #what Bio::EnsEMBL::BatchSubmission module to use

my $script_dir = `pwd`;
chomp $script_dir;

my $output_dir;
my $logic_name;
my $config_sanity = 1;
my $queue;
my $runnabledb_path;
my $write;
my $update;
my $script_verbosity;
my $resource;
my $sub_args;
my $input_id_type;
my $module;
my $insert_analysis;
my $ids_to_run;
my $make_input_ids;
my $slice_size;
my $slice_overlap;
my $coord_system;
my $coord_system_version;
my $slice;
my $file;
my $dir;
my $regex;
my $single;
my $translation_id;
my $name = 'genome';
my $seq_level;
my $top_level;
my $pre_exec;
my $open;
my $submission_interval = 5;
my $commandline;
my $pre_exec_command;

GetOptions(
           'dbhost=s'                  => \$dbhost,
           'dbname=s'                  => \$dbname,
           'dbuser=s'                  => \$dbuser,
           'dbpass=s'                  => \$dbpass,
           'dbport=s'                  => \$dbport,
           'help!'                     => \$help,
           'verbose!'                  => \$verbose,
           'perldoc!'                  => \$perldoc,
           'utils_verbosity=s'         => \$utils_verbosity,
           'queue_manager=s'           => \$queue_manager,
           'script_dir=s'              => \$script_dir,
           'output_dir=s'              => \$output_dir,
           'queue:s'                   => \$queue,
           'logic_name=s'              => \$logic_name,
           'config_sanity!'            => \$config_sanity,
           'runnabledb_path=s'         =>\$runnabledb_path,
           'update_input_id_analysis!' => \$update,
           'pre_exec!'                 => \$pre_exec,
           'script_verbosity!'         => \$script_verbosity,
           'resource_requirements:s'   => \$resource,
           'sub_args:s'                => \$sub_args,
           'input_id_file:s'           => \$ids_to_run,
           'make_input_ids!'           => \$make_input_ids,
           'coord_system:s'            => \$coord_system,
           'coord_system_version:s'    => \$coord_system_version,
           'slice!'                    => \$slice,
           'slice_size:s'              => \$slice_size,
           'slice_overlap:s'           => \$slice_overlap,
           'file!'                     => \$file,
           'dir:s'                     => \$dir,
           'file_regex:s'              => \$regex,
           'single!'                   => \$single,
           'single_name:s'             => \$name,
           'translation_ids!'          => \$translation_id,
           'seq_level!'                => \$seq_level,
           'top_level!'                => \$top_level,
           'open!'                     => \$open,
           'submission_interval:s'     => \$submission_interval,
           'write!'                    => \$write,
           'command:s'                 => \$commandline,
           'pre_exec_command:s'        => \$pre_exec_command,
          ) or useage(\@command_args);

perldoc() if $perldoc;

verbose($utils_verbosity);

unless ($dbhost && $dbname && $dbuser) {
  print STDERR "Must specify database with -dbhost, -dbname, -dbuser and -dbpass\n";
  print STDERR "Currently have -dbhost $dbhost -dbname $dbname ".
               "-dbuser $dbuser -dbpass $dbpass -dbport $dbport\n";
  $help = 1;
}

if (!$logic_name) {
  print STDERR "You must specify a logic_name for the analysis you want to run\n";
  $help = 1;
}

if ($insert_analysis && (!$input_id_type || !$module)) {
  print STDERR "If you want to insert your analysis object into your ".
               "you must specify input_id_type and module on the commandline too\n";
  $help = 1;
}

if ($commandline && $pre_exec && !$pre_exec_command) {
  print STDERR "If you want a pre_exec statement with your custom ".
               "commandline you must specify it with the -pre_exec_command ".
               "option\n";
  $help = 1;
}

if ($help) {
  useage(\@command_args);
}

if (!$queue_manager) {
  $queue_manager = $QUEUE_MANAGER; #found in BatchQueue.pm
}

my $batch_q_module = "Bio::EnsEMBL::Pipeline::BatchSubmission::$queue_manager";

my $batch_q_file = "$batch_q_module.pm";
$batch_q_file =~ s{::}{/}g;

eval {
  require "$batch_q_file";
};
if ($@) {
  throw("Can't find $batch_q_file [$@]");
}


my $db = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new
  (
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

if ($config_sanity) {
  $sanity->config_sanity_check;
}

my %hash = &set_up_queues;

my $batch_q_key = $logic_name;

if (!$hash{$batch_q_key}) {
  $batch_q_key = 'default';
}

if (!$queue) {
  $queue = $hash{$batch_q_key}{'queue'};
}

if (!$runnabledb_path) {
  $runnabledb_path = $hash{$batch_q_key}{'runnabledb_path'};
}

if (!$resource) {
  $resource = $hash{$batch_q_key}{'resource'};
}

if (!$sub_args) {
  $sub_args = $hash{$batch_q_key}{'sub_args'};
}

if (!$output_dir) {
  $output_dir = $hash{$batch_q_key}{'output_dir'};

  if (!$output_dir) {
    $output_dir = $DEFAULT_OUTPUT_DIR;
  }
}

if (!-d $output_dir) {
  throw("Can't use output_dir ".$output_dir." as its not a directory");
}

my $analysis;
if ($insert_analysis) {
  my $analysis = Bio::EnsEMBL::Pipeline::Analysis->new
    (
     -logic_name => $logic_name,
     -module => $module,
     -input_id_type => $input_id_type,
    );
} else {
  $analysis = $db->get_AnalysisAdaptor->fetch_by_logic_name($logic_name);
}

if (!$analysis) {
  throw("Can't run without an Analysis object");
}

#print STDERR "3Have commandline ".$commandline." update ".$update."\n";
my $input_ids = setup_input_ids($analysis, $db, $ids_to_run, 
                                $make_input_ids, $slice, $file, 
                                $translation_id, $single, $slice_size, 
                                $slice_overlap, $dir, $regex, $name, 
                                $seq_level, $top_level, $verbose, 
                                $logic_name);

my @batch_submission_objects;
#print STDERR "4Have commandline ".$commandline." update ".$update."\n";

foreach my $input_id (@$input_ids) {
  my ($stdout, $stderr) = make_filenames($output_dir, $input_id, $logic_name);

  if ($pre_exec && !$pre_exec_command) {
    $pre_exec_command = $script_dir."/test_RunnableDB -check ";
  }

  my $batch_job = $batch_q_module->new
      (
       -STDOUT     => $stdout,
       -STDERR     => $stderr,
       -PARAMETERS => $sub_args,
       -PRE_EXEC   => $pre_exec_command,
       -QUEUE      => $queue,
       -JOBNAME    => $db->dbname . ':' . $logic_name,
       -RESOURCE   => $resource,
      );

  my $command;

  if (!$commandline) {
    $command = $script_dir."/test_RunnableDB ";
    $command .= "-dbhost ".$db->host." -dbuser ".$db->username.
      " -dbname ".$db->dbname." -dbport ".$db->port." ";
    $command .= " -dbpass ".$db->password." " if($dbpass);
    $command .= " -input_id ".$input_id." -logic_name ".
      $logic_name." ";
    $command .= " -write " if($write);
    $command .= " -update_input_id_analysis " if($update);
    $command .= " -runnabledb_path " if($runnabledb_path);

    # $command .= " -utils_verbosity ".$utils_verbosity." " 
    #   if($utils_verbosity);

    $command .= " -input_id_type ".$input_id_type.
      " " if($input_id_type);
    $command .= " -module ".$module." " if($module);

    #print STDERR "Have commandline ".$command."\n";

  }else{
    $command = $commandline;
    $command =~ s/INPUT_ID/$input_id/;
  }

  $batch_job->construct_command_line($command);
  push(@batch_submission_objects, $batch_job);
}

print STDERR "YOU WILL NOT BE SUBMITTING ANY JOBS BECAUSE YOU DIDN'T ".
  "SPECIFY -open \n" if(!$open);

foreach my $batch_object (@batch_submission_objects) {
  print $batch_object->bsub."\n" if($verbose);

  if ($open) {
    eval {
      $batch_object->open_command_line;
    };
    ifi ($@ || !$batch_object->id) {
      throw("Failed to open ".$batch_object->bsub." $@");
    }

    print $batch_object->bsub."\n" if($verbose);
    sleep($submission_interval);
  } else {
    #print $batch_object->bsub."\n";
  }
}

print STDERR "YOU DIDN'T BE SUBMITTING ANY JOBS BECAUSE YOU DIDN'T SPECIFY ".
  "-open \n" if(!$open);

sub setup_input_ids {
  my ($analysis, $db, $ids_to_run, $make_input_ids,
      $slice, $file, $translation_id, $single, $slice_size, 
      $slice_overlap, $dir, $regex, $name, $seq_level, $top_level, 
      $verbose, $logic_name) = @_;

  if ($ids_to_run) {
    my $id_hash = read_id_file($ids_to_run);
    my @types = keys(%$id_hash);

    if (scalar(@types) != 1){
      throw("You have passed in a file with ".@types." input ".
            "id types something funny is going on");
    }

    if ($types[0] ne $analysis->input_id_type) {
      throw("If your input_ids aren't the same tupe are your analysis ".
            $types[0]." compared to ".$analysis->input_id_type." this ".
            "won't work");
    }

    my @ids = keys(%{$id_hash->{$analysis->input_id_type}});
    return \@ids;

  } elsif ($make_input_ids) {

    print STDERR "Making input ids\n" if($verbose);

    my $ids = make_input_ids($slice, $file, $translation_id, $single, 
                          $slice_size, $slice_overlap, $dir, 
                          $regex, $name, $seq_level, $top_level, 
                          $verbose, $logic_name,  $db);
    return $ids;
  } else {
    my @ids = @{$db->get_StateInfoContainer
                  ->list_input_ids_by_type($analysis->input_id_type)};
    print STDERR "Have ".@ids." ids\n";
    return \@ids;
  }
}


sub read_id_file{
  my ($file) = @_;

  my @ids;
  open IDS, "< $file" or throw("Can't open $file");
  while (<IDS>) {
    chomp;
    my($id) = split;
    push(@ids, $id);
  }
  return \@ids
}


sub make_filenames {
  my ($output_dir, $input_id, $logic_name) = @_;
  
  my $num = int(rand(10));
  
  my $dir = $output_dir . "/$num/";
  if ( ! -e $dir) {
    system( "mkdir $dir" );
  }

  my $stub = $input_id.".";
  $stub .= $logic_name.".";
  $stub .= int(rand(1000));
 
  my $stdout_file = ($dir.$stub.".out");
  my $stderr_file = ($dir.$stub.".err");
  return $stdout_file, $stderr_file;
}


sub make_input_ids{
  my ($slice, $file, $translation_id, $single, $slice_size, 
      $slice_overlap, $dir, $regex, $name, $seq_level, $top_level, 
      $verbose, $logic_name,  $db) = @_;

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
  return $ids;
}


sub set_up_queues {
  my %q;

  foreach my $queue (@$QUEUE_CONFIG) {
    my $ln = $queue->{'logic_name'};
    next unless $ln;
    delete $queue->{'logic_name'};

    while (my($k, $v) = each %$queue) {
      $q{$ln}{$k} = $v;
      $q{$ln}{'jobs'} = [];
      $q{$ln}{'last_flushed'} = undef;
      $q{$ln}{'batch_size'}      ||= $DEFAULT_BATCH_SIZE;
      $q{$ln}{'queue'}           ||= $DEFAULT_BATCH_QUEUE;
      $q{$ln}{'retries'}         ||= $DEFAULT_RETRIES;
      $q{$ln}{'cleanup'}         ||= $DEFAULT_CLEANUP;
      $q{$ln}{'runnabledb_path'} ||= $DEFAULT_RUNNABLEDB_PATH;
    }

    # a default queue for everything else
    unless (defined $q{'default'}) {
      $q{'default'}{'batch_size'} = $DEFAULT_BATCH_SIZE;
      $q{'default'}{'retries'} = $DEFAULT_RETRIES;
      $q{'default'}{'last_flushed'} = undef;
      $q{'default'}{'queue'} = $DEFAULT_BATCH_QUEUE;
      $q{'default'}{'jobs'} = [];
      $q{'default'}{'cleanup'} = $DEFAULT_CLEANUP;
      $q{'default'}{'runnabledb_path'} = $DEFAULT_RUNNABLEDB_PATH;
      $q{'default'}{'output_dir'} = $DEFAULT_OUTPUT_DIR
    }
  }
  return %q;
}


sub useage {
  my ($command_args) = @_;

  print "Your commandline was :\n".
    "lsf_submission.pl ".join("\t", @$command_args), "\n\n";

  print ("lsf_submission.pl is a script which will create and open ".
         "commandlines for the\nappropriate batch submission system ".
         "As standard it will produce command lines to run the\n".
         "test_RunnableDB script but it can be used to create".
         "submission statements\nfor other commands but see the perl".
         "docs for more information\n".
         "Everytime you run job_submission.pl you must pass in the ".
         "database options\n\n-dbhost The host where the pipeline ".
         "database is.\n-dbport The port.\n-dbuser The user to ".
         "connect as.\n-dbpass The password to use.\n".
         "-dbname The database name.\n\n".
         "This script also requires a logic_name of an analysis ".
         "specified with -logic_name.\n\n".
         "Other useful options include:\n\n".
         "-queue_manager allows you to change which ".
         "Bio::EnsEMBL::Pipeline::BatchSubmission module\nbeing ".
         "used\n".
         "-output_dir the directory to put the stderr and stdout".
         "of your analysis runs \n".
         "-script_dir the directory where the test_RunnableDB script".
         "you want to run lives otherwise\nit assumes the script is ".
         "in the current working directory\n\n".
         "-commandline is the option allows you to specify a ".
         "custom commandline which must\ninclude the word ".
         "INPUT_ID at the location in the commandline you wish the ".
         "input_id put\n".
         "-help will print out this".
         "information again\n-perldoc will print out the perl ".
         "documentation\n\n");
  exit(0);
}

j
sub perldoc{
  exec('perldoc', $0);
  exit(0);
}

=pod 

=head1 NAME

lsf_submission.pl

=head1 SYNOPSIS

lsf_submission.pl is a script which create and opens commandlines
appropriate to your batch submission system to run the test_RunnableDB
script

=head1 DESCRIPTION

This script creates submission lines to run the test_RunnableDB script
which will run RunnableDBs which fit the standard pipeline model which 
means they must take three arguments, Input_id, database adaptor and
analysis object and have the standard methods fetch_input, run, output
and write_output. The script does allow for the insertion of analysis
objects and creation of input_ids before commandline creation


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
  -insert_analysis this tells the script to insert the analysis object
  and requires both the follow options are also used
  -input_id_type this specified what type of input_id
  -module this is the name of the RunnableDB to be used

Submission system options

  -queue_manager this specifies which 
  Bio::EnsEMBL::Pipeline::BatchSubmission module is used
  -output_dir which directory to put the stderr and stdout files in
  This is the base directory the script creates 10 directorys from 0-9
  below that where the files are actually placed to ensure that no
  directory gets too many files in it
  -queue the submission system queue the jobs will be submitted to
  -pre_exec toggle if a pre exec command is to be used in the submission
  command
  -resource_requirements any resource requirements for the commandline
  in the LSF system this string will placed after the -R flag
  -sub_args any other arguments you want to pass to your submission system

  All of these options are discretionary and if not specified the
  settings will be taken from the BatchQueue config file

  -open this tells the script to actually open the commandlines to the
  submission system otherwise it will just print the commandlines
  -submission_interval the length of time to sleep for between opening
  submissions by default it is 5 seconds

test_RunnableDB options

  -script_dir the directory where you copy of test_RunnableDB lives
  if not specified the current working directory is used
  -runnabledb_path the perl path for the runnabledb to be used. Again
  this is discretionary and can be taken from BatchQueue if needed
  -write whether test_RunnableDB should write its output to the database
  -update_input_id_analysis whether test_RunnableDB should update the 
  input_id_analysis table on sucessful completion of a job
  -script_verbosity whether test_RunnableDB should have its verbose flag 
  set

other command options

  -command this option allows you to specify an alternative script to 
  test_RunnableDB but you must specify the full alternatve commandline
  on in quotes after this placing the word INPUT_ID in the position in
  the commandline you wish the input_id to appear

  -pre_exec_command if you wish a pre_exec_command to be used with a
  custom commandline you must also specify the pre_exec on the 
  commandline using this options

Input ID options


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
  -utils_verbosity, this affects the amount of chatter you recieve from
  the core module Bio::EnsEMBL::Utils::Exception. By default this is set
  to WARNING which means you see the prints from warnings and throws but
  not from deprecate and info calls. See the modules itself for more 
  information
  -verbose print out information about how the script is running

=head1 CONTACT

Post general queries to <ensembl-dev@ebi.ac.uk>

=head1 EXAMPLES

./lsf_submission.pl -dbhost myhost -dbuser user -dbpass password 
-dbport 3306 -dbname my_pipeline_database -logic_name RepeatMask -open

this will create and submit commandlines to your submission system to
run test_RunnableDB with RepeatMask anallysis

./lsf_submission.pl -dbhost myhost -dbuser user -dbpass password 
-dbport 3306 -dbname my_pipeline_database -logic_name RepeatMask -open
-insert_analysis -input_id_type CONTIG -module RepeatMasker  

this will create and submit the commandlines and it will also insert the 
analysis object called RepeatMask with the input_id_type CONTIG and module
RepeatMasker

./lsf_submission.pl -dbhost myhost -dbuser user -dbpass password 
-dbport 3306 -dbname my_pipeline_database -logic_name RepeatMask -open
-input_id_file input_ids

this will create and submit the commandlines based on the input_ids in
the file input_ids

./lsf_submission.pl -dbhost myhost -dbuser user -dbpass password 
-dbport 3306 -dbname my_pipeline_database -logic_name RepeatMask

this will just print the commandlines constucted to STDOUT

./lsf_submission.pl -dbhost myhost -dbuser user -dbpass password 
-dbport 3306 -dbname my_pipeline_database -logic_name RepeatMask
-open -commandline '/path/to/my_script.pl -dbhost myhost -dbuser user 
-dbpass password -dbport 3306 -dbname my_pipeline_database -input_id
INPUT_ID -analysis RepeatMask'

this would run the custom commandline

./lsf_submission.pl -dbhost myhost -dbuser user -dbpass password 
-dbport 3306 -dbname my_pipeline_database -logic_name RepeatMask
-open -commandline '/path/to/my_script.pl -dbhost myhost -dbuser user 
-dbpass password -dbport 3306 -dbname my_pipeline_database -input_id
INPUT_ID -analysis RepeatMask' -pre_exec 
-pre_exec_command '/path/to/my_script.pl -check'

this would run the custom commandline with a pre_exec check

=head1 SEE ALSO

  rulemanager.pl
  job_submission.pl and
  pipeline_sanity.pl all here in ensembl-pipeline/scripts

and also using_the_ensembl_pipeline.txt in the ensembl-docs cvs module

=cut

