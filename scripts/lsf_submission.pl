#!/usr/local/ensembl/bin/perl 

use strict;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use Bio::EnsEMBL::Pipeline::Utils::PipelineSanityChecks;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Utils::InputIDFactory;
use Bio::EnsEMBL::Pipeline::Config::General;
use Bio::EnsEMBL::Pipeline::Config::BatchQueue;
use Bio::EnsEMBL::Pipeline::Analysis;
use Bio::EnsEMBL::Pipeline::RuleManager;
$| = 1;



#database arguments
my $dbhost    = $ENV{'ENS_DBHOST'};
my $dbname    = $ENV{'ENS_DBNAME'};
my $dbuser    = $ENV{'ENS_DBUSER'};
my $dbpass    = $ENV{'ENS_DBPASS'};
my $dbport    = $ENV{'ENS_DBPORT'} || 3306;
my $help; #get docs about script
my $verbose; #print statements about the running script
my $perldoc;
my @command_args = @ARGV;
my $utils_verbosity = 'WARNING'; #how verbose do you want the 
#Bio::EnsEMBL::Utils::Exceptions module to be by default it is set to
#WARNING as this gives warning and throws but not deprecates or infos
my $queue_manager; #what Bio::EnsEMBL::BatchSubmission module to use
my $script_dir = `pwd`;
chomp $script_dir;
my $output_dir;
my $logic_name;
my $config_sanity = 1;
my $queue;
my $runnabledb_path;
my $write;
my $update;
my $pre-exec;
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
my $input_id_type;
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
GetOptions(
           'dbhost=s'      => \$dbhost,
           'dbname=s'      => \$dbname,
           'dbuser=s'      => \$dbuser,
           'dbpass=s'      => \$dbpass,
           'dbport=s'      => \$dbport,
           'help!' => \$help,
           'verbose!' => \$verbose,
           'perldoc!' => \$perldoc,
           'utils_verbosity=s' => \$utils_verbosity,
           'queue_manager=s' => \$queue_manager,
           'script_dir=s' => \$script_dir,
           'output_dir=s' => \$output_dir,
           'queue=s' => \$queue,
           'logic_name=s' => \$logic_name,
           'config_sanity!' => \$config_sanity,
           'runnabledb_path=s' =>\$runnabledb_path,
           'write!' => $write,
           'update_input_id_analysis!' => \$update,
           'pre_exec!' => \$pre_exec,
           'script_verbosity!' => \$script_verbosity,
           'resource_requirements:s' => \$resource,
           'sub_args:s' => \$sub_args,
           'input_id_file:s' => $ids_to_run,
           'make_input_ids!' => \$make_input_ids,
           'coord_system:s'       => \$coord_system,
           'coord_system_version:s' => \$coord_system_version,
           'slice!'        => \$slice,
           'slice_size:s' => \$slice_size,
           'slice_overlap:s' => \$slice_overlap,
           'logic_name:s' => \$logic_name,
           'file!'         => \$file,
           'dir:s'        => \$dir,
           'file_regex:s' => \$regex,
           'single!'       => \$single,
           'single_name:s'=> \$name,
           'translation_ids!' => \$translation_id,
           'seq_level!' => \$seq_level,
           'top_level!' => \$top_level,
           'open!' => \$open,
           'submission_interval:s' => \$submission_interval,
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
  print STDERR "You must specify a logic_name for the analysis you ".
    "want to run\n";
  $help = 1;
}
if($insert_analysis && (!$input_id_type || !$module)){
  print STDERR "If you want to insert your analysis object into your ".
    "you must specify input_id_type and module on the commandline too\n";
  $help = 1;
}
if($help){
  useage(\@command_args);
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
  throw("Can't find $file [$@]");
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
if($config_sanity){
  $sanity->config_sanity_check;
}
my %BATCH_QUEUES = &set_up_queues;
my $batch_q_key = $logic_name;
if(!$BATCH_QUEUES{$batch_q_key}){
  $batch_q_key = 'default';
}
if(!$queue){
  $queue = $BATCH_QUEUES{$batch_q_key}{'queue'};
}
if(!$runnabledb_path){
  $runnabledb_path = $BATCH_QUEUES{$batch_q_key}{'runnabledb_path'};
}
if(!$resource){
  $resource = $BATCH_QUEUES{$batch_q_key}{'resource'};
}
if(!$sub_args){
  $sub_args = $BATCH_QUEUES{$batch_q_key}{'sub_args'};
}
if(!$output_dir){
  $output_dir = $BATCH_QUEUES{$batch_q_key}{'output_dir'};
  if(!$output_dir){
    $output_dir = $DEFAULT_OUTPUT_DIR;
  }
}
if(!-d $output_dir){
  throw("Can't use output_dir ".$output_dir." as its not a directory");
}
my $analysis;
if($insert_analysis){
  my $analysis = Bio::EnsEMBL::Pipeline::Analysis->new
    (
     -logic_name => $logic_name,
     -module => $module,
     -input_id_type => $input_id_type,
    );
}else{
  $analysis = $db->get_AnalysisAdaptor->fetch_by_logic_name($logic_name);
}
if(!$analysis){
  throw("Can't run without an Analysis object");
}

my $input_ids = setup_input_ids($analysis, $db, $ids_to_run, 
                                $make_input_ids, $slice, $file, 
                                $translation_id, $single, $slice_size, 
                                $slice_overlap, $dir, $regex, $name, 
                                $seq_level, $top_level, $verbose, 
                                $logic_name);

my @batch_submission_objects;

foreach my $input_id(@$input_ids){
  my ($stdout, $stderr) = make_filenames($output_dir, $input_id, 
                                         $logic_name);
  my $pre_command;
  if($pre_exec){
    $pre_command = $script_dir."/test_RunnableDB -check ";
  }
  my $batch_job = $batch_q_module->new
      (
       -STDOUT     => $stdout,
       -STDERR     => $stderr,
       -PARAMETERS => $sub_args,
       -PRE_EXEC   => $pre_command,
       -QUEUE      => $queue,
       -JOBNAME    => $db->dbname . ':' . $logic_name,
       -RESOURCE   => $resource,
      );
  my $command = $script_dir."/test_RunnableDB ";
  $command .= "-dbhost ".$db->host." -dbuser ".$db->username." -dbname ".
    $db->dbname." -dbport ".$db->port." ";
  $command .= " -dbpass ".$db->password." " if($dbpass);
  $command .= " -input_id ".$input_id." -logic_name ".$logic_name." ";
  $command .= " -write " if($write);
  $command .= " -update_input_id_analysis " if($update);
  $command .= " -runnabledb_path " if($runnabledb_path);
  #$command .= " -utils_verbosity ".$utils_verbosity." " 
  #  if($utils_verbosity);
  $command .= " -input_id_type ".$input_id_type." " if($input_id_type);
  $command .= " -module ".$module." " if($module);
  $batch_job->construct_command_line($command);
  push(@batch_submission_objects, $batch_job);
}

foreach my $batch_object(@batch_submission_objects){
  print $batch_object->bsub."\n" if($verbose);
  if($open){
    eval{
      $batch_object->open_command_line;
    };
    if($@ || !$batch_object->id){
      throw("Failed to open ".$batch_object->bsub." $@");
    }
    sleep($submission_interval);
  }
}
sub setup_input_ids{
  my ($analysis, $db, $ids_to_run, $make_input_ids,
      $slice, $file, $translation_id, $single, $slice_size, 
      $slice_overlap, $dir, $regex, $name, $seq_level, $top_level, 
      $verbose, $logic_name) = @_;

  if($ids_to_run){
    my $id_hash = Bio::EnsEMBL::Pipeline::RuleManager->read_id_file
      ($ids_to_run);
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
    my @ids = keys(%{$id_hash->{$analysis->input_id_type}});
    return \@ids;
  }elsif($make_input_ids){
    print STDERR "Making input ids\n" if($verbose);
    my $ids = make_input_ids($slice, $file, $translation_id, $single, 
                          $slice_size, $slice_overlap, $dir, 
                          $regex, $name, $seq_level, $top_level, 
                          $verbose, $logic_name, $input_id_type);
    return $ids;
  }else{
    my @ids = @{$db->get_StateInfoContainer
                  ->list_input_ids_by_type($analysis->input_id_type)};
    return \@ids;
  }
}

sub make_filenames {
  my ($output_dir, $input_id, $logic_name) = @_;
  
  my $num = int(rand(10));
  
  my $dir = $output_dir . "/$num/";
  if( ! -e $dir ) {
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
      $q{$ln}{'batch_size'} ||= $DEFAULT_BATCH_SIZE;
      $q{$ln}{'queue'} ||= $DEFAULT_BATCH_QUEUE;
      $q{$ln}{'retries'} ||= $DEFAULT_RETRIES;
      $q{$ln}{'cleanup'} ||= $DEFAULT_CLEANUP;
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
sub useage{
	exec('perldoc', $0);
	exit(0);
}

