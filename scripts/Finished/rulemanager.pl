#!/usr/local/ensembl/bin/perl -w

use strict;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use Bio::EnsEMBL::Pipeline::Utils::PipelineSanityChecks;
use Bio::EnsEMBL::Pipeline::Finished::RuleManager;
use Bio::EnsEMBL::Pipeline::DBSQL::Finished::DBAdaptor;
use BlastableVersion;

$ENV{BLASTDB} = '/data/blastdb/Finished';

$| = 1;

my $term_sig = 0;

$SIG{TERM} = \&termhandler;
$SIG{INT}  = \&termhandler;

#database arguments
my $dbhost = $ENV{'ENS_DBHOST'};
my $dbname = $ENV{'ENS_DBNAME'};
my $dbuser = $ENV{'ENS_DBUSER'};
my $dbpass = $ENV{'ENS_DBPASS'};
my $dbport = $ENV{'ENS_DBPORT'} || 3306;

#command line args
my $help;             #get docs about script
my $verbose;          #print statements about the running script
my $queue_manager;    #what Bio::EnsEMBL::BatchSubmission module to use

#The rulemanger sleeps when a predefined number of jobs has been reached
my $max_job_time;   # the longest time the rulemanger should sleep for
my $min_job_time;   # the shortest time the rulemanger sleeps for
my $sleep_per_job;  # the length of time to sleep per job when limit is
                    # reached
my $runner;         # the runner script to use when running jobs (this will be
                    # over ridden by anything in config
my $output_dir;     # the output_dir to use when running jobs (this won't be
                    # over ridden when running jobs
my $job_limit;      # the maximum number of jobs of a perdefined status that are
    # in the system. The predefined statuses are set in the BatchQueue.pm config
    # module and currently refer to LSF
my $mark_awol       = 1;  # Flag as whether to mark jobs which have gone missing
                          # from the system
my $rename_on_retry = 1;  #Whether to rename jobs stdout/err which are
                          #being retried
my $once;              # Only run the loop once
my @starts_from;       # array of logic_names which are used to get the starting
                       # input_id hash
my @analyses_to_run;   # array of logic_names of analyses to run
my @analyses_to_skip;  # array of logic_names of analyses to skip
my @types_to_run;      # array of input_id_types to consider for the run
my @types_to_skip;     # array of input_id_types to not consider for the run
my $ids_to_run;        # filepath to file of input_ids to run
my $ids_to_skip;       # filepath to file of input_ids to skip
my $config_sanity      = 1;    #flag as to whether to check configuration sanity
my $accumulator_sanity = 1;    #as above for accumulators
my $rules_sanity       = 1;    #same for rules
my $db_sanity          = 1;    # same for the database tables
my $kill_jobs;      # whether to kill jobs which have been running for to long
my $killed_time;    # how long to let jobs have before being killed
my $kill_file;      # a path to a file to record what jobs have been killed
my $rerun_sleep = 200;    #how long to sleep for at the end of a loop if no
                          #jobs get submitted
my $utils_verbosity = 'WARNING';    #how verbose do you want the
     #Bio::EnsEMBL::Utils::Exceptions module to be by default it is set to
     #WARNING as this gives warning and throws but not deprecates or infos
my $unlock = 0
  ; # deletes the lock of a pipeline (the entry in the meta table 'pipeline.lock'

my $shuffle;
my $accumulators       = 1;
my $force_accumulators = 0;
my $reread_input_ids   = 0;    # toggle whether to reread input_id each time the
                               # script loops
my $reread_rules = 0;      # toggle whether to reread rules each time the script
                           # loops
my $reupdate_analysis = 1;
my $perldoc      = 0;
my @command_args = @ARGV;
my $submission_limit;
my $submission_number = 1000;
GetOptions(
	'dbhost=s'               => \$dbhost,
	'dbname=s'               => \$dbname,
	'dbuser=s'               => \$dbuser,
	'dbpass=s'               => \$dbpass,
	'dbport=s'               => \$dbport,
	'help!'                  => \$help,
	'verbose!'               => \$verbose,
	'queue_manager=s'        => \$queue_manager,
	'max_job_time=s'         => \$max_job_time,
	'min_job_job=s'          => \$min_job_time,
	'sleep_per_job=s'        => \$sleep_per_job,
	'runner=s'               => \$runner,
	'output_dir=s'           => \$output_dir,
	'job_limit=s'            => \$job_limit,
	'mark_awol!'             => \$mark_awol,
	'rename_on_retry'        => \$rename_on_retry,
	'starts_from=s@'         => \@starts_from,
	'analysis|logic_name=s@' => \@analyses_to_run,
	'skip_analysis=s@'       => \@analyses_to_skip,
	'input_id_type=s@'       => \@types_to_run,
	'skip_input_id_type=s@'  => \@types_to_skip,
	'input_id_file=s'        => \$ids_to_run,
	'skip_input_id_file=s'   => \$ids_to_skip,
	'config_sanity!'         => \$config_sanity,
	'accumulator_sanity!'    => \$accumulator_sanity,
	'db_sanity!'             => \$db_sanity,
	'rules_sanity!'          => \$rules_sanity,
	'kill_jobs!'             => \$kill_jobs,
	'killed_time=s'          => \$killed_time,
	'kill_file=s'            => \$kill_file,
	'rerun_sleep=s'          => \$rerun_sleep,
	'utils_verbosity=s'      => \$utils_verbosity,
	'shuffle!'               => \$shuffle,
	'accumulators!'          => \$accumulators,
	'force_accumulators!'    => \$force_accumulators,
	'reread_input_ids!'      => \$reread_input_ids,
	'reread_rules!'          => \$reread_rules,
	'reupdate_analysis!'	 => \$reupdate_analysis,
	'once!'                  => \$once,
	'perldoc!'               => \$perldoc,
	'submission_limit!'      => \$submission_limit,
	'submission_number=s'    => \$submission_number,
	'unlock|delete_lock'     => \$unlock,
  )
  or useage( \@command_args );

perldoc() if $perldoc;
verbose($utils_verbosity);

@analyses_to_run = map {split/,/} @analyses_to_run ;

unless ( $dbhost && $dbname && $dbuser ) {
	print STDERR
	  "Must specify database with -dbhost, -dbname, -dbuser and -dbpass\n";
	print STDERR "Currently have -dbhost $dbhost -dbname $dbname "
	  . "-dbuser $dbuser -dbpass $dbpass -dbport $dbport\n";
	$help = 1;
}
if ($help) {
	useage( \@command_args );
}

if ( !defined($mark_awol) ) {
	$mark_awol = 0;
}

my $db = Bio::EnsEMBL::Pipeline::DBSQL::Finished::DBAdaptor->new(
	-host   => $dbhost,
	-dbname => $dbname,
	-user   => $dbuser,
	-pass   => $dbpass,
	-port   => $dbport,
);

my $sanity =
  Bio::EnsEMBL::Pipeline::Utils::PipelineSanityChecks->new( -DB => $db, );
my $rulemanager = Bio::EnsEMBL::Pipeline::Finished::RuleManager->new(
	-DB              => $db,
	-QUEUE_MANAGER   => $queue_manager,
	-MARK_AWOL       => $mark_awol,
	-RENAME_ON_RETRY => $rename_on_retry,
	-MAX_JOB_SLEEP   => $max_job_time,
	-MIN_JOB_SLEEP   => $min_job_time,
	-SLEEP_PER_JOB   => $sleep_per_job,
	-VERBOSE         => $verbose,
	-RUNNER          => $runner,
	-OUTPUT_DIR      => $output_dir,
	-JOB_LIMIT       => $job_limit,
	-UNLOCK          => $unlock,
);

my $done;
my $reset;

if ($config_sanity) {
	$sanity->config_sanity_check;
}

if ($db_sanity) {
	$sanity->db_sanity_check;
}

my %analyses_to_run  = %{ $rulemanager->logic_name2dbID( \@analyses_to_run ) };
my %analyses_to_skip = %{ $rulemanager->logic_name2dbID( \@analyses_to_skip ) };

if ( $ids_to_run && !-e $ids_to_run ) {
	throw("Must be able to read $ids_to_run");
}

if ( $ids_to_skip && !-e $ids_to_skip ) {
	throw("Must be able to read $ids_to_skip");
}

if (   $ids_to_run
	|| @analyses_to_run
	|| @types_to_run
	|| @starts_from
	|| @analyses_to_skip
	|| @types_to_skip
	|| $ids_to_skip
	|| $submission_limit )
{
	print STDERR "You are running with options which may break "
	  . "accumulators they are being switched off. You can prevent this with "
	  . "the -force_accumulators option\n";
	$accumulators = 0;
}

update_analysis_dbversion();

my $all_rules = $rulemanager->rules;

if ( $accumulators && $accumulator_sanity ) {
	$accumulators =
	  $sanity->accumulator_sanity_check( $all_rules, $accumulators );
}

$accumulators = 1 if ($force_accumulators);

if ($rules_sanity) {
	$sanity->rule_type_sanity( $all_rules, $verbose );
}


my $space = $rulemanager->job_stats($job_limit);

my %always_incomplete_accumulators;
my %accumulator_analyses;

foreach my $rule (@$all_rules) {
	if ( $rule->goalAnalysis->input_id_type eq 'ACCUMULATOR' ) {
		$accumulator_analyses{ $rule->goalAnalysis->logic_name } =
		  $rule->goalAnalysis;
	}
}

setup_pipeline(
	\%analyses_to_run,                \%analyses_to_skip,
	$all_rules,                       \%accumulator_analyses,
	\%always_incomplete_accumulators, $ids_to_run,
	$ids_to_skip,                     \@types_to_run,
	\@types_to_skip,                  \@starts_from,
	$rulemanager
);


my %completed_accumulator_analyses;
my $submission_count = 0;
LOOP: while (1) {
	print "Reading IDs \n" if $verbose;
	my $submitted = 0;
	update_analysis_dbversion() if($reupdate_analysis);

	%completed_accumulator_analyses =
	  %{ $rulemanager->fetch_complete_accumulators };
	my %incomplete_accumulator_analyses = %always_incomplete_accumulators;

	if ($reread_input_ids) {
		$rulemanager->input_id_setup(
			$ids_to_run,     $ids_to_skip, \@types_to_run,
			\@types_to_skip, \@starts_from
		);
	}
	if ($reread_rules) {
		$rulemanager->empty_rules;
		$all_rules = $rulemanager->rules;
		$rulemanager->rules_setup( \%analyses_to_run, \%analyses_to_skip,
			$all_rules, \%accumulator_analyses,
			\%always_incomplete_accumulators );
	}

	my $id_hash = $rulemanager->input_ids;

  INPUT_ID_TYPE:
	foreach my $type ( keys(%$id_hash) ) {
		next INPUT_ID_TYPE if ( $type eq 'ACCUMULATOR' );

		my @id_list = keys( %{ $id_hash->{$type} } );
		@id_list = shuffle(@id_list) if $shuffle;

	  INPUT_ID:
		foreach my $input_id (@id_list) {

			if ($term_sig) {
				print "Got term signal\n" if ($verbose);
				$done = 1;
				last INPUT_ID_TYPE;
			}

			print "\n===\nInput id " . $input_id . "\n" if ($verbose);

			my %analHash;
			my @anals = @{
				$rulemanager->stateinfocontainer->fetch_analysis_by_input_id(
					$input_id)
			  };

			foreach my $rule ( @{ $rulemanager->rules } ) {
				my $anal =
				  $rule->check_for_analysis( \@anals, $type,
					\%completed_accumulator_analyses, $verbose );
				if ( UNIVERSAL::isa( $anal, 'Bio::EnsEMBL::Pipeline::Analysis' )
				  )
				{
					if ( $anal->input_id_type ne 'ACCUMULATOR' ) {
						$analHash{ $anal->dbID } = $anal;
					}
				}
				else {
					if (   $rule->goalAnalysis->input_id_type eq 'ACCUMULATOR'
						&& $rule->has_condition_of_input_id_type($type) )
					{
						$incomplete_accumulator_analyses{ $rule->goalAnalysis
							  ->logic_name } = 1;
					}
				}
			}
			my $current_jobs_hash =
			  $rulemanager->job_adaptor->fetch_hash_by_input_id($input_id);

			print "\n\n Ended up with "
			  . scalar( keys %analHash )
			  . " analyses which pass conditions\n"
			  if ($verbose);
			for my $anal ( values %analHash ) {
				print  "For " . $anal->logic_name . ":\n" if $verbose;;
				if (
					$rulemanager->can_job_run(
						$input_id, $anal, $current_jobs_hash
					)
				  )
				{
					$submitted++;
					$submission_count++;
					if (   $submission_limit
						&& $submission_count >= $submission_number )
					{
						$done = 1;
						last INPUT_ID_TYPE;
					}
				}
				else {
					print  "cannot run\n" if $verbose;;
				}
			}
		}

	}

	my $naccum_submitted = 0;
	if ( !$done && !$reset ) {
		if ($accumulators) {
			%completed_accumulator_analyses =
			  %{ $rulemanager->fetch_complete_accumulators };

			my $accumulator_jobs_hash =
			  $rulemanager->job_adaptor->fetch_hash_by_input_id('ACCUMULATOR');

			foreach my $logic_name ( keys %accumulator_analyses ) {
				if (   !exists( $incomplete_accumulator_analyses{$logic_name} )
					&& !exists( $completed_accumulator_analyses{$logic_name} ) )
				{

					if (
						$rulemanager->can_job_run(
							'ACCUMULATOR',
							$accumulator_analyses{$logic_name},
							$accumulator_jobs_hash
						)
					  )
					{
						$submitted++;
						$submission_count++;
						$naccum_submitted++;
					}
					elsif (
						exists( $incomplete_accumulator_analyses{$logic_name} )
					  )
					{
						print "Accumulator type analysis $logic_name "
						  . "conditions unsatisfied\n"
						  if $verbose;
					}
					else {
						print "Accumulator type analysis $logic_name "
						  . "already run\n"
						  if $verbose;
					}
				}
			}
		}
	}

	$space = $rulemanager->job_stats($job_limit);

	if ( !$done && !$naccum_submitted ) {
		if ( !$rulemanager->check_if_done ) {
			$done = 1;
		}
		else {
			$done = 0;
		}
	}

	if ( $done || $once ) {
		$rulemanager->db->pipeline_unlock;
		exit 0;
	}
	sleep($rerun_sleep)                if ($submitted == 0);
	print "Waking up and run again!\n" if $verbose;
}

sub setup_pipeline {
	my (
		$analyses_to_run,                $analyses_to_skip,
		$all_rules,                      $accumulator_analyses,
		$always_incomplete_accumulators, $ids_to_run,
		$ids_to_skip,                    $types_to_run,
		$types_to_skip,                  $starts_from,
		$rulemanager
	  )
	  = @_;

	$rulemanager->rules_setup( $analyses_to_run, $analyses_to_skip, $all_rules,
		$accumulator_analyses, $always_incomplete_accumulators );
	$rulemanager->input_id_setup(
		$ids_to_run,    $ids_to_skip, $types_to_run,
		$types_to_skip, $starts_from
	);
}

# Update the db version for Halfwise, Blast, EST and DepthFilter analysis in the analysis table
my @analysis;
sub update_analysis_dbversion {
	my $analysis_adaptor = $db->get_AnalysisAdaptor();
	@analysis         = @{$analysis_adaptor->fetch_all()} unless(@analysis);
	my $a_hash;
	my $uniprot_db_version;
	my $pfam_db_version;

	foreach my $ana (@analysis) {
		my $db_file;
		my $db_version;
		my $ln = $ana->logic_name;
		my $mod = $ana->module();
		my $ana_df;
		if ( $mod && ( $mod eq 'Blast' || $mod eq 'EST' || $mod =~ /Exonerate/ || $mod =~ /Halfwise/ ) ) {
			$a_hash->{$ln} = $ana;
			$db_file    = fetch_databases( $ana->db_file );
			$db_version = get_db_version( $db_file , 1);
			save_db_version($analysis_adaptor,$ana,$db_version) unless $mod =~ /Halfwise/;
			# save uniprot and pfam db versions
			$uniprot_db_version = $db_version if $ana->logic_name =~ /Uniprot_raw/;
			$pfam_db_version = $db_version if $mod =~ /Halfwise/;
		}
	}

	foreach my $ana (@analysis) {
		next unless $ana->module;
		# update the DepthFilter analysis
		if($ana->module =~ /DepthFilter/) {
			my ($ori_ana) = $ana->parameters =~ /ori_analysis\s+=>\s+(\w+)/;
			save_db_version($analysis_adaptor,$ana,$a_hash->{$ori_ana}->db_version);
		}
		# update Halfwise
		elsif ($ana->module =~ /Halfwise/) {
			my $new_db_version  = $uniprot_db_version.'_'.$pfam_db_version;
			save_db_version($analysis_adaptor,$ana,$new_db_version);
		}
	}
}

sub save_db_version {
	my ($a_ad,$ana,$db_version) = @_;
	if($db_version && ($db_version ne $ana->db_version))
	{
		print STDOUT "Analysis update: ".$ana->logic_name." old version ".$ana->db_version.
					 " new version $db_version\n" if $verbose;
		$ana->db_version($db_version);
		$a_ad->update($ana);
	}
}

my %db_filename = ();
sub fetch_databases {
	my ($db) = @_;
	my @databases;
	$db =~ s/\s//g;
	my $path;
	foreach my $dbname ( split( ",", $db ) ) {
		# fetch the update file db version
		$dbname =~ s/^embl_/emnew_/g;
		return $db_filename{$dbname} if($db_filename{$dbname});
		$path = $dbname;
		unless ( $dbname =~ m!^/! ) {
			$path = $ENV{BLASTDB} . "/" . $dbname;
		}
		if ( -f "${path}" ) {
			$db_filename{$dbname} = $path;
		} elsif( -f "${path}-1") {
			$db_filename{$dbname} = "${path}-1";
		}
		return $db_filename{$dbname} if($db_filename{$dbname});
	}
	throw( "No databases exist for " . $db );
}

sub get_db_version {
	my ($db,$force_dbi) = @_;
	my $ver = eval {
		my $blast_ver = BlastableVersion->new();
		$blast_ver->force_dbi($force_dbi);
		$blast_ver->get_version($db);
		$blast_ver;
	};

	my $dbv = '';
	if($@ =~ /not in blast tracking db/ ) {
		$dbv = get_db_version($db,0);
	} elsif($@) {
		throw("I failed to get a BlastableVersion for $db [$@]");
	} else {
		$dbv  = $ver->version();
	}

	return $dbv;
}

sub shuffle {
	my (@in) = @_;
	my @out;

	srand;

	push @out, splice( @in, rand @in, 1 ) while @in;

	return @out;
}

sub termhandler {
	$term_sig = 1;
}

sub useage {
	my ($command_args) = @_;
	print "Your commandline was :\n"
	  . "rulemanager.pl "
	  . join( "\t", @$command_args ), "\n\n";
	print "ensembl-pipeline/scripts/rulemanager.pl is the main script used "
	  . "to run the pipeline\n\n";
	print "Everytime you run the rulemanager you must pass in the database "
	  . "options\n\n";
	print "-dbhost     The host where the pipeline database is.\n"
	  . "-dbport     The port.\n"
	  . "-dbuser     The user to connect as.\n"
	  . "-dbpass     The password to use.\n"
	  . "-dbname   The database name.\n\n";
	print "Other options you may find useful are:\n\n"
	  . "-once, which means the RuleManager loop only executes once "
	  . "before exiting\n"
	  . "-analysis and -skip_analysis, where you can specific an "
	  . "individual analysis to run or skip\n"
	  . "-input_id_file and -skip_input_id_file, files containing a "
	  . "list of input ids to run in the format input_id "
	  . "input_id_type\n\n";
	print " -perldocs will print out the perl documentation of this module "
	  . "and -help will print out the help again \n";
	exit(0);
}

sub perldoc {
	exec( 'perldoc', $0 );
	exit(0);
}

=pod

=head1 NAME

rulemanager.pl, a script for running the ensembl-pipeline.
This is a script which replicates most of the behaviour previously
described by RuleManager3.pl

for more information about the ensembl-pipeline or using rulemanager.pl
read using_the_ensembl_pipeline.txt in ensembl-doc,

There is a genome research paper
The ensembl analysis pipeline.
Genome Res. 2004 May;14(5):934-41.
PMID: 15123589

or mail ensembl-dev@ebi.ac.uk

=head1 SYNOPSIS

This script takes a series of input_ids, first group by input_id_type
then establishes on the basis of rules described in the pipeline database
what analyses it can run with what input_ids

=head1 OPTIONS

DB Connection Details


   -dbhost     The host where the pipeline database is.
   -dbport     The port.
   -dbuser     The user to connect as.
   -dbpass     The password to use.
   -dbname     The database name.

These arguments are always required all other arguments are optional

Affecting the scripts behaviour

   -once this flag means the script will only run through its loop once
    this is useful for testing
   -reread_rules, this flag will force the rulemanager to reread the rule
    table every loop
   -reread_input_ids, this will force the rulemanager to reread the input
    ids each and every loop
   -reupdate_analysis, this flag will update the database version stored in
    the analysis table each and every loop, by default set to true.
   -accumulators, this is a flag to indicate if the accumulators are
    running. By default it is on but it is turned off if you specify any
    flag which affects the rules or input_ids sets. You can switch
    it off by specifying -noaccumulators
   -force_accumulators this forces accumulators on even if other conditions
    would switch it on
   -shuffle this reorders the array of input_ids each time it goes through
    through the loop. The purpose of this is to ensure the ids aren't
    always checked in the same order'
   -utils_verbosity, this affects the amount of chatter you recieve from
    the core module Bio::EnsEMBL::Utils::Exception. By default this is set
    to WARNING which means you see the prints from warnings and throws but
    not from deprecate and info calls. See the modules itself for more
    information
   -verbose, toggles whether some print statements are printed in this
    script and in the RuleManager
   -rerun_sleep, this is the amount of time the script will sleep for
    if no jobs have been submitted in the previous loop

Overridable Configurations options from Bio::EnsEMBL::Pipeline::Config
modules

  Options also defined in Bio::EnsEMBL::Pipeline::Config modules

  BatchQueue.pm

  -queue_manager this specifies which
   Bio::EnsEMBL::Pipeline::BatchSubmission module is used by
   Bio::EnsEMBL::Pipeline::RuleManager
  -job_limit the maximun number of jobs of specified status allowed in
   system
  -max_job_time the maximum time the rulemanager will sleep for when
   job limit is reached
  -min_job_time the minimum time the rulemanager will sleep for when
   job limit is reached
  -sleep_per_job the amount of time per job the rulemanager will sleep
   for if the job limit is reached
  -output_dir the path to an output directory when the jobs stderr and
   stdout will be redirected. This alwaysoverrides the values specified in
   BatchQueue
  -mark_awol toggle to specify whether to mark jobs as awol if lost from
   the submission system this can apply strain to the job table. It is on
   by default and can be switched off using the -nomark_awol flag

  General.pm
   -runner path to a default runner script. This will override what is
    set in General.pm but will be overidden by any analyses specific
    settings found in BatchQueue.pm
   -rename_on_retry a toggle to specify whether to rename stderr/stdout
    file when a job is retried as otherwise the submission system just
    cats them all together

Sanity Checking

   -config_sanity this is a test to check certain values are defined in
    your General.pm and BatchQueue.pm config files
   -db_sanity this checks various tables in the databases for missing
   values or bad foreign key relationships
   -accumulator_sanity this checks you accumulators are set up correctly
   -rules_sanity this checks you rules are setup correctly

   this are all on by default but can be switched off using -nooption_name

Altering what input_ids or analyses are considered for submission

  -starts_from this should be a logic_name of an analysis which has been
  run. The input_ids used for the central loop are retrived on the basis of
  these logic_names. This option can appear in the commandline multiple
  times
  -analysis a logic_name of an analysis you want run. If you specify
  this option you will only run this analysis or analyses as the option can
  appear on the command line multiple times
  -skip_analysis a logic_name of an analysis you don't want to run. If
  this option is specified these are the only analyses which won't be run
  note this option and -analysis are mutually exclusive
  -input_id_type a type of input_id you want to run with. If these is
  specified no other type of input will be considered. Again this can
  appear multiple times on the commandline
  -skip_input_id_type a type of input to no be considered. If this appears
  only these input_ids will not be run. Again this can appear multiple
  times on the commandline and is multially exclusive with -input_id_type
  -input_id_file path to a file of input_ids which should be in the format
  input_id input_id_type
  -skip_input_id_file path to a file of input_ids to not run. The file
  should be in the same format as for -input_id_file
  -submission_limit, only submit a defined number of jobs, by default 1000
  -submission_number, the number of jobs to submit

Currently unused options

  -kill_jobs, flag as whether to check hwo long jobs have been running for
   and kill them if they have run for too long
  -killed_time the length of time to let jobs run for before killing them
  -kill_file a file where to record the input_ids and logic_names of killed
   jobs

these options are currently unused as we aren't sure of the best way to
check job running time. They have been brought over from RuleManager3.pl
as the functionality may be reinstated at a later stage'

=head1 EXAMPLES

a standard commandline for running rulemanager

./rulemanager -dbhost myhost -dbuser user -dbpass password -dbport 3306
-dbname my_pipeline_database -shuffle

it is generally a good idea to redirect your output to a file and run the
process in the background

./rulemanager -dbhost myhost -dbuser user -dbpass password -dbport 3306
-dbname my_pipeline_database -shuffle >& rulemanager_output.txt&

if you still want to see the output on the screen you can do this by piping
the output though tee

./rulemanager -dbhost myhost -dbuser user -dbpass password -dbport 3306
-dbname my_pipeline_database -shuffle | tee rulemanager_output.txt

if you have just set up you database you may want to test everything is
okay. the input_id_file option is a good way to do this

./rulemanager -dbhost myhost -dbuser user -dbpass password -dbport 3306
-dbname my_pipeline_database -shuffle -input_id_file my_input_ids

as this way you can get the pipeline to considered a limited number of
input_ids and as such is causes less problems if things go wrong


=head1 SEE ALSO

  lsf_submission.pl
  job_submission.pl and
  pipeline_sanity.pl all here in ensembl-pipeline/scripts

and also using_the_ensembl_pipeline.txt in the ensembl-docs cvs module

=cut

