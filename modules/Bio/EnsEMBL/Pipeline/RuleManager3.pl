# Script for operating the analysis pipeline
#
# You may distribute this code under the same terms as perl itself


use strict;
use Getopt::Long;
use Sys::Hostname;
use Socket;

use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;


use Bio::EnsEMBL::Pipeline::Config::General;
use Bio::EnsEMBL::Pipeline::Config::BatchQueue;


unless (&config_sanity_check) {
    exit 1;
}


# Also, an alarm is set to go off every (say) 2 minutes (or not at all if
# $wakeup is false). When the script receives this it looks to see how many
# jobs are in the queue and, if there too many, goes to sleep for a while.
# You could also do this by counting the number of distinct(LSF_id) entries
# in the job table but you would have to also check for dead processes etc.
# It looks for failed jobs as well and restarts them.

# signal parameters
my $term_sig =  0;
my $rst_sig  =  0;
my $alarm    =  0;
my $wakeup   =  120;   # period to check batch queues; set to 0 to disable

# the signal handlers
$SIG{USR1} = \&sighandler;
$SIG{ALRM} = \&alarmhandler;
$SIG{TERM} = \&termhandler;
$SIG{INT} = \&termhandler;


# dynamically load appropriate queue manager (e.g. LSF)


my $max_time = $MAX_JOB_TIME;
my $killed_file = $KILLED_INPUT_IDS;
my $batch_q_module = "Bio::EnsEMBL::Pipeline::BatchSubmission::$QUEUE_MANAGER";

my $file = "$batch_q_module.pm";
$file =~ s{::}{/}g;
eval {
    require "$file";
};
if ($@) {
    print STDERR "Error trying to load $batch_q_module;\ncan't find $file\n";
    exit 1;
}

my $get_pend_jobs;
if ($batch_q_module->can("get_pending_jobs")) {
    my $f = $batch_q_module . "::get_pending_jobs";
    $get_pend_jobs = \&$f;
}


# command line options override 
# anything set in the environment variables.

my $dbhost    = $ENV{'ENS_DBHOST'};
my $dbname    = $ENV{'ENS_DBNAME'};
my $dbuser    = $ENV{'ENS_DBUSER'};
my $dbpass    = $ENV{'ENS_DBPASS'};
my $dbport    = $ENV{'ENS_DBPORT'} || 3306;

$| = 1;

my $local        = 0;       # Run failed jobs locally
my @analyses;               # Only run this analysis ids
my $submitted;
my $idlist_file;
my ($done, $once);
my $runner;
my $shuffle;  
my $output_dir;
my %analyses;
my $verbose;
my $rerun_sleep = 3600;
my $overload_sleep = 300;
my $db_sanity = 1; 
my $help;

GetOptions(
    'dbhost=s'      => \$dbhost,
    'dbname=s'      => \$dbname,
    'dbuser=s'      => \$dbuser,
    'dbpass=s'      => \$dbpass,
    'dbport=s'      => \$dbport,
    'local'         => \$local,
    'idlist_file=s' => \$idlist_file,
    'runner=s'      => \$runner,
    'output_dir=s'  => \$output_dir,
    'once!'         => \$once,
    'shuffle!'      => \$shuffle,
    'analysis=s@'   => \@analyses,
    'v!'            => \$verbose,
    'dbsanity!'     => \$db_sanity,
    'h|help'	    => \$help,	   
)
or useage();

if(!$dbhost || !$dbname || !$dbuser){
  print STDERR " you must provide a host and a database name and a user".
  " for you db connection\n";
  $help = 1;
}

useage() if $help;

unless ($dbhost && $dbname && $dbuser) {
    print STDERR "Must specify database with -dbhost, -dbname, -dbuser and -dbpass\n";
    exit 1;
}


my $db = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new(
    -host   => $dbhost,
    -dbname => $dbname,
    -user   => $dbuser,
    -pass   => $dbpass,
    -port   => $dbport,
);

if($db_sanity){
  &db_sanity_check($db);
}else{
  print STDERR "You are not checking db sanity are you sure you set it ".
    "up correctly?\n" if($verbose);
}

my $rule_adaptor = $db->get_RuleAdaptor;
my $job_adaptor  = $db->get_JobAdaptor;
my $ana_adaptor  = $db->get_AnalysisAdaptor;
my $sic          = $db->get_StateInfoContainer;

# analysis options are either logic names or analysis dbID's



%analyses   = logic_name2dbID($ana_adaptor, @analyses);


if ($idlist_file && ! -e $idlist_file) {
    die "Must be able to read $idlist_file";
}


# Lock to prevent >1 instances of the same pipeline running together
#
# Puts an entry in the meta table keyed on 'pipeline.lock', which
# gives the host running the script, user, PID and when it was
# started

if (my $lock_str = $db->pipeline_lock) {
    # Another pipeline is running: describe it
    my($user, $host, $pid, $started) = $lock_str =~ /(\w+)@(\w+):(\d+):(\d+)/;
    $started = scalar localtime $started;

    print STDERR <<EOF;
Error: this pipeline appears to be running!

    db       $dbname\@$dbhost
    pid      $pid on host $host
    started  $started

The process above must be terminated before this script can be run.
If the process does not exist, remove the lock by removing the lock
from the pipeline database:

    delete from meta where meta_key = 'pipeline.lock';

Thank you


EOF

    exit 1;
}

# create lock

my $host = &qualify_hostname(hostname());
my $user = scalar getpwuid($<);
my $lock_str = join ":", "$user\@$host", $$, time();
$db->pipeline_lock($lock_str);


# Fetch all the analysis rules.  These contain details of all the
# analyses we want to run and the dependences between them. e.g. the
# fact that we only want to run blast jobs after we've repeat masked etc.

my @rules       = $rule_adaptor->fetch_all;

my %accumulator_analyses;

foreach my $rule (@rules) {
  if ($rule->goalAnalysis->input_id_type eq 'ACCUMULATOR') {
    $accumulator_analyses{$rule->goalAnalysis->logic_name} = $rule->goalAnalysis;
  }
}

# Need here to strip rules which don't need to be run.

my @id_list;     # All the input ids to check

alarm $wakeup if $wakeup;  
                # Signal to the script to do something in the future,
                # Such as check for failed jobs, no of jobs in queue

my %completed_accumulator_analyses;

while (1) {
    $submitted = 0;

    # This loop reads input ids from the database a chunk at a time
    # until we have all the input ids.

    my $id_type_hash;

    if (defined $idlist_file) {

        # read list of id's from a file,
        # e.g. a list of contigs on the golden path

        open IDS, "< $idlist_file" or die "Can't open $idlist_file";
        while (<IDS>) {
            chomp;
            my($id) = split;
            ($id) or die "Invalid id $_";
            push @id_list, $id;
        }
        close IDS;
    }
    else {

        # This loop reads input ids from the database a chunk at a time
        # until we have all the input ids.
        # NB It's almost as much work to get one ID as the whole lot, so setting
        # the 'chunksize' variable to a small number doesn't really achieve much.

        print "Reading IDs ... ";

        $id_type_hash = $sic->get_all_input_id_analysis_sets;
    }

    my @anals = @{$sic->fetch_analysis_by_input_id('ACCUMULATOR')};

    foreach my $anal (@anals) {
        if ($anal->input_id_type eq 'ACCUMULATOR') {
            print "Adding completed accumulator for " . $anal->logic_name . "\n";

            $completed_accumulator_analyses{$anal->logic_name} = 1;
        } else {
            print STDERR "WARNING: Expected input_id_type to be ACCUMULATOR for input_id ACCUMULATOR\n"
        }
    }

    # Now we loop over all the input ids. We fetch all the analyses
    # from the database that have already run. We then check all the
    # rules to see which new analyses can be run.  e.g. if we've run
    # RepeatMasker we can run genscan. If we've run genscan we can run
    # blast jobs.
    #
    # All the analyses we're allowed to run are stored in a hash %analHash

    my %incomplete_accumulator_analyses;

    INPUT_ID_TYPE: foreach my $input_id_type (keys %$id_type_hash) {
	
        next INPUT_ID_TYPE if ($input_id_type eq 'ACCUMULATOR');
	
	if($batch_q_module->can('get_job_time')){
	  my @running_jobs = $job_adaptor->fetch_by_Status('RUNNING');
	  &job_time_check($batch_q_module, $verbose, \@running_jobs, 
			  $killed_file, $max_time);
	}

        @id_list = keys %{$id_type_hash->{$input_id_type}};

        @id_list = &shuffle(@id_list) if $shuffle;

        print "Checking $input_id_type ids\n";

        JOBID: foreach my $id (@id_list) {
    
            # handle signals. they are 'caught' in the handler subroutines but
            # it is only here they we do anything with them. so if the script
            # is doing something somewhere else (getting IDs at the start of
            # the while loop) or dozing, we have to wait until it gets here to
            # do anything...
    
            if ($alarm == 1) {
                $alarm = 0;
    
                # retry_failed_jobs($job_adaptor, $DEFAULT_RETRIES);
                while ($get_pend_jobs && !$term_sig && &$get_pend_jobs >= $MAX_PENDING_JOBS) {
                    sleep $overload_sleep;
                }
                alarm $wakeup;
            }
    
            if ($term_sig) {
                $done = 1;
                last INPUT_ID_TYPE;
            }
    
            if ($rst_sig) {
                $done = 0;
                @rules = $rule_adaptor->fetch_all;
                last INPUT_ID_TYPE;
            }
    
            my @anals = @{$sic->fetch_analysis_by_input_id($id)};
    
            my %analHash;
    
            # check all rules, which jobs can be started
    
            my @current_jobs = $job_adaptor->fetch_by_input_id($id);
    
            RULE: for my $rule (@rules)  {
                if (keys %analyses && ! defined $analyses{$rule->goalAnalysis->dbID}) {
                    if ($rule->goalAnalysis->input_id_type eq 'ACCUMULATOR') {
                        $incomplete_accumulator_analyses{$rule->goalAnalysis->logic_name} = 1;
                    }
                    next RULE;
                }
                print "Check ",$rule->goalAnalysis->logic_name, " - " . $id if $verbose;
    
                my $anal = $rule->check_for_analysis (\@anals, $input_id_type, \%completed_accumulator_analyses);
    
                if ($anal) {
                    print " fullfilled.\n" if $verbose;
                    if ($rule->goalAnalysis->input_id_type ne 'ACCUMULATOR') {
                      $analHash{$anal->dbID} = $anal;
                    }
                } else {
                    print " not fullfilled.\n" if $verbose;

                    if ($rule->goalAnalysis->input_id_type eq 'ACCUMULATOR' &&
                        $rule->has_condition_of_input_id_type($input_id_type) ) {

                        print " Makes ACCUMULATOR " . $rule->goalAnalysis->logic_name  . " incomplete\n" if($verbose);
                        $incomplete_accumulator_analyses{$rule->goalAnalysis->logic_name} = 1;
                    }
                }
            }
    
            # Now we loop over all the allowed analyses in the hash. We
            # first check the database to see if the job is already running.
            # If so we skip it.
            #
            # If all is ok we create a new job, store it in the database and
            # submit it to the batch runner. This will keep a check of the
            # number of jobs created. When $flushsize jobs have been stored
            # send to LSF.
    
            ANAL: for my $anal (values %analHash) {
                my $result_flag = run_if_new($id,
                                             $anal,
                                             \@current_jobs,
                                             $local,
                                             $verbose,
                                             $output_dir,
                                             $job_adaptor);

                if ($result_flag == -1) {
                    next JOBID;
                } elsif ($result_flag == 0) {
                    next ANAL;
                } else { 
                    $submitted++;
                }
            }
        }
    }

    if ( ! $done) {
        my @current_jobs = $job_adaptor->fetch_by_input_id('ACCUMULATOR');
        foreach my $accumulator_logic_name (keys %accumulator_analyses) {
            print "Checking accumulator type analysis $accumulator_logic_name\n" if $verbose;
            if (!exists($incomplete_accumulator_analyses{$accumulator_logic_name}) &&
                !exists($completed_accumulator_analyses{$accumulator_logic_name})) {
                my $result_flag = run_if_new('ACCUMULATOR',
                                             $accumulator_analyses{$accumulator_logic_name},
                                             \@current_jobs,
                                             $local,
                                             $verbose,
                                             $output_dir,
                                             $job_adaptor);
                if ($result_flag == 1 && $verbose) { print "Started accumulator type job for anal $accumulator_logic_name\n"; }
    
            } elsif (exists($incomplete_accumulator_analyses{$accumulator_logic_name})) {
                print "Accumulator type analysis $accumulator_logic_name conditions unsatisfied\n" if $verbose;
            } else {
                print "Accumulator type analysis $accumulator_logic_name already run\n" if $verbose;
            }
        }
    }

    &shut_down($db) if $done || $once;
    sleep($rerun_sleep) if $submitted == 0;
    @id_list = ();
    print "Waking up and run again!\n" if $verbose;
}


sub run_if_new {
    my ($id, $anal, $current_jobs, $local, $verbose, $output_dir, $job_adaptor) = @_;

    print "Checking analysis " . $anal->dbID . "\n\n" if $verbose;
    # Check whether it is already running in the current_status table?

    my $retFlag=0;
    eval {
        foreach my $cj (@$current_jobs) {

             print "Comparing to current_job " . $cj->input_id . " " .
                  $cj->analysis->dbID . " " .
                  $cj->current_status->status . " " .
                  $anal->dbID . "\n" if $verbose;

            if ($cj->analysis->dbID == $anal->dbID) {
                if ($cj->current_status->status eq 'FAILED' && $cj->retry_count <= $DEFAULT_RETRIES) {
                    $cj->batch_runRemote;
                    print "Retrying job\n";
                }
                else {
                    print "\nJob already in pipeline with status : " . $cj->current_status->status . "\n" if $verbose ;
                }
                $retFlag = 1;
            }
        }
    };


    if ($@) {
        print "ERROR: comparing to current jobs. Skipping analysis for " . $id . " [$@]\n";
        return -1;
    } elsif ($retFlag) {
        return 0;
    }

    my $job = Bio::EnsEMBL::Pipeline::Job->new(-input_id => $id,
                                               -analysis => $anal,
                                               -output_dir => $output_dir,
                                               -runner => $PIPELINE_RUNNER_SCRIPT);


    print "Store ", $id, " - ", $anal->logic_name, "\n" if $verbose;
    $job_adaptor->store($job);

    if ($local) {
        print "Running job locally\n" if $verbose;
        eval {
          $job->runLocally;
        };
        if ($@) {
           print STDERR "ERROR running job " . $job->dbID .  " "  . $job->stderr_file . "[$@]\n";
        }
    } else {
        eval {
            print "\tBatch running job\n" if $verbose;
            $job->batch_runRemote;
        };
        if ($@) {
            print STDERR "ERROR running job " . $job->dbID . " " . $job->stderr_file . " [$@]\n";
        }
    }
    return 1;
}

# remove 'lock'
sub shut_down {
    my ($db) = @_;

    my ($a_job) = $db->get_JobAdaptor->fetch_by_Status("CREATED");
    if ($a_job) {
        $a_job->flush_runs($db->get_JobAdaptor);
    }
    $db->pipeline_unlock;
    exit 0;
}


# handler for SIGTERM
sub termhandler {
    $term_sig = 1;
}


# handler for SIGUSR1
sub sighandler {

    $rst_sig = 1;
    $SIG{SIG1} = \&sighandler;
};


# handler for SIGALRM
sub alarmhandler {
    $alarm = 1;
    $SIG{ALRM} = \&alarmhandler;
}


# turn a name of the form 'ecs1a' into 'ecs1a.sanger.ac.uk'
sub qualify_hostname {
    my ($hostname) = @_;

    my $addr = gethostbyname($hostname);

    my $host = gethostbyaddr($addr, AF_INET);

    return $host;
}


sub retry_failed_jobs {
    my ($ja, $retry) = @_;

    my @failed_jobs = @{$ja->list_job_id_by_status('FAILED')};

    foreach my $jobId (@failed_jobs) {
        my $job = $ja->fetch_by_dbID($jobId);
        if ($job->retry_count <= $retry) {
            $job->batch_runRemote;

        }
    }
}


# NB Incomplete
# Delete RUNNING jobs which are not known to LSF.
sub check_lsf_status {
    my ($LSF_params) = @_;

    my $queue = $LSF_params->{'queue'};
    my (%lsf_id, %jobs);
    my $ja = $db->get_JobAdaptor;

    my $age = 300;   # 5 hours

    my $lsf = Bio::EnsEMBL::Pipeline::LSF->new(
        -queue => $queue,
        # -user  => $username,
    );

    foreach my $job ($lsf->get_all_jobs()) {
        my $stat = $job->status;
        my $id   = $job->id;
        $jobs{$stat}++;
        $lsf_id{$id} = 1;
    }

    foreach my $job ($db->fetch_by_Age($age)) {
        my $lsf = $job->submission_id;
        local $_ = $job->get_last_status->status;   # Are you local?
        STAT: {
            /RUNNING/ && do {
                #$job->remove unless defined $lsf_id{$lsf};
                last STAT;
            };
            /SUBMITTED/ && do {
                last STAT;
            };
        } # STAT
    }
}


sub shuffle {
    my (@in) = @_;
    my @out;

    srand;

    push @out, splice(@in, rand @in, 1) while @in;

    return @out;
}


sub config_sanity_check {
    my $ok = 1;
    no strict 'vars';
    print STDERR "checking config sanity\n";
    unless ($QUEUE_MANAGER) {
        print "Need to specify QUEUE_MANAGER in Config/BatchQueue.pm\n";
	$ok = 0;
    }
    unless ($LIB_DIR) {
        print "Need to specify LIB_DIR in Config/General.pm\n";
	$ok = 0;
    }
    unless ($DATA_DIR) {
        print "Need to specify DATA_DIR in Config/General.pm\n";
	$ok = 0;
    }
    unless ($BIN_DIR) {
        print "Need to specify BIN_DIR in Config/General.pm\n";
	$ok = 0;
    }

    return $ok;
}


sub logic_name2dbID {
    my ($ana_adaptor, @analyses) = @_;
    my %analyses;

    foreach my $ana (@analyses) {
        if ($ana =~ /^\d+$/) {
            $analyses{$ana} = 1;
        }
        else {
	    my $id = $ana_adaptor->fetch_by_logic_name($ana)->dbID;
	    if ($id) {
                $analyses{$id} = 1;
	    }
	    else {
	        print STDERR "Could not find analysis $ana\n";
	    }
        }
    }
    return %analyses;
}

sub db_sanity_check{
  my ($db) = @_;

  my ($sth, $query, $count);
  #check all rules in the rule_goal table have existing analyses
  $query = "SELECT count(distinct g.rule_id) "
          ."from rule_goal g "
          ."left join analysis a on g.goal = a.analysis_id "
	  ."where a.analysis_id is null";

  $sth = $db->prepare($query);
  $sth->execute;
  $count = $sth->fetchrow;
  if($count){
    die "Some of your goals in the rule_goal table don't seem ".
      "to have entries in the analysis table";
  }
  #check all rules in the rule_condition table have existing analyses
  $query = "select count(distinct c.rule_id) ".
    "from rule_conditions c ".
      "left join analysis a on c.condition = a.logic_name ".
	"where a.logic_name is null";

  $sth = $db->prepare($query);
  $sth->execute;
  $count = $sth->fetchrow;
  if($count){
    die "Some of your conditions in the rule_condition table don't ".
      "seem to have entries in the analysis table";
  }
  #check all the analyses have types
  $query = "select count(distinct(a.analysis_id)) ".
    "from analysis a left join input_id_type_analysis t ".
      "on a.analysis_id = t.analysis_id ".
	"where t.analysis_id is null";
  $sth = $db->prepare($query);
  $sth->execute;
  $count = $sth->fetchrow;
  if($count){
    die "Some of you analyses don't have entries in the ".
      " input_id_type_analysis table"; 
  }
  #check that all types which aren't accumulators have entries in
  #input__id_analysis table
  $query = "select distinct(t.input_id_type) ".
    "from input_id_type_analysis t ".
      "left join input_id_analysis i ".
	"on t.input_id_type = i.input_id_type ".
	  "where i.input_id_type is null ".
	    "and t.input_id_type != 'ACCUMULATOR'";
  
  $sth = $db->prepare($query);
  $sth->execute;
  $count = $sth->fetchrow;
  if($count){
    die "Some of you types  don't have entries in the ".
      " input_id_type_analysis table"; 
  }
}


sub job_time_check{
  my ($batch_q_module, $verbose, $running_jobs, $file, $time) = @_;
   open(KILLED, "+>>".$file) or die("couldn't open ".$file." $!");
  foreach my $job(@$running_jobs){
    my $time = $batch_q_module->get_job_time($job->submission_id);
    if($time >= $max_time){
      $batch_q_module->kill_job($job->submission_id);
      print KILLED $job->input_id." ".$job->analysis->logic_name." ".$job->analysis->module."\n";
      print STDERR $job->input_id." ".$job->analysis->logic_name." ".$job->analysis->module."\n" if($verbose);
      $job->set_status('KILLED');
    }
  }

}

sub useage{
	exec('perldoc', $0);
	exit;
}

=pod

=head1 NAME

monitor

=head1 SYNOPSIS

Pipeline Monitor Script

A Simple script using the Monitor.pm module to display information on the status of the pipeline.


=head1 OPTIONS

=head2 [DB Connection Details]

   -dbhost     The host where the pipeline database is.
   -dbport     The port.
   -dbuser     The user to connect as.
   -dbpass     The password to use.
   -dbname   The database name.

=head2 [Other Options]

   -local run the pipeline locally and not using LSF
   -idlist_file a path to a file containing a list of input ids to use
   -runner path to a runner script (if you want to overide the setting
				    in BatchQueue.pm)
   -output_dir path to a output dir (if you want to overide the 
				     setting in BatchQueue.pm)
   -once only run through the RuleManager loop once
   -shuffle before running though the loop shuffle the order of the 
    input ids
   -analysis only run with these analyses objects, can be logic_names or
    analysis ids
   -v verbose mode
   -dbsanity peform some db sanity checks, can be switched of with 
    -nodbsanity
   
-h or -help will print out the help again

=head1 EXAMPLES




=head1 SEE ALSO

  Bio::EnsEMBL::Pipeline::Job
  Bio::EnsEMBL::Pipeline::Config::General
  Bio::EnsEMBL::Pipeline::Config::BatchQueue

and also using_the_ensembl_pipeline.txt in the ensembl-docs cvs module

=cut
