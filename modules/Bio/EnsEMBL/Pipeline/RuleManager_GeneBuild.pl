# Script for operating the analysis pipeline
#
# You may distribute this code under the same terms as perl itself


use strict;
use Getopt::Long;
use Sys::Hostname;
use Socket;

use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;


use Bio::EnsEMBL::Pipeline::Config::GeneBuild::General;
use Bio::EnsEMBL::Pipeline::Config::General;
use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Scripts;
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
my $alarm    =  0;
my $wakeup   =  120;   # period to check batch queues; set to 0 to disable

# the signal handlers
# $SIG{USR1} = \&sighandler;
$SIG{ALRM} = \&alarmhandler;
$SIG{TERM} = \&termhandler;


# dynamically load appropriate queue manager (e.g. LSF)

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

my $dbhost    = 'ecs1c';
my $dbname    = 'pipeline_genebuild_test';
my $dbuser    = 'ensadmin';
my $dbpass    = 'ensembl';

$| = 1;

my $chunksize    = 1000000;   # How many Input_ids to fetch at one time
my $currentStart = 0;       # Running total of job ids
my $completeRead = 0;       # Have we got all the input ids yet?
my $local        = 0;       # Run failed jobs locally
my @analyses;               # Only run this analysis ids
my $submitted;
my $idlist;
my ($done, $once);
my $runner;
my $shuffle;  
my $output_dir;
my @start_from;
my %analyses;

GetOptions(
    'dbhost=s'        => \$dbhost,
    'dbname=s'      => \$dbname,
    'dbuser=s'      => \$dbuser,
    'dbpass=s'      => \$dbpass,
    'local'         => \$local,
    'idlist=s'      => \$idlist,
    'runner=s'      => \$runner,
    'output_dir=s'  => \$output_dir,
    'once!'         => \$once,
    'shuffle!'      => \$shuffle,
    'start_from=s@' => \@start_from,
    'analysis=s@'   => \@analyses
)
or die ("Couldn't get options");

#print STDERR "dbhost ".$dbhost." dbuse ".$dbuser." dbname ".$dbname."\n";

unless ($dbhost && $dbname && $dbuser) {
    die "Must specify database with -dbhost, -dbname, -dbuser and -dbpass $! \n";
    exit 1;
}


my $RUNNER_SCRIPT = $GB_RUNNER || $runner;




my $db = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new(
    -host   => $dbhost,
    -dbname => $dbname,
    -user   => $dbuser,
    -pass   => $dbpass,
);

my $rule_adaptor = $db->get_RuleAdaptor;
my $job_adaptor  = $db->get_JobAdaptor;
my $ana_adaptor  = $db->get_AnalysisAdaptor;
my $sic          = $db->get_StateInfoContainer;


# analysis options are either logic names or analysis dbID's

unless (@start_from) {
    print "Need to specify at least 1 analysis to start from with -start_from\n";
    exit 1;
}


%analyses   = logic_name2dbID($ana_adaptor, @analyses);
my %tmp = logic_name2dbID($ana_adaptor, @start_from);
@start_from = keys %tmp;

if ($idlist && ! -e $idlist) {
    die "Must be able to read $idlist";
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

# Need here to strip rules which don't need to be run.

my @id_list;     # All the input ids to check

alarm $wakeup if $wakeup;  
                # Signal to the script to do something in the future,
                # Such as check for failed jobs, no of jobs in queue

while (1) {
    $submitted = 0;

    # This loop reads input ids from the database a chunk at a time
    # until we have all the input ids.

    if($batch_q_module->can('get_job_time')){
      my @running_jobs = $job_adaptor->fetch_by_Status('RUNNING');
      open(KILLED, "+>>".$GB_KILLED_INPUT_IDS) or die("couldn't open ".$GB_KILLED_INPUT_IDS." $!");
      foreach my $job(@running_jobs){
	my $time = $batch_q_module->get_job_time($job->submission_id);
	if($time >= $GB_MAX_JOB_TIME){
	  $batch_q_module->kill_job($job->submission_id);
	  print KILLED $job->input_id." ".$job->analysis->logic_name." ".$job->analysis->module."\n";
	  $job->set_status('KILLED');
	}
      }
    }
    if (defined $idlist) {

        # read list of id's from a file,
        # e.g. a list of contigs on the golden path

        open IDS, "< $idlist" or die "Can't open $idlist";
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

        if (!$completeRead) {
            #print "Reading IDs ... ";

	    foreach my $a (@start_from) {
		push @id_list, @{$sic->list_input_id_by_Analysis($a)};
	    }

            $completeRead = 1;
        }
    }

    @id_list = &shuffle(@id_list) if $shuffle;

    # Now we loop over all the input ids. We fetch all the analyses
    # from the database that have already run. We then check all the
    # rules to see which new analyses can be run.  e.g. if we've run
    # RepeatMasker we can run genscan. If we've run genscan we can run

    # blast jobs.
    #
    # All the analyses we're allowed to run are stored in a hash %analHash

    JOBID: while (@id_list) {
        my $id = shift @id_list;

        # handle signals. they are 'caught' in the handler subroutines but
        # it is only here they we do anything with them. so if the script
        # is doing something somewhere else (getting IDs at the start of
        # the while loop) or dozing, we have to wait until it gets here to
        # do anything...

        if ($alarm == 1) {
            $alarm = 0;

            # retry_failed_jobs($job_adaptor, $DEFAULT_RETRIES);
            while ($get_pend_jobs && !$term_sig && &$get_pend_jobs >= $MAX_PENDING_JOBS) {
                sleep 600;
            }
            alarm $wakeup;
        }

	if ($term_sig) {
	    $done = 1;
	    last JOBID;
	}

        my @anals = @{$sic->fetch_analysis_by_input_id($id)};

        my %analHash;

        # check all rules, which jobs can be started

        my @current_jobs = $job_adaptor->fetch_by_input_id($id);

        RULE: for my $rule (@rules)  {
            if (keys %analyses && ! defined $analyses{$rule->goalAnalysis->dbID}) {
               next RULE;
            }
            #print "Check ",$rule->goalAnalysis->logic_name, " - " . $id;

            my $anal = $rule->check_for_analysis (@anals);

            if ($anal) {

             #   print " fullfilled.\n";
                $analHash{$anal->dbID} = $anal;
            } else {
             #   print " not fullfilled.\n";
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
            # print "Checking analysis " . $anal->dbID . "\n\n";
            # Check whether it is already running in the current_status table?

            eval {
                foreach my $cj (@current_jobs) {

                    # print "Comparing to current_job " . $cj->input_id . " " .
                           $cj->analysis->dbID . " " .
                           $cj->current_status->status . " " .
                           $anal->dbID . "\n";

                    if ($cj->analysis->dbID == $anal->dbID) {
                        if ($cj->current_status->status eq 'FAILED' && $cj->retry_count <= $DEFAULT_RETRIES) {
                            $cj->batch_runRemote;
                            # print "Retrying job\n";
                        }
                        else {
                            # print "\nJob already in pipeline with status : " . $cj->current_status->status . "\n";
                        }
                        next ANAL;
                    }
                }
            };


            if ($@) {
                print "ERROR: comparing to current jobs. Skipping analysis for " . $id . " [$@]\n";
                next JOBID;
              }

            my $job = Bio::EnsEMBL::Pipeline::Job->new(-input_id => $id,
						       -analysis => $anal,
						      -output_dir => $output_dir);


            #print "Store ", $id, " - ", $anal->logic_name, "\n";
            $submitted++;
            $job_adaptor->store($job);

            if ($local) {
                print "Running job locally\n";
                eval {
                  $job->runLocally;
                };
                if ($@) {
                   print STDERR "ERROR running job " . $job->dbID .  " "  . $job->stderr_file . "[$@]\n";
                }
            } else {
              eval {
                print "\tBatch running job\n";
                $job->batch_runRemote;
              };
              if ($@) {
                print STDERR "ERROR running job " . $job->dbID . " " . $job->stderr_file . " [$@]\n";
              }
	      sleep($GB_RULEMANAGER_SLEEP);
            }

        }

    }
    &shut_down($db) if $done || $once;
    if ($completeRead == 1) {
        sleep(3600) if $submitted == 0;
        $completeRead = 0;
        $currentStart = 0;
        @id_list = ();
        #print "Waking up and run again!\n";
    }

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
# sub sighandler {

    # SIG: {
    # }
    # $SIG{SIG1} = \&sighandler;
# };


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
