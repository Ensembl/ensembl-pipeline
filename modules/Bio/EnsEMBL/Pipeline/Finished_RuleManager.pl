#!/usr/local/ensembl/bin/perl -w
#
# Script for operating the analysis pipeline
#
# You may distribute this code under the same terms as perl itself

BEGIN{
    $ENV{'BLASTDB'}     = '/data/blastdb/Finished';
    $ENV{'BLASTMAT'}    = '/usr/local/ensembl/data/blastmat';
    $ENV{'BLASTFILTER'} = '/usr/local/ensembl/bin';
}

use strict;
use Getopt::Long;
use Sys::Hostname;
use Socket;
use POSIX qw(:signal_h setsid);
use FindBin ();
use File::Basename ();
use File::Spec::Functions;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Data::Dumper;
use Bio::EnsEMBL::Pipeline::Config::Blast;
use Bio::EnsEMBL::Pipeline::Config::General;
use Bio::EnsEMBL::Pipeline::Config::BatchQueue;
use IO::Socket;
use Data::Dumper; 

$| = 1; # localise?

my $SCRIPT                = File::Basename::basename($0);
my $SELF                  = catfile $FindBin::Bin, $SCRIPT;
my $DEFAULT_HOME          = '';
my $HOME                  = $DEFAULT_HOME;
my $HOSTNAME              = Sys::Hostname::hostname();
my $PID_FILE              = $DEFAULT_HOME . "$SELF.$HOSTNAME.pid";
my $DEBUG                 = 1;
my $NAMED                 = 'FINISHED_PIPELINE';
my $OUT                   = $ENV{$NAMED."_OUTFILE"};

my $OVERLOAD_SLEEP        = 10;  # THIS SHOULD BE SIG_BLOCKED [3600]
my $PAUSE                 = 6;     # THIS DOESN'T NEED TO BE    [60]
my $GET_PEND_JOBS         = undef; # subroutine to check number of pending jobs in queue

# Also, an alarm is set to go off every (say) 2 minutes (or not at all if
# $wakeup is false). When the script receives this it looks to see how many
# jobs are in the queue and, if there too many, goes to sleep for a while.
# You could also do this by counting the number of distinct(LSF_id) entries
# in the job table but you would have to also check for dead processes etc.
# It looks for failed jobs as well and restarts them.

# signal parameters
my $term_sig = 0;
my $rst_sig  = 0;
my $alarm    = 0;
my $wakeup   = 6;   # period to check batch queues; set to 0 to disable

# the signal handlers

# handler for SIGTERM
sub sigTERMhandler {
    print STDERR "Handling a TERM/INT SIGNAL\n";
    $term_sig = 1;
}
# handler for SIGUSR1
sub sigUSR1handler {
    print STDERR "Handling a USR1 SIGNAL\n";
    $rst_sig = 1;
    $SIG{SIG1} = \&sigUSR1handler;
};
# handler for SIGUSR2
sub sigUSR2handler {
    print STDERR "Handling a USR2 SIGNAL\n";
    &pause(1);
};
# handler for SIGALRM
sub sigALRMhandler {
    print STDERR "Handling a ALRM SIGNAL\n";
    &check_not_overloaded;
    alarm $wakeup if $wakeup;
    $SIG{ALRM} = \&sigALRMhandler;
}

# POSIX unmasks the sigprocmask properly
my $sigset1       = POSIX::SigSet->new();
my $sigHUPaction  = POSIX::SigAction->new('sigHUP_handler', $sigset1, &POSIX::SA_NODEFER);
POSIX::sigaction(&POSIX::SIGHUP, $sigHUPaction)   or warn "Couldn't $sigHUPaction";

my $sigset2       = POSIX::SigSet->new();
my $sigUSR1action = POSIX::SigAction->new('sigUSR1handler', $sigset2, &POSIX::SA_NODEFER);
POSIX::sigaction(&POSIX::SIGUSR1, $sigUSR1action) or warn "Couldn't $sigUSR1action";

my $sigset3       = POSIX::SigSet->new();
my $sigUSR2action = POSIX::SigAction->new('sigUSR2handler', $sigset3, &POSIX::SA_NODEFER);
POSIX::sigaction(&POSIX::SIGUSR2, $sigUSR2action) or warn "Couldn't $sigUSR2action";

my $sigset4       = POSIX::SigSet->new();
my $sigALRMaction = POSIX::SigAction->new('sigALRMhandler', $sigset4);#, &POSIX::SA_NODEFER);
POSIX::sigaction(&POSIX::SIGALRM, $sigALRMaction) or warn "Couldn't $sigALRMaction";

my $sigset5       = POSIX::SigSet->new();
my $sigTERMaction = POSIX::SigAction->new('sigTERMhandler', $sigset5);
POSIX::sigaction(&POSIX::SIGTERM, $sigTERMaction) or warn "Couldn't $sigTERMaction";

my $sigset6       = POSIX::SigSet->new();
my $sigINTaction  = POSIX::SigAction->new('sigTERMhandler', $sigset6);
POSIX::sigaction(&POSIX::SIGINT, $sigINTaction)   or warn "Couldn't $sigINTaction";

sigprocmask(SIG_BLOCK, $sigset6) or die "Can't block $sigset6"; # blocks TERM

unless (&config_sanity_check) {
    exit 1;
}

# command line options override 
# anything set in the environment variables.

my $dbhost         = $ENV{'ENS_DBHOST'};
my $dbname         = $ENV{'ENS_DBNAME'};
my $dbuser         = $ENV{'ENS_DBUSER'};
my $dbpass         = $ENV{'ENS_DBPASS'};
my $dbport         = $ENV{'ENS_DBPORT'} || 3306;

my $local          = 0;                # Run failed jobs locally
my @analyses       = ();               # Only run these analysis ids
my $submitted      = undef;
my $idlist         = undef;
my ($done, $once)  = (0) x 2;
my $no_run         = 0; # makes it only create jobs but doesn't submit to LSF
my $runner         = undef;
my $shuffle        = 0;  
my $output_dir     = undef;
my @start_from     = ();
my @assembly_types = ();
my $db_sanity      = 1;
my %analyses       = ();
my $priority       = 0;
my $verbose        = 1;
my $rerun_sleep    = 100;

my ($start, $stop, $restart, $status, $daemonFree, $pause) = (0) x 6;
my @controls = ('start'     => \$start,
                'stop'      => \$stop,
                'restart'   => \$restart,
                'status'    => \$status,
                'buffyMode' => \$daemonFree,
		'pause'     => \$pause);
restart_string(@ARGV);

GetOptions(
    'dbhost=s'      => \$dbhost,
    'dbname=s'      => \$dbname,
    'dbuser=s'      => \$dbuser,
    'dbpass=s'      => \$dbpass,
    'dbport=s'      => \$dbport,
    'local'         => \$local,
    'idlist=s'      => \$idlist,
    'runner=s'      => \$runner,
    'output_dir=s'  => \$output_dir,
    'once!'         => \$once,
    'shuffle!'      => \$shuffle,
    'norun'         => \$no_run,
    'start_from=s@' => \@start_from,
    'assembly=s@'   => \@assembly_types,
    'analysis=s@'   => \@analyses,
    'dbsanity!'     => \$db_sanity,
    'v!'            => \$verbose,
    'priority'      => \$priority,
    @controls
)
or die ("Couldn't get options");

unless ($dbhost && $dbname && $dbuser) {
    print STDERR "Must specify database with -dbhost, -dbname, -dbuser and -dbpass\n";
    exit 1;
}
# -------------------------------------
# Handle daemon control all sub routines called must exit!
# check uniqueness of control parameters
if($start + $stop + $status + $restart + $daemonFree + $pause > 1){ useage(); }
# decide which one
elsif($start)      { daemonize($OUT); }
elsif($stop)       { stop_daemon();   }
elsif($status)     { daemon_status(); }
elsif($restart)    { restart();       }
elsif($pause)      { pause();         }
elsif($daemonFree) { print "$NAMED : Odd I'm not running as a daemon. Hit Ctrl-C quick\n" unless eval{ get_pid() }; }
else               { useage();        } # Catch when no daemon control specified

sigprocmask(SIG_UNBLOCK, $sigset6);

# -------------------------------------
# this is now the daemon's code

# dynamically load appropriate queue manager (e.g. LSF)

my $batch_q_module = "Bio::EnsEMBL::Pipeline::BatchSubmission::$QUEUE_MANAGER";

my $file = "$batch_q_module.pm";
$file =~ s{::}{/}g;
eval {
    require "$file";
};
if ($@) {
    print STDERR "Error trying to load $batch_q_module;\ncan't find $file : $@\n";
    exit 1;
}

if ($batch_q_module->can("get_pending_jobs")) {
    my $f = $batch_q_module . "::get_pending_jobs";
    $GET_PEND_JOBS = \&$f;
}

sub check_not_overloaded{
    warn "$NAMED : Checking Queue Manager status\n";
    # retry_failed_jobs($job_adaptor, $DEFAULT_RETRIES);
    my $queue = {
        '-USER'  => 'rds', 
        '-QUEUE' => $DEFAULT_BATCH_QUEUE, 
        '-DEBUG' => 0
        };
    #while ((ref($GET_PEND_JOBS) eq 'CODE') && !$term_sig && (&$GET_PEND_JOBS($queue) >= $MAX_PENDING_JOBS)) {
    while ($GET_PEND_JOBS && !$term_sig && &$GET_PEND_JOBS($queue) >= $MAX_PENDING_JOBS) {
	&pause(1);
    }
}

my $DBADAPTOR = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new(
    -host   => $dbhost,
    -dbname => $dbname,
    -user   => $dbuser,
    -pass   => $dbpass,
    -port   => $dbport,
);
if($db_sanity){
    &db_sanity_check($DBADAPTOR);
}else{
    print STDERR "$NAMED : NO DB SANITY CHECK.  OK?\n" if($verbose);
}

################################################################################
##### SIG HUP HANDLER
my @restart = restart_array();

sub sigHUP_handler {
    print STDERR "\n";
    print STDERR "$NAMED : got SIGHUP @ ".scalar(localtime())."\n";
    print STDERR "$NAMED : exec'ing: $SELF @restart\n";
    $DBADAPTOR->pipeline_unlock() if $DBADAPTOR;
    exec($SELF, @restart) or die "Couldn't restart: $!\n";
}
################################################################################

my $rule_adaptor = $DBADAPTOR->get_RuleAdaptor;
my $job_adaptor  = $DBADAPTOR->get_JobAdaptor;
my $ana_adaptor  = $DBADAPTOR->get_AnalysisAdaptor;
my $sic          = $DBADAPTOR->get_StateInfoContainer;

# analysis options are either logic names or analysis dbID's
unless (@start_from) {
    print "Need to specify at least 1 analysis to start from with -start_from\n";
    exit 1;
}


%analyses   = logic_name2dbID($ana_adaptor, @analyses);
my %tmp     = logic_name2dbID($ana_adaptor, @start_from);
@start_from = keys %tmp;


die "Must be able to read $idlist"          if ($idlist && ! -e $idlist);

alarm $wakeup if $wakeup;
#my $i = 0;
#while($i < 10000){
#    print "$i \n";
#    &pause(1) unless($i % 100);
#    if($term_sig){
#	print STDERR "$NAMED : got SIGNAL\n";
#	exit;
#    }
#    sleep(1);
#    $i++;
#}

#exit();



# Lock to prevent >1 instances of the same pipeline running together
#
# Puts an entry in the meta table keyed on 'pipeline.lock', which
# gives the host running the script, user, PID and when it was
# started

if (my $lock_str = $DBADAPTOR->pipeline_lock && !$no_run) {
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
$DBADAPTOR->pipeline_lock($lock_str);


# Fetch all the analysis rules.  These contain details of all the
# analyses we want to run and the dependences between them. e.g. the
# fact that we only want to run blast jobs after we've repeat masked etc.

my @rules       = $rule_adaptor->fetch_all;

# Need here to strip rules which don't need to be run.

my @id_list;     # All the input ids to check

alarm $wakeup if $wakeup;  
                # Signal to the script to do something in the future,
                # Such as check for failed jobs, no of jobs in queue

{
    my $no_of_loops = {};
    my $max_loops   = 10;
    sub flush_created_if_they_might_never_get_done{
	my ($submitteds, $db) = @_;
	my $overall = 0;
	foreach my $k(keys(%$submitteds)){
	    my $a = $submitteds->{$k} || 0;
	    warn "Checking analysis $k, which has $a submitted this loop";
	    if($a == 0){
		$no_of_loops->{$k} ||= 1;
		warn $no_of_loops->{$k} . " of $max_loops loops with no submits passed for analysis $k";
		$no_of_loops->{$k}++;
		$overall++;
	    }
	    if($no_of_loops->{$k} > $max_loops){
		warn "max_loops reached (".$no_of_loops->{$k}.") for $k. flushing runs.";
		my $job_a   = $db->get_JobAdaptor();
		my ($a_job) = $job_a->fetch_by_Status("CREATED");
		if ($a_job) {
                    my $sigset2 = POSIX::SigSet->new(&POSIX::SIGALRM);
                    sigprocmask(SIG_BLOCK, $sigset2) or die "Can't block SIGALRM for flush runs: $!\n";
                    $a_job->flush_runs($job_a);
		    $overall--;
                    sigprocmask(SIG_UNBLOCK, $sigset2) or die "Can't unblock SIGALRM after flush runs: $!\n";
		}
		$no_of_loops->{$k} = 0;
	    }
	}
	return $overall > 0 ? 1 : 0;
    }
}

while (1) {
    
    $submitted = {};

    # This loop reads input ids from the database a chunk at a time
    # until we have all the input ids.

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

        print "Reading IDs ... ";

        foreach my $a (@start_from) {
	    push @id_list, @{$sic->list_input_id_by_Analysis_assembly_type_priority($a, \@assembly_types, $priority)};
	}
#	print scalar(@id_list);
	#&Bio::EnsEMBL::Pipeline::DBSQL::StateInfoContainer::increment_priority();
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
	    &check_not_overloaded();
            alarm $wakeup;
        }

	if ($term_sig) {
	    $done = 1;
	    last JOBID;
	}

	if ($rst_sig) {
	    $rst_sig = 0;
            @rules = $rule_adaptor->fetch_all;
	    last JOBID;
	}

        my @anals = @{$sic->fetch_analysis_by_input_id($id)};

        my %analHash;

        # check all rules, which jobs can be started

        my @current_jobs = $job_adaptor->fetch_by_input_id($id);

        RULE: for my $rule (@rules)  {
            # what does this check ????
            if (keys %analyses && ! defined $analyses{$rule->goalAnalysis->dbID}) {
               next RULE;
            }
            print "Check ",$rule->goalAnalysis->logic_name, " - " . $id if $verbose;
            # hard wiring makes devil happy.
            my $anal = $rule->check_for_analysis (\@anals, 'CONTIG');

            if (UNIVERSAL::isa($anal,'Bio::EnsEMBL::Pipeline::Analysis')) { # checks whether its a reference[obj] maybe should be UNIVERSAL::isa() check
                print " fullfilled. got a ".ref($anal)."\n" if $verbose;
                $analHash{$anal->dbID} = $anal;
            } elsif(!($anal & 1)) { # checks if its returned a value passing Input_Id_Type Check.
                # now do some more bit checking
                # only interesting to see if condition check[4] NOT failed && complete check[2] failed
                # translates to  !($anal & 4) && ($anal & 2)
                if(!($anal & 4) && ($anal & 2) && !$sic->check_is_current($rule->goalAnalysis->dbID, $id)){
                    print " fullfilled **** IS_CURRENT FALSE ****\n" if $verbose;
                    $analHash{$rule->goalAnalysis->dbID} = $rule->goalAnalysis();
		}elsif($anal & 4){
		    #&Bio::EnsEMBL::Pipeline::DBSQL::StateInfoContainer::reset_priority();
		    print " not fullfilled *** is_current true BUT CONDITION NOT MET ***\n" if $verbose;
                }else{
                    print " not fullfilled **** is_current true ****\n" if $verbose;
                }
            } else {
                print " not fullfilled.\n" if $verbose;
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
	    $submitted->{$anal->logic_name} ||= 0;
            my $submittable = 1;
            eval {
                foreach my $cj (@current_jobs) {

                    #print "Comparing to current_job " . $cj->input_id . " " .
                    #$cj->analysis->dbID . " " .
                    #$cj->current_status->status . " " .
                    #$anal->dbID . "\n";

                    if ($cj->analysis->dbID == $anal->dbID) {
                        if ($cj->current_status->status eq 'FAILED' && $cj->retry_count <= $DEFAULT_RETRIES) {
                            my $sigset2 = POSIX::SigSet->new(&POSIX::SIGALRM);
                            sigprocmask(SIG_BLOCK, $sigset2) or die "Can't block SIGALRM for batch_runRemote: $!\n";
                            $cj->batch_runRemote;
                            sigprocmask(SIG_UNBLOCK, $sigset2) or die "Can't unblock SIGALRM after batch_runRemote: $!\n";
                            # print "Retrying job\n";
                        }elsif($cj->current_status->status eq 'VOID'){
                            $submittable = 0;
                            my $input_id = $cj->input_id();
                            my $analysis = $cj->analysis();
                            my $dbID     = $cj->dbID();
                            eval{
                                $sic->store_input_id_analysis($input_id,
                                                              $analysis,
                                                              'localhost'
                                                              );
                            };
                            if($@){
                                print STDERR "Error updating 'VOID' job $dbID :\n[$@]\n";
                            }else{
                                print STDERR " *** Updated 'VOID' job <$dbID> *** \n";
                            }
                        }
                        else {
                            $submittable = 0;
                            # print "\nJob already in pipeline with status : " . $cj->current_status->status . "\n";
                        }
                        #next ANAL;
                    }
                }
            };
            next ANAL unless $submittable;

            if ($@) {
                print "ERROR: comparing to current jobs. Skipping analysis for " . $id . " [$@]\n";
                next JOBID;
              }

            my $job = Bio::EnsEMBL::Pipeline::Job->new(-input_id => $id,
						       -analysis => $anal,
						       -output_dir => $output_dir,
						       -runner => $PIPELINE_RUNNER_SCRIPT);


            print "Store ", $id, " - ", $anal->logic_name, "\n" if $verbose;
            $submitted->{$anal->logic_name}++;
            $job_adaptor->store($job);

            my $sigset2 = POSIX::SigSet->new(&POSIX::SIGALRM);
            sigprocmask(SIG_BLOCK, $sigset2) or die "Can't block SIGALRM for batch_runRemote: $!\n";
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
                $job->batch_runRemote unless $no_run;
		#&Bio::EnsEMBL::Pipeline::DBSQL::StateInfoContainer::reset_priority();
              };
              if ($@) {
                print STDERR "ERROR running job " . $job->dbID . " " . $job->stderr_file . " [$@]\n";
              }
            }
            sigprocmask(SIG_UNBLOCK, $sigset2) or die "Can't unblock SIGALRM after batch_runRemote: $!\n";
            
        }

    }

    &shut_down($DBADAPTOR) if $done || $once;
    my $this_time;# = flush_created_if_they_might_never_get_done($submitted,$DBADAPTOR);
    my $slept;
    if($this_time){ 
        # not much was submitted, probably because
	# there was too many unfulfilled dependancies
	# wait while the jobs finish
	local $SIG{ALRM} = 'IGNORE';
	warn "time to sleep for $rerun_sleep seconds... IGNORING \$SIG{ALRM}";
	$slept = sleep $rerun_sleep;
    }else{
	warn "time to sleep for $rerun_sleep seconds...";
	$slept = sleep $rerun_sleep;
    }
    @id_list = ();
    print STDERR "After sleeping for $slept . Waking up and run again!\n" if $verbose;
}

# --------------------------------------------------------
# Subroutines
# --------------------------------------------------------
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
# sub check_lsf_status {
#     my ($LSF_params) = @_;

#     my $queue = $LSF_params->{'queue'};
#     my (%lsf_id, %jobs);
#     my $ja = $db->get_JobAdaptor;

#     my $age = 300;   # 5 hours

#     my $lsf = Bio::EnsEMBL::Pipeline::LSF->new(
#         -queue => $queue,
#         # -user  => $username,
#     );

#     foreach my $job ($lsf->get_all_jobs()) {
#         my $stat = $job->status;
#         my $id   = $job->id;
#         $jobs{$stat}++;
#         $lsf_id{$id} = 1;
#     }

#     foreach my $job ($db->fetch_by_Age($age)) {
#         my $lsf = $job->submission_id;
#         local $_ = $job->get_last_status->status;   # Are you local?
#         STAT: {
#             /RUNNING/ && do {
#                 #$job->remove unless defined $lsf_id{$lsf};
#                 last STAT;
#             };
#             /SUBMITTED/ && do {
#                 last STAT;
#             };
#         } # STAT
#     }
# }


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
    print STDERR "$NAMED : Checking config sanity\n";
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
sub db_sanity_check{
  my ($db) = @_;

  my ($sth, $query, $msg);
  #check all rules in the rule_goal table have existing analyses
  $query = qq{SELECT COUNT(DISTINCT g.rule_id)
                FROM rule_goal g
                LEFT JOIN analysis a ON g.goal = a.analysis_id
	        WHERE a.analysis_id IS NULL};
  $msg = "Some of your goals in the rule_goal table don't seem".
         " to have entries in the analysis table";
  execute_sanity_check($db, $query, $msg);
  #check all rules in the rule_condition table have existing analyses
  $query = qq{SELECT COUNT(DISTINCT c.rule_id)
                FROM rule_conditions c
                LEFT JOIN analysis a ON c.condition = a.logic_name || c.condition IS NULL
	        WHERE a.logic_name IS NULL};
  $msg = "Some of your conditions in the rule_condition table don't" .
         " seem to have entries in the analysis table\n $query";
  execute_sanity_check($db, $query, $msg);
  #check all the analyses have types
  $query = qq{SELECT COUNT(DISTINCT(a.analysis_id))
                FROM analysis a
                LEFT JOIN input_id_type_analysis t ON a.analysis_id = t.analysis_id
	        WHERE t.analysis_id IS NULL};
  $msg = "Some of your analyses don't have entries in the".
         " input_id_type_analysis table"; 
  execute_sanity_check($db, $query, $msg);
  #check that all types which aren't accumulators have entries in
  #input__id_analysis table
  $query = qq{SELECT DISTINCT(t.input_id_type)
                FROM input_id_type_analysis t
                LEFT JOIN input_id_analysis i ON t.input_id_type = i.input_id_type
	        WHERE i.input_id_type IS NULL
                && t.input_id_type != 'ACCUMULATOR'};
  $msg = "Some of your types don't have entries in the".
         " input_id_type_analysis table";
  execute_sanity_check($db, $query, $msg);
}
sub execute_sanity_check{
    my ($db, $query, $msg) = @_;
    my $sth = $db->prepare($query);
    $sth->execute();
    die $msg if $sth->fetchrow();
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

# ---------------------------------------------
# for the daemon control
sub daemonize {
    my $output = shift || '/dev/null';
    if(! -e $PID_FILE){
        print STDERR "$NAMED : Starting up Rule Manager on $HOSTNAME....\n"    if $DEBUG;
        print STDERR "$NAMED : Using database configured from $HOME\n"         if $DEBUG;
        print STDERR "$NAMED : PID file at $PID_FILE\n"                        if $DEBUG;
        chdir '/'		    or die "Can't chdir to /: $!";
        open (STDIN, '/dev/null')   or die "Can't read /dev/null: $!";
        open (STDOUT, ">$output")   or die "Can't write to $output: $!";
        defined(my $pid = fork)	    or die "Can't fork: $!";
        if($pid){
            open(PID,">$PID_FILE")  or die "Can't write to $PID_FILE: $!";
            print PID $pid;
            close PID;
            exit if $pid;
        }
        setsid			    or die "Can't start a new session: $!";
        open (STDERR, '>&STDOUT')   or die "Can't dup stdout: $!";
    }else{
        print STDERR "$NAMED : $0 is already running as a daemon using pid file: <$PID_FILE>\n";
        exit;
    }
}  
sub stop_daemon{
    if(my $pid = eval{get_pid()}){
        kill 9, $pid;
        unlink $PID_FILE;
        my $db;
        eval{ 
            $db = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new(
                                                                -host   => $dbhost,
                                                                -dbname => $dbname,
                                                                -user   => $dbuser,
                                                                -pass   => $dbpass,
                                                                -port   => $dbport,
                                                            );
            
        };
        if($@){
            print STDERR $@;
            print STDERR "\nTo remove the lock  from the database I need host port name etc...\n";
        }else{
            if(my @lock = split("[\@:]", $db->pipeline_lock())){
                unless($pid ne $lock[2]){
                    $db->pipeline_unlock();
                    print STDERR "$NAMED : killed PID <$lock[2]> and removed pipeline.lock\n";
                    exit;
                }
                print STDERR "PID: <$pid> from file does not match PID: <$lock[2]> from db\n";
            }else{
                print STDERR "Couldn't find the lock. Does this process exist?\n";
            }
        }
    }else{
        print STDERR "$NAMED : Instance using pid file <$PID_FILE> is not currently running\n";
    }
    exit;
}
sub daemon_status{
    my $checked;
    if(my $pid = eval{get_pid()}){
	if($checked = check_pid($pid)){
	    print STDERR "$NAMED : Instance running under PID: $pid\n";
        }else{
	    print STDERR "$NAMED : Instance NOT running under PID: $pid on this host check below\n";
	}
        print STDERR "$NAMED : Checking the database...\n";
        my $db;
        eval{ 
            $db = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new(
                                                                -host   => $dbhost,
                                                                -dbname => $dbname,
                                                                -user   => $dbuser,
                                                                -pass   => $dbpass,
                                                                -port   => $dbport,
                                                            );
            
        };
        if($@){
            print STDERR $@;
            print STDERR "\nTo check the databases for a lock I need host port name etc...\n";
        }else{
            if(my @lock = split("[\@:]", $db->pipeline_lock())){
                $lock[3] = localtime($lock[5] = $lock[3]);
		print STDERR "The db suggests the pipeline is running! \n\n";
		print STDERR "    db       $dbname\@$dbhost              \n";
		print STDERR "    pid      $lock[2] on host $lock[1]     \n";
		print STDERR "    started  $lock[3] ($lock[5])         \n\n";
		if($checked){
		    print STDERR "The process above can be terminated by rerunning this script with\n";
		    print STDERR "-stop option rather than -status.\n";
		}else{
		    print STDERR "The process does not exist, remove the lock by removing the lock   \n";
		    print STDERR "from the pipeline database:                                      \n\n";
		    print STDERR "    echo \" delete from meta where meta_key = 'pipeline.lock'; \" ";
		    print STDERR "| mysql -h $dbhost -D $dbname -u $dbuser -P $dbport \n";
		    print STDERR "-p$dbpass" if $dbpass;
		    print STDERR "\n\n You probably need to remove $PID_FILE too.                  \n\n";
		    print STDERR "    rm -f $PID_FILE\n";
		}
            }else{
                print STDERR "$NAMED : No lock present in the database! This is wrong\n";
            }
        }
    }else{
        print STDERR "$NAMED : Instance using pid file <$PID_FILE> is not currently running\n";
    }
    exit;
}
sub get_pid{
    open(PID,"$PID_FILE") or die "Can't read $PID_FILE: $!";
    my $pid = <PID>;
    close PID;
    return $pid;
}
sub check_pid{
    my $pid = shift;
    my $found = 0;
    open(my $fh, "ps -p $pid |") || die "Can't open pipe with ps";
    while(<$fh>){
	my @F = split;
	$found = 1 if $F[0] eq $pid;
    }
    close $fh;
    return $found ? 1 : 0;
}

sub pause{
    my $script = shift;
    if($script & 1){ # called from within SIG HANDLER
	print STDERR "$NAMED : Going to sleep for $OVERLOAD_SLEEP secs [blocking SIGnals]\n";
	my $sigset = POSIX::SigSet->new(&POSIX::SIGINT);
	$sigset->addset(&POSIX::SIGALRM); # If we want to block ALARMS TOO.
	sigprocmask(SIG_BLOCK, $sigset)   or die "Can't Block SIGINT for sleep  : $!\n";
	sleep($OVERLOAD_SLEEP);
	sigprocmask(SIG_UNBLOCK, $sigset) or die "Can't Unblock SIGINT for sleep: $!\n";
    }elsif($script & 2){
	print STDERR "$NAMED : Going to sleep for $PAUSE secs\n";
	sleep($PAUSE);	
    }else{       # called from this script as a control
	my $SIG = 31;
	if(my $pid = eval{get_pid()}){
	    kill $SIG, $pid;
	}else{
	    print STDERR "$NAMED : COULDN'T DO IT\n";
	}
    }
}
sub check_db{
    
}
sub restart{
    print STDERR "$NAMED : ****************************************\n";
    print STDERR "$NAMED : Please use -stop then -start to restart\n";
    print STDERR "$NAMED : ****************************************\n";
    if(my $pid = eval{get_pid()}){
        kill 1, $pid;
    }else{
        print STDERR "$NAMED : \n";
    }
    exit;
}
{
    my $restart_string;
    sub restart_string{
        my (@passed) = @_;
        my %controls = @controls;
        my @control = keys %controls;
        if(@passed){
            foreach my $option(splice @passed){
                if($option =~ /^\-/){
                    $option = substr($option,1);
                    my $ok = 0;
                    map { $ok = 1 if $_ =~ /^$option(.+)?/ } @control;
                    $restart_string .= " -" . $option unless $ok;
                }else{
                    $restart_string .= " " . $option;
                }
            }
            print STDERR "$NAMED : RestartString <$SELF$restart_string -buffyMode>\n";
        }
        return $restart_string . " -buffyMode\n";
    }
    sub restart_array{
        unless($restart_string){
            return ();
        }
        return split(" ", $restart_string), '-buffyMode';
    }

}
sub useage{
    exit( exec('perldoc', $0) or die "Can't perldoc $0: $!" );
}

=pod

=head1

FinishedRuleManager [-options]

 -start
 -stop
 -restart
 -status
 -buffy
 -dbhost
 -dbport
 -dbuser
 -dbuser
 -dbname

=cut
