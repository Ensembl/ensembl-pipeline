#!/usr/local/bin/perl5.6.1
#
# Script for operating the analysis pipeline
#
# You may distribute this code under the same terms as perl itself

$ENV{'BLASTDB'}     = '/data/blastdb/Ensembl';
$ENV{'BLASTMAT'}    = '/usr/local/ensembl/data/blastmat';
$ENV{'BLASTFILTER'} = '/usr/local/ensembl/bin';


use strict;
use Getopt::Long;
use Sys::Hostname;
use Socket;
use POSIX 'setsid';
use FindBin ();
use File::Basename ();
use File::Spec::Functions;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Data::Dumper;
use Bio::EnsEMBL::Pipeline::Config::Blast;
use Bio::EnsEMBL::Pipeline::Config::General;
use Bio::EnsEMBL::Pipeline::Config::BatchQueue;

my $script = File::Basename::basename($0);
my $SELF = catfile $FindBin::Bin, $script;

my $DEFAULT_HOME          = '';
my $HOME                  = $DEFAULT_HOME;
my $HOSTNAME              = Sys::Hostname::hostname();
my $PID_FILE              = $DEFAULT_HOME . "$SELF.$HOSTNAME.pid";
my $DEBUG                 = 1;
my $NAMED                 = 'FINISHED_PIPELINE';
my $OUT                   = $ENV{$NAMED."_OUTFILE"};


# POSIX unmasks the sigprocmask properly
my $sigset = POSIX::SigSet->new();
my $action = POSIX::SigAction->new('sigHUP_handler',
                                $sigset,
                                &POSIX::SA_NODEFER);
POSIX::sigaction(&POSIX::SIGHUP, $action);
# remove daemon control args from ARGV copy



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
my $idlist;
my ($done, $once);
my $no_run; # makes it only create jobs but doesn't submit to LSF
my $runner;
my $shuffle;  
my $output_dir;
my @start_from;
my @assembly_types = ();
my $db_sanity      = 1;
my %analyses;
my $verbose        = 1;
my $rerun_sleep    = 3600;
my $overload_sleep = 300;
my ($start, $stop, $restart, $status, $daemonFree);
my @controls = ('start'     => \$start,
                'stop'      => \$stop,
                'restart'   => \$restart,
                'status'    => \$status,
                'buffyMode' => \$daemonFree);
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
if($start + $stop + $status + $restart + $daemonFree > 1){ useage(); }
# decide which one
elsif($start)      { daemonize($OUT); }
elsif($stop)       { stop_daemon();   }
elsif($status)     { daemon_status(); }
elsif($restart)    { restart();       }
elsif($daemonFree) { print "$NAMED : Odd I'm not running as a daemon. Hit Ctrl-C quick\n" unless eval{ get_pid() }; }
else               { useage();        } # Catch when no daemon control specified
# -------------------------------------
# this is now the daemon's code

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
    print STDERR "$NAMED : NO DB SANITY CHECK.  OK?\n" if($verbose);
}
my @restart = restart_array();

sub sigHUP_handler {
    print STDERR "\n";
    print STDERR "$NAMED : got SIGHUP @ ".localtime()."\n";
    print STDERR "$NAMED : execing: $SELF @restart\n";
    $db->pipeline_unlock() if $db;
    exec($SELF, @restart) or die "Couldn't restart: $!\n";
}

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
	    push @id_list, @{$sic->list_input_id_by_Analysis_assembly_type($a, \@assembly_types)};
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
                sleep $overload_sleep;
            }
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
						      -output_dir => $output_dir,
						      -runner => $PIPELINE_RUNNER_SCRIPT);


            print "Store ", $id, " - ", $anal->logic_name, "\n" if $verbose;
            $submitted++;
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
                $job->batch_runRemote unless $no_run;
              };
              if ($@) {
                print STDERR "ERROR running job " . $job->dbID . " " . $job->stderr_file . " [$@]\n";
              }
            }

        }

    }

    &shut_down($db) if $done || $once;
    sleep($rerun_sleep) if $submitted == 0;
    @id_list = ();
    print "Waking up and run again!\n" if $verbose;
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
                LEFT JOIN analysis a ON c.condition = a.logic_name
	        WHERE a.logic_name IS NULL};
  $msg = "Some of your conditions in the rule_condition table don't" .
         " seem to have entries in the analysis table";
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
    if(my $pid = eval{get_pid()}){
        print STDERR "$NAMED : Instance running under PID: $pid\n";
        
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
                print STDERR <<EOF;
This pipeline appears to be running!

    db       $dbname\@$dbhost
    pid      $lock[2] on host $lock[1]
    started  $lock[3] ($lock[5])

The process above can be terminated by rerunning this script with
-stop option rather than -status.
If the process does not exist, remove the lock by removing the lock
from the pipeline database:

    delete from meta where meta_key = 'pipeline.lock';

You probably need to remove $PID_FILE too.

Thank you.

EOF
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

=head1 Useage

FinishedRuleManager [-options]

(general)

 -help      (show this pod)
 -start_from __ANALYSIS_ID__ (analysis id to start from)
 -once      (only go through main loop once)
 -norun     (just create jobs in job & job_status tables, useful for debug)
 -idlist __IDLISTFILE__      (file of input_ids to use)
 -assembly __ASSEMBLY_TYPE__ (assembly type of input_ids to use)

(database stuff)

 -dbhost __HOSTNAME__
 -dbport __PORT__    
 -dbname __DBNAME__  
 -dbuser __USERNAME__
 -dbpass __PASSWORD__
   
(control stuff pick one)
 -start   (starts the script as a daemon)
 -stop    (stops the daemon previously started)
 -restart (tries to restart, not recommended, use -stop then -start)
 -status  (tries to work out the status of the daemon)
 -buffy   (just runs the script without any vampires)

=cut
