# Script for operating the analysis pipeline

#
# Creator: Arne Stabenau <stabenau@ebi.ac.uk>
# Date of creation: 05.09.2000
# Last modified : 15.06.2001 by Simon Potter
#
# Copyright EMBL-EBI 2000
#
# You may distribute this code under the same terms as perl itself


BEGIN {
    require "Bio/EnsEMBL/Pipeline/pipeConf.pl";
    # Can we have a way of reading a (local) pipeConf.pl as well?
    # e.g. if it exists in the current dir, use that one in preference
}

use strict;
use Getopt::Long;
use Sys::Hostname;
use Socket;

use Bio::EnsEMBL::Pipeline::DBSQL::RuleAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::JobAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::AnalysisAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::StateInfoContainer;
use Bio::EnsEMBL::Pipeline::DBSQL::Protein::DBAdaptor;
use Bio::EnsEMBL::Pipeline::LSF;


# Signals and events
#
# Rather than restarting this script all the time when something
# changes, send a signal (USR1) to the process. It looks in the
# ~/.ens-pipe.* directory and takes action depending on what's there
#
# Recommended way of shutting down the script.
#   touch ~/.ens-pipe../stop
#   kill -USR1 <process ID>
#
# Also, an alarm is set to go off every (say) 2 minutes (or not at all if
# $wakeup is false). When the script receives this it looks to see how many
# jobs are in the queue and, if there too many, goes to sleep for a while.
# You could also do this by counting the number of distinct(LSF_id) entries
# in the job table but you would have to also check for dead processes etc.
# It looks for failed jobs as well and restarts them.

# signal parameters
my $gotsig   =  0;
my $alarm    =  0;
my $wakeup   =  120;   # period to check LSF queues; set to 0 to disable
my $EXIT     = \1;
my $NEWRULES = \2;
my $NEWINPUT = \3;
my $SLEEP    = \4;

# the signal handlers
$SIG{USR1} = \&sighandler;
$SIG{ALRM} = \&alarmhandler;


# defaults: command line options override pipeConf variables,
# which override anything set in the environment variables.

my $dbhost    = $::pipeConf{'dbhost'} || $ENV{'ENS_DBHOST'};
my $dbname    = $::pipeConf{'dbname'} || $ENV{'ENS_DBNAME'};
my $dbuser    = $::pipeConf{'dbuser'} || $ENV{'ENS_DBUSER'};
my $dbpass    = $::pipeConf{'dbpass'} || $ENV{'ENS_DBPASS'};
my $queue     = $::pipeConf{'queue'}  || $ENV{'ENS_QUEUE'};
my $nodes     = $::pipeConf{'usenodes'};
my $workdir   = $::pipeConf{'nfstmp.dir'};
my $flushsize = $::pipeConf{'batchsize'};
my $retry     = $::pipeConf{'retry'} || 3;
my $max_jobs  = $::pipeConf{'maxjobs'} || 1000; # max number of (pend) jobs
my $jobname   = $::pipeConf{'jobname'};
                            # Meaningful name displayed by bjobs
                            # aka "bsub -J <name>"
                            # maybe this should be compulsory, as
                            # the default jobname really isn't any use
my $bsub      = $::pipeConf{'bsub_opt'};

$| = 1;

my $chunksize    = 1000000;   # How many InputIds to fetch at one time
my $currentStart = 0;       # Running total of job ids
my $completeRead = 0;       # Have we got all the input ids yet?
my $local        = 0;       # Run failed jobs locally
my @analysis;               # Only run this analysis ids
my %analysis;
my $submitted;
my $idlist;
my ($done, $once);

GetOptions(
    'host=s'      => \$dbhost,
    'dbname=s'    => \$dbname,
    'dbuser=s'    => \$dbuser,
    'dbpass=s'    => \$dbpass,
    'flushsize=i' => \$flushsize,
    'local'       => \$local,
    'idlist=s'    => \$idlist,
    'queue=s'     => \$queue,
    'jobname=s'   => \$jobname,
    'usenodes=s'  => \$nodes,
    'once!'       => \$once,
    'retry=i'     => \$retry,
    'analysis=s@' => \@analysis
)
or die ("Couldn't get options");

%analysis = map{$_, 1} @analysis;

my $db = Bio::EnsEMBL::Pipeline::DBSQL::Protein::DBAdaptor->new(
    -host   => $dbhost,
    -dbname => $dbname,
    -user   => $dbuser,
    -pass   => $dbpass,
);


if ($idlist && ! -e $idlist) {
    die "Must be able to read $idlist";
}

# Lock to prevent two RuleManager's connecting to the same DB.
# (i.e. same dbname and dbhost)
#
# Makes directory in $HOME and writes a DBM file in this directory
# which stores useful things like process id/host and the time it
# was started.
#
# This lock should really be in the form of a table in the database

$dbhost = &qualify_hostname($dbhost);
my $lock_dir = $ENV{"HOME"} . "/.ens-pipe.$dbhost.$dbname";
my $username = scalar getpwuid($<);



if (-e $lock_dir) {
    # Another pipeline is running: describe it
    my($subhost, $pid, $started) = &running_pipeline($lock_dir);
    $started = scalar localtime $started;

    print STDERR <<EOF;
Error: a pipeline appears to be running!

    db       $dbname\@$dbhost
    pid      $pid on host $subhost
    started  $started

You cannot have two RuleManagers connecting to the same database.
The process above must be terminated before this script can be run.
If the process does not exist, remove the lock by executing

    rm -r $lock_dir

Thankyou


EOF


    exit 1;
}


my $host = &qualify_hostname(hostname());
&create_lock($lock_dir, $host, $$);


my $ruleAdaptor = $db->get_RuleAdaptor;
my $jobAdaptor  = $db->get_JobAdaptor;
my $sic         = $db->get_StateInfoContainer;


# $LSF_params - send certain (LSF) parameters to Job. This hash contains
# things LSF wants to know, i.e. queue name, nodelist, jobname (things that
# go on the bsub command line), plus the queue flushsize. This hash is
# passed to batch_runRemote which passes them on to flush_runs.
#
# The idea is that you could have more than one of these hashes to suit
# different types of jobs, with different LSF options. You would then define
# a queue 'resolver' function. This would take the Job object and return the
# queue type, based on variables in the Job/underlying Analysis object.
#
# For example, you could put slow (e.g., blastx) jobs in a different queue,
# or on certain nodes, or simply label them with a different jobname.
# This should probably be defined in a conf file.

my $LSF_params = {};

$LSF_params->{'queue'}     = $queue     if defined $queue;
$LSF_params->{'nodes'}     = $nodes     if $nodes;
$LSF_params->{'flushsize'} = $flushsize if defined $flushsize;
$LSF_params->{'jobname'}   = $jobname   if defined $jobname;
$LSF_params->{'bsub'}      = $bsub      if defined $bsub;

# Fetch all the analysis rules.  These contain details of all the
# analyses we want to run and the dependences between them. e.g. the
# fact that we only want to run blast jobs after we've repeat masked etc.

my @rules       = $ruleAdaptor->fetch_all;

# Need here to strip rules which don't need to be run.

my @idList;     # All the input ids to check

alarm $wakeup if $wakeup;  
                # Signal to the script to do something in the future,
                # Such as check for failed jobs, no of jobs in queue

while (1) {
    $submitted = 0;

    # This loop reads input ids from the database a chunk at a time
    # until we have all the input ids.

    if (defined $idlist) {

        # read list of id / class pairs from a file,
        # e.g. a list of contigs on the golden path

        open IDS, "< $idlist" or die "Can't open $idlist";
        while (<IDS>) {
            chomp;
            my($id, $class) = split;
            ($id && $class) or die "Invalid id/class $_";
            push @idList, [ $id, $class ];
        }
        close IDS;
    }
    else {

        # This loop reads input ids from the database a chunk at a time
        # until we have all the input ids.
        # NB It's almost as much work to get one ID as the whole lot, so setting
        # the 'chunksize' variable to a small number doesn't really achieve much.

        if (!$completeRead) {
            print "Reading IDs ... ";
            my @tmp = $sic->list_inputId_class_by_start_count($currentStart, $chunksize);

            print "got ", scalar(@tmp), "\n";

            push(@idList,@tmp);

            $currentStart += scalar(@tmp);

            if (scalar(@tmp) < $chunksize) {
                $completeRead = 1;
            }
        }
    }

    @idList = shuffle(@idList);

    # Now we loop over all the input ids. We fetch all the analyses
    # from the database that have already run. We then check all the
    # rules to see which new analyses can be run.  e.g. if we've run
    # RepeatMasker we can run genscan. If we've run genscan we can run

    # blast jobs.
    #
    # All the analyses we're allowed to run are stored in a hash %analHash

    JOBID: while (@idList) {
        my $id = shift @idList;

        # handle signals. they are 'caught' in the handler subroutines but
        # it is only here they we do anything with them. so if the script
        # is doing something somewhere else (getting IDs at the start of
        # the while loop) or dozing, we have to wait until it gets here to
        # do anything...

        if ($alarm == 1) {
            $alarm = 0;

            # retry_failed_jobs($jobAdaptor, $retry, $LSF_params);
            while (&lsf_jobs($LSF_params) >= $max_jobs && $gotsig == 0) {
                sleep 300;
            }
            alarm $wakeup;
        }
        if ($gotsig == $EXIT) {
            print STDERR "Shutting down...\n";
            $done = 1;
            $gotsig = 0;
            last JOBID;
        }
        if ($gotsig == $NEWRULES) {
            @rules = $ruleAdaptor->fetch_all;
            $gotsig = 0;
        }
        if ($gotsig == $NEWINPUT) {
            $completeRead = 1;
            $gotsig = 0;
            last JOBID;
        }
        if ($gotsig == $SLEEP) {
            # this should be two signals: suspend and restart
            system("sleep 300"); # system call so that this sleep is
                                 # not "woken up" by other signals
            $gotsig = 0;
        }

        my @anals = $sic->fetch_analysis_by_inputId_class($id->[0], $id->[1]);

        my %analHash;

        # check all rules, which jobs can be started

        my @current_jobs = $jobAdaptor->fetch_by_inputId($id->[0]);

        RULE: for my $rule (@rules)  {
            if (keys %analysis && ! defined $analysis{$rule->goalAnalysis->dbID}) {
               next RULE;
            }
            print "Check ",$rule->goalAnalysis->logic_name, " - " . $id->[0];

            my $anal = $rule->check_for_analysis (@anals);

            if ($anal) {

                print " fullfilled.\n";
                $analHash{$anal->dbID} = $anal;
            } else {
                print " not fullfilled.\n";
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
                        if ($cj->current_status->status eq 'FAILED' && $cj->retry_count <= $retry) {
                            $cj->batch_runRemote($LSF_params);
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
                print "ERROR: comparing to current jobs. Skipping analysis for " . $id->[0] . " [$@]\n";
                next JOBID;
              }

            my $job = Bio::EnsEMBL::Pipeline::Job->create_by_analysis_inputId($anal, $id->[0], $id->[1]);


            print "Store ", $id->[0], " - ", $anal->logic_name, "\n";
            $submitted++;
            $jobAdaptor->store($job);

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
                $job->batch_runRemote($LSF_params);
              };
              if ($@) {
                print STDERR "ERROR running job " . $job->dbID . " " . $job->stderr_file . " [$@]\n";
              }
            }

        }

    }
    Bio::EnsEMBL::Pipeline::Job->flush_runs($jobAdaptor, $LSF_params);
    &close_down($lock_dir) if $done || $once;
    if ($completeRead == 1) {
        sleep(3600) if $submitted == 0;
        $completeRead = 0;
        $currentStart = 0;
        @idList = ();
        print "Waking up and run again!\n";
    }

}



# remove 'lock' file
sub close_down {
    my ($dir) = @_;

    unlink "$dir/db.pag";
    unlink "$dir/db.dir";
    rmdir $dir;
    exit 0;
}


# handler for SIGUSR1
sub sighandler {

    SIG: {
        if (-e "$lock_dir/stop") {
            $gotsig = $EXIT;
            unlink "$lock_dir/stop";
            last SIG;
        }
        if (-e "$lock_dir/newrules") {
            $gotsig = $NEWRULES;
            unlink "$lock_dir/newrules";
            last SIG;
        }
        if (-e "$lock_dir/newinput") {
            $gotsig = $NEWINPUT;
            unlink "$lock_dir/newinput";
            last SIG;
        }
        if (-e "$lock_dir/sleep") {
            $gotsig = $SLEEP;
            unlink "$lock_dir/sleep";
            last SIG;
        }
    }
    $SIG{USR1} = \&sighandler;
};



# handler for SIGALRM
sub alarmhandler {
    $alarm = 1;
    $SIG{ALRM} = \&alarmhandler;
}



# create lock file in home directory
sub create_lock {
    my ($dir, $host, $pid) = @_;
    my %db;

    mkdir $dir, 0775 or die "Can't make lock directory";

    dbmopen %db, "$dir/db", 0666;
    $db{'subhost'} = $host;
    $db{'pid'}     = $pid;
    $db{'started'} = time();
    dbmclose %db;
}



# running pipelines should have lock files in ~/ens-pipe.*
sub running_pipeline {
    my ($dir) = @_;
    my %db;

    dbmopen %db, "$dir/db", undef;
    my $host = $db{'subhost'};
    my $name = $db{'pid'};
    my $time = $db{'started'};
    dbmclose %db;

    return $host, $name, $time;
}



# turn a name of the form 'ecs1a' into 'ecs1a.sanger.ac.uk'
sub qualify_hostname {
    my ($hostname) = @_;

    my $addr = gethostbyname($hostname);

    my $host = gethostbyaddr($addr, AF_INET);

    return $host;
}



sub retry_failed_jobs {
    my ($ja, $retry, $LSF_params) = @_;

    my @failed_jobs = $ja->list_jobId_by_status('FAILED');
    $LSF_params->{'nodes'} = 'ecs1d';

    foreach my $jobId (@failed_jobs) {
        my $job = $ja->fetch_by_dbID($jobId);
        if ($job->retry_count <= $retry) {
            $job->batch_runRemote($LSF_params);

        }
    }
}



# Returns number of pending jobs.
sub lsf_jobs {
    my ($LSF_params) = @_;

    my $queue = $LSF_params->{'queue'};
    my $cmd;

    if (defined $queue) {
        $cmd = "bjobs -p -q $queue | grep PEND";
    }
    else {
        $cmd = "bjobs -p | grep PEND";
    }

    my $bjobs = `$cmd`;
    my $njobs = $bjobs =~ tr!\n!\n!; # no of newlines -> no of jobs

    return $njobs;

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
        -user  => $username,
    );

    foreach my $job ($lsf->get_all_jobs()) {
        my $stat = $job->status;
        my $id   = $job->id;
        $jobs{$stat}++;
        $lsf_id{$id} = 1;
    }

    foreach my $job ($db->fetch_by_age($age)) {
        my $lsf = $job->LSF_id;
        local $_ = $job->get_last_status->status;   # Are you local?
        STAT: {
            /RUNNING/ && do {
                $job->remove unless defined $lsf_id{$lsf};
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

    while (@in) {
        my $idx = int rand scalar @in;
        next if $idx > $#in;
        my $elem = splice @in, $idx, 1;
        push @out, $elem;
    }


    return @out;
}

