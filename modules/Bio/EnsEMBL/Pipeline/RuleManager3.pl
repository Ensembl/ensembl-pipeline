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
}

use strict;
use Getopt::Long;

use Bio::EnsEMBL::Pipeline::DBSQL::RuleAdaptor;;
use Bio::EnsEMBL::Pipeline::DBSQL::JobAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::AnalysisAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::StateInfoContainer;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;

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
my $jobname   = $::pipeConf{'jobname'};
my $retry     = $::pipeConf{'retry'} || 3;
my $sleep     = $::pipeConf{'sleep'} || 3600;

$| = 1;

my $chunksize    = 500000;  # How many InputIds to fetch at one time
my $currentStart = 0;       # Running total of job ids
my $completeRead = 0;       # Have we got all the input ids yet?
my $local        = 0;       # Run failed jobs locally
my $analysis;               # Only run this analysis ids
my $submitted;
my $jobname;                # Meaningful name displayed by bjobs
			    # aka "bsub -J <name>"
			    # maybe this should be compulsory, as
			    # the default jobname really isn't any use
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
    'analysis=s'  => \$analysis
)
or die ("Couldn't get options");

my $db = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new(
    -host   => $dbhost,
    -dbname => $dbname,
    -user   => $dbuser,
    -pass   => $dbpass,
);

my $ruleAdaptor = $db->get_RuleAdaptor;
my $jobAdaptor  = $db->get_JobAdaptor;
my $sic         = $db->get_StateInfoContainer;


# scp
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

my $LSF_params = {};
$LSF_params->{'queue'}     = $queue if defined $queue;
$LSF_params->{'nodes'}     = $nodes if $nodes;
$LSF_params->{'flushsize'} = $flushsize if defined $flushsize;
$LSF_params->{'jobname'}   = $jobname if defined $jobname;

# Fetch all the analysis rules.  These contain details of all the
# analyses we want to run and the dependences between them. e.g. the
# fact that we only want to run blast jobs after we've repeat masked etc.

my @rules       = $ruleAdaptor->fetch_all;

my @idList;     # All the input ids to check

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
	    print "Reading input IDs ...\n";
	    my @tmp = $sic->list_inputId_class_by_start_count($currentStart, $chunksize);

	    print "Read ", scalar(@tmp), " input ids.\n";

	    push(@idList,@tmp);

	    $currentStart += scalar(@tmp);

	    if (scalar(@tmp) < $chunksize) {
		$completeRead = 1;
	    }
	}
    }

    # Now we loop over all the input ids. We fetch all the analyses
    # from the database that have already run. We then check all the
    # rules to see which new analyses can be run.  e.g. if we've run
    # RepeatMasker we can run genscan. If we've run genscan we can run
    # blast jobs.
    #
    # All the analyses we're allowed to run are stored in a hash %analHash

    JOBID: while (@idList) {
	my $id = shift @idList;

	# A 'hack'. This enables you to halt the script if something
	# goes awry. Will introduce proper process control in future...

	if (-e '.rm.stop') {
	    print "Shutting down...\n";
	    $done = 1;
	    last JOBID;
	}

	my @anals = $sic->fetch_analysis_by_inputId_class($id->[0], $id->[1]);

	my %analHash;

	# check all rules, which jobs can be started

	my @current_jobs = $jobAdaptor->fetch_by_inputId($id->[0]);

	RULE: for my $rule (@rules)  {
            print STDERR $analysis . " " . $rule->goalAnalysis->dbID . "\n";
            if ($analysis && $analysis != $rule->goalAnalysis->dbID) {
               next RULE;
            }
	    print "\nChecking rule ",$rule->goalAnalysis->logic_name," for " . $id->[0] . "\n\n";

	    my $anal = $rule->check_for_analysis (@anals);

	    if ($anal) {
		print "\tfullfilled.\n";
		$analHash{$anal->dbID} = $anal;
	    } else {
		print "\tnot fullfilled.\n";
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
	    print "\nChecking analysis " . $anal->dbID . "\n\n";
	    # Check whether it is already running in the current_status table?

	    eval {
		foreach my $cj (@current_jobs) {

		    print "Comparing to current_job " . $cj->input_id . " " .
			   $cj->analysis->dbID . " " .
			   $cj->current_status->status . " " .
			   $anal->dbID . "\n";

		    if ($cj->analysis->dbID == $anal->dbID) {
			if ($cj->current_status->status eq 'FAILED' && $cj->retry_count < $retry) {
			    $cj->batch_runRemote($LSF_params);
			    print "Retrying job\n";
			}
			else {
			    print "\nJob already in pipeline with status : " . $cj->current_status->status . "\n";
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


	    # print "\n\tStoring job\n";
	    print "\n\tStoring job ", $id->[0], " ", $anal->logic_name, "\n";
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
    # FIXME: need a conditional here - only execute next line if jobs have
    # been created. This is needed if you halt the script with 'stop'
    # and need to flush created jobs.
    Bio::EnsEMBL::Pipeline::Job->flush_runs($jobAdaptor, $LSF_params);
    exit 0 if $done || $once;
    sleep($sleep) if $submitted == 0;
    $completeRead = 0;
    $currentStart = 0;
    @idList = ();
    print "Waking up and run again!\n";

}
