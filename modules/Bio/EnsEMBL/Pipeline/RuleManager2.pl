# Script for operating the analysis pipeline
#
# Creator: Arne Stabenau <stabenau@ebi.ac.uk>
# Date of creation: 05.09.2000
# Last modified : 05.09.2000 by Arne Stabenau
#
# Copyright EMBL-EBI 2000
#
# You may distribute this code under the same terms as perl itself




# outline of the algorithm

use strict;
use Getopt::Long;

use Bio::EnsEMBL::Pipeline::DBSQL::RuleAdaptor;;
use Bio::EnsEMBL::Pipeline::DBSQL::JobAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::AnalysisAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::StateInfoContainer;
use Bio::EnsEMBL::Pipeline::DBSQL::Obj;

my $dbhost = $::pipeConf{'dbhost'};
my $dbname = $::pipeConf{'dbname'};
my $dbuser = $::pipeConf{'dbuser'};

$| = 1;

my $chunksize    = 5000;                # How many InputIds to fetch at one time
my $currentStart = 0;                   # Running total of job ids
my $completeRead = 0;                   # Have we got all the input ids yet?
my $flushsize    = 50;                  # How many jobs to store up before we submit to LSF
my $local        = 0;                   # Run failed jobs locally


GetOptions(     'host=s'     => \$dbhost,
                'dbname=s'   => \$dbname,
                'dbuser=s'   => \$dbuser,
                'local'      => \$local,
                ) or die ("Couldn't get options");


my $db = Bio::EnsEMBL::Pipeline::DBSQL::Obj->new ( -host   => $dbhost,
						   -dbname => $dbname,
						   -user   => $dbuser,
						   );

my $ruleAdaptor = $db->get_RuleAdaptor;
my $jobAdaptor  = $db->get_JobAdaptor;
my $sic         = $db->get_StateInfoContainer;


# Fetch all the analysis rules.  These contain details of all the analyses we
# want to run and the dependences between them. e.g. the fact that we only
# want to run blast jobs after we've repeat masked etc.

my @rules       = $ruleAdaptor->fetch_all;

my $jobcount;                          # Running total of jobs - needed so we know
                                       # when to flush jobs to the LSF queue

my @idList;                            # All the input ids to check

while( 1 ) {
  
    # This loop reads input ids from the database a chunk at a time
    # until we have all the input ids.

    if( !$completeRead ) {
	print ( "Read a chunk\n" );
	my @tmp = $sic->list_inputId_class_by_start_count($currentStart, $chunksize );
	
	print ("Read ",scalar( @idList )," input ids.\n" );
	
	push(@idList,@tmp);
	
	$currentStart += scalar( @tmp );
	
	if( scalar( @tmp ) < $chunksize ) {
	    $completeRead = 1;
	}
    }

    # Now we loop over all the input ids. We fetch all the analyses from the database
    # that have already run.  We then check all the rules to see which new
    # analyses can be run.  e.g. if we've run RepeatMasker we can run genscan.  If
    # we've run genscan we can run blast jobs.
    #
    # All the analyses we're allowed to run are stored in a hash %analHash
    
    JOBID: while ( @idList ) {
	my $id = shift( @idList );

	my @anals = $sic->fetch_analysis_by_inputId_class( $id->[0], $id->[1] );

	my %analHash;
      
	# check all rules, which jobs can be started

	my @current_jobs = $jobAdaptor->fetch_by_inputId($id->[0]);

	for my $rule ( @rules )  {
	    print( "\nChecking rule ",$rule->goalAnalysis->logic_name," for " . $id->[0] . "\n\n" );
	    
	    my $anal = $rule->check_for_analysis( @anals );
	    
	    if( $anal ) {
		print( "\tfullfilled.\n" );
		$analHash{$anal->dbID} = $anal;
	    } else {
		print( "\tnot fullfilled.\n" );
	    }
	}
	
	# Now we loop over all the allowed analyses in the hash.
	# We first check the database to see if the job is already running.
	# If so we skip it.
	#
	# If all is ok we create a new job, store it in the database and submit it to the 
	# batch runner.
	# 
	# At the end we check what our current job count is. If greater than $flushsize we
	# send all jobs to LSF
	
	for my $anal ( values %analHash ) {
	    print ("\nChecking analysis " . $anal->dbID . "\n\n");
	    # Check whether it is already running in the current_status table?

	    eval {
		foreach my $cj (@current_jobs) {
		    
		    print ("Comparing to current_job " . $cj->input_id . " " . 
			   $cj->analysis->dbID . " " . 
			   $cj->current_status->status . " " . 
			   $anal->dbID . "\n");
		    
		    if ($cj->analysis->dbID == $anal->dbID && $cj->current_status->status ne "FAILED") {
			print ("\nJob already in pipeline with status : " . $cj->current_status->status . "\n");
			next JOBID;
		    }
		}
	    };

	    if ($@) {
		print ("ERROR: comparing to current jobs. Skipping analysis for " . $id->[0] . " [$@]\n");
		next JOBID;
	    }
	    $jobcount++;    
	    my $job = Bio::EnsEMBL::Pipeline::Job->create_by_analysis_inputId( $anal, $id->[0] ); 

	
	    print "\n\tStoring job\n";
	    $jobAdaptor->store( $job );

            if ($local) {
                print "Running job locally\n";
                eval {
                  $job->runLocally;
                };
                if ($@) {
                   print STDERR "ERROR running job " . $job->dbID .  " "  . $job->stderr_file . "[$@]\n";
                }
            } else {
	        print "\tBatch running job\n";
	        $job->batch_runRemote('ultra_blast_farm');
            }

	}
	
	
	if ( $jobcount == $flushsize ) {
	    Bio::EnsEMBL::Pipeline::Job->flush_runs( $jobAdaptor );
	    $jobcount = 0;
	}
	
    }
    sleep( 30 );
    $completeRead = 0;
    $currentStart = 0;
    @idList = ();
    print( "Waking up and run again!\n" );
	
}




  
__END__



fetch_all_rules ..

fetch_object_attributes
? How to get the initial hot id list??



for hotJobs (jobs which have been made and are running)
  check successful state,
  add the analysis to the inputIds attributes
  add the inputId to the hotIds list

for hotIds
  for each of the ids get the full attribute covering
  check all rule which could apply
  create all jobs and add them to hotJob

for all failed hot jobs
  reissue if appropriate
  mail if mailcount is not exceeded
  die if mailcount is exceeded

