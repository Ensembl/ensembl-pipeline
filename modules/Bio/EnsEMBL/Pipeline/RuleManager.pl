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

use Bio::EnsEMBL::Pipeline::DBSQL::RuleAdaptor;;
use Bio::EnsEMBL::Pipeline::DBSQL::JobAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::AnalysisAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::StateInfoContainer;
use Bio::EnsEMBL::Pipeline::DBSQL::Obj;


my $mailReceiver = "scp\@sanger.ac.uk";
my $maxMails = 50000;
my $currentMail = 0;
$| = 1;
my %stats = ( );
my $stopfile = '/nfs/acari/enspipe/.rmstop_m';

# how many second before job becomes old?
my $oldJob = 36000;

my $chunksize = 50000;
# my $currentStart = 300000;
my $currentStart = 0;
my $completeRead = 0;

my $pid = fork();
die( "Cant fork hoover" ) unless defined($pid);

my $db = Bio::EnsEMBL::Pipeline::DBSQL::Obj->new
  ( -host => 'ecs1a.sanger.ac.uk',
    -dbname => 'simonmouse',
    -user => 'ensadmin',
  );

my $jobAdaptor = $db->get_JobAdaptor;

if( $pid == 0 ) {
  # Im the hoover
  print STDERR ("Hoover started.\n" );
  my ($laststart,$sleeptime);
  my @jobIds;
  
  while( 1 ) {
    $laststart = time();
    # remove everything which is at least 5 minutes old
    @jobIds = $jobAdaptor->list_jobId_by_status_age( 'HOOVER', 5 );
    print STDERR ("Hoover ",scalar( @jobIds ), " jobs.\n" );
    foreach my $jobId ( @jobIds ) {
      my $job = $jobAdaptor->fetch_by_dbID( $jobId );
      $job->remove();
    }
    print STDERR ("Files deleted.\n" );
    $jobAdaptor->remove_by_dbID( @jobIds );
    # min 5 minutes between two rounds
    $sleeptime = 300 - time() + $laststart;
    print STDERR ( "Hoover done, sleep $sleeptime seconds.\n" );
    if( $sleeptime > 0 ) {
      sleep( $sleeptime );
    }
  }
}    

my $ruleAdaptor = $db->get_RuleAdaptor;
my @rules = $ruleAdaptor->fetch_all;
my %analHash;

my $sic = $db->get_StateInfoContainer;

my @hotJobs = ();
my @hotIds = ();
my %workOn;

my $mailCount = 0;

&intialize_hotJobs_workon;

while( 1 ) {

  $completeRead = 1 if -e $stopfile;
  
  if( !$completeRead ) {
    print ( "Read a chunk\n" );
    my @idList = $sic->list_inputId_class_by_start_count
      ( $currentStart, $chunksize );
    push( @hotIds, @idList );
    print ("Read ",scalar( @idList )," input ids.\n" );
    $currentStart += scalar( @idList );
    $stats{'sizeHotIds'} = scalar( @hotIds );
    if( scalar( @idList ) < $chunksize ) {
      $completeRead = 1;
    }
  }


  my $jobsStarted = $stats{'jobsStarted'};
  while ( @hotIds ) {
    my $hotId = shift( @hotIds );
    $stats{'sizeHotIds'}--; 
    my @anals = $sic->fetch_analysis_by_inputId_class
      ( $hotId->[0], $hotId->[1] );
    %analHash = ();
      
    # check all rules, which jobs can be started
    for my $rule ( @rules )  {
#      print( "Checking rule ",$rule->goalAnalysis->logic_name,".\n" );
#      print( " for ",$hotId->[0]," is " );
      	
      my $anal = $rule->check_for_analysis( @anals );
      if( $anal ) {
#       print( "fullfilled.\n" );
	$analHash{$anal->dbID} = $anal;
      } else {
#       print( "not fullfilled.\n" );
      }
    }
      
    # start all jobs ... put the to hotJobs
    for my $anal ( values %analHash ) {
      if( defined $workOn{ $hotId->[0]."-".$anal->dbID } ) {
      	next;
      }
#      print "Building new job for analysis $anal...";
      my $job = Bio::EnsEMBL::Pipeline::Job->create_by_analysis_inputId
	( $anal, $hotId->[0] ); 
      $workOn{ $hotId->[0]."-".$anal->dbID } = 1;
#     print "...storing...";
      $jobAdaptor->store( $job );
      $job->batch_runRemote( resolve_queue( $anal ));
#      $job->runLocally;
      $stats{'jobsStarted'}++;
#      print "Running",$anal->logic_name," on ",$hotId->[0],".\n";
      push( @hotJobs, [ $job, time ]  );
    }
  }
  if( $jobsStarted == $stats{'jobsStarted' } ) {
    Bio::EnsEMBL::Pipeline::Job->flush_runs( $jobAdaptor );
  }
  
  if( scalar(@hotJobs) ==0 &&
      scalar(@hotIds) == 0 &&
      $completeRead ) {
    print ( "Nothing left to do. Sleep 5 minutes.\n" );
    &print_stats;
    $completeRead = 0;
    $currentStart = 0;
    sleep( 300 );
    print( "Waking up and run again!\n" );
  }
  # now try to find about the jobs in
  # the hotlist
#  sleep( 3 );

  print "going to check.. success\n";
  &check_jobs_against_success;
  print "going to check .. failed\n";
  &check_jobs_against_failed;
  print "going to check old\n";
  &check_jobs_for_old;
  print "back for a new set\n";
  &print_stats;
}


# works on global @hotJobs, @hotIds
sub check_jobs_against_success {

  my @success = $jobAdaptor->list_jobId_by_status( "SUCCESSFUL" );
  my %successHash = map { ( $_,1) } @success;
  my @newHot;
  print ( "check ",scalar( @success )," successful jobs\n" );
  
  while( @hotJobs ) {
    my( $job, $time ) = @{shift( @hotJobs )};
    if( defined $successHash{$job->dbID} ) {
      # yep, successful job
      my $class = resolve_class( $job->analysis );
      $sic->store_inputId_class_analysis
	( $job->input_id, $class, 
	  $job->analysis );
      delete $workOn{$job->input_id."-".$job->analysis->dbID};
      # $job->remove; very costly ...
      $job->set_status( "HOOVER" );
      $stats{'jobsSuccessful'}++;
      push( @hotIds, [ $job->input_id, $class ] );
      $stats{'sizeHotIds'}++;
    } else {
      push( @newHot, [ $job, $time ] );
    }
  }
  @hotJobs = @newHot;
}


# works on global @hotJobs
sub check_jobs_against_failed {
  my @failed = $jobAdaptor->list_jobId_by_status( 'FAILED' );
  my %failedHash = map { ( $_, 1 ) } @failed;
  
  my @newHot;
  print "Failed jobs: ", scalar( @failed ), ".\n";
  $stats{'jobFailed'} += scalar( @failed );
  
  while( @hotJobs ) {
    my ( $job, $time ) = @{shift( @hotJobs )};
    if( defined $failedHash{$job->dbID} ) {
      # yep, failed job, get an updated version of it ...
      $job = $jobAdaptor->fetch_by_dbID( $job->dbID );
      if( $job->retry_count < 3 ) {
        print( "Job restart ",$job->dbID,"\n" );
	# $job->runRemote( resolve_queue( $job->analysis ));
	$job->batch_runRemote( resolve_queue( $job->analysis ));
	$stats{'jobRestarts'}++;
	push( @newHot, [ $job, time ] );
      } else {
	# dont retry, send mail ....
	# mail_job_problem( $job ); # sod this...
	$stats{'jobsFailedOften'}++;
      }
    } else {
      push( @newHot, [ $job, $time ] );
    }
  }
  @hotJobs = @newHot;
}


# what happend to the jobs which didnt anounce 
# there whereabouts in the db?
# works on @hotJobs
sub check_jobs_for_old {
  my $now = time;
  my $count = 0;

  print "Check for old jobs ";

  for my $jobTime ( @hotJobs ) {
    my ( $job, $time ) = @{$jobTime};
    if( $now-$time > $oldJob ) {
      $stats{'jobsOld'}++;
      $count++;
      if( -s $job->stdout_file ) {
	# job finished but didnt announce in the db
	# mail_job_problem( $job ); # osod this
	$job->set_status( "FAILED" );
	$stats{'jobsDied'}++;
      } else {
	# job not finished, kill it ...
	system( "bkill ".$job->LSF_id );
	$stats{'jobsLost'}++;
	$job->set_status( "FAILED" );
      }
    }
  }
  print "... found $count.\n";
}


sub mail_job_problem {
  my $job = shift;

  print ( "Problem with job mailed\n" );
  open( PIPE, "| Mail -s \"Pipeline problem [mouse]\" $mailReceiver" ) or 
    die( "Cant contact babysitter..." );
  print PIPE ( "Tried ",$job->analysis->module," 3 times and failed\n" );
  print PIPE ( "on input id ", $job->input_id,"\n" );
  
  if( open( FILE, $job->stdout_file )) {
    print PIPE "----begin stdout----\n";
    while( <FILE> ) {
      print PIPE $_;
    }
    close( FILE );
    print PIPE "----end stdout----\n";
  }
  
  if( open( FILE, $job->stderr_file )) {
    print PIPE "----begin stderr----\n";
    while( <FILE> ) {
      print PIPE $_;
    }
    close( FILE );
    print PIPE "----end stderr----\n";
	  }
  close PIPE;
  print ( "Mailed a problem.\n" );
  $mailCount++;
  if( $mailCount >= $maxMails ) {
    &to_many_error_exit;
  }
}

sub intialize_hotJobs_workon {
  # go through existing jobs and make them hotJobs
  # so they are observed.
  # as well make them workOn entries so nothing is
  # started twice.
}


sub to_many_error_exit {
  open( PIPE, "|Mail -s \"PIPELINE STOPPED [mouse]\" $mailReceiver" ) or
    die( "Cant contact babysitter..." );
  print PIPE ( "RuleManager reached maximum error report number and stopped.\n" );
  close PIPE;
  print( "Die off too many Mails.\n" );
  exit 1;
}

sub resolve_queue {
  my $anal = shift;
  # return "slow_blast_farm";
  return "acarichunky";
}

sub resolve_class {
  my $anal = shift;
  return "contig";
}


sub print_stats {
  my ( $key, $val );
  print "-- Statistics --\n";
  while( ($key,$val) = each %stats ) {
    print "   $key = $val\n";
  }
  print "--\n";
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

