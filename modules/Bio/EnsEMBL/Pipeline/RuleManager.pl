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

$mailReceiver = "stabenau@ebi.ac.uk";
$maxMails = 3;
$currentMail = 0;

%stats = ( );

# how many second before job becomes old?
$oldJob = 60;

# $statusDir = '/tmp/ruleManager';

$chunksize = 500;
$currentStart = 0;
$completeRead = 0;

my $ruleAdaptor = $db->get_RuleAdaptor;
my @rules = $ruleAdaptor->fetch_all;
my $jobAdaptor = $db->get_JobAdaptor;
my %analHash;

my $sic = $db->get_StateInfoContainer;

my @hotJobs = ();
my @hotIds = ();

my $mailCount = 0;

while( 1 ) {

  if( !completeRead ) {
    my @idList = $sic->list_inputId_class_by_start_count
      ( $currentStart, $chunksize );
    push( @hotIds, @idList );
    $currentStart += scalar( @idList );
    $stats{'sizeHotIds'} = scalar( @hotIds );
    if( scalar( @idList ) < $chunksize ) {
      $completeRead = 1;
    }
  }

  print ( "Read a chunk\n" );

  if( ! check_stop_condition ) {
    while ( @hotIds ) {
      my $hotId = shift( @hotIds );
      
      my @anals = $sic->fetch_analysis_by_inputId_class
	( $hotId->[0], $hotId->[1] );
      %analHash = undef;
      
      # check all rules, which jobs can be started
      for my $rule ( @rules ) {
	$anal = $rule->check_for_analysis( @anals );
	if( $anal ) {
	  $analHash{$anal->dbID} = $anal;
	}
      }
      
      # start all jobs ... put the to hotJobs
      for my $anal ( values %analHash ) {
	my $job = Bio::EnsEMBL::Pipeline::Job->create_by_analysis_inputId
	  ( $anal, $hotId->[0] );
	$jobAdaptor->store( $job );
	$queue = resolve_queue( $anal );
	$job->runRemote( $queue );
	$stats{'jobsStarted'}++;
	push( @hotJobs, [ $job, time ]  );
      }
    }
  }
  if( scalar(@hotJobs) ==0 &&
      scalar(@hotIds) == 0 &&
      $completeRead ) {
    print ( "Nothing left to do. Sleep 5 minutes.\n" );
    print_stats;
    $completeRead = 0;
    $currentStart = 0;
    sleep( 300 );
    print( "Waking up and run again!\n" );

  # now try to find about the jobs in
  # the hotlist
  check_jobs_against_success;
  check_jobs_against_failed;
  check_jobs_startup_problem;
  check_jobs_for_old;
  print_stats;
}


# works on global @hotJobs, @hotIds
sub check_jobs_against_success {

  my @success = $jobAdaptor->list_jobId_by_status( "SUCCESSFUL" );
  my %successHash = map { ( $_,1) } @success;
  my @newHot;
  while( @hotJobs ) {
    ( $job, $time ) = @{shift( @hotJobs )};
    if( defined $successHash{$job->dbID} ) {
      # yep, successful job
      my $class = resolve_class( $job->analysis );
      $sic->store_inputId_class_analysis
	( $job->input_id, $class, 
	  $job->analysis );
      $job->remove;
      $stats{'jobsSuccessful'}++;
      push( @hotId, [ $job->input_id, $class ] );
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

  while( @hotJobs ) {
    ( $job, $time ) = @{shift( @hotJobs )};
    if( defined $failedHash{$job->dbID} ) {
      # yep, failed job, get an updated version of it ...
      $job = $jobAdaptor->fetch_by_dbID( $job->dbID );
      if( $job->retry_count < 1 ) {
	$queue = resolve_queue( $job->analysis );
	$job->runRemote( $queue );
	$stats{'jobRestarts'}++;
	push( @newHot, [ $job, time ] );
      } else {
	# dont retry, send mail ....
	mail_job_problem( $job );
	$stats{'jobsFailedOften'}++;
      }
    } else {
      push( @newHot, [ $job, $time ] );
    }
  }
  @hotJobs = @newJobs;
}


# what happend to the jobs which didnt anounce 
# there whereabouts in the db?
# works on @hotJobs
sub check_jobs_for_old {
  $now = time;

  my @newHot;
  for my $jobTime ( @hotJobs ) {
    my ( $job, $time ) = @{$jobTime};
    if( $now-$time > $oldJob ) {
      $stats{'jobsOld'}++;
      if( -s $job->stdout_file ) {
	# job finished but didnt announce in the db
	mail_job_problem( $job );
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
}


sub mail_job_problem {
  my $job = shift;

  open( PIPE, "| Mail -s \"Pipeline problem\" $mailReceiver" ) or 
    die( "Cant contact babysitter..." );
  print PIPE ( "Tried ",$job->analysis->module," 3 times and faild\n" );
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
    close( FILE )
      print PIPE "----end stderr----\n";
  }
  close PIPE;
  print ( "Mailed a problem.\n" );
  $mailCount++;
  if( $mailCount >= $maxMails ) {
    to_many_error_exit;
  }
}



sub to_many_error_exit {
  open( PIPE, "|Mail -s \"PIPELINE STOPPED\" $mailReceiver" ) or
    die( "Cant contact babysitter..." );
  print PIPE ( "RuleManager reached maximum error report number and stopped.\n" );
  close PIPE;
  print( "Die off too many Mails.\n" );
  exit 1;
}

sub resolve_queue {
  my $anal = shift;
  return "blastfarm";
}

sub resolve_class {
  my $anal = shift;
  return "contig";
}


sub print_stats {
  my ( $key, $val );
  print "-- Statistics --\n";
  for ($key,$val) ( each %stats ) {
    print "   $key=$val\n";
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

