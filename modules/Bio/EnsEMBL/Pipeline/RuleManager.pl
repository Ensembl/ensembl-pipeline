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

$mailReport = "stabenau@ebi.ac.uk";
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
    if( scalar( @idList ) < $chunksize ) {
      $completeRead = 1;
    }
  }

  if( ! check_stop_condition ) {
    while ( @hotIds ) {
      my $hotId = shift( @hotIds );
      
      my @anals = $sic->fetch_analysis_by_inputId_class
	( $hotId->[0], $hotId->[1] );
      %analHash = ();
      
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

	push( @hotJobs, [ $job, time ]  );
      }
    }
  }

  # now try to find about the jobs in
  # the hotlist
  
  # note successful jobs
  my @success = $jobAdaptor->list_jobId_by_status( "SUCCESSFUL" );
  my %successHash = map { ( $_,1) } @success;
  my @newHot;
  while( @hotJobs ) {
    $job = shift( @hotJobs );
    if( defined $successHash{$job->dbID} ) {
      # yep, successful job
      $sic->store_inputId_class_analysis
	( $job->input_id, resolve_class( $job->analysis ), 
	  $job->analysis );
      $job->remove;
    } else {
      push( @newHot, $job );
    }
  }

  @hotJobs = @newHot;
  
  
  


}




sub resolve_queue {
  my $anal = shift;
  return "blastfarm";
}

sub resolve_class {
  my $anal = shift;
  return "contig";
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

