#!/usr/local/ensembl/bin/perl

=pod

=head1 NAME

dequeuer.pl

=head1 DESCRIPTION

a script that dequeue and submit a limited number of jobs
(depending on JOB_LIMIT in BatchQueue) from a MySQL based priority queue
(QUEUE_HOST and QUEUE_NAME) into the LSF farm queue. The jobs are progressively
added in the MySQL queue by the pipeline rulemanager scripts.

=head1 SYNOPSIS

The MySQL queue is a single table that orders the jobs by job's priority
and by job's creation date. A job object can be recovered from the
pipeline database with the parameters 'job_id', 'host' and 'pipeline'.
Note that, the pipeline database and queue connexion parameters
(login, password and port) are fetched from the ~/.netrc file.
See the Net::Netrc module for more details.

=head1 OPTIONS

	-help|h		displays this documentation with PERLDOC
	-verbose
	-sleep		this is the amount of time the script will sleep for
			after each loop and wait for free slots (default: 180s)
	-flush		the number of loop after which the jobs batch queues should be flushed (default: 30)
	-fetch_number	number of jobs fetched from the queue (default: 100)
	-analysis 	a logic_name of an analysis you want to dequeue and submit.
			If you specify this option you will only submit this analysis
			or analyses as the option can appear on the command line
			multiple times
  	-skip_analysis	a logic_name of an analysis you don't want to dequeue and
  			submit. If this option is specified these are the only
  			analyses which won't be submit
  	-pipeline	only dequeue jobs stored in this pipeline(s)
  	-skip_pipeline	don't dequeue jobs stored in this pipeline(s)

These arguments are overridable configurations
options from Bio::EnsEMBL::Pipeline::Config::BatchQueue.pm

	-queue_manager		this specifies which
				Bio::EnsEMBL::Pipeline::BatchSubmission module is used
	-job_limit		the maximun number of jobs of the specified status allowed in the
	 			system
	-queue_name		database name of the MySQL based priority queue
	-queue_host		host name of the queue

=head1 SEE ALSO

rulemanager.pl in ensembl-pipeline/scripts/Finished

=head1 CONTACT

Mustapha Larbaoui B<email> ml6@sanger.ac.uk

=cut

use strict;
use Bio::EnsEMBL::Pipeline::DBSQL::Finished::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Config::General;
use Bio::EnsEMBL::Pipeline::Config::BatchQueue;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use DBI;
use Net::Netrc;
use Sys::Hostname;
use Getopt::Long;

my $job_id;
my $queue_manager;
my $verbose = 0;
my $job_limit;
my $queue_host;
my $queue_name;
my $sleep = 180;
my $flush = 30; # batch queue flush frequency: n => flush the batch queue once every n loops
my $fetch_number = 100;
my @analyses_to_run;
my @analyses_to_skip;
my @pipeline_to_run;
my @pipeline_to_skip;
my $db_adaptors;

my $usage = sub { exec( 'perldoc', $0 ); };

$SIG{TERM} = \&termhandler;
$SIG{INT}  = \&termhandler;

GetOptions(
	'verbose!'        => \$verbose,
	'job_limit=s'     => \$job_limit,
	'queue_manager=s' => \$queue_manager,
	'queue_name=s'    => \$queue_name,
	'queue_host=s'    => \$queue_host,
	'sleep=s'		  => \$sleep,
	'flush=s'		  => \$flush,
	'fetch_number=s'  => \$fetch_number,
	'analysis|logic_name=s@' => \@analyses_to_run,
	'skip_analysis=s@'       => \@analyses_to_skip,
	'pipeline=s@'		=> \@pipeline_to_run,
	'skip_pipeline=s@'       => \@pipeline_to_skip,
	'h|help!'         => $usage

  )
  or die("Couldn't get options");

$job_limit     = $JOB_LIMIT     unless ($job_limit);
$queue_manager = $QUEUE_MANAGER unless ($queue_manager);
$queue_host    = $QUEUE_HOST    unless ($queue_host);
$queue_name    = $QUEUE_NAME    unless ($queue_name);

@analyses_to_run = map {split/,/} @analyses_to_run ;
@analyses_to_skip = map {split/,/} @analyses_to_skip ;
@pipeline_to_run = map {split/,/} @pipeline_to_run ;
@pipeline_to_skip = map {split/,/} @pipeline_to_skip ;

# Job fetch statement handle
my $sql_fetch = qq(SELECT id, created, priority, job_id, host, pipeline, analysis, is_update FROM queue);
my @where = ();
if(@analyses_to_run) {
  my $ana_string = join("', '", @analyses_to_run);
  push @where, "analysis IN ('$ana_string')";
}
if(@analyses_to_skip) {
  my $skip_ana_string = join("', '", @analyses_to_skip);
  push @where, "analysis NOT IN ('$skip_ana_string')";
}
if(@pipeline_to_run) {
  my $pipe_string = join("', '", @pipeline_to_run);
  push @where, "pipeline IN ('$pipe_string')";
}
if(@pipeline_to_skip) {
  my $skip_pipe_string = join("', '", @pipeline_to_skip);
  push @where, "pipeline NOT IN ('$skip_pipe_string')";
}
if(scalar(@where)) {
    $sql_fetch .= ' WHERE '.join(' AND ', @where);
}
$sql_fetch .= " ORDER BY priority DESC, CREATED ASC LIMIT ? ";

my $fetch = &get_dbi( $queue_name, $queue_host )->prepare($sql_fetch);

# Job delete statement handle
my $delete = &get_dbi( $queue_name, $queue_host )->prepare(
	qq{
		DELETE FROM queue
		WHERE id = ?
	}
);

# Load the BatchSubmission module (LSF)
my $batch_q_module = "Bio::EnsEMBL::Pipeline::BatchSubmission::$queue_manager";
my $file           = "$batch_q_module.pm";
$file =~ s{::}{/}g;
eval { require "$file"; };
if ($@) {
	print STDERR "Error trying to load $batch_q_module;\ncan't find $file\n";
	exit 1;
}
my $batch_q_object = $batch_q_module->new();

my $loop = 1;
# The main loop
while ($loop) {
	my $free_slots = &job_stats;
	&flush_queue($free_slots);
	&flush_batch() if( $loop%$flush ==0);
	sleep($sleep);
	print "Waking up and run again (loop $loop) !\n" if $verbose;
	$loop++;
}

## methods

sub termhandler {
	&flush_batch();
	print "Exit DQdequeuer ....\n" if $verbose;
	exit 0;
}

# flush the batch jobs
sub flush_batch {
	print "Flushing the batch queues ....\n" if $verbose;
	foreach my $host ( keys %$db_adaptors ) {
		foreach my $dbname ( keys %{ $db_adaptors->{$host} } ) {
			print "\tpipeline $host -> $dbname\n" if $verbose;
			my $job_adaptor =
			  $db_adaptors->{$host}->{$dbname}->get_JobAdaptor();
			my ($a_job) = $job_adaptor->fetch_by_Status("CREATED",1,1);
			($a_job) = $job_adaptor->fetch_by_Status("CREATED") unless($a_job) ;
			if ($a_job) {
				$a_job->flush_runs($job_adaptor,'',$verbose);
			}
		}
	}
}

# remove a certain number of jobs from the queue
# and submit them into the farm
sub flush_queue {
	my ($slots) = @_;

	SLOT:while($slots) {
		$fetch->execute($fetch_number);
		if($fetch->rows == 0){
			print "No job in queue\n";
			last SLOT;
		}
		JOB:while ( my @row = $fetch->fetchrow_array ) {
			my ( $id, $created, $priority, $job_id, $host, $pipe_name, $analysis, $update ) =
			  @row[ 0 .. 7 ];
			my $submitted = 0;
			my $job = &get_job_from_db( $job_id, $pipe_name, $host );
			if ($job) {
				$job->priority($priority);
				$job->update($update);

				eval {
					$submitted = $job->batch_runRemote;
					my $location = $submitted ? 'LSF':'BATCH';
					print "\t$location\tsubmitted job " . $job_id
					  . " ${host}/${pipe_name} ".$job->analysis->logic_name."\tpriority "
					  . $priority . "\n"
					  if $verbose;
				};
				if ($@) {
					warn(   "ERROR running job "
						  . $job->dbID . " "
						  . $job->analysis->logic_name . " "
						  . $job->stderr_file
						  . " [$@]" );
				}
				else {
					&delete_job($id);
				}
			}
			else {
				warn(   "Job " . $job_id
					  . " not in database "
					  . $host . "/"
					  . $pipe_name );
				&delete_job($id);
			}
			$slots-- if $submitted;
			last SLOT unless $slots;
		}
	}
}

# Get some stats about farm jobs
sub job_stats {

	# Do job_stats call before getting jobs
	if ( !$batch_q_object->can('job_stats') ) {
		throw( $batch_q_object . " doesn't have the job_stats method" );
	}
	my %statuses_to_count = map { $_, 1 } @{$JOB_STATUSES_TO_COUNT};   #found in
	       #BatchQueue.pm
	my %job_stats = %{ $batch_q_object->job_stats };

	my $global_job_count = 0;    # job count for all pipelines

  GLOBAL: foreach my $sub_id ( keys %job_stats ) {
		if ( $statuses_to_count{ $job_stats{$sub_id} } ) {
			$global_job_count++;
		}
	}

	print "$global_job_count / $job_limit Pending jobs in the farm\n"
	  if ($verbose);

	my $free_slots = $job_limit - $global_job_count; # number of free farm slots
	$free_slots =
	    $free_slots > 0
	  ? $free_slots
	  : 0;    # total nb. of jobs must not exceeds job limit
	print "$free_slots slots available\n" if $verbose;
	return $free_slots;
}

sub get_job_from_db {
	my ( $job_id, $pipe_name, $host ) = @_;
	my $job_adaptor = &get_db_adaptor( $pipe_name, $host )->get_JobAdaptor();
	return $job_adaptor->fetch_by_dbID($job_id);
}

sub get_db_adaptor {
	my ( $dbname, $dbhost ) = @_;
	if ( $db_adaptors->{$dbhost}->{$dbname} ) {
		return $db_adaptors->{$dbhost}->{$dbname};
	}
	my ( $dbuser, $dbpass, $dbport ) = &get_db_param( $dbhost );

	my $db = Bio::EnsEMBL::Pipeline::DBSQL::Finished::DBAdaptor->new(
		-host   => $dbhost,
		-user   => $dbuser,
		-dbname => $dbname,
		-pass   => $dbpass,
		-port   => $dbport
	  )
	  or die(
"Failed to create Bio::EnsEMBL::Pipeline::DBSQL::Finished::DBAdaptor to db $dbname \n"
	  );

	$db_adaptors->{$dbhost}->{$dbname} = $db;

	return $db;
}

sub get_dbi {
	my ( $dbname, $dbhost ) = @_;
	my ( $dbuser, $dbpass, $dbport ) = &get_db_param( $dbhost );
	my $dsn = "DBI:mysql:host=$dbhost;dbname=$dbname;port=$dbport";

	return DBI->connect( $dsn, $dbuser, $dbpass );
}

sub get_db_param {
	my ( $dbhost ) = @_;
	my ( $dbuser, $dbpass, $dbport );

	my $ref = Net::Netrc->lookup($dbhost);
	throw("$dbhost entry is missing from ~/.netrc") unless ($ref);
	$dbuser = $ref->login;
	$dbpass = $ref->password;
	$dbport = $ref->account;
	throw(
		"Missing parameter in the ~/.netrc file:\n
			machine " .  ( $dbhost || 'missing' ) . "\n
			login " .    ( $dbuser || 'missing' ) . "\n
			password " . ( $dbpass || 'missing' ) . "\n
			account "
		  . ( $dbport || 'missing' )
		  . " (should be used to set the port number)"
	  )
	  unless ( $dbuser && $dbpass && $dbport );

	return ( $dbuser, $dbpass, $dbport );
}

sub delete_job {
	my ($id) = @_;
	return $delete->execute($id);
}

