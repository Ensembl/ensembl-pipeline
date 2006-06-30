#!/usr/local/ensembl/bin/perl 

# this script is run as a background process that
# dequeue and submit a limited number of jobs
# (depending on JOB_LIMIT in BatchQueue) from
# a directory based priority queue, DirQueue.
# Jobs are added in this queue by the pipeline rule_manager
# scripts.

use strict;
use Bio::EnsEMBL::Pipeline::DBSQL::Finished::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Config::General;
use Bio::EnsEMBL::Pipeline::Config::BatchQueue;
use IPC::DirQueue;
use Net::Netrc;
use Sys::Hostname;
use Getopt::Long;

my $job_id;
my $queue_manager;
my $verbose = 0;
my $job_limit;    # the maximum number of jobs of a perdefined status that are
    # in the system. The predefined statuses are set in the BatchQueue.pm config
    # module and currently refer to LSF

my $dir_queue = $JOB_DIR_QUEUE;    # path to the DirQueue

my $db_adaptors;

$SIG{TERM} = \&termhandler;
$SIG{INT}  = \&termhandler;

GetOptions(
	'verbose!'        => \$verbose,
	'job_limit=s'     => \$job_limit,
	'queue_manager=s' => \$queue_manager,
	'dir_queue=s'     => \$dir_queue,

  )
  or die("Couldn't get options");

$job_limit     = $JOB_LIMIT     unless ($job_limit);
$queue_manager = $QUEUE_MANAGER unless ($queue_manager);
my $options = { dir => $dir_queue, };

my $dq = IPC::DirQueue->new($options);

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

while (1) {
	my $free_slots = &job_stats;
	&flush_queue($free_slots);
	&flush_batch();
	sleep(180);
	print "Waking up and run again!\n" if $verbose;
}

sub termhandler {
	&flush_batch();
	print "Exit DQdequeuer ....\n" if $verbose;
	exit 0;
}

sub flush_batch {
	foreach my $host ( keys %$db_adaptors ) {
		foreach my $dbname ( keys %{ $db_adaptors->{$host} } ) {
			my $job_adaptor =
			  $db_adaptors->{$host}->{$dbname}->get_JobAdaptor();
			my ($a_job) = $job_adaptor->fetch_by_Status("CREATED");
			if ($a_job) {
				$a_job->flush_runs($job_adaptor);
			}
		}
	}
}

sub flush_queue {
	my ($slots) = @_;
	for ( ; $slots > 0 ; $slots-- ) {
		my $queued_job = $dq->pickup_queued_job();
		if ($queued_job) {
			my $submitted = 1;
			my $md        = $queued_job->{metadata};
			my $job_id    = $md->{job_id};
			my $pipe_name = $md->{pipeline};
			my $host      = $md->{host};
			my $priority  = $md->{priority};
			my $job       = &get_job_from_db( $job_id, $pipe_name, $host );
			if ($job) {
				$job->priority($priority);
				eval {
					    print "\tBatch running job " . $job_id
					  . " priority "
					  . $priority . "\n"
					  if $verbose;
					$submitted = $job->batch_runRemote;
				};
				if ($@) {

					#$queued_job->return_to_queue();
					$submitted = 0;
					warn(   "ERROR running job "
						  . $job->dbID . " "
						  . $job->analysis->logic_name . " "
						  . $job->stderr_file
						  . " [$@]" );
				}
				else {
					$queued_job->finish();
				}
			}
			else {

				#$queued_job->return_to_queue();
				$submitted = 0;
				warn(   "Job " . $job_id
					  . " not in database "
					  . $host . "/"
					  . $pipe_name );
			}
			$slots++ unless $submitted;
		}
		else {
			print "No job in DirQueue\n" if $verbose;
			last;
		}
	}
}

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
	my ( $dbuser, $dbpass, $dbport );
	if ( $db_adaptors->{$dbhost}->{$dbname} ) {
		return $db_adaptors->{$dbhost}->{$dbname};
	}
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

