# Mar 9, 2006 10:20:44 AM
#
# Created by Mustapha Larbaoui <ml6@sanger.ac.uk>

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::Finished::Job

=head1 SYNOPSIS

=head1 DESCRIPTION

Run a Finished analysis job.
Allow to save the database version and the runtime info in input_id_analysis.

=head1 FEEDBACK

=head1 AUTHOR - Mustapha Larbaoui

Mustapha Larbaoui E<lt>ml6@sanger.ac.ukE<gt>

=head1 CONTACT

Post general queries to B<anacode@sanger.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _


=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::Finished::Job;

use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Pipeline::Job;
use Bio::EnsEMBL::Pipeline::Config::BatchQueue;
use Bio::EnsEMBL::Pipeline::Config::General;
use Bio::EnsEMBL::Analysis::Tools::Logger;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

@ISA = qw(Bio::EnsEMBL::Pipeline::Job);

# dynamically load appropriate queue manager (e.g. LSF)

my $batch_q_module = "Bio::EnsEMBL::Pipeline::BatchSubmission::$QUEUE_MANAGER";

my $file = "$batch_q_module.pm";
$file =~ s{::}{/}g;
require "$file";

# BATCH_QUEUES is a package variable which stores available
# 'queues'. this allows different analysis types to be sent
# to different nodes etc. it is keyed on analysis logic_name
# and has different parameters such as resource (e.g. node set),
# number of jobs to be batched together etc.
#
# this may be better encapsulated as a separate class, but
# i can't think at the moment how best to do this.

my %BATCH_QUEUES = &set_up_queues;

=head2 run_module

Run the job

=cut

sub priority {
	my ($self,$priority) = @_;
	if ($priority) {
    	$self->{priority} = $priority;
  	}
	if ( !$self->{priority} ) {
		my $ln       = $self->analysis->logic_name;
		my $p = $BATCH_QUEUES{$ln}{priority};
		throw("Priority for $ln not set in BatchQueue file\n") unless $p;
		if ( scalar(@$p) == 1 ) {
			$self->{priority} = $p->[0];
		}
		else {
			$self->set_update_value;
			$self->{priority} = $self->update ? $p->[1] : $p->[0];
		}
	}

	return $self->{priority};
}

# update values: 0 no update, new clone analysis; 1 update only patch file; 2 update patch and release files; 
sub set_update_value {
	my ($self) = @_;
	my $update_value = 0;
	my $sic  = $self->adaptor->db->get_StateInfoContainer;
	my $db_version_saved = $sic->fetch_db_version($self->input_id, $self->analysis);
	my $db_version_current = $self->analysis->db_version;
	if($db_version_saved) {
		$update_value = 2;
		# split the embl blast db version "12-Mar-06 (85)" to
		# patch version "12-Mar-06" and release version "85"
		my ($patch_sv,$release_sv) = $db_version_saved =~ /^(\S+)\s+\((\d+)\)$/;
		my ($patch_cv,$release_cv) = $db_version_current =~ /^(\S+)\s+\((\d+)\)$/;
		if($release_sv && ($release_sv eq $release_cv)){
			$update_value = 1;
		}	
	}
	$self->{update} = $update_value;
}   

sub update {
	my ($self,$update) = @_;
	if ($update) {
    	$self->{update} = $update;
  	}
  	
	return $self->{update} || 0;
}

sub run_module {
	my $self               = shift;
	my $is_dbversion_saved = 0;

	# start timer
	my $start = time;

	my $module   = $self->analysis->module;
	my $hash_key = $self->analysis->logic_name;
	my $rdb;
	my ( $err, $res );

	print "Running " . $module . " with " . $self . "\n";

	if ( !exists( $BATCH_QUEUES{$hash_key} ) ) {
		$hash_key = 'default';
	}

	my $runnable_db_path = $BATCH_QUEUES{$hash_key}{runnabledb_path};
	my $verbosity        = $BATCH_QUEUES{$hash_key}{verbosity};

	my $perl_path;

	#print STDERR "Getting ".$hash_key." batchqueue value\n";

	if ( $module =~ /::/ ) {

		#print STDERR "Module contains path info already\n";
		$module =~ s/::/\//g;
		$perl_path = $module;
	}
	elsif ($runnable_db_path) {
		$perl_path = $runnable_db_path . "/" . $module;
	}
	else {
		$perl_path = $module;
	}

	my $current_verbosity = logger_verbosity;
	verbose($verbosity);
	logger_verbosity($verbosity);

	#print STDERR "have perlpath ".$perl_path."\n";
  STATUS:
	{
		eval {
			require $perl_path . ".pm";
			$perl_path =~ s/\//::/g;
			$rdb = $perl_path->new(
				-analysis  => $self->analysis,
				-input_id  => $self->input_id,
				-db        => $self->adaptor->db,
				-verbosity => $verbosity
			);
		};

		if ( $err = $@ ) {
			print( STDERR "CREATE: Lost the will to live Error\n" );
			$self->set_status("FAILED");
			throw(  "Problems creating runnable $module for "
				  . $self->input_id
				  . " [$err]\n" );
		}

		# "READING"
		print "READING\n";
		eval {
			$self->set_status("READING");
			$res = $rdb->fetch_input;
		};
		if ( $err = $@ ) {
			$self->set_status("FAILED");
			print( STDERR "READING: Lost the will to live Error\n" );
			throw(  "Problems with $module fetching input for "
				  . $self->input_id
				  . " [$err]\n" );
		}

		if ( $rdb->input_is_void ) {
			$self->set_status("VOID");
		}
		else {
			print "RUNNING\n";

			# "RUNNING"
			eval {
				$self->set_status("RUNNING");
				$rdb->db->dbc->disconnect_when_inactive(1);
				$rdb->run;
				$rdb->db->dbc->disconnect_when_inactive(0);
			};
			if ( $err = $@ ) {
				print STDERR $@ . "\n";

				if ( my $err_state = $rdb->failing_job_status ) {
					$self->set_status($err_state);
					if ( $err_state eq 'VOID' ) {
						last STATUS;
					}
				}
				else {
					$self->set_status("FAILED");    # default to just failed
					                                #these jobs get retried
				}

				print( STDERR "RUNNING: Lost the will to live Error\n" );
				throw(  "Problems running $module for "
					  . $self->input_id
					  . " [$err]\n" );
			}

			# "WRITING"
			print "WRITING\n";
			eval {
				$self->set_status("WRITING");
				$rdb->write_output;

				if ( $rdb->can('db_version_searched') ) {
					my $new_db_version = $rdb->db_version_searched();
					my $analysis       = $self->analysis();
					my $old_db_version = $analysis->db_version();

					$analysis->db_version($new_db_version);

					#$self->adaptor->db->get_AnalysisAdaptor->update($analysis);
					$is_dbversion_saved = 1;
				}

				$self->set_status("SUCCESSFUL");
			};
			if ( $err = $@ ) {
				$self->set_status("FAILED");
				print( STDERR "WRITING: Lost the will to live Error\n" );
				throw(  "Problems for $module writing output for "
					  . $self->input_id
					  . " [$err] [ $?]" );
			}
		}
	}

	# end timer
	my $end = time;

	# Run time in seconds
	my $runtime = $end - $start;

	logger_verbosity($current_verbosity);

	# update job in StateInfoContainer
	eval {
		my $sic = $self->adaptor->db->get_StateInfoContainer;
		$sic->store_input_id_analysis( $self->input_id, $self->analysis,
			$self->execution_host, $is_dbversion_saved, $runtime );
	};
	if ( $err = $@ ) {
		my $error_msg =
		    "Job finished successfully, but could not be "
		  . "recorded as finished.  Job : ["
		  . $self->input_id
		  . "]\n[$err]";
		eval { $self->set_status("FAIL_NO_RETRY"); };
		$error_msg .= (
"(And furthermore) Encountered an error in updating the job to status failed_no_retry.\n[$@]"
		  )
		  if $@;
		throw($error_msg);
	}
	else {
		print STDERR "Updated successful job " . $self->dbID . "\n";
	}
}

=head2 batch_runRemote

  Title   : batch_runRemote
  Usage   : $job->batch_runRemote
  Function: see parent class
  Returns : 
  Args    : Is static, private function, dont call with arrow notation.

=cut

sub batch_runRemote {
	my ($self) = @_;
	my $submitted = 1;
	my $queue;
	my $dbname = $self->adaptor->db->dbc->dbname();
	my $host = $self->adaptor->db->dbc->host();

	if ( !exists( $BATCH_QUEUES{ $self->analysis->logic_name } ) ) {
		$queue = 'default';
	}
	else {
		$queue = $self->analysis->logic_name;
	}
	# add job to batch jobs array
	my $batch_jobs = $BATCH_QUEUES{$queue}{'jobs'};
	$batch_jobs->{$host}->{$dbname} = [] unless $batch_jobs->{$host}->{$dbname};
	push @{ $batch_jobs->{$host}->{$dbname} }, $self->dbID;
	# maximum batch jobs size
	my $batch_size = $BATCH_QUEUES{$queue}{'batch_size'};
	if(scalar(@$batch_size) > 1){
		$batch_size = ($self->update == 1) ? @$batch_size[1] : @$batch_size[0];
	} else {
		$batch_size = @$batch_size[0];
	}
	
	if (
		scalar( @{ $batch_jobs->{$host}->{$dbname} } ) >= $batch_size )
	{
		$self->flush_runs( $self->adaptor, $queue );
	} else {
		$submitted = 0;
	}
	
	return $submitted;
}

=head2 flush_runs

  Title   : flush_runs
  Usage   : $job->flush_runs( jobadaptor, [queue] );
  Function: Methode extended to handle the failed 'out of memory' Jobs and use the big memory queue. 
  Returns : 
  Args    : 

=cut

sub flush_runs {
	my ( $self, $adaptor, $queue, $verbose ) = @_;
	# flush_runs is optionally sent a queue to deal with
	# @analyses is a list of logic_names (strings)

	my @analyses = ($queue) || ( keys %BATCH_QUEUES );

	if ( !defined $adaptor ) {
		throw("Cannot run remote without db connection");
	}

	local *FILE;

	my $dbc      = $adaptor->db->dbc;
	my $host     = $dbc->host;
	my $username = $dbc->username;
	my $dbname   = $dbc->dbname;
	my $pass     = $dbc->password;
	my $port     = $dbc->port;

		# runner.pl: first look at value set in RuleManager ($RUNNER_SCRIPT)
	# then in same directory as Job.pm,
	# and fail if not found

		my $runner = $self->runner;

		if ( !$runner || !-x $runner ) {
		$runner = __FILE__;
		$runner =~ s:/[^/]*$:/runner.pl:;
		my $caller = caller(0);
		throw(  "runner " . $runner
			  . " not found - needs to be set in "
			  . "$caller\n" )
		  unless -x $runner;
	}
	  ANAL:
	for my $anal (@analyses) {
		my $queue = $BATCH_QUEUES{$anal};
		my @job_ids;
		@job_ids = @{ $queue->{'jobs'}->{$host}->{$dbname} } 
				if ($queue->{'jobs'}->{$host}->{$dbname});
		if ( !@job_ids ) {
			next ANAL;
		}
		print "\t\t$anal\t".scalar(@job_ids)." jobs\n" if $verbose;
			my $this_runner = $queue->{'runner'};
		$this_runner = ( -x $this_runner ) ? $this_runner : $runner;

		my $lastjob = $adaptor->fetch_by_dbID( $job_ids[-1] );
		
		while( !$lastjob && @job_ids) {
			pop @job_ids;
			$lastjob = $adaptor->fetch_by_dbID( $job_ids[-1] ) if($job_ids[-1]);
		}
		
		if ( !$lastjob ) {
			next ANAL;
		}

			my $pre_exec =
		  $this_runner . " -check -output_dir " . $self->output_dir;

		my $farm_queue    = $queue->{'queue'};
		my $farm_resource = $queue->{'resource'};
		my $param = $queue->{'sub_args'}.' -sp '.$self->priority.' ';
		
		if ( $self->priority == $BIG_MEM_PRIORITY ) {
			$farm_queue    = $BIG_MEM_QUEUE;
			$farm_resource = $BIG_MEM_RESOURCE;
			$param 		  .= $BIG_MEM_PARAM; 
		}
		
		if ( $self->priority == $LONG_JOB_PRIORITY ) {
			$farm_queue    = $LONG_JOB_QUEUE;
		}

			my $batch_job = $batch_q_module->new(
			-STDOUT     => $lastjob->stdout_file,
			-PARAMETERS => $param,
			-PRE_EXEC   => $pre_exec,
			-QUEUE      => $farm_queue,
			-JOBNAME    => $dbname . ':' . $anal,
			-NODES      => $queue->{'nodes'},
			-RESOURCE   => $farm_resource
		);

			my $cmd;

			if ( !$self->cleanup ) {
			$batch_job->stderr_file( $lastjob->stderr_file );
		}

			# check if the password has been defined, and write the
		# "connect" command line accordingly otherwise -pass gets the
		# first job id as password, instead of remaining undef

			if ($pass) {
			$cmd = $runner
			  . " -dbhost $host -dbuser $username -dbname $dbname -dbpass $pass -dbport $port";
		}
		else {
			$cmd = $runner
			  . " -dbhost $host -dbuser $username -dbname $dbname -dbport $port";
		}
		$cmd .= " -output_dir " . $self->output_dir;
		$cmd .= " -queue_manager $QUEUE_MANAGER  ";
		if ( $self->cleanup ) {
			$cmd .= " -cleanup ";
		}
		$cmd .= " @job_ids";

			$batch_job->construct_command_line($cmd);

			eval {

				# SMJS LSF Specific for debugging
			#print "Submitting: ", $batch_job->bsub, "\n";
			$batch_job->open_command_line();
		};

			if ($@) {
			print STDERR "Couldnt batch submit @job_ids \n[$@]\n";
			print STDERR "Using " . $batch_job->bsub . "\n";
			foreach my $job_id (@job_ids) {
				my $job = $adaptor->fetch_by_dbID($job_id);
				$job->set_status("FAILED");
			}
		}
		else {
			my @jobs = $adaptor->fetch_by_dbID_list(@job_ids);
			foreach my $job (@jobs) {

					if ( $job->retry_count > 0 ) {
					for ( $job->stdout_file, $job->stderr_file ) {
						open( FILE, ">" . $_ );
						close(FILE);
					}
				}

					if ( $batch_job->id ) {
					$job->submission_id( $batch_job->id );
				}
				else {

						# submission seems to have succeeded, but we didnt
					# get a job ID. Safest NOT to raise an error here,
					# (a warning would have already issued) but flag
					print STDERR
"Job: Null submission ID for the following, but continuing: @job_ids\n";
					$job->submission_id(0);
				}
				$job->retry_count( $job->retry_count + 1 );
				$job->set_status("SUBMITTED");
				$job->stdout_file( $lastjob->stdout_file );
				$job->stderr_file( $lastjob->stderr_file );
			}
			$adaptor->update(@jobs);
		}
		$queue->{'jobs'}->{$host}->{$dbname}         = [];
		$queue->{'last_flushed'}->{$host}->{$dbname} = time;
	}
}

sub set_up_queues {
	my %q;

	foreach my $queue (@$QUEUE_CONFIG) {
		my $ln = $queue->{logic_name};

		next unless $ln;

		while ( my ( $k, $v ) = each %$queue ) {
			$q{$ln}{$k} = $v;
		}
		$q{$ln}{jobs}         = {};
		$q{$ln}{last_flushed} = undef;
		$q{$ln}{batch_size}      ||= $DEFAULT_BATCH_SIZE;
		$q{$ln}{queue}           ||= $DEFAULT_BATCH_QUEUE;
		$q{$ln}{retries}         ||= $DEFAULT_RETRIES;
		$q{$ln}{cleanup}         ||= $DEFAULT_CLEANUP;
		$q{$ln}{runnabledb_path} ||= $DEFAULT_RUNNABLEDB_PATH;
		$q{$ln}{output_dir}      ||= $DEFAULT_OUTPUT_DIR;
		$q{$ln}{runner}          ||= $DEFAULT_RUNNER;
		$q{$ln}{verbosity}       ||= $DEFAULT_VERBOSITY;
	}

	# a default queue for everything else
	if ( !exists( $q{default} ) ) {
		$q{default}{jobs}         = {};
		$q{default}{last_flushed} = undef;
	}

	# Need these set, so do the ||= thing
	$q{default}{batch_size}      ||= $DEFAULT_BATCH_SIZE;
	$q{default}{queue}           ||= $DEFAULT_BATCH_QUEUE;
	$q{default}{retries}         ||= $DEFAULT_RETRIES;
	$q{default}{cleanup}         ||= $DEFAULT_CLEANUP;
	$q{default}{runnabledb_path} ||= $DEFAULT_RUNNABLEDB_PATH;
	$q{default}{output_dir}      ||= $DEFAULT_OUTPUT_DIR;
	$q{default}{runner}          ||= $DEFAULT_RUNNER;
	$q{default}{verbosity}       ||= $DEFAULT_VERBOSITY;
	return %q;
}

1;
