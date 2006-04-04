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
	if ( !$runnable_db_path ) {
		$runnable_db_path = $DEFAULT_RUNNABLEDB_PATH;
	}

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
					$self->adaptor->db->get_AnalysisAdaptor->update($analysis);
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

=head2 set_up_queues

=cut

# scp: set up package variable for queue config
# i'm not sure this is the best way of doing this
# should ideally have this stuff in object(s)

sub set_up_queues {
  my %q;

  foreach my $queue (@$QUEUE_CONFIG) {
    my $ln = $queue->{logic_name};

    next unless $ln;

    #delete $queue->{logic_name};

    while (my($k, $v) = each %$queue) {
      $q{$ln}{$k}     = $v;
    }
    $q{$ln}{jobs} = [];
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
  if ( ! exists($q{default})) {
    $q{default}{jobs} = [];
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
