# May 16, 2006 3:59:20 PM
#
# Created by Mustapha Larbaoui <ml6@sanger.ac.uk>

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::Finished::RuleManager.pm

=head1 SYNOPSIS


=head1 DESCRIPTION

Finished group specific RuleManager module.
Allow a better job submission handling through a priority queue
and limit the number of submitted jobs. Note that a ~/.netrc
file is necessary to fetch the priority queue connexion parameters.
See Net::Netrc module for more details.


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

package Bio::EnsEMBL::Pipeline::Finished::RuleManager;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Pipeline::Config::General;
use Bio::EnsEMBL::Pipeline::Config::BatchQueue;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Pipeline::RuleManager;
use Bio::EnsEMBL::Pipeline::Finished::Job;
use File::stat;
use Net::Netrc;
use DBI;

@ISA = qw(Bio::EnsEMBL::Pipeline::RuleManager);


=head2 _job_db_queue

  Function  : Gets the database connection to the RuleManager queue.
  Exceptions: none

=cut

sub _job_db_queue {
	my ($self) = @_;
	if ( !$self->{'db_queue'} ) {
		$self->{'db_queue'} = &get_dbi( $QUEUE_NAME, $QUEUE_HOST );
	}

	return $self->{'db_queue'};
}

=head2 push_job
  Function  : add a job in the queue
  Exceptions: none

=cut

sub push_job {
	my ( $self, $job, $priority ) = @_;

	my $dbq    = $self->_job_db_queue;
	my $insert = $dbq->prepare(
		qq {
			INSERT INTO queue (created, priority, job_id, host, pipeline, analysis, is_update)
			VALUES ( NOW() , ? , ? , ? , ? , ? , ?)
		}
	);
	my $job_id = $job->dbID;
	my $dbc    = $job->adaptor->db->dbc();
	my $dbname = $dbc->dbname;
	my $host   = $dbc->host;
	my $analysis = $job->analysis->logic_name;
	my $job_priority = $job->priority;
	$job_priority = $priority if $priority;
	$job_priority = $URGENT_JOB_PRIORITY
	  if ( $self->urgent_input_id->{ $job->input_id } );
	my $update = $job->update;

	return $insert->execute( $job_priority, $job_id, $host, $dbname, $analysis, $update );
}

=head2 can_run_job

  Arg [1]   : string, input_id
  Arg [2]   : Bio::EnsEMBL::Pipeline::Analysis
  Arg [3]   : string, directory path  (optional)
  Arg [4]   : string, runner script path (optional)
  Arg [5]   : int, for a boolean flag to mark verbosity (optional)
  Function  : Check if a job can be created for an input_id and analysis
              If a job already exists check if it needs to be retried,
              then push the job in a priority queue.
  Returntype: int
  Exceptions: throws if not passed an input_id or an analysis object and
              if fails to submit the job
  Example   : $rulemanager->can_run_job('filename', $analysis
                                        'path/to/dir',
                                        'path/to/runner', 1);

=cut

sub can_job_run {
	my ( $self, $input_id, $analysis, $current_jobs ) = @_;

	if ( !$input_id || !$analysis ) {
		throw(  "Can't create job without an input_id $input_id or analysis "
			  . "$analysis" );
	}

	my $job;
	my $status;

	if ( $current_jobs->{ $analysis->dbID } ) {
		my $cj = $current_jobs->{ $analysis->dbID };
		$status = $cj->{_status}->status;

		if (
			(
				   $status eq 'FAILED'
				|| $status eq 'AWOL'
				|| $status eq 'BUS_ERROR'
				|| $status eq 'OUT_OF_MEMORY'
				|| $status eq 'RUNTIME_LIMIT'
			)
			&& $cj->can_retry
		  )
		{
			print "Retrying job " . $cj->dbID . " with status $status\n"
			  if $self->be_verbose;

			if ( $self->rename_on_retry ) {
				$self->rename_files($cj);
			}
			$cj->set_status('CREATED');
			$job = $cj;
			$job->output_dir($self->output_dir) if $self->output_dir;
		}
	}
	else {
		$job = $self->create_and_store_job( $input_id, $analysis );
		print "Creating job " . $job->dbID . "\n" if $self->be_verbose;
	}

	if ($job) {
		my $priority = 0;
		if($status) {
			$priority = $BIG_MEM_PRIORITY if $status eq 'OUT_OF_MEMORY';
			$priority = $LONG_JOB_PRIORITY if $status eq 'RUNTIME_LIMIT';
		}

		$self->push_job( $job, $priority );

		return 1;
	}

	return 0;
}

sub create_and_store_job {
	my ( $self, $input_id, $analysis ) = @_;

	if ( !$input_id || !$analysis ) {
		throw(  "Can't create job without an input_id $input_id or analysis "
			  . "$analysis" );
	}
	my $job = Bio::EnsEMBL::Pipeline::Finished::Job->new(
		-input_id   => $input_id,
		-analysis   => $analysis,
		-output_dir => $self->output_dir,
		-runner     => $self->runner,
	);

	eval { $self->job_adaptor->store($job); };

	if ($@) {
		throw(  "Failed to store job "
			  . $job->input_id . " "
			  . $job->analysis->logic_name . " "
			  . $@ );
	}

	return $job;
}

=head2 urgent_input_id/read_input_file

  Arg : None
  Function  : gets a hash reference of the input_id that need to be completed at short notice.
  			  The path to the file that contains the list of urgent input_ids
  			  is set in the BatchQueue configuration file.
  			  (see variable URGENT_INPUTID_FILE)
  Returntype: Hash reference

=cut

my $input_list      = {};
my $s_last_modified = '';

sub urgent_input_id {
	my ($self) = @_;
	my $file = $URGENT_INPUTID_FILE;
	if ( -e $file ) {
		my $c_last_modified = stat($file)->[9];
		if ( $s_last_modified ne $c_last_modified ) {
			$s_last_modified = $c_last_modified;
			$input_list      = $self->read_input_file($file);
		}
	}
	else {
		$input_list = {};
	}

	return $input_list;
}

sub read_input_file {
	my ( $self, $file ) = @_;
	my $list = {};
	open( my $IN, "<$file" ) || throw("Unable to open input id file $file");
	while (<$IN>) {
		chomp;
		s/\s//g;
		$list->{$_} = 1;
	}
	close($IN);

	return $list;
}


sub get_dbi {
	my ( $dbname, $dbhost ) = @_;
	my ($dbuser, $dbpass, $dbport) = &get_db_param( $dbname, $dbhost );
	my $dsn = "DBI:mysql:host=$dbhost;dbname=$dbname;port=$dbport";

	return DBI->connect( $dsn, $dbuser, $dbpass );
}

sub get_db_param {
	my ( $dbname, $dbhost ) = @_;
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

	  return ($dbuser, $dbpass, $dbport);
}


sub check_if_done {
  my ($self) = @_;
  my @jobs = $self->job_adaptor->fetch_all;
  my $continue;

 JOB:
  foreach my $job (@jobs) {
    my $status = $job->current_status->status;

    if ($status eq 'KILLED' || $status eq 'SUCCESSFUL') {
      next JOB;
    } elsif ($status eq 'FAILED' || $status eq 'AWOL' || $status eq 'OUT_OF_MEMORY' || $status eq 'RUNTIME_LIMIT') {
      if (!$job->can_retry) {
        next JOB;
      } else {
        return 1;
      }
    } else {
      return 1;
    }
  }

  return 0;
}

1;
