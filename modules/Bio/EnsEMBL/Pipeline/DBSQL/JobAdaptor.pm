# Object for storing an idlist and calculating diffences between lists
#
# Cared for by ensembl (ensembl-dev@ebi.ac.uk)
#
# Copyright ensembl
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod

=head1 NAME

Bio::EnsEMBL::Pipeline::DBSQL::JobAdaptor

=head1 SYNOPSIS


=head1 DESCRIPTION

adaptor to access and store entries in the job and job_status tables

=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::DBSQL::JobAdaptor;

use vars qw(@ISA);
use strict;
use warnings;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Pipeline::Job;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


=head2 fetch_by_dbID

  Arg [1]   : database id
  Function  : to fetch an particular entry from the database based
  on dbID
  Returntype: Bio::EnsEMBL::Pipeline::Job
  Exceptions: none
  Caller    :
  Example   : my $job = $jobadaptor->fetch_by_dbID(1);

=cut

sub fetch_by_dbID{
  my ($self, $dbID) = @_;

	my $cols = $self->_cols();

  my $query = qq {
    SELECT $cols
    FROM job
    WHERE job_id = ? };
	
  my $sth = $self->prepare($query);
  $sth->execute($dbID);

  my $jobs = $self->_jobs_from_sth($sth);

  $sth->finish();

	return undef if(!@$jobs);

	return $jobs->[0];
}


=head2 fetch_by_dbID_list

  Arg [1]   : listref of ints $dbIDs
  Function  : Retrieves a list of jobs from the database
  Returntype: listref of Bio::EnsEMBL::Pipeline::Jobs
  Exceptions: none
  Caller    : general
  Example   : my $job = $jobadaptor->fetch_by_dbID(1);

=cut

sub fetch_all_by_dbID_list {
  my ($self, $dbIDs) = @_;

	if(@$dbIDs == 0) {
		return [];
	}

	my $id_list = join(', ', @$dbIDs);

	my $cols = $self->_cols;

  my $query = qq {
    SELECT $cols
    FROM job
    WHERE job_id IN ($id_list) };
	
  print STDERR "\n\nSQL:\n$query\n\n";

  my $sth = $self->prepare($query);
  $sth->execute;
	
  my $result = $self->_jobs_from_sth($sth);

  $sth->finish();

	return $result;
}



=head2 fetch_all_by_job_name

  Arg [1]    : string $job_name
  Arg [2]    : (optional) int $array_index
  Example    : Retrieves a list of jobs with a specified job_name and 
               optionally an array index.  If an array_index is specified
               then the expected return value is a listref with a single job
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub fetch_all_by_job_name {
  my ($self, $job_name, $index) = @_;

  $self->throw('job_name argument expected') if(!$job_name);

  my $columns = $self->_cols();

  my $q = "SELECT $columns
           FROM job
           WHERE job_name = ?";
  my @bindvals = ($job_name);

  if($index) {
    $q .= " AND array_index = ?";
    push @bindvals, $index;
  }

  my $sth = $self->prepare($q);
  $sth->execute(@bindvals);

  my $result = $self->_jobs_from_sth($sth);

  $sth->finish();

  return $result;
}



sub fetch_all_by_status{
  my ($self, $status) = @_;
  my @jobs;
  my $q = qq {
  SELECT j.job_id,
         j.taskname,
         j.input_id,
         j.submission_id,
         j.job_name,
         j.array_index,
         j.parameters,
         j.module,
         j.stderr_file,
	 j.stdout_file,
         j.retry_count,
  MAX(CONCAT(LPAD(js.sequence_num, 10, '0'),':',
             UNIX_TIMESTAMP(js.time), ':', js.status)) AS max_status
  FROM   job_status js,
         job j
  WHERE  js.job_id = j.job_id
  AND    js.status = ?  
  GROUP BY js.job_id };

  my $sth = $self->prepare( $q );
  $sth->execute($status);

  my ($job_id, $taskname, $input_id, $submission_id, $job_name, 
      $array_index, $parameters, $module, $stderr_file, $stdout_file,
      $retry_count, $max_status);

  $sth->bind_columns(\$job_id, \$taskname, \$input_id,
		     \$submission_id, \$job_name, \$array_index,
		     \$parameters, \$module, \$stderr_file, \$stdout_file,
		     \$retry_count, \$max_status);

  while($sth->fetch) {
    push @jobs, Bio::EnsEMBL::Pipeline::Job->new
      (-input_id => $input_id,
       -taskname => $taskname,
       -module => $module,
       -dbID => $job_id,
       -submission_id => $submission_id,
       -job_name => $job_name,
       -array_index => $array_index,
       -parameters => $parameters,
       -module => $module,
       -stderr_file => $stderr_file,
       -stdout_file => $stdout_file,
       -retry_count => $retry_count,
       -adaptor => $self);
  }

  $sth->finish();

  return \@jobs;
}

#
# private method, just returns the columns in the job table
# the implementation of _jobs_from_sth depends on this method
#
sub _cols {
	my $self = shift;
	
	return  'job_id
   , taskname
	 , input_id
	 , submission_id
	 , job_name
	 , array_index
	 , parameters
	 , module
	 , stderr_file
	 , stdout_file
	 , retry_count';
}


=head2 fetch_all_by_taskname

  Arg [1]   : taskname e.g RepeatMasker
  Function  : to get all the jobs with a particular taskname
  Returntype: listref
  Exceptions: none
  Caller    :
  Example   : @jobs = @{$jobadptr->fetch_all_by_taskname('RepeatMasker_task')}

=cut

sub fetch_all_by_taskname{
  my ($self, $taskname) = @_;

	my $cols = $self->_cols();

   my $query = qq {
    SELECT $cols
    FROM job
    WHERE taskname = ? };

  my $sth = $self->prepare($query);
  $sth->execute("$taskname");

	my $jobs = $self->_jobs_from_sth($sth);

  $sth->finish();

	foreach my $job (@$jobs) {
    my $status = $self->fetch_current_status($job);
    $job->status($status);
  }

  return $jobs;
}


=head2 _jobs_from_sth

  Arg [1]   : $sth an executed statement handle
  Function  : creates a list of job objects from a statment handle
  Returntype: Bio::EnsEMBL::Pipeline::Job
  Exceptions: none
  Caller    : 
  Example   : $jobs = @{$self->_jobs_from_sth($sth);

=cut

sub _jobs_from_sth{
  my ($self, $sth) = @_;

  my ($job_id, $taskname, $input_id, $submission_id, $job_name, 
      $array_index, $parameters, $module, $stderr_file, $stdout_file,
      $retry_count);

  $sth->bind_columns(\$job_id, \$taskname, \$input_id, 
                     \$submission_id, \$job_name, \$array_index,
                     \$parameters, \$module,
										 \$stderr_file, \$stdout_file, \$retry_count);

	my @jobs;

	while($sth->fetch) {
		push @jobs, Bio::EnsEMBL::Pipeline::Job->new
			(-input_id => $input_id,
	     -taskname => $taskname,
	     -module => $module,
		   -dbID => $job_id,
			 -submission_id => $submission_id,
			 -job_name => $job_name,
			 -array_index => $array_index,
			 -parameters => $parameters,
			 -module => $module,
			 -stderr_file => $stderr_file,
			 -stdout_file => $stdout_file,
			 -retry_count => $retry_count,
			 -adaptor => $self);
  }

  $sth->finish();

  return \@jobs;
}


=head2 fetch_current_status

  Arg [1]   : Bio::EnsEMBL::Pipeline::Job
  Function  : gets the current status of a particular job
  Returntype: string (status)
  Exceptions: none
  Caller    : 
  Example   : my $status = $jobadaptor->fetch_current_status($job);

=cut



sub fetch_current_status{
  my ($self, $job) = @_;

  my $query = qq {
    SELECT status
    FROM job_status
    WHERE job_id = ?
    ORDER by sequence_num
    DESC
    LIMIT 1};

  my $sth = $self->prepare($query);

  $sth->execute($job->dbID);

  my ($status) = $sth->fetchrow;

  $sth->finish();

  return $status;
}

=head2 list_current_status

  Arg [1]    : none
  Example    : none
  Description: list of lists with job_id, taskname, status and unixtime
  Returntype : listref of listrefs
  Exceptions : none
  Caller     : PipelineManager.pm

=cut

sub list_current_status {
  my $self = shift;

  my $q = qq {
  SELECT j.job_id,
         j.input_id,
         j.taskname,
  MAX(CONCAT(LPAD(js.sequence_num, 10, '0'),':',
             UNIX_TIMESTAMP(js.time), ':', js.status)) AS max_status
  FROM   job_status js,
         job j
  WHERE  js.job_id = j.job_id
  GROUP BY js.job_id };

  my $sth = $self->prepare( $q );
  $sth->execute();

  my $result = [];
  my ( $job_id, $input_id, $taskname, $max_status );
	$sth->bind_columns(\$job_id, \$input_id, \$taskname, \$max_status);
  while($sth->fetch() ) {
		my ($seqnum, $timestamp, $status) = split(':', $max_status);
    push( @$result,
	    [ $job_id, $taskname, $input_id,  $status, $timestamp ]);
  }

  $sth->finish();

  return $result;
}

=head2 update_status

  Arg [1]    : Bio::EnsEMBL::Pipeline::Job $job
  Example    : $job_adaptor->update_status($job, 'SUBMITTED');
  Description: Updates the status of a single job by adding a row to the
               job_status table.
  Returntype : none
  Exceptions : thrown if the job does not have a dbID, or the status is not
               set
  Caller     : Job

=cut

sub update_status {
  my $self = shift;
  my $job = shift;
  my $status = shift;

	$self->throw('status argument is required') if(!$status);
	$self->throw('job must be stored in database') if(!$job->dbID);

  my $query =
    "INSERT INTO job_status (job_id, time, status)
     VALUES (?, NOW(), ?)";

  my $sth = $self->prepare( $query );
  $sth->execute( $job->dbID(), $status );

  $sth->finish();
}



=head2 update

  Arg [1]    : Bio::EnsEMBL::Pipeline::Job
  Example    : $job_adaptor->update($job);
  Description: updates a job that has already been stored
  Returntype : none
  Exceptions : thrown if the job argument does not already have a dbID
  Caller     : general

=cut

sub update {
  my $self = shift;
	my $job  = shift;

  $self->throw('Cannot update job that does not have dbID') if(!$job->dbID);

	my $sth = $self->prepare(
		'UPDATE job SET taskname = ?,
                    input_id = ?,
	                  submission_id = ?,
	                  job_name = ?,
	                  array_index = ?,
	                  parameters = ?,
	                  module = ?,
	                  stderr_file = ?,
	                  stdout_file = ?,
                    retry_count = ?
    WHERE job_id = ?');

	$sth->execute($job->taskname,
								$job->input_id,
								$job->submission_id,
								$job->job_name,
								$job->array_index,
								$job->parameters,
								$job->module,
								$job->stderr_file,
								$job->stdout_file,
								$job->retry_count,
							$job->dbID);

  $sth->finish();
}



=head2 store

  Arg [1]   : Bio::EnsEMBL::Pipeline::Job
  Function  : stores a job in the job table
  Returntype: string, dbID
  Exceptions: none
  Caller    : 
  Example   : $jobadaptor->store($job);

=cut



sub store{
  my ($self, $job) = @_;

  my $query = qq {
    INSERT INTO job ( taskname,
		      input_id,
		      submission_id,
		      job_name,
		      array_index,
		      parameters,
		      module,
		      stdout_file,
		      stderr_file,
					retry_count)
             VALUES (  ?, ?, ?, ?, ?, ?, ?, ?, ?, ? ) };

  my $sth = $self->prepare($query);

  $sth->execute( $job->taskname,
		 $job->input_id,
		 $job->submission_id,
		 $job->job_name,
		 $job->array_index,
		 $job->parameters,
		 $job->module,
		 $job->stdout_file,
		 $job->stderr_file,
		 $job->retry_count);

  $sth->finish();

  $job->adaptor($self);
  $job->dbID($sth->{'mysql_insertid'});
  $self->update_status($job, 'CREATED');
  return $job->dbID;
}


=head2 remove

  Arg [1]   : Bio::EnsEMBL::Pipeline::Job
  Function  : removes a job from the job tables
  Returntype: none
  Exceptions: none
  Caller    : 
  Example   : $jobadaptor->remove($job);

=cut

sub remove{
  my ($self, $job) = @_;

  my $job_id = $job->dbID;

  my $query = qq {
    DELETE from job
    WHERE job_id = ? };

  my $sth = $self->prepare($query);
  $sth->execute($job_id);

  my $second = qq {
    DELETE from job_status
    WHERE job_id = ? };

  $sth = $self->prepare($second);
  $sth->execute($job_id);

  $sth->finish();
}

1;
