use strict;
use warnings;


package Bio::EnsEMBL::Pipeline::DBSQL::JobAdaptor;




sub fetch_by_dbID {
}



=head2 update_status

  Arg [1]    : Bio::EnsEMBL::Pipeline::Job $job
  Example    : $job_adaptor->update_status($job, 'SUBMITTED');
  Description: Updates the status of a single job by adding a row to the
               job_status table.
  Returntype : none
  Exceptions : none
  Caller     : Job

=cut

sub update_status {
  my $self = shift;
  my $job = shift;
  my $status = shift;

  my $current_time = time();
  my $query = 
    "INSERT INTO job_status (job_id, time, status)
     VALUES (?, FROM_UNIXTIME(?), ?)";

  my $sth = $self->prepare( $query );
  $sth->execute( $job->dbID(), $current_time, $status );
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
  my $q = "
  SELECT j.job_id, j.input_id,j.taskname,
         MAX(
           CONCAT(LPAD(UNIX_TIMESTAMP(js.time), 11,'0'),
	          ':', js.status)) AS max_status,
    FROM job_status js, job j
   WHERE js.job_id = j.job_id
   GROUP BY js.job_id ";
  my $sth = $self->prepare( $q );
  $sth->execute();

  my $result = [];
  my ( $job_id, $input_id, $taskname, $max_status );
	$sth->bind_columns(\$job_id, \$input_id, \$taskname, \$max_status);
  while($sth->fetch() ) {
    push( @$result,
	    [ $job_id, $taskname, $input_id, reverse split(':', $max_status) ]);
  }
  return $result;
}


=head2 store

  Arg [1]    : Bio::EnsEMBL::Pipeline::Job
  Example    : $job_adaptor->store($job);
  Description: Stores a single job into the pipeline database
  Returntype : Bio::EnsEMBL::Job
  Exceptions : none
  Caller     : general

=cut

sub store {
  my $self = shift;
  my $job = shift;

  my $query =
   "INSERT INTO job (taskname, input_id, parameter, module)
    VALUES (?, ?, ?, ?)";

  my $sth = $self->prepare( $query );
  $sth->execute( $job->taskname(), $job->input_id,
                 $job->parameter, $job->module );

  $job->adaptor( $self );
  $job->dbID( $sth->{'mysql_insertid'} );

  return $job->dbID();
}

1;
