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

Bio::EnsEMBL::Pipeline::JobAdaptor

=head1 SYNOPSIS


=head1 DESCRIPTION

adaptor to access and store entries in the job and job_status tables

=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::JobAdaptor;

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

  my $query = qq {
    SELECT job_id
         , taskname
	 , input_id
	 , submission_id
	 , array_index  
	 , parameters 
	 , module
	 , stderr_file
	 , stdour_file
    FROM job
    WHERE job_id = $dbID };  
	   
  my $sth = $self->prepare($query);
  my $hashRef;
  my $job;
  if($hashRef = $sth->fetchrow_hashref()){
    $job = $self->_job_from_hashref($hashRef);
  }

  return $job;
}


=head2 fetch_all_by_taskname

  Arg [1]   : taskname e.g RepeatMasker
  Function  : to get all the jobs with a particular taskname
  Returntype: listref
  Exceptions: none
  Caller    : 
  Example   : my @jobs = @{$jobadaptor->fetch_all_by_taskname('RepeatMasker')}

=cut



sub fetch_all_by_taskname{
  my ($self, $taskname) = @_;


   my $query = qq {
    SELECT job_id
         , taskname
	 , input_id
	 , submission_id
	 , array_index  
	 , parameters 
	 , module
	 , stderr_file
	 , stdour_file
    FROM job
    WHERE taskname = $taskname };  
	   
  my $sth = $self->prepare($query);
  my $hashRef;
  my @jobs;
  while($hashRef = $sth->fetchrow_hashref()){
    my $job = $self->_job_from_hashref($hashRef);
    my $status = $self->fetch_current_status($job);
    $job->status($status);
    push(@jobs, $job);
  }

  return \@jobs;
}


=head2 _job_from_hashref

  Arg [1]   : hashred
  Function  : creates a job object from hashref
  Returntype: Bio::EnsEMBL::Pipeline::Job
  Exceptions: none
  Caller    : 
  Example   : my $job = $self->_job_from_hashref($hashref);

=cut



sub _job_from_hashref{
  my ($self, $hashref) = @_;
  
  my $job = Bio::EnsEMBL::Pipeline::Job->new(
					     -input_id => $hashref->{'input_id'},
					     -taskname => $hashref->{'taskname'},
					     -module => $hashref->{'module'},
					    );

  $job->dbID($hashref->{'job_id'});
  $job->submission_id($hashref->{'submission_id'});
  $job->array_index($hashref->{'array_index'});
  $job->parameters($hashref->{'parameters'});
  $job->stderr_file($hashref->{'stderr_file'});
  $job->stdout_file($hashref->{'stdout_file'});
  $job->adaptor($self);

  return $job;
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
    ORDER by time  
    DESC 
    LIMIT 1};

  my $sth = $self->prepare($query);

  $sth->execute($job->dbID);

  my ($status) = $sth->fetchrow;

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
  MAX(CONCAT(LPAD(UNIX_TIMESTAMP(js.time), 11,'0'), ':', js.status)) AS max_status,
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
    push( @$result,
	    [ $job_id, $taskname, $input_id, reverse split(':', $max_status) ]);
  }
  return $result;
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
    INSERT INTO JOB ( taskname,
		      input_id,
		      submission_id,
		      array_index,
		      parameters,
		      module,
		      stdout_file,
		      stderr_file )
             VALUES (  ?,
		       ?,
		       ?,
		       ?,
		       ?,
		       ?,
		       ?,
		       ? ) };

  my $sth = $self->prepare($query);

  $sth->execute( $job->taskname,
		 $job->input_id,
		 $job->submission_id,
		 $job->array_index,
		 $job->parameters,
		 $job->module,
		 $job->stdout_file,
		 $job->stderr_file );

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
  
}

1;
