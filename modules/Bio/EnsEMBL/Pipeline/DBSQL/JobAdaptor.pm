# Perl module for Bio::EnsEMBL::Pipeline::DBSQL::JobAdaptor
#
# Creator: Arne Stabenau <stabenau@ebi.ac.uk>
# Based on Job from Michele Clamp
#
# Date of creation: 15.08.2000
# Last modified : 15.08.2000 by Arne Stabenau
#
# Copyright EMBL-EBI 2000
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::DBSQL::JobAdaptor

=head1 SYNOPSIS

  $jobAdaptor = $dbobj->get_JobAdaptor;
  $jobAdaptor = $jobobj->adaptor;


=head1 DESCRIPTION

  Module to encapsulate all db access for persistent class Job.
  There should be just one per application and database connection.


=head1 CONTACT

  Contact Arne Stabenau on implemetation/design detail: stabenau@ebi.ac.uk
  Contact Ewan Birney on EnsEMBL in general: birney@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::Pipeline::DBSQL::JobAdaptor;

use Bio::EnsEMBL::Pipeline::Job;
use Bio::EnsEMBL::Pipeline::Status;
use Bio::EnsEMBL::Utils::Exception qw(stack_trace_dump 
                                      verbose throw warning);

use vars qw(@ISA);
use strict;

@ISA = qw();

sub new {
  my ($class,$dbobj) = @_;
  my $self = bless {}, $class;

  $self->db( $dbobj );
  return $self;
}

=head2 fetch_by_dbID

  Title   : fetch_by_dbID
  Usage   : my $job = $adaptor->fetch_by_dbID
  Function: Retrieves a job from database by internal id
  Returns : throws exception when something goes wrong.
            undef if the id is not in the db.
  Args    :

=cut

sub fetch_by_dbID {
  my $self = shift;
  my $id = shift;

  my $sth = $self->prepare(q{
    SELECT job_id, input_id, analysis_id, submission_id,
           stdout_file, stderr_file, retry_count, temp_dir, exec_host
    FROM   job
    WHERE  job_id = ?
  });

  $sth->execute($id);
  my $rowHashRef = $sth->fetchrow_hashref;
  if( ! defined $rowHashRef ) {
    return undef;
  }
  my $job = $self->_objFromHashref( $rowHashRef );
  $sth->finish;
  return $job;
}


sub fetch_by_submission_id {
  my $self = shift;
  my $submissionid = shift;
  my @result;

  my $sth = $self->prepare(q{
    SELECT job_id, input_id, analysis_id, submission_id,
           stdout_file, stderr_file, retry_count, temp_dir, exec_host
    FROM   job
    WHERE  submission_id = ?
  });

  $sth->execute($submissionid);
  while( my $rowHashRef = $sth->fetchrow_hashref ) {
    push( @result, $self->_objFromHashref( $rowHashRef ));
  }
  $sth->finish;
 
  return @result;
}

=head2 fetch_by_dbID_list

  Title   : fetch_by_dbID_list
  Usage   : my $job = $adaptor->fetch_by_dbID_list
  Function: Retrieves jobs from database by internal id
  Returns : throws exception when something goes wrong.
            undef if the id is not in the db.
  Args    :

=cut

sub fetch_by_dbID_list {
  my ($self, @id) = @_;

  return undef unless @id;
  my @jobs;
  local $" = ',';   # are you local?

  my $sth = $self->prepare( qq{
    SELECT job_id, input_id, analysis_id, submission_id,
           stdout_file, stderr_file, retry_count, temp_dir, exec_host
    FROM job
    WHERE job_id in (@id) } );

  $sth->execute();
  while (my $row = $sth->fetchrow_hashref) {
    my $job = $self->_objFromHashref($row);
    push(@jobs,$job);
  }
  $sth->finish;
  return @jobs or undef;
}


=head2 fetch_by_Status_Analysis {

  Title   : fetch_by_Status_Analysis
  Usage   : my @jobs = $adaptor->fetch_by_Status_Analysis($id, $status)
  Function: Retrieves all jobs in the database matching status and
            an analysis id
  Returns : @Bio::EnsEMBL::Pipeline::Job
  Args    : Analysis obj, string status, and optional start and end limits

=cut

sub fetch_by_Status_Analysis {
  my ($self,$status, $analysis, $start, $end) = @_;

  throw("Require status and analysis id for fetch_by_Status_Analysis")
                            unless ($analysis && $status);
  if( ! defined $analysis->dbID ){
    throw( "Analysis needs to be in database" );
  }
  my $analysisId = $analysis->dbID;

  my $query = q{
        SELECT   j.job_id, j.input_id, j.analysis_id, j.submission_id,
                 j.stdout_file, j.stderr_file, j.retry_count,
           j.temp_dir, j.exec_host
        FROM     job j, job_status js
        WHERE    j.job_id = js.job_id
        AND      j.analysis_id = ?
        AND      js.status = ?
        AND      js.is_current = 'y'
        ORDER BY time desc
  };
    
  $query .= " LIMIT $start, $end" if ($start && $end);

  my $sth = $self->prepare($query);
  my $res = $sth->execute($analysisId, $status);

  my @jobs;

  while (my $row = $sth->fetchrow_hashref)
  {
    my $job = $self->_objFromHashref($row);
    push(@jobs,$job);
  }
  $sth->finish;
  return @jobs;
}

sub fetch_by_Status {
  my ($self, $status, $start, $end) = @_;

  throw("Require status for fetch_by_Status")
                            unless ($status);
    

  my $query = q{
        SELECT   j.job_id, j.input_id, j.analysis_id, j.submission_id,
                 j.stdout_file, j.stderr_file, j.retry_count, j.temp_dir, 
           j.exec_host
        FROM     job j, job_status js
        WHERE    j.job_id = js.job_id
        AND      js.status = ?
        AND      js.is_current = 'y'
        ORDER BY time desc
  };
    
  $query .= " LIMIT $start, $end" if ($start && $end);

  my $sth = $self->prepare($query);
  my $res = $sth->execute($status);

  my @jobs;

  while (my $row = $sth->fetchrow_hashref)
  {
    my $job = $self->_objFromHashref($row);
    push(@jobs,$job);
  }

  $sth->finish;

  return @jobs;
}



sub list_dbIDs{
  my ($self) = @_;

  my $query = q{
	SELECT   j.job_id
	FROM     job j
    };

  my $sth = $self->prepare($query);
  $sth->execute();

  my @ids;

  while( my ($id) = $sth->fetchrow){
    push(@ids, $id);
  }

  $sth->finish;

  return \@ids;
}

sub fetch_all{
  my ($self) = @_;

  my $query = q{
	SELECT   j.job_id, j.input_id, j.analysis_id, j.submission_id,
	         j.stdout_file, j.stderr_file, j.retry_count, j.temp_dir, 
           j.exec_host
	FROM     job j
    };

  my $sth = $self->prepare($query);
  my $res = $sth->execute();

  my @jobs;

  while (my $row = $sth->fetchrow_hashref) {
     my $job = $self->_objFromHashref($row);
     push(@jobs,$job);
  }

  $sth->finish;
  return @jobs;
}

=head2 fetch_by_age {

  Title   : fetch_by_age
  Usage   : my @jobs = $db->fetch_by_age($duration)
  Function: Retrieves all jobs in the database
            that are older than than a certain duration given in minutes.
  Returns : @Bio::EnsEMBL::Pipeline::Job
  Args    : int

=cut

sub fetch_by_age {
  my ($self,$age) = @_;

  throw("No input status for get_JobsByAge")
        unless defined($age);
    #convert age from minutes to seconds

  my $sth = $self->prepare(qq{
        SELECT j.job_id, j.input_id, j.analysis_id, j.submission_id,
         j.stdout_file, j.stderr_file, j.retry_count, j.temp_dir, 
         j.exec_host
        FROM   job j, job_status js
        WHERE  j.job_id = js.job_id
        AND    is_current = 'y'
        AND    js.time < DATE_SUB(NOW(), INTERVAL $age MINUTE)
  });
    
  my $res = $sth->execute();

  my @jobs;

  while (my $row = $sth->fetchrow_hashref) {
    my $job = $self->_objFromHashref($row);
    push(@jobs,$job);
  }

  $sth->finish;
  return @jobs;
}


=head2 fetch_by_input_id

  Title   : fetch_by_input_id
  Usage   : my @job = $adaptor->fetch_by_input_id
  Function: Retrieves all jobs from adaptor with certain input id
  Returns : list of job objects
            throws exception when something goes wrong.
  Args    :

=cut

sub fetch_by_input_id {
  my $self = shift;
  my $inputid = shift;
  my @result;

  my $sth = $self->prepare(q{
    SELECT job_id, input_id, analysis_id, submission_id,
           stdout_file, stderr_file, retry_count, temp_dir, exec_host
    FROM   job
    WHERE  input_id = ?
  });

  $sth->execute($inputid);
  while( my $rowHashRef = $sth->fetchrow_hashref ) {
    push( @result, $self->_objFromHashref( $rowHashRef ));
  }
  $sth->finish;

  return @result;
}


sub fetch_hash_by_input_id{
  my $self = shift;
  my $inputid = shift;

  my $sth = $self->prepare(q{
    SELECT j.job_id, j.input_id, j.analysis_id, j.submission_id,
           j.stdout_file, j.stderr_file, j.retry_count, j.temp_dir, j.exec_host,
           js.status, js.time, js.is_current
    FROM   job j, job_status js
    WHERE  j.input_id = ? 
      AND  j.job_id = js.job_id
      AND  js.is_current = 'y'
  });

  my @results;

  $sth->execute($inputid);
  while( my $rowHashRef = $sth->fetchrow_hashref ) {
    my $job = $self->_objFromHashref($rowHashRef);
    push(@results, $job);

    my $status = Bio::EnsEMBL::Pipeline::Status->new
          (
           '-jobid'   => $rowHashRef->{job_id},
           '-status'  => $rowHashRef->{status},
           '-created' => $rowHashRef->{time},
          );
        
    $self->current_status($job, $status);
  }

  my %hash;
  foreach my $result (@results) {
    $hash{$result->analysis->dbID} = $result;
  }
  $sth->finish;

  return \%hash;
}

=head2 store

  Title   : store
  Usage   : $job->store
  Function: puts a job in the db and gives it an internal id
            expects analysis to be already in db.
  Returns : throws exception when something goes wrong.
  Args    :

=cut

sub store {
  my $self = shift;
  my $job = shift;

  if( ! defined( $job->analysis->dbID )) {
    throw( "Need to store analysis first" );
  }

  my $sth = $self->prepare(q{
    INSERT into job (input_id, analysis_id,
                     submission_id, stdout_file, stderr_file, 
                     retry_count, temp_dir, exec_host)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?)
  });

  $sth->execute( 
                $job->input_id,
                $job->analysis->dbID,
                $job->submission_id,
                $job->stdout_file,
                $job->stderr_file,
                $job->retry_count,
                $job->temp_dir,
                $job->execution_host,
               );

 
  my $dbId = $sth->{'mysql_insertid'};
  $job->dbID( $dbId );
  $job->adaptor( $self );

  $sth->finish;

  $self->set_status( $job, "CREATED" );
}

=head2 remove

  Title   : remove
  Usage   : $jobadaptor->remove( $job )
  Function: deletes entries for job from database tables.
            deletes also history of status.
  Returns : throws exception when something goes wrong.
  Args    :

=cut

sub remove {
  my $self = shift;
  my $job = shift;

  if( ! $job->dbID ) {
    throw( "Cant remove job without dbID" );
  }
  my $dbID = $job->dbID;

  my $sth = $self->prepare(qq{
    DELETE FROM job
    WHERE  job_id = $dbID
  });
  $sth->execute;

  $sth = $self->prepare( qq{
    DELETE FROM job_status
    WHERE  job_id = $dbID
  });
  $sth->execute;
  $sth->finish;
}


=head2 remove_by_dbID

  Title   : remove_by_dbID
  Usage   : $jobadaptor->remove_by_dbID( $dbID )
  Function: deletes entries for job from database tables.
            deletes also history of status. Can take a list of ids.
  Returns : throws exception when something goes wrong.
  Args    :

=cut

sub remove_by_dbID {
  my $self = shift;
  my @dbIDs = @_;

  if( $#dbIDs == -1 ) { return }

  my $inExpr = "(".join(",", @dbIDs).")";

  my $sth = $self->prepare(qq{
    DELETE FROM job
    WHERE       job_id IN $inExpr
  });
  $sth->execute;
  $sth->finish;

  $sth = $self->prepare(qq{
    DELETE FROM job_status
    WHERE       job_id IN $inExpr
  });
  $sth->execute;
  $sth->finish;
}


=head2 update

  Title   : update
  Usage   : $job->update;
            $jobAdaptor->update($job);
            $jobAdaptor->update(@jobs)
  Function: a job which is already in db can update its contents
            it only updates stdout_file stderr_file, retry_count, 
            temp_dir, exec_host.
            and submission_id
  Returns : throws exception when something goes wrong.
  Args    : an array of Pipeline::Job objects

=cut

sub update {
  my ($self, @jobs) = @_;

  # only stdout, stderr, retry, submission_id and status are likely to be updated

  my $sth = $self->prepare(q{
                             UPDATE job
                             SET    stdout_file   = ?,
                             stderr_file   = ?,
                             retry_count   = ?,
                             submission_id = ?,
                             exec_host = ?,
                             temp_dir = ?
                             WHERE  job_id = ?
                            });

  foreach my $job (@jobs) {
    
    $sth->execute( $job->stdout_file,
                   $job->stderr_file,
                   $job->retry_count,
                   $job->submission_id,
                   $job->execution_host,
                   $job->temp_dir,
                   $job->dbID );
  }
  $sth->finish;
}



=head2 _objFromHashref

  Title   : _objFromHashref
  Usage   : my $job = $self->objFromHashref( $queryResult )
  Function: Creates a Job object from given hash reference.
            The hash contains column names and content of the column.
  Returns : the object or undef if that wasnt possible
  Args    : a hash reference

=cut

sub _objFromHashref {
  # create the appropriate job object

  my $self = shift;
  my $hashref = shift;
  my $job;
  my $analysis;
 
  
  $analysis =
    $self->db->get_AnalysisAdaptor->
      fetch_by_dbID( $hashref->{analysis_id} );
  $job = Bio::EnsEMBL::Pipeline::Job->new
    (
     '-dbobj'     => $self->db,
     '-adaptor'   => $self,
     '-id'        => $hashref->{'job_id'},
     '-submission_id' => $hashref->{'submission_id'},
     '-input_id'  => $hashref->{'input_id'},
     '-stdout'    => $hashref->{'stdout_file'},
     '-stderr'    => $hashref->{'stderr_file'},
     '-analysis'  => $analysis,
     '-retry_count' => $hashref->{'retry_count'},
     '-exec_host' => $hashref->{'exec_host'},
     '-temp_dir' => $hashref->{'temp_dir'},
  );
  return $job;
}


=head2 set_status

  Title   : set_status
  Usage   : my $status = $job->set_status
  Function: Sets the job status
  Returns : nothing
  Args    : Pipeline::Job Bio::EnsEMBL::Pipeline::Status

=cut

sub set_status {
    my ($self, $job, $stat_str) = @_;
    my $status;
    my $jobId;

    if( ! defined ($jobId = $job->dbID)) {
      throw( "Job has to be in database" );
    }


    eval {	
        my ($sth, $sth_upd, $sth_ins, $res);

        $sth_upd = $self->prepare(q{
                                UPDATE job_status
                                SET    is_current = 'n'
                                WHERE  job_id = ?
                                  AND  is_current = 'y'
                               });
        
        $sth_ins = $self->prepare(q{
                                INSERT into job_status
                                (job_id, status, time, is_current)
                                VALUES (?, ?, NOW(), 'y')
                               });

        $sth_upd->execute($jobId);
        $sth_upd->finish;

        $sth_ins->execute($jobId, $stat_str);
        $sth_ins->finish;
        
        $sth = $self->prepare("SELECT NOW()");
        $sth->execute();
        
        my $time = ($sth->fetchrow_arrayref())->[0];
        
        $status = Bio::EnsEMBL::Pipeline::Status->new
          (
           '-jobid'   => $jobId,
           '-status'  => $stat_str,
           '-created' => $time,
          );
        
        $self->current_status($job, $status);
        $sth->finish;
      };
    
    if ($@) {
      print( " $@ " );
      
      throw("Error setting status to $stat_str");
    } else {
      return $status;
    }
}


=head2 current_status

  Title   : current_status
  Usage   : my $status = $job->current_status
  Function: Get/set method for the current status
  Returns : Bio::EnsEMBL::Pipeline::Status
  Args    : Bio::EnsEMBL::Pipeline::Status

=cut

sub current_status {
    my ($self, $job, $arg) = @_;

    if (defined($arg))
    {
      throw("[$arg] is not a Bio::EnsEMBL::Pipeline::Status object")
        unless $arg->isa("Bio::EnsEMBL::Pipeline::Status");
      $job->{'_status'} = $arg;
    }
    else
      {
        throw("Can't get status if id not defined")
          unless defined($job->dbID);
        my $id =$job->dbID;
        my $sth = $self->prepare(q{
                                   SELECT status
                                   FROM   job_status
                                   WHERE  job_id = ?
                                   AND    is_current = 'y'
                                  });
        my $res = $sth->execute($id);
        my $status;
        while (my  $rowhash = $sth->fetchrow_hashref() ) {
          $status = $rowhash->{'status'};
        }
        $sth->finish;

        $sth = $self->prepare("SELECT NOW()");
        $res = $sth->execute();
        my $time;
        while (my $rowhash = $sth->fetchrow_arrayref()) {
          $time    = $rowhash->[0];
        }
        if(!$status){
          my ($p, $f, $l) = caller;
          warning("Have found no status for ".$job->dbID." ".
                      $job->input_id." ".$job->analysis->dbID.
                      " assuming is sucessful $f:$l\n");
          my $std = stack_trace_dump();
          print STDERR "$std\n";
          $status = 'SUCCESSFUL';
        }
        my $statusobj = Bio::EnsEMBL::Pipeline::Status->new(
                                                            '-jobid'   => $id,
                                                            '-status'  => $status,
                                                            '-created' => $time,
                                                           );
        $job->{'_status'} = $statusobj;
        $sth->finish;
      }
    return $job->{'_status'};
  }

=head2 get_all_status

  Title   : get_all_status
  Usage   : my @status = $job->get_all_status
  Function: Get all status objects associated with this job
  Returns : @Bio::EnsEMBL::Pipeline::Status
  Args    : Bio::EnsEMBL::Pipeline::Job

=cut

sub get_all_status {
  my ($self, $job) = @_;
  my @status;

  throw("Can't get status if id not defined")
    unless defined($job->dbID);

  my $sth = $self->prepare(q{
    SELECT   job_id, status, UNIX_TIMESTAMP(time)
    FROM     job_status
    WHERE    id = ?
    ORDER BY time desc
  });

  my $res = $sth->execute($job->dbID);

  while (my $rowhash = $sth->fetchrow_hashref() ) {
    my $time      = $rowhash->{'UNIX_TIMESTAMP(time)'};
    my $status    = $rowhash->{'status'};
    my $statusobj = Bio::EnsEMBL::Pipeline::Status->new(
      '-jobid'   => $job->dbID,
      '-status'  => $status,
      '-created' => $time,
    );

    push(@status,$statusobj);
  }

  $sth->finish;

  return @status;
}

=head2 get_last_status

  Title   : get_last_status
  Usage   : my @status = $job->get_all_status
  Function: Get most recent status object associated with this job
  Returns : Bio::EnsEMBL::Pipeline::Status
  Args    : Bio::EnsEMBL::Pipeline::Job, status string

=cut

sub get_last_status {
  my ($self, $job) = @_;

  throw("Can't get status if id not defined")
    unless defined($job->dbID);

  my $sth = $self->prepare (qq{
    SELECT job_id, status, UNIX_TIMESTAMP(time)
    FROM   job_status
    WHERE  is_current = 'y'
    AND    job_id = ?
  });

  my $res = $sth->execute($job->dbID);
  my $rowHashRef = $sth->fetchrow_hashref();
  if( ! defined $rowHashRef ) {
    return undef;
  }    

  my $time      = $rowHashRef->{'UNIX_TIMESTAMP(time)'};
  my $status    = $rowHashRef->{'status'};
  my $statusobj = Bio::EnsEMBL::Pipeline::Status->new(
    '-jobid'   => $job->dbID,
    '-status'  => $status,
    '-created' => $time,
  );

  $sth->finish;

  return $statusobj;
}

sub list_job_id_by_status {
  my ($self,$status) = @_;
  my @result;
  my @row;

  my $sth = $self->prepare(qq{
    SELECT   j.job_id
      FROM   job j, job_status js
     WHERE   j.job_id = js.job_id
       AND   js.status = '$status'
       AND   is_current = 'y'
    ORDER BY job_id
  });
  $sth->execute;

  while( @row = $sth->fetchrow_array ) {
    push( @result, $row[0] );
  }

  $sth->finish;
  return \@result;
}


sub list_job_id_by_status_age {
  my ($self,$status,$age) = @_;

  my @result;
  my @row;
  my $sth = $self->prepare(qq{
    SELECT   job_id
      FROM   job_status
       AND   status = '$status'
       AND   time < DATE_SUB(NOW(), INTERVAL $age MINUTE)
    ORDER BY job_id
  });
  $sth->execute;

  while( @row = $sth->fetchrow_array ) {
    push( @result, $row[0] );
  }

  $sth->finish;
  return \@result;
}


sub db {
  my ( $self, $arg )  = @_;
  if(  defined $arg ) {
      $self->{'_db'} = $arg;
  }
  $self->{'_db'};
}

sub prepare {
  my ( $self, $query ) = @_;
  $self->db->prepare( $query );
}

sub deleteObj {
  my ($self) = @_;
  my @dummy = values %{$self};
  foreach my $key ( keys %$self ) {
    delete $self->{$key};
  }
  foreach my $obj ( @dummy ) {
    eval {
      $obj->deleteObj;
    }
  }
}


sub lock_tables{
  my ($self) = @_;
  
  my $sql = "LOCK TABLES job READ, job_status READ";

  my $sth = $self->db->prepare($sql);

  $sth->execute;
  $sth->finish;

}

sub unlock_tables{
  my ($self) = @_;
  
  my $sql = "UNLOCK TABLES";

  my $sth = $self->db->prepare($sql);

  $sth->execute;
  $sth->finish;
}


1;
