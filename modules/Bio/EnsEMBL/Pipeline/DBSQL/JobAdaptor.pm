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
use Bio::Root::RootI;

use vars qw(@ISA);
use strict;

@ISA = qw( Bio::Root::RootI );

sub new {
  my ($class,$dbobj) = @_;
  my $self = $class->SUPER::new();
  
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

  my $sth = $self->prepare( q{
    SELECT jobId, input_id, class, analysisId, LSF_id, object_file,
      stdout_file, stderr_file, retry_count
    FROM job
    WHERE jobId = ? } );
  
  $sth->execute( $id );
  my $rowHashRef = $sth->fetchrow_hashref;
  if( ! defined $rowHashRef ) {
    return undef;
  }

  return $self->_objFromHashref( $rowHashRef );
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

    $self->throw("Require status and analysis id for fetch_by_Status_Analysis") 
                            unless ($analysis && $status); 
    if( ! defined $analysis->dbID ){
       $self->throw( "Analysis needs to be in database" );
    }
    my $analysisId = $analysis->dbID;

    my $query = "select j.jobId, j.input_id, j.class, j.analysisId, j.LSF_id," .
	                   "j.stdout_file, j.stderr_file,".
                       "j.object_file, j.retry_count".
                       "j.status_file " . 
	            "from job as j, current_status as cs, jobstatus as js ". 
                "where j.jobId = cs.jobId and js.jobId = cs.jobId and ". 
                       "cs.status = js.status and ".
                       "j.analysis = $analysisId and ".
                       "cs.status = \'$status\' ".
                       "order by js.time DESC";
    $query .= " limit $start, $end" if ($start && $end);
                 
    my $sth = $self->prepare($query);
    my $res = $sth->execute();
    
    my @jobs;
    
    while (my $row = $sth->fetchrow_hashref) 
    {
	    my $job = $self->_objFromHashref($row);
	    push(@jobs,$job);
    }
    return @jobs;
}


=head2 fetch_by_Age {

  Title   : fetch_by_Age
  Usage   : my @jobs = $db->fetch_by_Age($duration)
  Function: Retrieves all jobs in the database
            that are older than than a certain duration given in minutes.
  Returns : @Bio::EnsEMBL::Pipeline::Job
  Args    : int

=cut

sub fetch_by_Age {
    my ($self,$age) = @_;

    $self->throw("No input status for get_JobsByAge") 
        unless defined($age);
    #convert age from minutes to seconds

    my $query = 'SELECT j.jobId, j.input_id, j.class, j.analysisId, j.LSF_id, '
                .'j.stdout_file, j.stderr_file, j.object_file, '
                .'j.retry_count '     
                .'FROM job as j, jobstatus as js, current_status as cs ' 
                .'WHERE cs.jobId = js.jobId '
                    .'AND cs.status = js.status '
                    .'AND cs.jobId = j.jobId '
                    ."AND js.time < DATE_SUB( NOW(), INTERVAL $age MINUTE )";
            
    my $sth = $self->prepare($query);
    my $res = $sth->execute();
    
    my @jobs;

    while (my $row = $sth->fetchrow_hashref) {
	my $job = $self->_objFromHashref($row);
	push(@jobs,$job);
    }

    return @jobs;
}


=head2 fetch_by_inputId

  Title   : fetch_by_inputId
  Usage   : my @job = $adaptor->fetch_by_inputId
  Function: Retrieves all jobs from adaptor with certain input id
  Returns : list of job objects
            throws exception when something goes wrong.
  Args    : 

=cut

sub fetch_by_inputId {
  my $self = shift;
  my $inputid = shift;
  my @result;

  my $sth = $self->prepare( q{
    SELECT jobId, input_id, class, analysisId, LSF_id, object_file,
      stdout_file, stderr_file, retry_count
    FROM job
    WHERE input_id = ? } );
  
  $sth->execute( $inputid );
  while( my $rowHashRef = $sth->fetchrow_hashref ) {
    push( @result, $self->_objFromHashref( $rowHashRef ));
  }

  return @result;
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
    $self->throw( "Need to store analysis first" );
  }

  my $sth = $self->prepare( q{
    INSERT into job( input_id, class, analysisId,
      LSF_id, stdout_file, stderr_file, object_file,
      retry_count ) 
    VALUES ( ?, ?, ?, ?, ?, ?, ?, ? ) } );

  $sth->execute( $job->input_id,
                 $job->class,
                 $job->analysis->dbID,
                 $job->LSF_id,
                 $job->stdout_file,
                 $job->stderr_file,
                 $job->input_object_file,
                 $job->retry_count );

  $sth = $self->prepare( "SELECT LAST_INSERT_ID()" );
  $sth->execute;

  my $dbId = ($sth->fetchrow_arrayref)->[0];
  $job->dbID( $dbId );
  $job->adaptor( $self );

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

  if( ! defined $job->dbID ) {
    $self->throw( "Cant remove job without dbID" );
  }
  my $dbID = $job->dbID;

  my $sth = $self->prepare( qq{
    DELETE FROM job
     WHERE jobId = $dbID } );
  $sth->execute;

  $sth = $self->prepare( qq{
    DELETE FROM current_status
     WHERE jobId=$dbID } );
  $sth->execute;

  $sth = $self->prepare( qq{
    DELETE FROM jobstatus
     WHERE jobId = $dbID } );
  $sth->execute;
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
  
  my $inExpr = "(".join( ",",@dbIDs ).")";
  
  my $sth = $self->prepare( qq{
    DELETE FROM job
     WHERE jobId IN $inExpr } );
  $sth->execute;

  $sth = $self->prepare( qq{
    DELETE FROM current_status
     WHERE jobId IN $inExpr } );
  $sth->execute;
  $sth = $self->prepare( qq{
    DELETE FROM jobstatus
     WHERE jobId IN $inExpr } );
  $sth->execute;
}


=head2 update

  Title   : update
  Usage   : $job->update; $jobAdaptor->update( $job )
  Function: a job which is already in db can update its contents
            it only updates stdout_file, stderr_file, retry_count
            and LSF_id
  Returns : throws exception when something goes wrong.
  Args    : 

=cut

sub update {
  my $self = shift;
  my $job = shift;
  
  # only stdout, stderr, retry, LSF_id and status are likely to be updated

  my $sth = $self->prepare( q{
    UPDATE job
       SET stdout_file = ?,
           stderr_file = ?,
           object_file = ?,
           retry_count = ?,
           LSF_id = ?
     WHERE jobId = ? } );

  $sth->execute( $job->stdout_file,
		 $job->stderr_file,
		 $job->input_object_file,
		 $job->retry_count,
		 $job->LSF_id,
		 $job->dbID );
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
      fetch_by_dbID( $hashref->{analysisId} );

  $job = Bio::EnsEMBL::Pipeline::Job->new
  (
   '-dbobj'    => $self->db,
   '-adaptor'  => $self,
   '-id'       => $hashref->{'jobId'},
   '-lsf_id'   => $hashref->{'LSF_id'},
   '-input_id' => $hashref->{'input_id'},
   '-class'    => $hashref->{'class'},
   '-stdout'   => $hashref->{'stdout_file'},
   '-stderr'   => $hashref->{'stderr_file'},
   '-input_object_file' => $hashref->{'object_file'},
   '-analysis' => $analysis,
   '-retry_count' => $hashref->{'retry_count'}
  );

  return $job;
}


# provide a hashref
# each value in it used for a combined query, if the described object is in
# returns a job object, if it is in, else undef

sub exists {
  my $self = shift;
  my $hashref = shift;

  $self->throw( "Not implemented yet" );
}

# Code directly from Michele

=head2 set_status

  Title   : set_status
  Usage   : my $status = $job->set_status
  Function: Sets the job status
  Returns : nothing
  Args    : Pipeline::Job Bio::EnsEMBL::Pipeline::Status

=cut

sub set_status {
    my ($self,$job,$arg) = @_;

    if( ! defined $job->dbID ) {
      $self->throw( "Job has to be in database" );
    }

    my $status;

    eval {	
	my $sth = $self->prepare("insert delayed into jobstatus(jobId,status,time) values (" .
					 $job->dbID . ",\"" .
					 $arg      . "\"," .
					 "now())");
	my $res = $sth->execute();

	$sth = $self->prepare("replace into current_status(jobId,status) values (" .
				      $job->dbID . ",\"" .
				      $arg      . "\")");

	$res = $sth->execute();
	
	$sth = $self->prepare("select now()" );
	
	$res = $sth->execute();
	
	my $rowhash = $sth->fetchrow_arrayref();
	my $time    = $rowhash->[0];

	$status = Bio::EnsEMBL::Pipeline::Status->new
	  (  '-jobid'   => $job->dbID,
	     '-status'  => $arg,
	     '-created' => $time,
	  );
	
	$self->current_status($job, $status);
	
#	print("Status for job [" . $job->dbID . "] set to " . $status->status . "\n");
    };

    if ($@) {
#      print( " $@ " );

	$self->throw("Error setting status to $arg");
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
	$self->throw("[$arg] is not a Bio::EnsEMBL::Pipeline::Status object") 
	    unless $arg->isa("Bio::EnsEMBL::Pipeline::Status");
	$job->{'_status'} = $arg;
    }
    else 
    {
	$self->throw("Can't get status if id not defined") 
	    unless defined($job->dbID);
	my $id =$job->dbID;
	my $sth = $self->prepare
	    ("select status from current_status where jobId=$id");
	my $res = $sth->execute();
	my $status;
	while (my  $rowhash = $sth->fetchrow_hashref() ) {
	    $status = $rowhash->{'status'};
	}

	$sth = $self->prepare("select now()");
	$res = $sth->execute();
	my $time;
	while (my  $rowhash = $sth->fetchrow_hashref() ) {
	    $time    = $rowhash->{'now()'};
	}
	my $statusobj = new Bio::EnsEMBL::Pipeline::Status
	    ('-jobid'   => $id,
	     '-status'  => $status,
	     '-created' => $time,
	     );
	$job->{'_status'} = $statusobj;
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
  
  $self->throw("Can't get status if id not defined") 
    unless defined($job->dbID);

  my $sth = $self->prepare
    ("select jobId,status, UNIX_TIMESTAMP(time) from  jobstatus " . 
     "where id = \"" . $job->dbID . "\" order by time desc");
  
  my $res = $sth->execute();
  
  my @status;
  while (my $rowhash = $sth->fetchrow_hashref() ) {
    my $time      = $rowhash->{'UNIX_TIMESTAMP(time)'};#$rowhash->{'time'};
    my $status    = $rowhash->{'status'};
    my $statusobj = new Bio::EnsEMBL::Pipeline::Status(-jobid   => $job->dbID,
						       -status  => $status,
						       -created => $time,
						      );
                               
    push(@status,$statusobj);
    
  }
  
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

  $self->throw("Can't get status if id not defined")
    unless defined($job->dbID);

  my $sth = $self->prepare (qq{
    SELECT js.jobId, cs.status, UNIX_TIMESTAMP(time)
      FROM jobstatus js, current_status cs
     WHERE js.jobId = cs.jobId
       AND js.status = cs.status
       AND js.jobId = ?} );

  my $res = $sth->execute($job->dbID);
  my $rowHashRef = $sth->fetchrow_hashref();
  if( ! defined $rowHashRef ) {
    return undef;
  }

  my $time      = $rowHashRef->{'UNIX_TIMESTAMP(time)'};#$rowhash->{'time'};
  my $status    = $rowHashRef->{'status'};
  my $statusobj = new Bio::EnsEMBL::Pipeline::Status(-jobid   => $job->dbID,
						     -status  => $status,
						     -created => $time,
						     );
  return $statusobj;
}

sub list_jobId_by_status {
  my ($self,$status) = @_;
  my @result;
  my @row;

  my $sth = $self->prepare( qq{
    SELECT j.jobId
      FROM job j, current_status c
     WHERE j.jobId = c.jobId
       AND c.status = '$status'
     ORDER BY jobId } );
  $sth->execute;
  
  while( @row = $sth->fetchrow_array ) {
    push( @result, $row[0] );
  }
  
  return @result;
}


sub list_jobId_by_status_age {
  my ($self,$status,$age) = @_;
  
  my @result;
  my @row;
  my $sth = $self->prepare( qq{
    SELECT js.jobId
      FROM current_status c, jobstatus js
     WHERE js.jobId = c.jobId
       AND c.status = '$status'
       AND js.status = '$status'
       AND js.time < DATE_SUB( NOW(), INTERVAL $age MINUTE )
     ORDER BY jobId } );
  $sth->execute;
  
  while( @row = $sth->fetchrow_array ) {
    push( @result, $row[0] );
  }
  
  return @result;
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

# creates all tables for this adaptor - job, jobstatus and current_status
# if they exist they are emptied and newly created
sub create_tables {
  my $self = shift;
  my $sth;

  $sth = $self->prepare("drop table if exists job");
  $sth->execute();

  $sth = $self->prepare(qq{
    CREATE TABLE job (
    job_id        int(10) unsigned  default 0  not null auto_increment,
    input_id      varchar(40)       default '' not null,
    class         enum("clone","contig","vc","gene") not null,
    analysis_id   int(10) unsigned  default 0  not null,
    LSF_id        int(10) unsigned  default 0,
    stdout_file   varchar(100)      default '' not null,
    stderr_file   varchar(100)      default '' not null,
    object_file   varchar(100)      default '' not null,
    retry_count   int               default 0,

    PRIMARY KEY   (job_id),
    KEY input     (input_id),
    KEY analysis  (analysis_id)
    );
  });
  $sth->execute();

  $sth = $self->prepare("drop table if exists jobstatus");
  $sth->execute();

  $sth = $self->prepare(qq{
    CREATE TABLE jobstatus (
    job_id     int(10) unsigned  default 0 not null,
    status     varchar(40)       default 'CREATED' not null,
    time       datetime          default '0000-00-00 00:00:00' not null,

    KEY job    (job_id),
    KEY status (status)
    );
  });
  $sth->execute();

  $sth = $self->prepare("drop table if exists current_status");
  $sth->execute();

  $sth = $self->prepare(qq{
    CREATE TABLE current_status (
    job_id  int(10) unsigned  default 0 not null,
    status  varchar(40)       default '' not null,

    PRIMARY KEY (job_id),
    KEY status  (status)
    );
  });
  $sth->execute();
  $sth->finish;
}

1;
