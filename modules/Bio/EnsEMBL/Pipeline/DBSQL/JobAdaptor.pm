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
use Bio::EnsEMBL::Pipeline::Status;

use vars qw(@ISA);
use strict;

@ISA = qw( Bio::Root::RootI );


=head2 new

  Arg [1]   : dbadpator object 
  Function  : makes a new Jobadaptor object
  Returntype: Bio::EnsEMBL::Pipeline::DBSQL::JobAdaptor;
  Exceptions: none 
  Caller    : Bio::EnsEMBL::Pipeline::DBSQL::Jobadaptor
  Example   : Bio::EnsEMBL::Pipeline::DBSQL::Jobadaptor->new($db);

=cut


sub new {
  my ($class,$dbobj) = @_;
  my $self = $class->SUPER::new();
  
  $self->db( $dbobj );
  return $self;
}

###############
#fetching jobs#
###############

=head2 fetch_by_dbID

  Arg [1]   : dbId 
  Function  : gets all the information about a particualr job from the job table
  Returntype: hash ref
  Exceptions: none
  Caller    : Bio::EnsEMBL::Pipeline::DBSQL::JobAdaptor
  Example   : 

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


=head2 fetch_by_Status_Analysis

  Arg [1]   : scalar of status
  Function  : returns an array of hashes about each jobs with a particular status and analysis type
  Returntype: array of Bio::EnsEMBL::Pipeline::Job
  Exceptions: throws if no analysis object or status given
  Caller    : Bio::EnsEMBL::Pipeline::DBSQL::JobAdaptor
  Example   : 

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


=head2 fetch_by_Age

  Arg [1]   : scalar of age 
  Function  : selects jobs from database with specified age
 Returntype: array of Bio::EnsEMBL::Pipeline::Job
  Exceptions: thows if no ages is defined
  Caller    : Bio::EnsEMBL::Pipeline::DBSQL::JobAdaptor
  Example   : 

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



################
#SETTING STATUS#
################


=head2 setting statuses

  Arg [1]   : job object  
  Function  : the following methods set the status of single jobs or return the status of many jobs depending on the arguments passed in
 Returntype: single or array of Bio::EnsEMBL::Pipeline::Status object
  Exceptions: they throw errors if no job object or no status object are passed in
  Caller    : Bio::EnsEMBL::Pipeline::JobAdaptor
  Example   : $self->adaptor->set_status($self, 'submitted'); (from job.pm)

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
        
#       print("Status for job [" . $job->dbID . "] set to " . $status->status . "\n");
    };

    if ($@) {
      print( " $@ " );

        $self->throw("Error setting status to $arg");
    } else {
        return $status;
    }
}




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

#these both return a list of jobids depending on the status and for the 2nd method the age of the job

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

###################
#database editting#
###################


=head2 store

  Arg [1]   : Job object 
  Function  : stores a job object to the database
  Returntype: none
  Exceptions: throws if job doesn't have an analysis object or analysis object doesn't have an id
  Caller    : Bio::EnsEMBL::Pipeline:DBSQL::JobAdaptor
  Example   : 

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

  Arg [1]   : Job object 
  Function  : removes given job from db
  Returntype: none
  Exceptions: throws if job doesn't have a dbid '
  Caller    : Bio:EnsEMBL::Pipeline::DBSQL::JobAdaptor
  Example   : 

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

  Arg [1]   : array of dbIds
  Function  : delete jobs of given dbId from db
  Returntype: none
  Exceptions: none
  Caller    : Bio::EnsEMBL::Pipeline::DBSQL::JobAdaptor
  Example   : 

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

  Arg [1]   : Job object 
  Function  : updates the stdout, stderr, retry, LSF_id of a job in the job table
  Returntype: none
  Exceptions: none
  Caller    : Bio::EnsEMBL::Pipeline::DBSQL::JobAdaptor
  Example   : 

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

#######################
#MISCELLANEOUS METHODS#
#######################


=head2 create_tables

  Arg [1]   :none 
  Function  :drops existing job tables and creates new ones
  Returntype: none
  Exceptions: none
  Caller    : Bio::EnsEMBL::Pipeline::DBSQL::JobAdaptor
  Example   : 

=cut


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

#accessor method which defines the db adaptor

sub db {
  my ( $self, $arg )  = @_;
  if(  defined $arg ) {
      $self->{'_db'} = $arg;
  }
  $self->{'_db'};
}

#calls the prepare method of the adaptor
sub prepare {
  my ( $self, $query ) = @_;
  $self->db->prepare( $query );
}

#was used when object files were being used is no longer required

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

#presumable was to check if a job existed but hasn't been implemented

sub exists {
  my $self = shift;
  my $hashref = shift;

  $self->throw( "Not implemented yet" );
}

#creates a job object given an hashref from the job table

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

1;
