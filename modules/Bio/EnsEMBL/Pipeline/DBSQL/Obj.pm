#
# Object for storing the connection to the analysis database
#
# Cared for by Michele Clamp  <michele@sanger.ac.uk>
#
# Copyright Michele Clamp
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::DB::Obj

=head1 SYNOPSIS

=head1 DESCRIPTION

Interface for the connection to the analysis database

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::DBSQL::Obj;


use vars qw(@ISA);
use strict;
use DBI;

use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::Pipeline::DBSQL::RuleAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::AnalysisAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::JobAdaptor;

use Bio::Root::Object;

# use FreezeThaw qw(freeze thaw);

# Inherits from the base bioperl object

@ISA = qw(Bio::EnsEMBL::DBSQL::Obj Bio::Root::Object);

sub _initialize {
    my ($self,@args) = @_;

    my $make = $self->SUPER::_initialize(@args);

    return $make; # success - we hope!

}

=head2 get_JobAdaptor

 Title   : get_JobAdaptor
 Usage   : $db->get_JobAdaptor
 Function: The Adaptor for Job objects in this db
 Example : 
 Returns : Bio::EnsEMBL::Pipeline::DBSQL::JobAdaptor
 Args    : 


=cut

sub get_JobAdaptor {
  my ($self) = @_;

  if( ! defined $self->{_JobAdaptor} ) {
    $self->{_JobAdaptor} = Bio::EnsEMBL::Pipeline::DBSQL::JobAdaptor->new
      ( $self );
  }

  return $self->{_JobAdaptor};
}

=head2 get_AnalysisAdaptor

 Title   : get_AnalysisAdaptor
 Usage   : $db->get_AnalysisAdaptor
 Function: The Adaptor for Analysis objects in this db
 Example : 
 Returns : Bio::EnsEMBL::Pipeline::DBSQL::AnalysisAdaptor
 Args    : 


=cut

sub get_AnalysisAdaptor {
  my ($self) = @_;

  if( ! defined $self->{_AnalysisAdaptor} ) {
    $self->{_AnalysisAdaptor} = Bio::EnsEMBL::Pipeline::DBSQL::AnalysisAdaptor->new
      ( $self );
  }

  return $self->{_AnalysisAdaptor};
}


=head2 get_RuleAdaptor

 Title   : get_RuleAdaptor
 Usage   : $db->get_RuleAdaptor
 Function: The Adaptor for Rule objects in this db
 Example : 
 Returns : Bio::EnsEMBL::Pipeline::DBSQL::RuleAdaptor
 Args    : 


=cut

sub get_RuleAdaptor {
  my ($self) = @_;

  if( ! defined $self->{_RuleAdaptor} ) {
    $self->{_RuleAdaptor} = Bio::EnsEMBL::Pipeline::DBSQL::RuleAdaptor->new
      ( $self );
  }

  return $self->{_RuleAdaptor};
}

sub get_JobsByCurrentStatus {
    my ($self,$status) = @_;


    $self->throw("No input status for get_JobsByCurrentStatus") unless defined($status);

    my $query = "select j.id,j.input_id,j.analysis,j.LSF_id,j.machine,j.queue," .
   	        " j.stdout_file,j.stderr_file,j.input_object_file,j.output_file,j.status_file,cs.status " . 
		" from job as j, current_status as cs where cs.id = j.id and cs.status = \"" . 
		    $status . "\""; 
    
    my $sth = $self->prepare($query);
    my $res = $sth->execute();

    my @jobs;

    while (my $row = $sth->fetchrow_hashref) {
	my $job = $self->_parseJob($row);
	push(@jobs,$job);
    }

    return @jobs;
}

=head2 get_JobsByAge {

  Title   : get_JobsByAge
  Usage   : my @jobs = $db->get_JobsByAge($duration)
  Function: Retrieves all jobs in the database
            that are older than than a certain duration given in minutes.
  Returns : @Bio::EnsEMBL::Pipeline::DB::JobI
  Args    : int

=cut

sub get_JobsByAge {
    my ($self,$age) = @_;

    $self->throw("No input status for get_JobsByAge") 
        unless defined($age);
    #convert age from minutes to seconds
    my $ageinseconds = $age * 60;
    my $query = 'SELECT cs.id, j.input_id, j.analysis, j.LSF_id, j.machine, '
                    .'j.queue, j.stdout_file, j.stderr_file, j.input_object_file, '
                    .'j.output_file, j.status_file '
                .'FROM job as j, jobstatus as js, current_status as cs ' 
                .'WHERE cs.id = js.id '
                    .'AND cs.status = js.status '
                    .'AND cs.id = j.id '
                    .'AND Unix_TIMESTAMP(js.time) > Unix_TIMESTAMP()-'.$ageinseconds;
            
    my $sth = $self->prepare($query);
    my $res = $sth->execute();
    
    my @jobs;

    while (my $row = $sth->fetchrow_hashref) {
	my $job = $self->_parseJob($row);
	push(@jobs,$job);
    }

    return @jobs;
}

sub get_JobsByAnalysis {
    my ($self,$analysis) = @_;

    
    $self->throw("No input analysis for get_JobsByCurrentAnalysis") unless defined($analysis);
    
    my $query = "select id,input_id,analysis,LSF_id,machine,queue," .
	            "stdout_file,stderr_file,input_object_file,output_file,status_file " . 
	            " from job where analysis = $analysis"; 
    my $sth = $self->prepare($query);
    my $res = $sth->execute();
    
    my @jobs;
    
    while (my $row = $sth->fetchrow_hashref) {
	my $job = $self->_parseJob($row);
	push(@jobs,$job);
    }

    return @jobs;
}

=head2 get_JobsbyStatus_and_Analysis {

  Title   : get_JobsbyStatus_and_Analysis
  Usage   : my @jobs = $db->get_JobsbyStatus_and_Analysis($id, $status)
  Function: Retrieves all jobs in the database matching status and 
            an analysis id
  Returns : @Bio::EnsEMBL::Pipeline::DB::JobI
  Args    : int analysis id, string status, and optional start and end limits

=cut

sub get_JobsbyStatus_and_Analysis {
    my ($self,$status, $analysis, $start, $end) = @_;

    $self->throw("Require status and analysis id for get_JobsbyStatus_and_Analysis") 
                            unless ($analysis && $status); 

    my $query = "select j.id, j.input_id, j.analysis, j.LSF_id," .
	                   "j.machine,queue, j.stdout_file, j.stderr_file,".
                       "j.input_object_file, j.output_file,".
                       "j.status_file " . 
	            "from job as j, current_status as cs, jobstatus as js ". 
                "where j.id = cs.id and js.id = cs.id and ". 
                       "cs.status = js.status and ".
                       "j.analysis = $analysis and ".
                       "cs.status = \'$status\' ".
                       "order by js.time DESC";
    $query .= " limit $start, $end" if ($start && $end);
                 
    my $sth = $self->prepare($query);
    my $res = $sth->execute();
    
    my @jobs;
    
    while (my $row = $sth->fetchrow_hashref) 
    {
	    my $job = $self->_parseJob($row);
	    push(@jobs,$job);
    }
    return @jobs;
}

sub _parseJob {
    my ($self,$row) = @_;

    $self->throw("No row object input") unless defined($row);

    my $jobid             = $row->{id};
    my $input_id          = $row->{input_id};
    my $analysis_id       = $row->{analysis};
    my $LSF_id            = $row->{LSF_id};
    my $machine           = $row->{machine};
    my $object            = $row->{object};
    my $queue             = $row->{queue};
    my $stdout            = $row->{stdout_file};
    my $stderr            = $row->{stderr_file};
    my $input_object_file = $row->{input_object_file};
    my $output_file       = $row->{output_file};
    my $status_file       = $row->{status_file};
    
    my $analysis          = $self->get_Analysis($analysis_id);
    
       $analysis->id($analysis_id);
    
    my $job = new Bio::EnsEMBL::Pipeline::DBSQL::Job(-id       => $jobid,
						     -input_id => $input_id,
						     -analysis => $analysis,
						     -LSF_id   => $LSF_id,
						     -machine  => $machine,
						     -object   => $object,
						     -queue    => $queue,
						     -dbobj    => $self,
						     -stdout   => $stdout,
						     -stderr   => $stderr,
						     -input_object_file => $input_object_file,
						     -output_file => $output_file,
                                                     -status_file => $status_file,        
						     );
    
    return $job;
}


sub delete_Job {
    my ($self,$id) = @_;

    $self->throw("No job id for delete_Job") unless defined($id);
    my $job = $self->get_Job($id);
    my $query = "delete from job where id = $id";
    my $sth   = $self->prepare($query);
    my $res   = $sth ->execute;

    $query = "delete from jobstatus where id = $id";
    $sth   = $self->prepare($query);
    $res   = $sth->execute;

    $query = "delete from current_status where id = $id";
    $sth  = $self->prepare($query);
    $res  = $sth->execute;
    
    unlink $job->status_file;
    unlink $job->stdout_file;
    unlink $job->stderr_file;
    unlink $job->input_object_file;
}

=head2 get_all_Status {

  Title   : get_all_Status
  Usage   : my @status_list = $db->get_all_Status()
  Function: Retrieves list of all status present in jobstatus
  Returns : array of status present in DB
  Args    : none

=cut

sub get_all_Status {
    my ($self) = @_;

    my $query = 'SELECT status from jobstatus group by status ';
            
    my $sth = $self->prepare($query);
    my $res = $sth->execute();
    
    my @statuslist;
    
    while (my $row = $sth->fetchrow_hashref) 
    {
	    push(@statuslist,$row->{'status'});
    }
    return @statuslist;
}


sub get_InputIdsByAnalysis {
    my ($self,$analid) = @_;
    
    $self->throw("No input analysis id  for get_InputIdsByAnalysis") unless defined($analid);

    my $sth = $self->prepare("select distinct input_id  from job where analysis = $analid");
    my $res = $sth->execute();
    my $row = $sth->fetchrow_hashref;

    my @input_ids;

    while ($row = $sth->fetchrow_hashref) {
	my $input_id       = $row->{input_id};
	print("id is $input_id\n");
	push(@input_ids,$input_id);
    }
    
    return @input_ids;
}




=head2 prepare

 Title   : prepare
 Usage   : $sth = $dbobj->prepare("select * from job where id = $id");
 Function: prepares a SQL statement on the DBI handle
 Example :
 Returns : A DBI statement handle object
 Args    : a SQL string


=cut

sub prepare{
   my ($self,$string) = @_;

   if( ! $string ) {
       $self->throw("Attempting to prepare an empty SQL query!");
   }
   
   return $self->_db_handle->prepare($string);
}

=head2 _db_handle

 Title   : _db_handle
 Usage   : $sth = $dbobj->_db_handle($dbh);
 Function: Get/set method for the database handle
 Example :
 Returns : A database handle object
 Args    : A database handle object


=cut

sub _db_handle {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_db_handle} = $arg;
    }
    return $self->{_db_handle};
}



=head2 _lock_tables

 Title   : _lock_tables
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _lock_tables{
   my ($self,@tables) = @_;
   
   my $state;
   foreach my $table ( @tables ) {
       if( $self->{'_lock_table_hash'}->{$table} == 1 ) {
	   $self->warn("$table already locked. Relock request ignored");
       } else {
	   if( $state ) { $state .= ","; } 
	   $state .= "$table write";
	   $self->{'_lock_table_hash'}->{$table} = 1;
       }
   }

   my $sth = $self->prepare("lock tables $state");
   my $rv = $sth->execute();
   $self->throw("Failed to lock tables $state") unless $rv;

}

=head2 _unlock_tables

 Title   : _unlock_tables
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut

sub _unlock_tables{
   my ($self,@tables) = @_;

   my $sth = $self->prepare("unlock tables");
   my $rv  = $sth->execute();
   $self->throw("Failed to unlock tables") unless $rv;
   %{$self->{'_lock_table_hash'}} = ();
}


=head2
 Title   : DESTROY
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut


sub DESTROY{
   my ($obj) = @_;

   $obj->_unlock_tables();

   if( $obj->{'_db_handle'} ) {
       $obj->{'_db_handle'}->disconnect;
       $obj->{'_db_handle'} = undef;
   }
}

sub deleteObj {
  my $self = shift;
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




