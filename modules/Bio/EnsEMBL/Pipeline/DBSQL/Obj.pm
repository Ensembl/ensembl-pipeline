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

use Bio::EnsEMBL::Pipeline::DB::ObjI;
use Bio::EnsEMBL::Pipeline::DBSQL::Job;
use FreezeThaw qw(freeze thaw);

# Inherits from the base bioperl object
@ISA = qw(Bio::EnsEMBL::Pipeline::DB::ObjI Bio::Root::Object);

sub _initialize {
    my ($self,@args) = @_;

    my $make = $self->SUPER::_initialize;

    my ($db,$host,$driver,$user,$password,$debug) = 
	$self->_rearrange([qw(DBNAME
			      HOST
			      DRIVER
			      USER
			      PASS
			      DEBUG
			      )],@args);
    
    $db   || $self->throw("Database object must have a database name");
    $user || $self->throw("Database object must have a user");
    
    if( ! $driver ) {
	$driver = 'mysql';
    }
    
    if( ! $host ) {
	$host = 'localhost';
    }
    
    my $dsn = "DBI:$driver:database=$db;host=$host";
    my $dbh = DBI->connect("$dsn","$user",$password);

    $dbh || $self->throw("Could not connect to database $db user $user using [$dsn] as a locator");
      
    $self->_db_handle($dbh);
    
    return $make; # success - we hope!

}



=head2 get_Job {

  Title   : get_Job
  Usage   : my $job = $db->get_Job($id)
  Function: Retrieves a job from the database
            when given its id
  Returns : Bio::EnsEMBL::Pipeline::DB::JobI
  Args    : int

=cut


sub get_Job {
    my ($self,$id) = @_;

    $self->throw("No input id for get_Job") unless defined($id);

    my $sth = $self->prepare("select id,input_id,analysis,LSF_id,machine,object,queue " . 
			     "from job " . 
			     "where id = $id");
    my $res = $sth->execute();
    my $row = $sth->fetchrow_hashref;

    my $input_id    = $row->{input_id};
    my $analysis_id = $row->{analysis};
    my $LSF_id      = $row->{LSF_id};
    my $machine     = $row->{machine};
    my $object      = $row->{object};
    my $queue       = $row->{queue};

    my $analysis    = $self->get_Analysis($analysis_id);

    my $job;

    if (defined($object)) {
	print(STDERR "Recreating job object from stored version\n");
	
	($job)  = FreezeThaw::thaw($object);

	if (! $job->isa("Bio::EnsEMBL::Pipeline::DB::JobI")) {
	    $self->throw("Object string didn't return a Bio::EnsEMBL::Pipeline::DB::JobI object [" . ref($job) ."]");
	}

	$job->_dbobj($self);

    } else {
	print(STDERR "Creating new job object - no stored version\n");
	
	$job = new Bio::EnsEMBL::Pipeline::DBSQL::Job(-id => $id,
						      -input_id => $input_id,
						      -analysis => $analysis,
						      -LSF_id   => $LSF_id,
						      -machine  => $machine,
						      -queue    => $queue,
						      );
	$job->_dbobj($self);
    }
    
    return $job;
}


=head2 get_JobsByInputId {

  Title   : get_JobsByInputId
  Usage   : my @jobs = $db->get_JobsByInputId($id)
  Function: Retrieves all jobs in the database
            that have a certain input id.
	    Input id will usually be associated with
            a sequence in the main ensembl database.
  Returns : @Bio::EnsEMBL::Pipeline::DB::JobI
  Args    : int

=cut


sub get_JobsByInputId {
    my ($self,$id) = @_;

    $self->throw("No input input_id for get_JobsByInputId") unless defined($id);

    my $sth = $self->prepare("select id,input_id,analysis,LSF_id,machine,object,queue " . 
			     "from job " . 
			     "where input_id = $id");
    my $res = $sth->execute();
    my $row = $sth->fetchrow_hashref;

    my @jobs;

    while ($row = $sth->fetchrow_hashref) {
	my $input_id    = $row->{input_id};
	my $analysis_id = $row->{analysis};
	my $LSF_id      = $row->{LSF_id};
	my $machine     = $row->{machine};
	my $object      = $row->{object};
	my $queue       = $row->{queue};
	
	my $analysis    = $self->get_Analysis($analysis_id);
	
	my $job = new Bio::EnsEMBL::Pipeline::DBSQL::Job(-id => $id,
							 -input_id => $input_id,
							 -analysis => $analysis,
							 -LSF_id   => $LSF_id,
							 -machine  => $machine,
							 -object   => $object,
							 -queue    => $queue,
							 );
	
	push(@jobs,$job);
    }

    return @jobs;
}


=head2 write_Analysis {

  Title   : write_Analysis
  Usage   : $db->write_Analysis($analysis)
  Function: Write an analysis object into the database
            with a check as to whether it exists.
  Returns : nothing
  Args    : int

=cut


sub write_Analysis {
    my ($self,$analysis) = @_;

    $self    ->throw("Input must be Bio::EnsEMBL::Pipeline::Analysis [$analysis]") unless
    $analysis->isa("Bio::EnsEMBL::Pipeline::Analysis");


    my ($id,$created) = $self->exists_Analysis($analysis);

    if (defined($id)) {
	$analysis->id     ($id);
	$analysis->created($created);
	return $analysis;
    }
    my $query = 
	                     "insert into analysis(id,created,db,db_version,db_file," .
                	     "program,program_version,program_file," .
			     "parameters,module,module_version," .
			     "gff_source,gff_feature) " .
			     " values (NULL," . 
			     "now()"                    .   ",\"" .
			     $analysis->db              . "\",\"" .
			     $analysis->db_version      . "\",\"" .
			     $analysis->db_file         . "\",\"" .
			     $analysis->program         . "\",\"" .
			     $analysis->program_version . "\",\"" .
			     $analysis->program_file    . "\",\"" .
			     $analysis->parameters      . "\",\"" .
			     $analysis->module          . "\",\"" .
			     $analysis->module_version  . "\",\"" .
			     $analysis->gff_source      . "\",\"" .
			     $analysis->gff_feature     . "\")";

    my $sth = $self->prepare($query);
    my $res = $sth->execute();
    
#    $sth = $self->prepare("select LAST_INSERT_ID()");
#    $res = $sth->execute();

#    my $id = $sth->fetchrow_hashref->{'LAST_INSERT_ID()'};
    return $analysis;
}

=head2 exists_Analysis

 Title   : exists_Analysis
 Usage   : $obj->exists_Analysis($anal)
 Function: Tests whether this feature already exists in the database
 Example :
 Returns : Analysis id if the entry exists
 Args    : Bio::EnsEMBL::Analysis::Analysis

=cut

sub exists_Analysis {
    my ($self,$anal) = @_;

    $self->throw("Object is not a Bio::EnsEMBL::Pipeline::Analysis") unless $anal->isa("Bio::EnsEMBL::Pipeline::Analysis");

    my $query = "select id,created from analysis where " .
	" program =    \""      . $anal->program         . "\" and " .
        " program_version = \"" . $anal->program_version . "\" and " .
	" program_file = \""    . $anal->program_file    . "\" and " .
        " parameters = \""      . $anal->parameters      . "\" and " .
        " module = \""          . $anal->module          . "\" and " .
        " module_version = \""    . $anal->module_version  . "\" and " .
	" gff_source = \""      . $anal->gff_source      . "\" and" .
	" gff_feature = \""     . $anal->gff_feature     . "\"";

    my $sth = $self->prepare($query);
    my $rv  = $sth->execute();

    if ($rv && $sth->rows > 0) {
	my $rowhash = $sth->fetchrow_hashref;
	my $analid  = $rowhash->{'id'}; 
	my $created = $rowhash->{'created'};

	return ($analid,$created);
    } else {
	return;
    }
}

=head2 get_Analysis {

  Title   : get_Analysis
  Usage   : my $analyis = $db->get_Analysis($id)
  Function: Retrieves an analysis object with
            id $id
  Returns : Bio::EnsEMBL::Pipeline::Analysis
  Args    : int

=cut


sub get_Analysis {
    my ($self,$id) = @_;

    $self->throw("No analysis id input") unless defined($id);

    my $query = "select created, " .
	        "program,program_version,program_file," .
	        "db,db_version,db_file," .
		"module,module_version," .
		"gff_source,gff_feature,".
		"parameters " . 
		"from analysis where " . 
                "id = $id";

    my $sth = $self->prepare($query);
    my $rv  = $sth->execute();

    if ($rv && $sth->rows > 0) {
	my $rowhash         = $sth->fetchrow_hashref;
	my $analid          = $rowhash->{'id'}; 
	my $created         = $rowhash->{'created'};
	my $program         = $rowhash->{'program'};
	my $program_version = $rowhash->{'program_version'};
	my $program_file    = $rowhash->{'program_file'};
	my $db              = $rowhash->{'db'};
	my $db_version      = $rowhash->{'db_version'};
	my $db_file         = $rowhash->{'db_file'};
	my $module          = $rowhash->{'module'};
	my $module_version  = $rowhash->{'module_version'};
	my $parameters      = $rowhash->{'parameters'};

	my $analysis = new Bio::EnsEMBL::Pipeline::Analysis(-id => $analid,
							    -created         => $created,
							    -program         => $program,
							    -program_version => $program_version,
							    -program_file    => $program_file,
							    -db              => $db,
							    -db_version      => $db_version,
							    -db_file         => $db_file,
							    -module          => $module,
							    -module_version  => $module_version,
							    -parameters      => $parameters,
							    );


	return ($analysis);
    } else {
	$self->throw("Couldn't find analysis object with id [$id]");
    }

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







