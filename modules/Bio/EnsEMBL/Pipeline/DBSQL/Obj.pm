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

use Bio::EnsEMBL::Pipeline::Analysis;
use Bio::EnsEMBL::Pipeline::ExonPair;
use Bio::EnsEMBL::Pipeline::DB::ObjI;
use Bio::EnsEMBL::Pipeline::DBSQL::Job;
use Bio::Root::Object;

use FreezeThaw qw(freeze thaw);

# Inherits from the base bioperl object

@ISA = qw(Bio::EnsEMBL::Pipeline::DB::ObjI  Bio::EnsEMBL::DBSQL::Obj Bio::Root::Object);

sub _initialize {
    my ($self,@args) = @_;

    my $make = $self->SUPER::_initialize(@args);

    return $make; # success - we hope!

}

=head2 get_TimClone

 Title   : get_Clone
 Usage   : $db->get_Clone($disk_id)
 Function: Gets a Clone object with disk_id as its disk_id
 Example : $db->get_Clone($disk_id)
 Returns : Bio::EnsEMBL::Pipeline::DBSQL::Clone object 
 Args    : disk_id


=cut

sub get_TimClone{
   my ($self,$disk_id) = @_;

   $disk_id || $self->throw("Trying to delete a clone without a disk id!");

   my $sth = $self->prepare("select disk_id from clone where disk_id = \"$disk_id\";");
   my $res = $sth ->execute();
   my $rv  = $sth ->rows;

   if( ! $rv ) {
       # make sure we deallocate sth - keeps DBI happy!
       $sth = undef;
       $self->throw("Clone $disk_id does not seem to occur in the database!");
   }

   my $clone = new Bio::EnsEMBL::Pipeline::DBSQL::Clone( -disk_id    => $disk_id,
							 -dbobj      => $self );

   return $clone;
}

=head2 create_Clone

 Title   : create_Clone
 Usage   : $db->create_Clone($disk_id,$clone_group,$chromosome)
 Function: writes a new clone in the database
 Example : $db->create_Clone('dummy','SU','22')
 Returns : nothing
 Args    : disk_id,clone group,chromosome


=cut

sub create_Clone {
   my ($self,$disk_id,$clone_group,$chromosome) = @_;

   $disk_id || $self->throw("Trying to create a clone without a disk id!\n");
   $clone_group || $self->throw("Trying to create a clone without a clone group!");
   $chromosome || $self->throw("Trying to create a clone without a chromosome id!");
   
   my @sql;

   push(@sql,"lock tables clone write");
   push(@sql,"insert into clone(disk_id,clone_group,chromosome,last_check,created,modified) values('$disk_id','$clone_group','$chromosome',now(),now(),now())");
   push(@sql,"unlock tables");   

   foreach my $sql (@sql) {
     my $sth =  $self->prepare($sql);
     #print STDERR "Executing $sql\n";
     my $rv  =  $sth->execute();
     $self->throw("Failed to insert clone $disk_id") unless $rv;
   }
}

=head2 get_Job 

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
    my $query = "select id,input_id,analysis,LSF_id,machine,object,queue " . 
	",stdout_file,stderr_file,input_object_file,output_file,status_file " .
        "from job " . 
	"where input_id = \"$id\"";
    
    my $sth = $self->prepare($query);
    my $res = $sth->execute();
    
    my @jobs;
    while (my $row = $sth->fetchrow_hashref) {
	my $job = $self->_parseJob($row);
	push(@jobs,$job);
    }
    return @jobs;
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
  Args    : int analysis id, string status

=cut

sub get_JobsbyStatus_and_Analysis {
    my ($self,$status, $analysis) = @_;

    $self->throw("Require status and analysis id for get_JobsbyStatus_and_Analysis") 
                            unless ($analysis && $status);

    my $query = "select j.id, j.input_id, j.analysis, j.LSF_id," .
	                   "j.machine,queue, j.stdout_file, j.stderr_file,".
                       "j.input_object_file, j.output_file,".
                       "j.status_file " . 
	            "from job as j, current_status as cs ". 
                "where j.id = cs.id and j.analysis = $analysis ".
                       "and cs.status = \'$status\'";
                 
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
	                     "insert into analysisprocess(id,created,db,db_version,db_file," .
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
    
    $sth = $self->prepare("select LAST_INSERT_ID()");
    $res = $sth->execute();

    my $id = $sth->fetchrow_hashref->{'LAST_INSERT_ID()'};

    $analysis->id($id);

    # Should fetch the created stamp as well now.

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

    my $query = "select id,created from analysisprocess where " .
	" program =    \""      . $anal->program         . "\" and " .
        " program_version = \"" . $anal->program_version . "\" and " .
	" program_file = \""    . $anal->program_file    . "\" and " .
        " parameters = \""      . $anal->parameters      . "\" and " .
        " module = \""          . $anal->module          . "\" and " .
        " module_version = \""  . $anal->module_version  . "\" and " .
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
		        "from analysisprocess where " . 
                "id = $id";

    my $sth = $self->prepare($query);
    my $rv  = $sth->execute();

    if ($rv && $sth->rows > 0) {
	my $rowhash         = $sth->fetchrow_hashref;
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

	my $analysis = new Bio::EnsEMBL::Pipeline::Analysis(-id => $id,
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

=head2 get_all_Analysis 

  Title   : get_all_Analysis
  Usage   : my @analyses = $db->get_all_Analysis()
  Function: Retrieves all analysis objects
  Returns : a list of Bio::EnsEMBL::Pipeline::Analysis
  Args    : none

=cut

sub get_all_Analysis {
    my ($self) =@_;
    
    my $query = "select id, " .
	            "program,program_version," .
	            "db,db_version," .
		        "gff_source,gff_feature, ".
                "module, module_version,parameters ".
		        "from analysisprocess ";
    
    my $sth = $self->prepare($query);
    my $rv  = $sth->execute();
    my @analyses;
    
    if ($rv && $sth->rows > 0)
    {
        while (my $row = $sth->fetchrow_hashref)
        {
            my $id              = $row->{'id'};
            my $program         = $row->{'program'};
            my $program_version = $row->{'program_version'};
            my $db              = $row->{'db'};
            my $db_version      = $row->{'db_version'};
            my $gff_source      = $row->{'gff_source'};
            my $gff_feature     = $row->{'gff_feature'};
            my $module          = $row->{'module'};
	        my $module_version  = $row->{'module_version'};
	        my $parameters      = $row->{'parameters'};
            
                        
            my $analysis = new Bio::EnsEMBL::Pipeline::Analysis(
                                -id              => $id,
							    -program         => $program,
							    -program_version => $program_version,
							    -db              => $db,
							    -db_version      => $db_version,
							    -gff_source      => $gff_source,
                                -gff_feature     => $gff_feature,
                                -module          => $module,
							    -module_version  => $module_version,
							    -parameters      => $parameters,
                                );
            push (@analyses, $analysis);
        }
    }
    else
    {
        $self->throw("No data in analysis table\n");
    }
    return (@analyses);
}

=head2 get_AnalysisSummary 

  Title   : get_AnalysisSummary
  Usage   : my $analyis = $db->get_AnalysisSummary($id)
  Function: Retrieves summary of analyses objects matching analysis id
  Returns : a hash containing summary of analyses
  Args    : int

=cut

sub get_AnalysisSummary {
    my ($self, $id) = @_;
    
    my $query = "select count(*), cs.status " .
	            "from job, current_status as cs " .
	            "where job.id = cs.id " .
		        "and job.analysis = $id ".
                "group by cs.status";
    
    my $sth = $self->prepare($query);
    my $rv  = $sth->execute();
    my %summary;
    
    while (my $row = $sth->fetchrow_hashref)
    {
            $summary{$row->{'status'}}  = $row->{'count(*)'};        
    }
    return %summary;
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


=head2 get_all_ExonPairs

 Title   : get_all_ExonPairs(@contigs)
 Usage   : my @pairs = get_all_ExonPairs(@contigs);
 Function: Returns all the exon pairs for exons that lie on 
           the input contigs
 Example : 
 Returns : @Bio::EnsEMBL::Analysis::ExonPair
 Args    : @Bio::EnsEMBL::DB::ContigI

=cut


sub get_all_ExonPairs {
    my ($self,@contigids) = @_;

    my @pairs;
    my $contigstr = "";

    foreach my $contig (@contigids) {
#	$self->throw("Not a Bio::EnsEMBL::DB::ContigI") unless $contig->isa("Bio::EnsEMBL::DB::ContigI");
	$contigstr .= "'" . $contig . "',";
    }

    $contigstr =~ s/(.*)\,$/$1/;

    my ($exon1_id,$exon2_id,$exon1_version,$exon2_version,$created);

    # This could be changed so it gets all the exons as well
    
    my $query = "select distinct ep.exon1_id,ep.exon2_id,ep.created,ep.type,exon1_version,exon2_version " . 
	"from exon_pair as ep,exon as e where " .
        "e.contig in (" . $contigstr . ") and "             .
	"((ep.exon1_id = e.id and ep.exon1_version = e.version) or " .
	"( ep.exon2_id = e.id and ep.exon2_version = e.version))";

    my $sth = $self->prepare($query);
    my $res = $sth->execute();

    while (my $row = $sth->fetchrow_arrayref) {
	my $ex1id    = $row->[0];
	my $ex2id    = $row->[1];
	my $created  = $row->[2];
	my $type     = $row->[3];
	my $version1 = $row->[4];
	my $version2 = $row->[5];
	my $contigid = $row->[6];

	my $exon1 = $self->get_Exon($ex1id);
	my $exon2 = $self->get_Exon($ex2id);

	my $pair = new Bio::EnsEMBL::Pipeline::ExonPair(-exon1 => $exon1,
							-exon2 => $exon2,
							-type  => $type);
	
	push(@pairs,$pair);
    }
 
    return (@pairs);

}


=head2 write_ExonPairs

 Title   : write_ExonPairs(@pairs)
 Usage   : my @pairs = write_ExonPairs(@pairs);
 Function: Returns all the exon pairs for exons that lie on 
           the input contigs
 Example : 
 Returns : @Bio::EnsEMBL::Pipeline::ExonPair
 Args    : @Bio::EnsEMBL::DB::ContigI

=cut

sub write_ExonPairs {
    my ($self,@pairs) = @_;

    my $sth = $self->prepare("insert into exon_pair(exon1_id,exon2_id,type,created,exon1_version,exon2_version) " .
			     "values (?,?,?,?,?,?)");

    foreach my $pair (@pairs) {
	$self->throw("Pair $pair is not an ExonPair") unless $pair->isa("Bio::EnsEMBL::Pipeline::ExonPair");


	my $sth2 = $self->prepare("select * from exon_pair where " .
				  "exon1_id = '"      . $pair->exon1->id      . "' and " .
				  "exon2_id = '"      . $pair->exon2->id      . "' and " .
				  "exon1_version = '" . $pair->exon1->version . "' and " .
				  "exon2_version = '" . $pair->exon2->version . "'");

	my $res = $sth2->execute;

	if ($sth2->rows == 0) {
	    $sth->execute($pair->exon1->id,
			  $pair->exon2->id,
			  $pair->type,
			  "now()",
			  $pair->exon1->version,
			  $pair->exon2->version);
	} else {
	    $self->warn("Exon pair [" . $pair->exon1->id . "][" . $pair->exon2->id . 
			" already exists. Not writing to database");
	}
    }
}

=head2 delete_ExonPairs

 Title   : delete_ExonPairs(@pairs)
 Usage   : $obj->delete_ExonPairs(@pairs)
 Function: Deletes all input exon pairs from the database
 Example : 
 Returns : Nothing
 Args    : @Bio::EnsEMBL::Pipeline::ExonPair

=cut

sub delete_ExonPairs {
    my ($self,@pairs) = @_;

    foreach my $pair (@pairs) {
	$self->throw("Not an Bio::EnsEMBL::Pipeline::ExonPair") unless $pair->isa("Bio::EnsEMBL::Pipeline::ExonPair");

	my $query = "delete from exon_pair where exon1_id = \""        . $pair->exon1->id      . "\" and " .
	                                        "exon2_id = \""        . $pair->exon2->id      . "\" and " .
						"exon1_version = \""   . $pair->exon1->version . "\" and " .
     				                "exon2_version = \""   . $pair->exon2->version . "\"";

	my $sth = $self->prepare($query);
	my $res = $sth->execute;

    }
}

=head2 write_Exon

 Title   : write_Exon
 Usage   : $obj->write_Exon($exon)
 Function: writes a particular exon into the database
 Example :
 Returns : 
 Args    :

=cut

sub write_Exon{
   my ($self,$exon) = @_;
   my $old_exon;
   
   if( ! $exon->isa('Bio::EnsEMBL::Exon') ) {
       $self->throw("$exon is not a EnsEMBL exon - not dumping!");
   }

   eval {
       $old_exon = $self->get_Exon($exon->id);

       $self->warn("Exon [" . $exon->id . "] already exists in database. Not writing");

       
   };
   
   if  ( $@ || ($exon->version > $old_exon->version)) {
       my $exonst = "insert into exon (id,version,contig,created,modified," .
	   "seq_start,seq_end,strand,phase,stored,end_phase) values ('" .
	       $exon->id()         . "'," .
               $exon->version()    . ",'".
	       $exon->contig_id()  . "', FROM_UNIXTIME(" .
	       $exon->created()    . "), FROM_UNIXTIME(" .
	       $exon->modified()   . ")," .
	       $exon->start        . ",".
	       $exon->end          . ",".
	       $exon->strand       . ",".
	       $exon->phase        . ",now(),".
	       $exon->end_phase    . ")";

       my $sth = $self->prepare($exonst);
       $sth->execute();
   }
}


=head2 delete_Exon

 Title   : delete_Exon
 Usage   : $obj->delete_Exon($exon_id)
 Function: Deletes exon, including exon_transcript rows
 Example : $obj->delete_Exon(ENSE000034)
 Returns : nothing
 Args    : $exon_id


=cut

sub delete_Exon{
    my ($self,$exon_id) = @_;

    $exon_id || $self->throw ("Trying to delete an exon without an exon_id\n");
    
    #Delete exon_transcript rows
#    my $sth = $self->prepare("delete from exon_transcript where transcript = '".$exon_id."'");
#    my $res = $sth ->execute;

    #Delete exon rows
    my $sth = $self->prepare("delete from exon where id = '".$exon_id."'");
    my $res = $sth->execute;

 #   $self->delete_Supporting_Evidence($exon_id);

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


=head2 get_Exon

 Title   : get_Exon
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut

sub get_Exon{
   my ($self,$exonid) = @_;

   my $exon = Bio::EnsEMBL::Exon->new();

   my $sth     = $self->prepare("select id,version,contig,UNIX_TIMESTAMP(created)," . 
				"UNIX_TIMESTAMP(modified),seq_start,seq_end,strand," . 
				"phase from exon where id = '$exonid'");
   my $res     = $sth->execute;
   my $rowhash = $sth->fetchrow_hashref;

   if( ! defined $rowhash ) {
       $self->throw("No exon of this id $exonid");
   }

   $exon->contig_id($rowhash->{'contig'});
   $exon->version  ($rowhash->{'version'});

   my $contig_id = $exon->contig_id();

   # we have to make another trip to the database to get out the contig to clone mapping.
#   my $sth2     = $self->prepare("select clone from contig where id = '$contig_id'");
#   my $res2     = $sth2->execute;
#   my $rowhash2 = $sth2->fetchrow_hashref;

#   $exon->clone_id($rowhash2->{'clone'});

   # rest of the attributes
   $exon->id      ($rowhash->{'id'});
   $exon->created ($rowhash->{'UNIX_TIMESTAMP(created)'});
   $exon->modified($rowhash->{'UNIX_TIMESTAMP(modified)'});
   $exon->start   ($rowhash->{'seq_start'});
   $exon->end     ($rowhash->{'seq_end'});
   $exon->strand  ($rowhash->{'strand'});
   $exon->phase   ($rowhash->{'phase'});
   
   # we need to attach this to a sequence. For the moment, do it the stupid
   # way perhaps?

#   my $seq;

#   if( $self->_contig_seq_cache($exon->contig_id) ) {
#       $seq = $self->_contig_seq_cache($exon->contig_id);
#   } else {
#       my $contig = $self->get_Contig($exon->contig_id());
#       $seq = $contig->seq();
#       $self->_contig_seq_cache($exon->contig_id,$seq);
#   }

#   $exon->attach_seq($seq);

   return $exon;
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




