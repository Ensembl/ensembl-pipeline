#
# Object for storing details of a simple analysis job
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

Bio::EnsEMBL::Pipeline::SimpleJob

=head1 SYNOPSIS
    my $db       = new Bio::EnsEMBL::Pipeline::DBSQL::Obj(...);
    my $analysis = new Bio::EnsEMBL::Pipeline::Analysis  (...);

    my $job      = new Bio::EnsEMBL::Pipeline::DBSQL::Job(-dbobj    => $db,
							  -input_id => 1,
							  -analysis => $analysis,
							  -queue    => 'fast_blast_farm',
							  -create   => 1,
							  );

    #### This is the module you change for different analysis types.

    my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::ProcessList();

    my $runjob   = new Bio::EnsEMBL::Pipeline::SimpleJob(-jobobj   => $job,
							 -runnable => $runjob);

    my $stat   = $runjob->run;
    my $output = $runjob->output;



=head1 DESCRIPTION

Stores run and status details of a simple analysis job.  This module is a 
container bringing together the two parts of the pipeline - the database access 
part (job ids,status etc) and the doing part which actually runs the analysis and 
parses the results.

SimpleJob -> JobI       - stores job status and frozen objects in the database
          -> RunnableI  - runs programs and parses output.

In this way we can use our Runnable objects as standalone modules for 
one off analysis and can also change our job tracking code without interfering
with the analysis part.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::SimpleJob;

use vars qw(@ISA);
use strict;

use FreezeThaw qw(freeze thaw);

# Object preamble - inherits from Bio::Root::Object;

use Bio::EnsEMBL::Pipeline::DB::JobI;
use Bio::EnsEMBL::Pipeline::RunnableDBI;
#use Bio::EnsEMBL::Pipeline::RunnableDB::Est2Genome;
use Bio::EnsEMBL::Pipeline::RunnableDB::Clone_Vert_Est2Genome;
use Bio::EnsEMBL::Pipeline::RunnableDB::Clone_MiniGenewise;
use Bio::EnsEMBL::Pipeline::RunnableDB::FPC_BlastMiniGenewise;
use Bio::EnsEMBL::Pipeline::RunnableDB::Gene_Builder;
use Bio::EnsEMBL::Pipeline::RunnableDB::BigGene_Builder;
use Bio::EnsEMBL::Pipeline::RunnableDB::Genewise_Builder;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDBI Bio::EnsEMBL::Pipeline::DB::JobI Bio::Root::Object);


sub _initialize {
    my ($self,@args) = @_;

    my $make = $self->SUPER::_initialize;
    my ($jobobj,$runnable) = $self->_rearrange([qw(JOBOBJ
						   RUNNABLE
						   )],@args);

    $jobobj  ->isa("Bio::EnsEMBL::Pipeline::DB::JobI")  || $self->throw("No input job object");
    $runnable->isa("Bio::EnsEMBL::Pipeline::RunnableDBI") || $self->throw("No input RunnableDBI object");

    $self->jobobj  ($jobobj);
    $self->runnable($runnable);

    return $make; # success - we hope!
}


=head2 jobobj
  Title   : jobojb
  Usage   : $self->jobobj
  Function: Get/set method for the job object
            that reads and writes to the job tracking 
            database
  Returns : Bio::EnsEMBL::Pipeline::DB::JobI
  Args    : Bio::EnsEMBL::Pipeline::DB::JobI

=cut

sub jobobj {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	if ($arg->isa("Bio::EnsEMBL::Pipeline::DB::JobI")) {
	    $self->{_jobobj} = $arg;
	} else {
	    $self->throw("Input [$arg] is not a Bio::EnsEMBL::Pipeline::DB::JobI");
	}
    }
    return $self->{_jobobj};
}

=head2 runnable
  Title   : runnable
  Usage   : $self->runnable
  Function: Get/set method for the runnable object
            which contains code to run the analysis
  Returns : Bio::EnsEMBL::Pipeline::RunnableI
  Args    : Bio::EnsEMBL::Pipeline::RunnableI

=cut


sub runnable {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	if ($arg->isa("Bio::EnsEMBL::Pipeline::RunnableI")) {
	    $self->{_runnable} = $arg;
	} else {
	    $self->throw("[$arg] is not a Bio::EnsEMBL::Pipeline::RunnableI");
	}
    }
    return $self->{_runnable};
}

### These are the RunnableI methods

=head2 run

  Title   : run
  Usage   : $self->run
  Function: Runs the analysis on the object
  Returns : nothing
  Args    : none

=cut



sub run {
    my ($self) = @_;
    print(STDERR "Running in simplejob\n");
    $self->runnable->run;
    print ("done\n");
}

=head2 output

  Title   : output
  Usage   : my $output = $self->output;
  Function: Returns the parsed output from the job
  Returns : ref
  Args    : none

=cut

sub output {
    my ($self) = @_;

    return $self->runnable->output;
}

=head2 dbobj

  Title   : dbobj
  Usage   : my $db = $self->dbobj
  Function: The database handle used to fetch input
  Returns : ref
  Args    : none

=cut

sub dbobj {
    my ($self) = @_;

    return $self->runnable->dbobj;
}

=head2 fetch_input

  Title   : fetch_input
  Usage   : $self->fetch_input
  Function: Fetches any input data needed for the job from the database
  Returns : ref
  Args    : none

=cut

sub fetch_input {
    my ($self) = @_;

    return $self->runnable->fetch_input
}

### These are all the JobI methods

=head2 id

  Title   : id
  Usage   : $self->id($id)
  Function: Get/set method for the id of the job itself
            This will usually be generated by the
            back end database the jobs are stored in
  Returns : int
  Args    : int

=cut


sub id {
    my ($self,$arg) = @_;

    return $self->jobobj->id($arg);

}


=head2 get_id

  Title   : get_id
  Usage   : my $newid = $self->get_id
  Function: Creates a new job entry in the database
            and returns the new id.
  Returns : int
  Args    : 

=cut


sub get_id {
    my ($self) = @_;

    return $self->jobobj->get_id();
}

=head2 input_id

  Title   : input_id
  Usage   : $self->input_id($id)
  Function: Get/set method for the id of the input to the job
  Returns : int
  Args    : int

=cut


sub input_id {
    my ($self,$arg) = @_;

    return $self->jobobj->input_id($arg);
}

=head2 analysis

  Title   : analysis
  Usage   : $self->analysis($anal);
  Function: Get/set method for the analysis object of the job
  Returns : Bio::EnsEMBL::Pipeline::Analysis
  Args    : bio::EnsEMBL::Pipeline::Analysis



=cut


sub analysis {
    my ($self,$arg) = @_;
    return $self->jobobj->analysis($arg);
}

=head2 LSF_id

  Title   : LSF_id
  Usage   : $self->LSF_id($id)
  Function: Get/set method for the LSF id of the job
  Returns : int
  Args    : int

=cut


sub LSF_id {
    my ($self,$arg) = @_;

    return $self->jobobj->LSF_id($arg);

}

=head2 queue

  Title   : queue
  Usage   : $self->queue
  Function: Get/set method for the LSF queue name
  Returns : String
  Args    : String

=cut

sub queue {
    my ($self,$arg) = @_;

    return $self->jobobj->queue($arg);

}


=head2 machine

  Title   : machine
  Usage   : $self->machine($machine)
  Function: Get/set method for the machine the job is running on
  Returns : string
  Args    : string

=cut

sub machine {
    my ($self,$arg) = @_;

    return $self->jobobj->machine();
}

sub useDB {
    my ($self,$arg) = @_;

    return $self->jobobj->useDB($arg);

}

sub runlocally {
    my ($self,$arg) = @_;

    return $self->jobobj->runlocally($arg);

}

=head2 submit

  Title   : submit
  Usage   : $self->submit
  Function: Submits the job to the specified LSF queue
  Returns : 
  Args    : 

=cut

sub submit {
    my ($self) = @_;

    $self->jobobj->submit($self);

}

=head2 store

  Title   : store
  Usage   : $self->store
  Function: Stores the object as a string in the database
  Returns : Bio::EnsEMBL::Pipeline::Status
  Args    : none

=cut

sub store {
    my ($self) = @_;

    return $self->jobobj->store($self);

}


=head2 submission_checks

  Title   : submission_checks
  Usage   : $self->submission_checks
  Function: After submission to the LSF queue when 
            the wrapper script is run - these are
            the checks to run (on binaries,databases etc)
            before the job is run.
  Returns : String
  Args    : None

=cut

sub submission_checks {
    my ($self) = @_;

    $self->throw("Method submission_checks not implemented");

}

=head2 set_status

  Title   : set_status
  Usage   : my $status = $job->set_status
  Function: Sets the job status
  Returns : nothing
  Args    : Bio::EnsEMBL::Pipeline::Status

=cut

sub set_status {
    my ($self,$arg) = @_;

    return $self->jobobj->set_status($arg);
}


=head2 current_status

  Title   : current_status
  Usage   : my $status = $job->current_status
  Function: Get/set method for the current status
  Returns : Bio::EnsEMBL::Pipeline::Status
  Args    : Bio::EnsEMBL::Pipeline::Status

=cut

sub current_status {
    my ($self,$arg) = @_;

    return $self->jobobj($arg);
}

=head2 get_all_status

  Title   : get_all_status
  Usage   : my @status = $job->get_all_status
  Function: Get all status objects associated with this job
  Returns : @Bio::EnsEMBL::Pipeline::Status
  Args    : @Bio::EnsEMBL::Pipeline::Status

=cut

sub get_all_status {
    my ($self) = @_;

    $self->throw("Method get_all_status not implemented");
}


=head2 _dbobj

  Title   : _dbobj
  Usage   : my $db = $self->_dbobj
  Function: Get/set method for the database handle
  Returns : @Bio::EnsEMBL::Pipeline::DBSQL::Obj
  Args    : @Bio::EnsEMBL::Pipeline::DBSQL::Obj,none

=cut

sub _dbobj {
    my ($self,$arg) = @_;
 
    if (defined($arg)) {
     $self->runnable->dbobj($arg);
    } 
    return $self->jobobj->_dbobj($arg);
}




=head2 stdout_file

  Title   : stdout_file
  Usage   : my $file = $self->stdout_file
  Function: Get/set method for stdout.
  Returns : string
  Args    : string

=cut

sub stdout_file {
    my ($self,$arg) = @_;

    return $self->jobobj->stdout_file($arg);
}

=head2 stderr_file

  Title   : stderr_file
  Usage   : my $file = $self->stderr_file
  Function: Get/set method for stderr.
  Returns : string
  Args    : string

=cut

sub stderr_file {
    my ($self,$arg) = @_;

    return $self->jobobj->stderr_file($arg);
}

=head2 output_file

  Title   : output_file
  Usage   : my $file = $self->output_file
  Function: Get/set method for output
  Returns : string
  Args    : string

=cut

sub output_file {
    my ($self,$arg) = @_;

    return $self->jobobj->output_file($arg);

}


=head2 input_object_file

  Title   : intput_object_file
  Usage   : my $file = $self->input_object_file
  Function: Get/set method for the input object file
  Returns : string
  Args    : string

=cut

sub input_object_file {
    my ($self,$arg) = @_;

    return $self->jobobj->input_object_file($arg);
}

sub status_file {
    my ($self,$arg) = @_;

    return $self->jobobj->status_file($arg);
}
sub disconnect {
   my ($self) = @_;

   $self->jobobj->disconnect;

}
1;












