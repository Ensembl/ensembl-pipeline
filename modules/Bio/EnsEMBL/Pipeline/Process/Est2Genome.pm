#
# Object for storing details of an analysis job
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

Bio::EnsEMBL::Pipeline::Process::Est2Genome

=head1 SYNOPSIS

=head1 DESCRIPTION

Implementation of AnalysisProcess to create Est2Genome jobs

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::Process:Est2Genome;

use vars qw(@ISA);
use strict;


# Object preamble - inherits from Bio::Root::Object;


@ISA = qw(Bio::EnsEMBL::Pipeline::DBSQL::AnalysisProcess Bio::Root::Object);

sub _initialize {
    my ($self,@args) = @_;

    my $make = $self->SUPER::_initialize;
    my ($dbobj,$input_ids)  $self->_rearrange([qw(DBOBJ
						  INPUT_IDS
						  )],@args);

    $input_ids  || $self->throw("Can't create an AnalysisProcess object without input_ids");
    $dbobj      || $self->throw("Can't create an AnalysisProcess object without a database handle");

    $dbobj->isa("Bio::EnsEMBL::Pipeline::DBSQL::Obj") || 
	$self->throw("Database object [$dbobj] is not a Bio::EnsEMBL::Pipeline::DBSQL::Obj");

    $self->input_ids  ($input_ids);
    $self->_dbobj     ($dbobj);

    return $make; # success - we hope!
}

=head2 check_dependencies

  Title   : check_dependencies
  Usage   : $self->check_dependencies
  Function: Checks in the database whether
            any necessary previous jobs have completed
            before this process can start
  Returns : 1,undef
  Args    : 1,undef

=cut

sub check_dependencies {
    my ($self) = @_;

    $self->throw("check_dependencies hasn't been implemented yet");
}


=head2 make_jobs

  Title   : make_jobs
  Usage   : my @ids = $self->make_jobs
  Function: Creates all necessary jobs and 
            writes them into the database
            Returns an array of the job ids
  Returns : @int
  Args    : none

=cut


sub make_jobs {
    my ($self) = @_;

    $self->throw("make_jobs hasn't been implemented yet");
}

;









