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

Bio::EnsEMBL::Pipeline::DB::ObjI

=head1 SYNOPSIS

=head1 DESCRIPTION

Interface for the connection to the analysis database

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::ObjI;

use vars qw(@ISA);
use strict;


use Bio::Root::RootI;
use Bio::EnsEMBL::DB::ObjI;

@ISA = qw(Bio::EnsEMBL::DB::ObjI Bio::Root::RootI);

=head1 ABSTRACT METHODS

These methods should all be implemented in an
object implementing C<Bio::EnsEMBL::Pipeline::ObjI>.

=head2 get_Job {

  Title   : get_Job
  Usage   : my $job = $db->get_Job($id)
  Function: Retrieves a job from the database
            when given its id
  Returns : Bio::EnsEMBL::Pipeline::DB::JobI
  Args    : int

=cut


sub get_Job {
    my ($self) = @_;

    $self->throw("Method get_Job not implemented");
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
    my ($self) = @_;

    $self->throw("Method get_JobsByInputId not implemented");
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
    my ($self) = @_;
    $self->throw("Method get_JobsByAge not implemented");

}

=head2 get_all_Analysis 

  Title   : get_all_Analysis
  Usage   : my @analyses = $db->get_all_Analysis()
  Function: Retrieves all analysis objects
  Returns : a list of Bio::EnsEMBL::Pipeline::Analysis
  Args    : none

=cut

sub get_all_Analysis {
    my ($self) = @_;
    $self->throw("Method get_all_Analysis not implemented");

}


=head2 get_AnalysisSummary 

  Title   : get_AnalysisSummary
  Usage   : my $analyis = $db->get_AnalysisSummary($id)
  Function: Retrieves summary of analyses objects matching analysis id
  Returns : a hash containing summary of analyses
  Args    : int

=cut

sub get_AnalysisSummary {
    my ($self) = @_;
    $self->throw("Method get_AnalysisSummary not implemented");

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
    $self->throw("Method get_all_Status not implemented");
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
    my ($self) = @_;

    $self->throw("Method write_Analysis not implemented");
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
    my ($self) = @_;

    $self->throw("Method get_Analysis not implemented");
}


