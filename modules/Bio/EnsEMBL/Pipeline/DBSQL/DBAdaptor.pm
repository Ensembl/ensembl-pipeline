#
# Object for storing the connection to the analysis database
#
# Written by Simon Potter <scp@sanger.ac.uk>
# Based on Michele Clamp's Bio::EnsEMBL::Pipeline::DBSQL::Obj
#
# Copyright GRL/EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code


=pod

=head1 NAME

Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor -
adapter class for EnsEMBL Pipeline DB

=head1 SYNOPSIS

    my $dbobj = new Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
    $dbobj->do_funky_db_stuff;

=head1 DESCRIPTION

Interface for the connection to the analysis database

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;


use vars qw(@ISA);
use strict;
use DBI;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::DB::ObjI;
use Bio::EnsEMBL::Pipeline::DBSQL::RuleAdaptor;
use Bio::EnsEMBL::DBSQL::AnalysisAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::JobAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::StateInfoContainer;

use Bio::Root::RootI;

# Inherits from the base bioperl object

@ISA = qw(Bio::EnsEMBL::DBSQL::DBAdaptor);


# sub new {
    # my ($class,@args) = @_;
    # my $self = $class->SUPER::new(@args);
    # return $self;
# }

# new() inherited from Bio::EnsEMBL::DBSQL::BaseAdaptor


=head2 get_JobAdaptor

 Title   : get_JobAdaptor
 Usage   : $db->get_JobAdaptor
 Function: The Adaptor for Job objects in this db
 Example :
 Returns : Bio::EnsEMBL::Pipeline::DBSQL::JobAdaptor
 Args    : nothing


=cut

sub get_JobAdaptor {
  my ($self) = @_;

  if( ! defined $self->{_JobAdaptor} ) {
    $self->{_JobAdaptor} = Bio::EnsEMBL::Pipeline::DBSQL::JobAdaptor->new
      ( $self );
  }

  return $self->{_JobAdaptor};
}


# DECRUFTME - do we need this?

sub get_OldAnalysis {
        my ($self,$id) = @_;

        return $self->SUPER::get_Analysis($id);
}


=head2 get_AnalysisAdaptor

 Title   : get_AnalysisAdaptor
 Usage   : $db->get_AnalysisAdaptor
 Function: The Adaptor for Analysis objects in this db
 Example :
 Returns : Bio::EnsEMBL::Pipeline::DBSQL::AnalysisAdaptor
 Args    : nothing

=cut

sub get_AnalysisAdaptor {
  my ($self) = @_;

  if( ! defined $self->{_AnalysisAdaptor} ) {
    require Bio::EnsEMBL::DBSQL::AnalysisAdaptor;
    $self->{_AnalysisAdaptor} = Bio::EnsEMBL::DBSQL::AnalysisAdaptor->new
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
 Args    : nothing

=cut

sub get_RuleAdaptor {
  my ($self) = @_;

  if( ! defined $self->{_RuleAdaptor} ) {
    $self->{_RuleAdaptor} = Bio::EnsEMBL::Pipeline::DBSQL::RuleAdaptor->new
      ( $self );
  }

  return $self->{_RuleAdaptor};
}


=head2 get_StateInfoContainer

 Title   : get_StateInfoContainer
 Usage   : $db->get_StateInfoContainer
 Function:
 Example :
 Returns : Bio::EnsEMBL::Pipeline::DBSQL::StateInfoContainer
 Args    : nothing

=cut

sub get_StateInfoContainer {
  my ($self) = @_;

  if( ! defined $self->{_StateInfoContainer} ) {
    $self->{_StateInfoContainer} = Bio::EnsEMBL::Pipeline::DBSQL::StateInfoContainer->new
      ( $self );
  }

  return $self->{_StateInfoContainer};
}


# CHECKME - this is probably old hat as well...

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

    $self->warn(q/You really should use "$job->remove" :)/);

    $self->get_JobAdaptor->fetch_by_dbID($id)->remove
     or $self->warn("Can't recreate job with ID $id");
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

sub DESTROY {
   my ($obj) = @_;

   $obj->_unlock_tables();

   if( $obj->{'_db_handle'} ) {
       $obj->{'_db_handle'}->disconnect;
       $obj->{'_db_handle'} = undef;
   }
}
