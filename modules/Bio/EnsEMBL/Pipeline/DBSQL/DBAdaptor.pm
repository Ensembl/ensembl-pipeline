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

use Bio::EnsEMBL::Root;

# Inherits from the base bioperl object

@ISA = qw(Bio::EnsEMBL::DBSQL::DBAdaptor);


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
    require Bio::EnsEMBL::Pipeline::DBSQL::JobAdaptor;
    $self->{_JobAdaptor} = Bio::EnsEMBL::Pipeline::DBSQL::JobAdaptor->new
      ( $self );
  }

  return $self->{_JobAdaptor};
}


sub get_PmatchFeatureAdaptor{
  my ($self) = @_;
  print STDERR "getting a pmatch feature adaptor\n";
  if( ! defined $self->{_PmatchFeatureAdaptor} ) {
    require Bio::EnsEMBL::Pipeline::DBSQL::PmatchFeatureAdaptor;
    $self->{_PmatchFeatureAdaptor} = Bio::EnsEMBL::Pipeline::DBSQL::PmatchFeatureAdaptor->new
      ( $self );
  }

  return $self->{_PmatchFeatureAdaptor};
  
}

=head2 get_AnalysisAdaptor

 Title   : get_AnalysisAdaptor
 Usage   : $db->get_AnalysisAdaptor
 Function: The Adaptor for Analysis objects in this db
 Example :
 Returns : Bio::EnsEMBL::DBSQL::AnalysisAdaptor
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
    require Bio::EnsEMBL::Pipeline::DBSQL::RuleAdaptor;
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
    require Bio::EnsEMBL::Pipeline::DBSQL::StateInfoContainer;
    $self->{_StateInfoContainer} = Bio::EnsEMBL::Pipeline::DBSQL::StateInfoContainer->new
      ( $self );
  }

  return $self->{_StateInfoContainer};
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


sub pipeline_lock {
    my ($self, $string) = @_;

    my $sth;

    if ($string) {
	$sth = $self->prepare(qq{
	    INSERT into meta (meta_key, meta_value)
	    VALUES ('pipeline.lock', ?)
	});
	$sth->execute($string);
    }
    else {
	$sth = $self->prepare(qq{
	    SELECT meta_value
	    FROM   meta
	    WHERE  meta_key = 'pipeline.lock'
	});

        $sth->execute;
        my $row = $sth->fetchrow_arrayref;

        if ($row) {
	    return $row->[0];
        }
        else {
	    return undef;
        }
    }
    $sth->finish;
}


sub pipeline_unlock {
    my ($self) = @_;

    my $sth;

    $sth = $self->prepare(qq{
	DELETE
	FROM   meta
	WHERE  meta_key = 'pipeline.lock'
    });

    $sth->execute;
    $sth->finish;
}

1;
