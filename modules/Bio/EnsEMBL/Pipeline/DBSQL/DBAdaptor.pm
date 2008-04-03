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

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;


use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Root;

# Inherits from the base bioperl object

@ISA = qw(Bio::EnsEMBL::DBSQL::DBAdaptor);


# new() inherited from Bio::EnsEMBL::DBSQL::BaseAdaptor


sub get_available_adaptors {
  my ($self) = @_;

  my $pairs = $self->SUPER::get_available_adaptors();

  $pairs->{'Analysis'}                  = 'Bio::EnsEMBL::Pipeline::DBSQL::AnalysisAdaptor';
  $pairs->{'Job'}                       = 'Bio::EnsEMBL::Pipeline::DBSQL::JobAdaptor';
  $pairs->{'PmatchFeature'}             = 'Bio::EnsEMBL::Pipeline::DBSQL::PmatchFeatureAdaptor';
  $pairs->{'Rule'}                      = 'Bio::EnsEMBL::Pipeline::DBSQL::RuleAdaptor';
  $pairs->{'StateInfoContainer'}        = 'Bio::EnsEMBL::Pipeline::DBSQL::StateInfoContainer';
  $pairs->{'CompressedDnaAlignFeature'} = 'Bio::EnsEMBL::Pipeline::DBSQL::CompressedDnaAlignFeatureAdaptor';

  return $pairs; 
}

=head2 get_JobAdaptor

 Title   : get_JobAdaptor
 Usage   : $db->get_JobAdaptor
 Function: The Adaptor for Job objects in this db
 Example :
 Returns : Bio::EnsEMBL::Pipeline::DBSQL::JobAdaptor
 Args    : nothing


=cut

#sub get_JobAdaptor {
#  my ($self) = @_;
#
#  if( ! defined $self->{_JobAdaptor} ) {
#    require Bio::EnsEMBL::Pipeline::DBSQL::JobAdaptor;
#    $self->{_JobAdaptor} = Bio::EnsEMBL::Pipeline::DBSQL::JobAdaptor->new
#      ( $self );
#  }
#
#  return $self->{_JobAdaptor};
#}


#sub get_PmatchFeatureAdaptor{
#  my ($self) = @_;
#  #print STDERR "getting a pmatch feature adaptor\n";
#  if( ! defined $self->{_PmatchFeatureAdaptor} ) {
#    require Bio::EnsEMBL::Pipeline::DBSQL::PmatchFeatureAdaptor;
#    $self->{_PmatchFeatureAdaptor} = Bio::EnsEMBL::Pipeline::DBSQL::PmatchFeatureAdaptor->new
#      ( $self );
#  }
#
#  return $self->{_PmatchFeatureAdaptor};
#  
#}


=head2 get_RuleAdaptor

 Title   : get_RuleAdaptor
 Usage   : $db->get_RuleAdaptor
 Function: The Adaptor for Rule objects in this db
 Example :
 Returns : Bio::EnsEMBL::Pipeline::DBSQL::RuleAdaptor
 Args    : nothing

=cut

#sub get_RuleAdaptor {
#  my ($self) = @_;
#
#  if( ! defined $self->{_RuleAdaptor} ) {
#    require Bio::EnsEMBL::Pipeline::DBSQL::RuleAdaptor;
#    $self->{_RuleAdaptor} = Bio::EnsEMBL::Pipeline::DBSQL::RuleAdaptor->new
#      ( $self );
#  }
#
#  return $self->{_RuleAdaptor};
#}


=head2 get_StateInfoContainer

 Title   : get_StateInfoContainer
 Usage   : $db->get_StateInfoContainer
 Function:
 Example :
 Returns : Bio::EnsEMBL::Pipeline::DBSQL::StateInfoContainer
 Args    : nothing

=cut

#sub get_StateInfoContainer {
#  my ($self) = @_;
#
#  if( ! defined $self->{_StateInfoContainer} ) {
#    require Bio::EnsEMBL::Pipeline::DBSQL::StateInfoContainer;
#    $self->{_StateInfoContainer} = Bio::EnsEMBL::Pipeline::DBSQL::StateInfoContainer->new
#      ( $self );
#  }
#
#  return $self->{_StateInfoContainer};
#}

#sub get_AnalysisAdaptor {
#  my ($self) = @_;
#
#  if( ! defined $self->{_AnalysisAdaptor} ) {
#    require Bio::EnsEMBL::Pipeline::DBSQL::AnalysisAdaptor;
#    $self->{_AnalysisAdaptor} = Bio::EnsEMBL::Pipeline::DBSQL::AnalysisAdaptor->new
#      ( $self );
#  }
#
#  return $self->{_AnalysisAdaptor};
#}

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


sub prepare {
  my ($self, @args) = @_;

  $self->dbc->prepare(@args);
}


=head2 DESTROY

 Title   : DESTROY
 Usage   :
 Function:
 Example :
 Returns :
 Args    :

=cut

sub DESTROY {
   my ($obj) = @_;

   if( $obj->{'_db_handle'} ) {
       $obj->{'_db_handle'}->disconnect;
       $obj->{'_db_handle'} = undef;
   }
}

=head2 Some Utility Stuff

 Access to the meta table of the schema including
 
 pipeline_lock             - write a lock to the meta table
 pipeline_unlock           - remove the lock
 get_meta_value_by_key     - retrieve value by key
 store_meta_key_value      - write key, value pair
 remove_meta_key           - delete by key name
 make_meta_value_from_hash - flatten a hash to a string (uses dbi->quote to escape)
 make_hash_from_meta_value - returns a hash from a previously flattened string
 
    $lock_string = 'whatever to lock';
    $dbA->pipeline_lock($lock_string);
    my $lock_string = $dbA->pipeline_lock();
    
    $dbA->pipeline_unlock();
 
    $dbA->store_meta_key_value('my_key', 'the value');

    my $value = $dbA->get_meta_value_by_key('my_key');
    
    $dbA->remove_meta_key('my_key');

    my %hash = ('-host' => 'pfam', '-port' => '3306'...);
    my $flat = $dbA->make_meta_value_from_hash(\%hash);
    my %retrieved = %{$dbA->make_hash_from_meta_value($flat)};

=cut

sub pipeline_lock {
    my ($self, $string) = @_;
    my $lock = 'pipeline.lock';

    return $string ? $self->store_meta_key_value($lock, $string) : $self->get_meta_value_by_key($lock);
}

sub pipeline_unlock {
    my ($self) = @_;
    $self->remove_meta_key('pipeline.lock');
}

sub get_meta_value_by_key{
    my ($self, $key, $value) = @_;
    $self->throw("No key to get supplied") unless $key;
    my $sth = $self->prepare(qq{
        SELECT meta_value
        FROM   meta
        WHERE  meta_key = ? LIMIT 1
    }); # ONLY RETRIEVES FIRST ENTRY,
    # SHOULD THERE BE A UNIQUE KEY ON meta_key??
    $sth->execute($key);
    my $row = $sth->fetchrow_arrayref();
    $value = $row->[0] if $row;
    $sth->finish();
    return $value;
}
sub store_meta_key_value{
    my ($self, $key, $value) = @_;
    $self->throw("No key|value to insert supplied") unless $key && $value;
    my $sth = $self->prepare(qq{
        INSERT INTO meta (meta_key, meta_value) VALUES (?, ?)
    });
    $sth->execute($key, $value);
    $sth->finish();
    return undef;
}
sub remove_meta_key{
    my ($self, $key) = @_;
    $self->throw("No key to remove supplied") unless $key;
    my $sth = $self->prepare(qq{
	DELETE
	FROM   meta
	WHERE  meta_key = ?
    });
    $sth->execute($key);
    $sth->finish;
    return undef;
}
sub make_meta_value_from_hash{
    my ($self, $hash) = @_;
    my $dbh = $self->_db_handle();
    return join(",\n", map{ $dbh->quote($_)." => ".$dbh->quote($hash->{$_}) } keys(%$hash));
}
sub make_hash_from_meta_value{
    my ($self,$string) = @_;
    if($string){
        my $hash = { eval $string };
        $@ ? die "error evaluating $string" : return $hash || {};
    }
    return {};
}
1;
