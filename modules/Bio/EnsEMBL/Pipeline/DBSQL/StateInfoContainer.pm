# Perl module for Bio::EnsEMBL::Pipeline::DBSQL::StateInfoContainer
#
# Creator: Arne Stabenau <stabenau@ebi.ac.uk>
#
# Date of creation: 15.09.2000
# Last modified : 20.09.2000 by Arne Stabenau
#
# Copyright EMBL-EBI 2000
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::DBSQL::StateInfoContainer

=head1 SYNOPSIS

  $sic = $dbobj->get_StateInfoContainer;

=head1 DESCRIPTION

Module which encapsulates state request for objects in the database.
Starts of with a table input_id_analysis, providing which analysis was
done to inputIds but every state information access should go via
this object.

Deliberatly NOT called an adaptor, as it does not serve objects.

=head1 CONTACT

=over 4

=item Arne Stabenau on implemetation/design detail: stabenau@ebi.ac.uk

=item Ewan Birney on EnsEMBL in general: birney@ebi.ac.uk

=back

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Pipeline::DBSQL::StateInfoContainer;

use Bio::Root::RootI;
use vars qw(@ISA);
use strict;

@ISA = qw( Bio::Root::RootI );


=head2 new

Object constructor. Takes a
Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor as argument.

Can also be constructed by calling
get_StateInfoContainer() on a
Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor.

=cut

sub new {
  my ($class, $dbobj) = @_;
  my $self = $class->SUPER::new();

  $self->db( $dbobj );
  return $self;
}


=head2 fetch_analysis_by_inputId_class

Fetches all analyses
which have been completed on a specific input ID
and class. Takes two strings - input and class
and returns a list of Bio::EnsEMBL::Analysis.

=cut

sub fetch_analysis_by_inputId_class {
  my ($self,$inputId, $class) = @_;

  my @result;
  my @row;

  my $anaAd = $self->db->get_AnalysisAdaptor();

  my $sth = $self->prepare( q {
    SELECT analysis_id
      FROM input_id_analysis
     WHERE input_id = ?
       AND class = ? } );
  $sth->execute( $inputId, $class );

  while( my @row = $sth->fetchrow_array ) {
    my $analysis = $anaAd->fetch_by_dbID( $row[0] );
    if( defined $analysis ) {
      push( @result, $analysis );
    }
  }

  return @result;
}


=head2 store_inputId_class_analysis

Stores an input ID, class and analysis.
Takes an input ID, class (as strings)
and a Bio::EnsEMBL::Analysis object. Throws an
exception if any of the inputs are invalid.

=cut

sub store_inputId_class_analysis {
  my ($self, $inputId, $class, $analysis ) = @_;

  $self->throw("[$analysis] is not a Bio::EnsEMBL::Analysis object")
   unless $analysis->isa("Bio::EnsEMBL::Analysis");

  $self->throw("Invalid inputId [$inputId]")
   unless $inputId;

  # do we want to use a default class here?

  $self->throw("Invalid class [$class]")
   unless $class;

  my $sth = $self->prepare(qq{
    INSERT INTO input_id_analysis
                (input_id, class, analysis_id, created)
                values (?, ?, ?, now())
  });
  $sth->execute($inputId, $class, $analysis->dbID);
}


=head2 list_inputId_by_analysis

Takes a Bio::EnsEMBL::Analysis object and
returns a list of input IDs which have had
that analysis performed on them.

=cut

sub list_inputId_by_analysis {
  my $self = shift;
  my $analysis = shift;
  my @result;
  my @row;

  my $sth = $self->prepare( q {
    SELECT input_id
      FROM input_id_analysis
     WHERE analysis_id = ? } );
  $sth->execute( $analysis->dbID );

  while( @row = $sth->fetchrow_array ) {
    push( @result, $row[0] );
  }

  return @result;
}


=head2 list_inputId_created_by_analysis

As list_inputId_by_analysis() but
also returns the time the analysis was
created. Returns an list of (anonymous)
arrays.

  foreach my $id_time ($sic->list_inputId_created_by_analysis) {
    my $id      = $id_time->[0];
    my $created = $id_time->[1];
  }

=cut

sub list_inputId_created_by_analysis {
  my $self = shift;
  my $analysis = shift;
  my @result;
  my @row;

  my $sth = $self->prepare( q {
    SELECT input_id, unix_timestamp(created)
      FROM input_id_analysis
     WHERE analysis_id = ? } );
  $sth->execute( $analysis->dbID );

  while( @row = $sth->fetchrow_array ) {
    push( @result, [$row[0], $row[1]] );
  }

  return @result;
}


=head2 list_inputId_class_by_start_count

Returns a list of all inputId and class,
with an optional start and end limit.

  # get 1st 100 entries from list

  foreach my $idlist ($sic->list_inputId_class_by_start_count(0, 100)) {
    my $id    = $idlist->[0];
    my $class = $idlist->[1];
  }


=cut

sub list_inputId_class_by_start_count {
  my $self = shift;
  my ($start,$count) = @_;
  my @result;
  my @row;

  my $query = qq{
    SELECT input_id, clasS
      FROM input_id_analysis
     GROUP by input_id, class };

  if( defined $start && defined $count ) {
    $query .= "LIMIT $start,$count";
  }
  my $sth = $self->prepare( $query );
  $sth->execute;

  while( @row = $sth->fetchrow_array ) {
    push( @result, [ $row[0], $row[1] ] );
  }

  return @result;
}


=head2 delete_inputId_class

Takes an input ID and class (as strings) and
removes all matching entries from the pipeline
database. See also
delete_inputId() and delete_inputId_analysis().

=cut

sub delete_inputId_class {
  my $self = shift;
  my ($inputId, $class) = @_;

  my $sth = $self->prepare( qq{
    DELETE FROM input_id_analysis
    WHERE  input_id = ?
    AND    class = ?} );
  $sth->execute($inputId, $class);
}


=head2 delete_inputId

Takes an input ID (as a string) and removes
all matching entries from the pipeline
database. See also delete_inputId_class()
and delete_inputId_analysis().

=cut

sub delete_inputId {
  my $self = shift;
  my ($inputId) = shift;

  my $sth = $self->prepare( qq{
    DELETE FROM input_id_analysis
    WHERE  input_id = ?} );
  $sth->execute($inputId);
}


=head2 delete_inputId_analysis

Takes an input ID (as a string) and an
analysis internal ID or Bio::EnsEMBL::Analysis
and removes
all matching entries from the pipeline
database. See also delete_inputId()
and delete_inputId_class().

=cut

sub delete_inputId_analysis {
  my ($self, $inputId, $analysis) = @_;
  my $analysisId;

  if (ref $analysis && $analysis->isa('Bio::EnsEMBL::Analysis')) {
    $analysisId = $analysis->dbID;
  }
  else {
    $analysisId = $analysis;
  }

  my $sth = $self->prepare( qq{
    DELETE FROM input_id_analysis
    WHERE  input_id    = ?
    AND    analysis_id = ?} );
  $sth->execute($inputId, $analysisId);
}


=head2 db

Method to get/set the DB object. Takes
and returns a Bio::EnsEMBL::DBSQL::DBAdaptor.

  my $job_adaptor = $sic->db->get_JobAdaptor;

=cut

sub db {
  my ( $self, $arg )  = @_;

  ( defined $arg ) &&
    ($self->{'_db'} = $arg);
  $self->{'_db'};
}


=head2 prepare

Convenience DBI prepare method.

  my $sth = $sic->prepare(q{
    SELECT distinct input_id from input_id_analysis
  });
  $sth->execute;

=cut

sub prepare {
  my ( $self, $query ) = @_;
  $self->db->prepare( $query );
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


=head2 create_tables

Create table needed by
Bio::EnsEMBL::Pipeline::DBSQL::StateInfoContainer

=cut

sub create_tables {
  my $self = shift;
  my $sth;

  $sth = $self->prepare("drop table if exists input_id_analysis");
  $sth->execute();

  $sth = $self->prepare(qq{
    CREATE TABLE input_id_analysis (
    input_id     varchar(40) not null,
    class        enum("clone","contig","vc","gene") not null,
    analysis_id  int not null,
    created      datetime not null,

    PRIMARY KEY (analysis_id, input_id, class),
    KEY input_id_created (input_id, created),
    KEY input_id_class   (input_id, class)
    );
  });
  $sth->execute();
}

1;
