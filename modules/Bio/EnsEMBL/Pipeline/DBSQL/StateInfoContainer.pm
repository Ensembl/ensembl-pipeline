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

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::Analysis;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Sys::Hostname;
use vars qw(@ISA);
use strict;

@ISA = qw( Bio::EnsEMBL::Root );


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


=head2 fetch_analysis_by_input_id

Fetches all analyses
which have been completed on a specific input ID
Takes one string - input ID
and returns a ref to a list of Bio::EnsEMBL::Analysis.

=cut

sub fetch_analysis_by_input_id {
  my ($self, $inputId) = @_;

  my @result;
  my @row;

  my $anaAd = $self->db->get_AnalysisAdaptor();

  my $sth = $self->prepare( q {
    SELECT analysis_id
      FROM input_id_analysis
     WHERE input_id = ? } );
  $sth->execute($inputId);

  while( my $row = $sth->fetchrow_arrayref ) {
    my $analysis = $anaAd->fetch_by_dbID( $row->[0] );
    if($analysis ) {
      push( @result, $analysis );
    }
  }
  
  return \@result;
}


=head2 store_input_id_analysis

Stores an input ID, type and analysis [optionally runtime info].
Takes an input ID (as string), Bio::EnsEMBL::Analysis object, 
and [optionally save_runtime_info Boolean].
Throws an exception if any of the inputs are invalid or if 
the analysis object does not have a type.
 called by: Job->run_module

=cut

sub store_input_id_analysis {
  my ($self, $inputId, $analysis, $host, $save_runtime_info ) = @_;

  throw("[$analysis] is not a Bio::EnsEMBL::Pipeline::Analysis object")
   unless $analysis->isa("Bio::EnsEMBL::Pipeline::Analysis");

  throw("Invalid inputId [$inputId]")
   unless $inputId;

  throw("No type defined in analysis obj") 
   unless defined($analysis->input_id_type);

  # do we want to use a default class here?
# -----------------------------------------------------------------
  if($save_runtime_info){
  	print "Saving runtime info\n";
  	print join("\t",($inputId, $analysis->dbID, $host, $analysis->db_version)) . "\n";
      my $sth = $self->prepare(qq{
      	REPLACE INTO input_id_analysis
      	(input_id, input_id_type, analysis_id, created, runhost, db_version)
      	values (?, ?, ?, now(), ?, ?)
      	});
    $sth->execute($inputId, $analysis->input_id_type, $analysis->dbID, $host, $analysis->db_version || 'ERROR');
  }else{
      my $sth = $self->prepare(qq{
      	INSERT INTO input_id_analysis
      	(input_id, input_id_type, analysis_id, created, runhost)
      	values (?, ?, ?, now(), ?)
      	});
      $sth->execute($inputId, $analysis->input_id_type, $analysis->dbID, $host);
  }
#-----------------------------------------------------------------
}


=head2 list_input_id_created_by_analysis

As list_input_id_by_analysis() but
also returns the time the analysis was
created. Returns a ref to a list of (anonymous)
arrays.

  foreach my $id_time ($sic->list_input_id_created_by_analysis) {
    my $id      = $id_time->[0];
    my $created = $id_time->[1];
  }

=cut

sub list_input_id_created_by_analysis {
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

  return \@result;
}


=head2 list_input_id_by_Analysis

Takes a Bio::EnsEMBL::Analysis object and
returns a ref to a list of input IDs which have had
that analysis performed on them.

  foreach my $idlist ($sic->list_input_id_class_by_Analysis('SubmitSlice')) {
    my $id    = $idlist->[0];
    my $class = $idlist->[1];
  }


=cut

sub list_input_id_by_Analysis {
  my ($self, $analysis) = @_;

  if (ref $analysis && $analysis->isa("Bio::EnsEMBL::Analysis")) {
    $analysis = $analysis->dbID;
  }

  my @result;

  my $sth = $self->prepare(qq{
    SELECT input_id
      FROM input_id_analysis
     WHERE analysis_id = ?
  });

  $sth->execute($analysis);

  my $table = $sth->fetchall_arrayref;

  foreach my $row (@{$table}) {
    push @result, $row->[0];
  }

  return \@result;
}

=head2 list_input_id_by_Analysis_assembly_type

 
 As for list_input_id_by_Analysis but takes an optional list ref
of assembly types to return for. 

 $obj->list_input_id_by_Analysis_assembly_type($AnalysisObj, ['test', 'other']);

 This is mainly for the Finished Analysis Pipeline.

=cut

sub list_input_id_by_Analysis_assembly_type {
  my ($self, $analysis, $assembly) = @_;
  
  if (ref $analysis && $analysis->isa("Bio::EnsEMBL::Analysis")) {
    $analysis = $analysis->dbID;
  }
  $assembly ||= [];
  my @result;
  my $query = qq{
      SELECT input_id
        FROM input_id_analysis i, 
             contig c 
   LEFT JOIN assembly a ON c.contig_id = a.contig_id 
       WHERE i.input_id = c.name 
          && i.analysis_id = ?
      };
  if (scalar(@$assembly) > 0){
      $query .= "&& ( " . join("||", map{ " a.type = ? " } @$assembly ) . " )";
  }else{
      $query .= "&& a.type IS NULL";
  }

  my $sth = $self->prepare($query);
  
  $sth->execute($analysis, @$assembly);

  my $table = $sth->fetchall_arrayref;

  foreach my $row (@{$table}) {
    push @result, $row->[0];
  }

  return \@result;
}

=head2 list_input_id_class_by_start_count

Returns a ref to a list of all input_id and class,
with an optional start and end limit.

  # get 1st 100 entries from list

  foreach my $idlist ($sic->list_input_id_class_by_start_count(0, 100)) {
    my $id    = $idlist->[0];
    my $class = $idlist->[1];
  }


=cut

sub list_input_id_class_by_start_count {
  my $self = shift;
  my ($start,$count) = @_;
  my @result;
  my @row;

  my $query = qq{
    SELECT input_id, class
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

  return \@result;
}


sub get_all_input_id_analysis_sets {
  my $self = shift;
  my %id_type_hash;
  my @row;
  my @types;

  my $sth = $self->prepare( qq{
    SELECT DISTINCT input_id_type FROM input_id_type_analysis } );
  $sth->execute();

  while( @row = $sth->fetchrow_array ) {
    push @types,$row[0];
  }

  foreach my $type (@types) {
    my $ids = $self->list_input_ids_by_type($type);
    foreach my $id  (@$ids) {
      $id_type_hash{$type}{$id} = 1;
    }
  }
  return \%id_type_hash;
}

sub list_input_ids_by_type {
  my $self = shift;
  my $type = shift;
  my @ids;
  my @row;

  my $sth = $self->prepare( qq{
    SELECT distinct input_id FROM input_id_analysis where input_id_type=?
    });

  $sth->execute($type);

  while( @row = $sth->fetchrow_array ) {
    push @ids,$row[0];
  }

  return \@ids;
}

sub list_input_ids_by_analysis {
  my $self = shift;
  my $analysis_id = shift;
  my @ids;
  my @row;
  if(!$analysis_id){
    throw("need an analysis id list from");
  }
  my $sth = $self->prepare( qq{
    SELECT distinct input_id FROM input_id_analysis where analysis_id=?
    });

  $sth->execute($analysis_id);

  while( @row = $sth->fetchrow_array ) {
    push @ids,$row[0];
  }

  return \@ids;
}



=head2 delete_input_id

Takes an input ID (as a string) and removes
all matching entries from the pipeline
database. See also delete_input_id_class()
and delete_input_id_analysis().

=cut

sub delete_input_id {
  my $self    = shift;
  my $inputId = shift;

  my $sth = $self->prepare( qq{
    DELETE FROM input_id_analysis
    WHERE  input_id = ?} );
  $sth->execute($inputId);
}


=head2 delete_input_id_analysis

Takes an input ID (as a string) and an
analysis internal ID or Bio::EnsEMBL::Analysis
and removes
all matching entries from the pipeline
database. See also delete_input_id()
and delete_input_id_class().

=cut

sub delete_input_id_analysis {
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

  if($arg){ 
    $self->{'_db'} = $arg;
  }
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

1;
