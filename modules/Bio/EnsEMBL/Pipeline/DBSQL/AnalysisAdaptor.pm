# Perl module for Bio::EnsEMBL::Pipeline::DBSQL::AnalysisAdaptor
#
# Creator: Arne Stabenau <stabenau@ebi.ac.uk>
# Date of creation: 05.09.2000
# Last modified : 05.09.2000 by Arne Stabenau
#
# Copyright EMBL-EBI 2000
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::DBSQL::AnalysisAdaptor

=head1 SYNOPSIS

  $analysisAdaptor = $dbobj->getAnalysisAdaptor;
  $analysisAdaptor = $analysisobj->getAnalysisAdaptor;


=head1 DESCRIPTION

  Module to encapsulate all db access for persistent class Analysis.
  There should be just one per application and database connection.


=head1 CONTACT

    Contact Arne Stabenau on implemetation/design detail: stabenau@ebi.ac.uk
    Contact Ewan Birney on EnsEMBL in general: birney@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::DBSQL::AnalysisAdaptor;

use Bio::EnsEMBL::Analysis;
use Bio::Root::RootI;

use vars qw(@ISA);
use strict;

@ISA = qw( Bio::Root::RootI );

sub new {
  my ($class, $dbobj) = @_;
  my $self = $class->SUPER::new();

  $self->db( $dbobj );
  return $self;
}

=head2 fetch_all

  Title   : fetch_all
  Usage   : @analyses = $self->fetch_all;
  Function: retrieves all analyses from db;
  Returns : List of Bio::EnsEMBL::Pipeline::Analysis
  Args    : -

=cut

sub fetch_all {
  my ($self) = @_;
  my %analyses;
  my ( $analysis, $dbID );
  my $rowHashRef;

  my $sth = $self->prepare( q {
    SELECT analysisId, logic_name,
           program,program_version,program_file,
           db,db_version,db_file,
           module,module_version,
           gff_source,gff_feature,
           created, parameters
    FROM analysisprocess } );
  $sth->execute;

  while( $rowHashRef = $sth->fetchrow_hashref ) {
    my $analysis = $self->_objFromHashref( $rowHashRef  );
    $analyses{$analysis->dbID} = $analysis;
  }

  return values %analyses;
}

=head2 fetch_by_dbID

  Title   : fetch_by_dbID
  Usage   : my $analysis = $adaptor->fetch_by_dbID
  Function: Retrieves an analysis from database by internal id
  Returns : throws exception when something goes wrong.
            undef if the id is not in the db.
  Args    :

=cut

sub fetch_by_dbID {
  my ($self,$id) = @_;

  if( defined $self->{'_cache'}->{$id} ) {
    return $self->{'_cache'}->{$id};
  }

  my $sth = $self->prepare( q{
    SELECT analysisId, logic_name,
           program,program_version,program_file,
           db,db_version,db_file,
           module,module_version,
           gff_source,gff_feature,
           created, parameters
    FROM analysisprocess
    WHERE analysisId = ? } );

  $sth->execute( $id );
  my $rowHashRef = $sth->fetchrow_hashref;
  if( ! defined $rowHashRef ) {
    return undef;
  }

  my $anal = $self->_objFromHashref( $rowHashRef );
  $self->{'_cache'}->{$anal->dbID} = $anal;
  return $anal;
}

sub fetch_by_newest_logic_name {
  my ($self,$logic_name) = @_;

  my $sth = $self->prepare( q{
    SELECT analysisId, logic_name,
           program,program_version,program_file,
           db,db_version,db_file,
           module,module_version,
           gff_source,gff_feature,
           created, parameters
    FROM analysisprocess
    WHERE logic_name = ?
    ORDER BY created DESC } );

  $sth->execute( $logic_name );
  my $rowHashRef = $sth->fetchrow_hashref;
  if( ! defined $rowHashRef ) {
    return undef;
  }

  return $self->_objFromHashref( $rowHashRef );
}


sub fetch_by_logic_name {
  my ($self,$logic_name) = @_;
  my @result;
  my $analysis;
  my $rowHash;

  my $sth = $self->prepare( q{
    SELECT analysisId, logic_name,
           program,program_version,program_file,
           db,db_version,db_file,
           module,module_version,
           gff_source,gff_feature,
           created, parameters
    FROM analysisprocess
    WHERE logic_name = ?
    ORDER BY created DESC } );

  $sth->execute( $logic_name );
  my $rowHashRef;
  while( $rowHashRef = $sth->fetchrow_hashref ) {
       $analysis = $self->_objFromHashref( $rowHashRef );
    if( defined $analysis ) {
      push( @result, $analysis );
    }
  }
  return @result;
}

# store makes dbID for analysis object
# sets the creation time in created if it wasnt set before
sub store {
  my ($self,$analysis) = @_;

  if( defined $analysis->created ) {
    my $sth = $self->prepare( q{
      INSERT INTO analysisprocess
      SET created = ?,
          logic_name = ?,
	  db = ?,
	  db_version = ?,
          db_file = ?,
          program = ?,
          program_version = ?,
          program_file = ?,
	  parameters = ?,
          module = ?,
          module_version = ?,
          gff_source = ?,
          gff_feature = ? } );
    $sth->execute
      ( $analysis->created,
	$analysis->logic_name,
	$analysis->db,
	$analysis->db_version,
	$analysis->db_file,
	$analysis->program,
	$analysis->program_version,
	$analysis->program_file,
	$analysis->parameters,
	$analysis->module,
	$analysis->module_version,
	$analysis->gff_source,
	$analysis->gff_feature
      );
    $sth = $self->prepare( q{
      SELECT last_insert_id() ;
    } );
    $sth->execute;
    my $dbID = ($sth->fetchrow_array)[0];
    $analysis->dbID( $dbID );
  } else {
    my $sth = $self->prepare( q{

      INSERT INTO analysisprocess
      SET created = now(),
          logic_name = ?,
	  db = ?,
	  db_version = ?,
          db_file = ?,
          program = ?,
          program_version = ?,
          program_file = ?,
	  parameters = ?,
          module = ?,
          module_version = ?,
          gff_source = ?,
          gff_feature = ? } );

    $sth->execute
      ( $analysis->logic_name,
	$analysis->db,
	$analysis->db_version,
	$analysis->db_file,
	$analysis->program,
	$analysis->program_version,
	$analysis->program_file,
	$analysis->parameters,
	$analysis->module,
	$analysis->module_version,
	$analysis->gff_source,
	$analysis->gff_feature
      );

    $sth = $self->prepare( q{
      SELECT last_insert_id()
    } );
    $sth->execute;

    my $dbID = ($sth->fetchrow_array)[0];
    $analysis->dbID( $dbID );
    if( defined $dbID ) {
      $sth = $self->prepare( q{
	SELECT created
	FROM analysisprocess
	WHERE analysisId = ? } );
      $sth->execute( $dbID );
      $analysis->created( ($sth->fetchrow_array)[0] );
    }
  }
}

=head2 exists

 Title   : exists
 Usage   : $adaptor->exists($anal)
 Function: Tests whether this Analysis already exists in the database
 Example :
 Returns : newest Analysis object which has all what given analysis have.
 Args    : Bio::EnsEMBL::Analysis

=cut

sub exists {
    my ($self,$anal) = @_;
    my $resultAnalysis;

    $self->throw("Object is not a Bio::EnsEMBL::Analysis") unless $anal->isa("Bio::EnsEMBL::Analysis");

    my $query;
    my @conditions;
    push( @conditions, "program=\"".$anal->program."\"" ), if( defined  $anal->program );
    push( @conditions, "program_version=\"".$anal->program_version."\""), if( defined  $anal->program_version );
    push( @conditions, "program_file=\"".$anal->program_file."\""), if( defined  $anal->program_file );
    push( @conditions, "db=\"".$anal->db."\""), if( defined  $anal->db );
    push( @conditions, "db_version=\"".$anal->db_version."\""), if( defined  $anal->db_version );
    push( @conditions, "db_file=\"".$anal->db_file."\""), if( defined  $anal->db_file );
    push( @conditions, "gff_source=\"".$anal->gff_source."\""), if( defined  $anal->gff_source );
    push( @conditions, "gff_feature=\"".$anal->gff_feature."\""), if( defined  $anal->gff_feature );
    push( @conditions, "module=\"".$anal->module."\""), if( defined  $anal->module );
    push( @conditions, "module_version=\"".$anal->module_version."\""), if( defined  $anal->module_version );
    push( @conditions, "parameters=\"".$anal->parameters."\""), if( defined  $anal->parameters );
    push( @conditions, "logic_name=\"".$anal->logic_name."\""), if( defined  $anal->logic_name );

    $query = qq { SELECT analysisId, logic_name,
		     program,program_version,program_file,
		     db,db_version,db_file,
		     module,module_version,
		     gff_source,gff_feature,
		     created, parameters
		     FROM analysisprocess
		     WHERE }.
		       join( " and ", @conditions )." order by created DESC";
    my $sth = $self->prepare($query);
    my $rv  = $sth->execute();

    if ($rv && $sth->rows > 0) {
      my $rowHash = $sth->fetchrow_hashref;
      $resultAnalysis = _objFromHashref( $rowHash );
      return $resultAnalysis;
    } else {
      return undef;
    }
}


sub _objFromHashref {
  my ($self,$rowHash) = @_;

  my $analysis = Bio::EnsEMBL::Analysis->new
    ( -id => $rowHash->{analysisId},
      -db => $rowHash->{db},
      -db_file => $rowHash->{db_file},
      -program => $rowHash->{program},
      -program_version => $rowHash->{program_version},
      -program_file => $rowHash->{program_file},
      -gff_source => $rowHash->{gff_source},
      -gff_feature => $rowHash->{gff_feature},
      -module => $rowHash->{module},
      -module_version => $rowHash->{module_version},
      -parameters => $rowHash->{parameters},
      -created => $rowHash->{created},
      -logic_name => $rowHash->{logic_name}
    );

  return $analysis;
}

# fixme: needs renaming to removeInputId_class_analysis?
sub removeInputId_analysis {
  my ($self,$inputid,$class,$analysis) = @_;

  if (!defined($inputid)) {
    $self->throw("No input id defined");
  }
  if (!defined($analysis)) {
    $self->throw("Analysis not defined");
  }

  my $query = "delete from InputIdAnalysis where inputId = '$inputid' and analysisId = " . $analysis->dbID . " and class = '" . $class . "'";

  my $sth = $self->prepare($query);
  my $rv  = $sth->execute();
}

# fixme: needs renaming to removeInputId_class?
sub removeInputId {
  my ($self,$inputid,$class) = @_;

  if (!defined($inputid)) {
    $self->throw("No input id defined");
  }

  my $query = "delete from InputIdAnalysis where inputId = '$inputid' and class = '" . $class . "'";

  my $sth = $self->prepare($query);
  my $rv  = $sth->execute();
}

sub submitInputId {
  my ($self,$inputid,$class,$analysis) = @_;

  if (!defined($inputid)) {
    $self->throw("No input id defined");
  }
  if (!defined($class)) {
    $self->throw("No class defined");
  }
  if (!defined($analysis)) {
    $self->throw("Analysis not defined");
  }

  my $query = "insert into InputIdAnalysis (inputId,class,analysisId,created) values(\'$inputid\',\'$class\', ". $analysis->dbID . ",now())";

  my $sth = $self->prepare($query);
  my $rv  = $sth->execute();
}

sub db {
  my ( $self, $arg )  = @_;
  ( defined $arg ) &&
    ($self->{'_db'} = $arg);
  $self->{'_db'};
}

sub prepare {
  my ( $self, $query ) = @_;
  $self->db->prepare( $query );
}


sub deleteObj {
  my ($self) = @_;
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
