# Perl module for Bio::EnsEMBL::DBSQL::AnalysisAdaptor
#
# Creator: Arne Stabenau <stabenau@ebi.ac.uk>
# Date of creation: 25.01.2001
# Last modified : 25.01.2001 by Arne Stabenau
#
# Copyright EMBL-EBI 2000
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::Pipeline::AnalysisAdaptor 

=head1 SYNOPSIS

  $analysisAdaptor = $db_adaptor->getAnalysisAdaptor;
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

use Bio::EnsEMBL::Pipeline::Analysis;
use Bio::EnsEMBL::DBSQL::AnalysisAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use vars qw(@ISA);
use strict;

@ISA = qw( Bio::EnsEMBL::DBSQL::AnalysisAdaptor);



=head2 new

  Args       : Bio::EnsEMBL::DBSQL::DBAdaptor
  Example    : my $aa = new Bio::EnsEMBL::Pipeline::DBSQL::AnalysisAdaptor();
  Description: Creates a new Bio::EnsEMBL::Pipeline::DBSQL::AnalysisAdaptor object and
               internally loads and caches all the Analysis objects from the 
               database.
  Returntype : Bio::EnsEMBL::Pipeline::DBSQL::AnalysisAdaptor
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBSQL::DBAdaptor

=cut

sub new {
  my ($class, $db) = @_;
 
  my $self = $class->SUPER::new($db);
 
  #load and cache all of the Analysis objects
  $self->fetch_all;

  return $self;
}






=head2 store

  Arg [1]    : Bio:EnsEMBL::Analysis $analysis 
  Example    : $analysis_adaptor->store($analysis);
  Description: stores $analysis in db. Does not if already equiped with dbID.
               Sets created date if not already set. Sets dbID and adaptor
               inside $analysis. Returns dbID.
  Returntype : int dbID of stored analysis
  Exceptions : thrown if analysis argument does not have a logic name
  Caller     : ?

=cut

sub store {

  my $self = shift;
  my $analysis = shift;
  
  if( !$analysis || !($analysis->isa('Bio::EnsEMBL::Pipeline::Analysis'))) {
    throw("called store on Pipeline::AnalysisAdaptor with a [$analysis]");
  }

  $analysis->dbID && $analysis->adaptor && ( $analysis->adaptor() == $self ) && 
    return $analysis->dbID;


  my $dbID;
 
  
  if( $dbID = $self->exists( $analysis )) {
    $analysis->adaptor( $self );
    $analysis->dbID( $dbID );
    return $dbID;
  }
 
  if( !$analysis->logic_name ) {
    throw("Must have a logic name on the analysis object");
  }
  if (!$analysis->input_id_type) {
    throw("Must have an input_id_type");
  }
 
  if($analysis->created ) {
    my $sth = $self->prepare( q{
      INSERT IGNORE INTO analysis
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
    $dbID = $sth->{'mysql_insertid'};
  } else {
    my $sth = $self->prepare( q{

      INSERT IGNORE INTO analysis
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

    $dbID = $sth->{'mysql_insertid'};

    if( $dbID ) {
      $sth = $self->prepare( q{
	SELECT created 
	FROM   analysis
	WHERE  analysis_id = ? } );
      $sth->execute( $dbID );
      $analysis->created( ($sth->fetchrow_array)[0] );
    }
  }
  
  $self->{_cache}->{$dbID} = $analysis;
  $self->{_logic_name_cache}{lc($analysis->logic_name)} = $analysis;

  $analysis->adaptor( $self );
  $analysis->dbID( $dbID );
 
  if($analysis->input_id_type){
    my $sql = 'insert into input_id_type_analysis(input_id_type, analysis_id) values(?, ?)';
    my $sth = $self->prepare($sql);
    $sth->execute($analysis->input_id_type, $analysis->dbID);
  }
  return $dbID;
}





=head2 _objFromHashref

  Arg [1]    : hashref $rowHash
  Description: Private helper function generates an Analysis object from a 
               mysql row hash reference.
  Returntype : Bio::EnsEMBL::Analysis
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBSQL::AnalsisAdaptor::fetch_* methods

=cut
  
sub _objFromHashref {
  my $self = shift;
  my $rowHash = shift;
  
  my $analysis = Bio::EnsEMBL::Pipeline::Analysis->new
    (
     -id              => $rowHash->{analysis_id},
     -adaptor         => $self,
     -db              => $rowHash->{db},
     -db_file         => $rowHash->{db_file},
     -db_version      => $rowHash->{db_version},
     -program         => $rowHash->{program},
     -program_version => $rowHash->{program_version},
     -program_file    => $rowHash->{program_file},
     -gff_source      => $rowHash->{gff_source},
     -gff_feature     => $rowHash->{gff_feature},
     -module          => $rowHash->{module},
     -module_version  => $rowHash->{module_version},
     -parameters      => $rowHash->{parameters},
     -created         => $rowHash->{created},
     -logic_name      => $rowHash->{logic_name},
    );
  
  my $type = $self->fetch_analysis_input_id_type($analysis);
  $analysis->input_id_type($type);
  return $analysis;
}

sub fetch_analysis_input_id_type{
  my ($self, $analysis) = @_;
  
  my $sql = "select input_id_type from input_id_type_analysis".
    " where analysis_id = ?";
 
  my $sth = $self->prepare($sql);
  $sth->execute($analysis->dbID);
  my ($type) = $sth->fetchrow;
  return $type;
}

=head2 db

  Arg [1]    : (optional) Bio::EnsEMBL::DBSQL::DBAdaptor $db
               the database used by this adaptor.
  Example    : my $db = $analysis_adaptor->db()
  Description: Getter/Setter for the database this adaptor uses internally
               to fetch and store database objects.
  Returntype : Bio::EnsEMBL::DBSQL::DBAdaptor
  Exceptions : none
  Caller     : BaseAdaptor::new, general

=cut

sub db {
  my ( $self, $arg )  = @_;
  ( defined $arg ) &&
    ($self->{_db} = $arg);
  $self->{_db};;
}


sub remove {

  my $self = shift;
  my $analysis = shift;
  
  if( !$analysis || !($analysis->isa('Bio::EnsEMBL::Pipeline::Analysis'))) {
    throw("called store on Pipeline::AnalysisAdaptor with a [$analysis]");
  }
  my $sql = 'delete from analysis where analysis_id = ?';
  my $sth = $self->prepare($sql);
  $sth->execute($analysis->dbID);

  if($analysis->input_id_type){
    my $sql = 'delete from input_id_type_analysis where analysis_id = ?';
    my $sth = $self->prepare($sql);
    $sth->execute($analysis->dbID);
  }

}


sub fetch_by_input_id_type{
  my ($self, $input_id_type) = @_;
  
  my $query = q{ SELECT analysis_id
                 FROM input_id_type_analysis
                 WHERE input_id_type = ? };
  my $sth = $self->prepare($query);
  $sth->execute($input_id_type);
  my @analysis_ids;
  while(my ($id) = $sth->fetchrow){
    push(@analysis_ids, $id);
  }
  my @analyses;
  foreach my $id(@analysis_ids){
    my $analysis = $self->fetch_by_dbID($id);
    push(@analyses, $analysis);
  }
  return(\@analyses);
}

1;
