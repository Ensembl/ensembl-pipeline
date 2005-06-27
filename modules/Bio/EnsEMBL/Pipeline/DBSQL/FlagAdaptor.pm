# Perl module for Bio::EnsEMBL::Pipeline::DBSQL::FlagAdaptor
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::DBSQL::FlagAdaptor

=head1 SYNOPSIS

 $Adaptor = $dbobj->getFlagAdaptor();


=head1 DESCRIPTION

Module to encapsulate all db access to flag table


=head1 CONTACT

Post general queries to ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::DBSQL::FlagAdaptor;

use Bio::EnsEMBL::Pipeline::Flag;
use Bio::EnsEMBL::Root;
use vars qw(@ISA);
use strict;

use Carp;
@ISA = qw( Bio::EnsEMBL::Root );

=head2 Constructor

  Title   : new
  Usage   : $dbobj->get_FlagAdaptor
  Function: Standard Adaptor Constructor
  Returns : Bio::EnsEMBL::Pipeline::DBSQL::FlagAdaptor
  Args    : Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor

=cut


sub new {
  my ($class,$dbobj) = @_;
  my $self = $class->SUPER::new();

  $self->db( $dbobj );
  return $self;
}

=head2 store

  Title   : store
  Usage   : $self->store( $flag );
  Function: Stores a flag in db
            Sets adaptor and dbID in Flag
  Returns : -
  Args    : Bio::EnsEMBL::Pipeline::Flag

=cut

sub store {
  my ( $self, $flag ) = @_;

  $self->check_flag($flag);
  my $sql = "INSERT INTO flag SET ensembl_id = '".
    $flag->ensembl_id."' , analysis_id = '".
    $flag->goalAnalysis->dbID."' , table_name = '".
      $flag->type."';";
  my $sth = $self->prepare($sql);
  $sth->execute;
  $sth = $self->prepare("SELECT last_insert_id()");
  $sth->execute;
  my $dbID = ($sth->fetchrow_array)[0];
  $flag->dbID( $dbID );
  $flag->adaptor( $self );
  $sth->finish;
  return 1;
}

=head2 remove

  Title   : remove
  Usage   : $self->remove( $flag );
  Function: removes given object from database.
  Returns : dbID of object that was removed
  Args    : Bio::EnsEMBL::Pipeline::Flag

=cut

sub remove {
  my ( $self, $flag ) = @_;

  my $dbID = $flag->dbID;
  if( !defined $dbID ) {
    $self->throw( "FlagAdaptor->remove called with non persistent Flag" );
  }

  my $sth = $self->prepare("
    DELETE FROM flag
    WHERE  flag_id = $dbID" );
  $sth->execute;
  $sth->finish;
  return $dbID;
}


=head2 fetch_all

  Title   : fetch_all
  Usage   : @flags = $self->fetch_all;
  Function: retrieves list ref  of all flags from db;
  Returns : List ref of Bio::EnsEMBL::Pipeline::Flag
  Args    : -

=cut

sub fetch_all {
  my $self = shift;
  my $anaAdaptor = $self->db->get_AnalysisAdaptor;
  my %flags;
  my ( $analysis,$id,$type, $dbID );


  my $sth = $self->prepare("
    SELECT flag_id, ensembl_id, table_name, analysis_id
    FROM   flag " );
  $sth->execute;

 FLAG: while( my($flag_id, $ensembl_id, $table_name, $analysis_id  ) = $sth->fetchrow_array ) {
    $analysis = $anaAdaptor->fetch_by_dbID($analysis_id) or do {
      $self->warn("Couldn't find analysis to match dbID $analysis_id");
      next FLAG;
    };
    $dbID = $flag_id;
    $id   = $ensembl_id;
    $type = $table_name;

  my  $flag = Bio::EnsEMBL::Pipeline::Flag->new
      ( '-dbid'    => $dbID,
	'-type'    => $type,
	'-ensembl_id'      => $id,
	'-goalAnalysis'    => $analysis,
        '-adaptor' => $self );
    $flags{$dbID} = $flag;
  }
  my @array = values %flags;
  $sth->finish;
  return \@array;
}

=head2 fetch_by_dbID

  Title   : fetch_by_dbID
  Usage   : $self->fetch_by_dbID
  Function: Fetches object by its db identifier
  Returns : Bio::EnsEMBL::Pipeline::Flag
  Args    : Scalar

=cut

sub fetch_by_dbID {
  my ($self, $dbID) = @_;
  my $anaAdaptor = $self->db->get_AnalysisAdaptor;
  my ( $analysis,$id,$type,$flag );
  my $queryResult;

  my $sth = $self->prepare("
    SELECT flag_id, ensembl_id, table_name, analysis_id
    FROM   flag
    WHERE  flag_id = $dbID");
  $sth->execute;

  my ($flag_id, $ensembl_id, $table_name, $analysis_id  ) = $sth->fetchrow;
  if( !defined $flag_id ) {
    return undef;
  }
  $analysis = $anaAdaptor->fetch_by_dbID($analysis_id)
    or $self->throw("Can't find analysis with dbID $analysis_id\n");

  $flag = Bio::EnsEMBL::Pipeline::Flag->new
      ( '-dbid'    => $dbID,
	'-type'    => $table_name,
	'-ensembl_id'      => $ensembl_id,
	'-goalAnalysis'    => $analysis,
        '-adaptor' => $self );  
  $sth->finish;
  return $flag;
}

=head2 fetch_by_analysis

  Title   : fetch_by_analysis
  Usage   : $self->fetch_by_analysis( $analysis );
  Function: fetches flag objects based on analysis object
  Returns : Array ref of Flag objects
  Args    : Bio::EnsEMBL::Analysis

=cut


sub fetch_by_analysis{
  my ($self, $goal_analysis) = @_;
  my @flags;
  if(!$goal_analysis || 
     !$goal_analysis->isa("Bio::EnsEMBL::Analysis")){
    $self->throw("analysis ".$goal_analysis." must be a ".
          "Bio:EnsEMBL::Analysis object");
  }
  my $sql = "SELECT flag_id
     FROM flag
     WHERE analysis_id  = ".$goal_analysis->dbID;
  my $sth = $self->prepare($sql);
  $sth->execute;
 FLAG: while( my($flag_id) = $sth->fetchrow_array ) {
   push @flags, $self->fetch_by_dbID($flag_id);
 }
  $sth->finish;
  return \@flags;
}

=head2 fetch_by_ensembl_id

  Title   : fetch_by_ensembl_id
  Usage   : $self->fetch_by_ensembl_id( $id );
  Function: fetches all flag objects with the specified  ensembl identifier 
            ie: a transcript dbid or gene dbid, not stable identifiers
  Returns : Array ref of Flag objects
  Args    : Scalar

=cut


sub fetch_by_ensembl_id{
  my ($self, $id) = @_;
  my @flags;
  my $sql = "SELECT flag_id FROM flag WHERE ensembl_id  = $id";
  my $sth = $self->prepare($sql);
  $sth->execute;  
  while ( my($result) = $sth->fetchrow_array){
    push @flags, $self->fetch_by_dbID($result);
  }
  $sth->finish;
  return \@flags;
}



sub check_flag{
  my ($self,$flag)=@_;
  unless ($flag->goalAnalysis->isa("Bio::EnsEMBL::Analysis")){
    $self->throw("analysis ".$flag->goalAnalysis." must be a ".
          "Bio:EnsEMBL::Analysis object");
  }

  my $sql = "show tables;";
  my $sth = $self->prepare($sql);
  $sth->execute;
  while ( my($result) = $sth->fetchrow_array){
    if ($result eq $flag->type){
      $sth->finish;
      return 1;
    }
  }
  $self->throw("Cannot find table that corresponds to flag type ".$flag->type."\n");
  $sth->finish;
  return 1;
}

=head2 db

  Title   : db
  Usage   : $self->db;
  Function: gets the DBSQL::DBAdaptor for the Adaptor. Set is private.
  Returns : Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
  Args    : -

=cut


sub db {
  my ($self,$db) = @_;
  ( defined $db ) &&
    ( $self->{'_db'} = $db );
  $self->{'_db'};
}


# Convenience prepare function
sub prepare {
  my ($self,$query)  = @_;
  $self->db->dbc->prepare( $query );
}


1;
