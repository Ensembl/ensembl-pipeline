# Perl module for Bio::EnsEMBL::Pipeline::DBSQL::RuleAdaptor
#
# Creator: Arne Stabenau <stabenau@ebi.ac.uk>
# Date of creation: 10.09.2000
# Last modified : 10.09.2000 by Arne Stabenau
#
# Copyright EMBL-EBI 2000
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::DBSQL::RuleAdaptor

=head1 SYNOPSIS

  $jobAdaptor = $dbobj->getRuleAdaptor;
  $jobAdaptor = $jobobj->getRuleAdaptor;


=head1 DESCRIPTION

  Module to encapsulate all db access for persistent class Rule.
  There should be just one per application and database connection.


=head1 CONTACT

    Contact Arne Stabenau on implemetation/design detail: stabenau@ebi.ac.uk
    Contact Ewan Birney on EnsEMBL in general: birney@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Pipeline::DBSQL::RuleAdaptor;

use Bio::EnsEMBL::Pipeline::Rule;
use Bio::Root::RootI;
use vars qw(@ISA);
use strict;


@ISA = qw( Bio::Root::RootI );

=head2 Constructor

  Title   : new
  Usage   : $$dbobj->get_RuleAdaptor
  Function: Standard Adaptor Constructor
  Returns : Bio::EnsEMBL::Pipeline::DBSQL::RuleAdaptor
  Args    : Bio::EnsEMBL::Pipeline::DBSQL::Obj

=cut


sub new {
  my ($class,$dbobj) = @_;
  my $self = $class->SUPER::new();

  $self->db( $dbobj );
  return $self;
}

=head2 store

  Title   : store
  Usage   : $self->store( $rule );
  Function: Stores a rule in db
            Sets adaptor and dbID in Rule
  Returns : -
  Args    : Bio::EnsEMBL::Pipeline::Rule

=cut

sub store {
  my ( $self, $rule ) = @_;
  my $sth = $self->prepare( q{
    INSERT INTO RuleGoal
       SET goalAnalysisId = ? } );
  $sth->execute( $rule->goalAnalysis->dbID );
  $sth = $self->prepare( q {
    SELECT last_insert_id() } );
  $sth->execute;
  my $dbID = ($sth->fetchrow_array)[0];
  my @literals = $rule->list_conditions;
  for my $literal ( @literals ) {
    $sth = $self->prepare( qq{
      INSERT INTO RuleConditions
         SET ruleId=$dbID,
             conditionLiteral='$literal' } );
    $sth->execute;
  }
  $rule->dbID( $dbID );
  $rule->adaptor( $self );
}

=head2 remove

  Title   : remove
  Usage   : $self->remove( $rule );
  Function: removes given object from database.
  Returns : -
  Args    : Bio::EnsEMBL::Pipeline::Rule which must be persistent.
            ( dbID set )

=cut

sub remove {
  my ( $self, $rule ) = @_;

  my $dbID = $rule->dbID;
  if( !defined $dbID ) {
    $self->throw( "RuleAdaptor->remove called with non persistent Rule" );
  }

  my $sth = $self->prepare( qq {
    DELETE FROM RuleGoal
    WHERE ruleId = $dbID } );
  $sth->execute;
  $sth = $self->prepare( qq {
    DELETE FROM RuleConditions
     WHERE ruleID = $dbID } );
  $sth->execute;
}


=head2 fetch_all

  Title   : fetch_all
  Usage   : @rules = $self->fetch_all;
  Function: retrieves all rules from db;
  Returns : List of Bio::EnsEMBL::Pipeline::Rule
  Args    : -

=cut

sub fetch_all {
  my $self = shift;
  my $anaAdaptor = $self->db->get_AnalysisAdaptor;
  my %rules;
  my ( $analysis, $rule, $dbID );
  my @queryResult;

  my $sth = $self->prepare( q {
    SELECT ruleId,goalAnalysisId
      FROM RuleGoal } );
  $sth->execute;

  while( @queryResult = $sth->fetchrow_array ) {
    $analysis = $anaAdaptor->fetch_by_dbID( $queryResult[1] );
    $dbID = $queryResult[0];

    $rule = Bio::EnsEMBL::Pipeline::Rule->new
      ( '-dbid'    => $dbID,
	'-goal'    => $analysis,
        '-adaptor' => $self );
    # print STDERR "Setting $dbID rule\n";
    $rules{$dbID} = $rule;
  }

  $sth= $self->prepare( q{
    SELECT ruleId, conditionLiteral
      FROM RuleConditions } );
  $sth->execute;

  while( @queryResult = $sth->fetchrow_array ) {
      # print STDERR "@queryResult\n";
      $rules{$queryResult[0]}->add_condition( $queryResult[1] );
  }
  # print STDERR "Found @{[scalar keys %rules]} rules\n";
  return values %rules;
}

=head2 fetch_by_dbID

  Title   : fetch_by_dbID
  Usage   : $self->fetch_by_dbID
  Function: Standard fetch used by linked to objects
  Returns : Bio::EnsEMBL::Pipeline::Rule
  Args    : -

=cut

sub fetch_by_dbID {
  my $self = shift;
  my $dbID = shift;
  
  my $anaAdaptor = $self->db->get_AnalysisAdaptor;
  my ( $analysis, $rule );
  my $queryResult;

  my $sth = $self->prepare( q {
    SELECT ruleId,goalAnalysisId
      FROM RuleGoal 
      WHERE ruleId = ? } );
  $sth->execute( $dbID  );

  $queryResult = $sth->fetchrow_hashref;
  if( ! defined $queryResult ) {
    return undef;
  }
  
  $analysis = $anaAdaptor->fetch_by_dbID( $queryResult->{goalAnalysisId} );
      
  $rule = Bio::EnsEMBL::Pipeline::Rule->new
    ( '-dbid'    => $dbID,
      '-goal'    => $analysis,
      '-adaptor' => $self );

  $sth= $self->prepare( q{
    SELECT ruleId, conditionLiteral
      FROM RuleConditions 
      WHERE ruleId = ?} );
  $sth->execute( $dbID );
  
  while( $queryResult = $sth->fetchrow_hashref ) {
    $rule->add_condition( $queryResult->{conditionLiteral} );
  }
  return $rule;
}

=head2 db

  Title   : db
  Usage   : $self->db;
Function: gets the DBSQL::Obj for the Adaptor. Set is private.
  Returns : Bio::EnsEMBL::Pipeline::DBSQL::Obj;
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

sub create_tables{
  my ($self) = @_;
  my $sth;

  $sth = $self->prepare("drop table if exists RuleGoal");
  $sth->execute();

  $sth = $self->prepare(qq{
    CREATE TABLE RuleGoal (
    rule_id           int unsigned default 0 not null auto_increment,
    goal_analysis_id  int unsigned,

    PRIMARY KEY (rule_id)
    );
  });
  $sth->execute();

  $sth = $self->prepare("drop table if exists RuleConditions");
  $sth->execute();

  $sth = $self->prepare(qq{
    CREATE TABLE RuleConditions (
    rule_id            int not null,
    condition_literal  varchar(40)
    );
  });
  $sth->execute();
}

1;
