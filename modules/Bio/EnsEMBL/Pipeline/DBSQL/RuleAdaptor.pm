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

 $jobAdaptor = $dbobj->getRuleAdaptor();
 $jobAdaptor = $jobobj->getRuleAdaptor();

=head1 DESCRIPTION

Module to encapsulate all db access for persistent class Rule.
There should be just one per application and database connection.

=head1 CONTACT

Contact Arne Stabenau on implemetation/design detail: stabenau@ebi.ac.uk
Contact Ewan Birney on EnsEMBL in general: birney@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Pipeline::DBSQL::RuleAdaptor;

use Bio::EnsEMBL::Pipeline::Rule;
use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Utils::Exception qw (throw warning) ; 
use vars qw(@ISA);
use strict;

use Carp;
@ISA = qw( Bio::EnsEMBL::Root );

=head2 Constructor

  Title   : new
  Usage   : $dbobj->get_RuleAdaptor
  Function: Standard Adaptor Constructor
  Returns : Bio::EnsEMBL::Pipeline::DBSQL::RuleAdaptor
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
  Usage   : $self->store( $rule );
  Function: Stores a rule in db
            Sets adaptor and dbID in Rule
  Returns : -
  Args    : Bio::EnsEMBL::Pipeline::Rule

=cut

sub store {
  my ( $self, $rule ) = @_;

  my $sth = $self->prepare( q{
    INSERT INTO rule_goal
       SET goal = ? } ); 

  throw("Can't store a rule_goal without a dbID") 
  unless ( $rule->goalAnalysis->dbID ) ; 

  $sth->execute( $rule->goalAnalysis->dbID );

  $sth = $self->prepare( q {
    SELECT last_insert_id() } );
  $sth->execute;
  my $dbID = ($sth->fetchrow_array)[0];

  my @literals = @{$rule->list_conditions};
  for my $literal ( @literals ) {
    $sth = $self->prepare( qq{
      INSERT INTO rule_conditions
         SET rule_id = $dbID,
             rule_condition = '$literal' } );
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
    DELETE FROM rule_goal
    WHERE  rule_id = $dbID } );
  $sth->execute;

  $sth = $self->prepare( qq {
    DELETE FROM rule_conditions
    WHERE  rule_id = $dbID } );
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
    SELECT rule_id, goal
    FROM   rule_goal } );
  $sth->execute;

  RULE: while( my($ruleId, $goal) = $sth->fetchrow_array ) {
    if ($goal =~ /\D/) {
      $analysis = $anaAdaptor->fetch_by_newest_logic_name($goal) or do {
        $self->warn("Couldn't find analysis to match logic_name $goal");
	next RULE;
      };
    }
    else {
      $analysis = $anaAdaptor->fetch_by_dbID($goal) or do {
        $self->warn("Couldn't find analysis to match dbID $goal");
	next RULE;
      };
    }
    $dbID = $ruleId;


    $rule = Bio::EnsEMBL::Pipeline::Rule->new
      ( '-dbid'    => $dbID,
	'-goalAnalysis'    => $analysis,
        '-adaptor' => $self );
    $rules{$dbID} = $rule;
  }



  $sth= $self->prepare( q{
    SELECT rule_id, rule_condition
    FROM   rule_conditions } );
  $sth->execute;

  while( my($ruleId, $cond) = $sth->fetchrow_array ) {
      $rules{$ruleId}->add_condition($cond) if defined $rules{$ruleId};
  }
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
    SELECT rule_id, goal
    FROM   rule_goal
    WHERE  rule_id = ? } );
  $sth->execute( $dbID  );

  $queryResult = $sth->fetchrow_hashref;
  if( ! defined $queryResult ) {
    return undef;
  }

  if ((my $goal = $queryResult->{goal}) =~ /\D/) {
    $analysis = $anaAdaptor->fetch_by_newest_logic_name($goal)
     or $self->throw("Can't find analysis with logic_name $goal");
  }
  else {
    $analysis = $anaAdaptor->fetch_by_dbID($goal)
     or $self->throw("Can't find analysis with dbID $goal");
  }


  $rule = Bio::EnsEMBL::Pipeline::Rule->new
    ( '-dbid'    => $dbID,
      '-goalAnalysis'    => $analysis,
      '-adaptor' => $self );

  $sth= $self->prepare( q{
     SELECT rule_id, rule_condition
     FROM   rule_conditions
     WHERE  rule_id = ?} );
  $sth->execute( $dbID );

  while( $queryResult = $sth->fetchrow_hashref ) {
    $rule->add_condition( $queryResult->{rule_condition} );
  }
  return $rule;
}

sub fetch_by_goal{
  my ($self, $goal_analysis) = @_;
  
  if(!$goal_analysis || 
     !$goal_analysis->isa("Bio::EnsEMBL::Pipeline::Analysis")){
    throw("analysis ".$goal_analysis." must be a ".
          "Bio:EnsEMBL::Pipeline::Analysis object");
  } 

  my ( $sql, $rule, $sth ) ;

  unless ( $goal_analysis->dbID ) {  
    my $anaAdaptor = $self->db->get_AnalysisAdaptor;
    my $ta= $anaAdaptor->fetch_by_logic_name($goal_analysis->logic_name) ; 
    $goal_analysis->dbID($ta->dbID) ;  
  } 

  if ( $goal_analysis->dbID ) { 
    $sql = q{ SELECT rule_id
                 FROM rule_goal
                 WHERE goal = ?
               };
    $sth = $self->prepare($sql);
    $sth->execute($goal_analysis->dbID);
    my ($dbID) = $sth->fetchrow;
    $rule = $self->fetch_by_dbID($dbID);

  } else { 
    throw(" can't find analysis in database\n" ) ; 
  } 

  warning("Can't get the dbID for the analysis with logic_name ".
        $goal_analysis->logic_name . "in database " . 
        $self->db->dbname ) unless $rule ;

  return $rule;
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
