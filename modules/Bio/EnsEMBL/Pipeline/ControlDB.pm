# Copyright EMBL-EBI 2000
# Author: Arne Stabenau
# Creation: 11.07.2000


# module to handle object_state table entries
# contains convenience functions
# not yet objects for handeling the states.

package Bio::EnsEMBL::Pipeline::ControlDB;

use strict;
use DBI;
use vars qw( @ISA );

# maybe import the configuration 
BEGIN {
  require "Bio/EnsEMBL/Pipeline/pipeConf.pl";
}

@ISA = qw( Bio::Root::RootI );

sub new {
  # without the pipeConf file nothing goes, so all parameters are there
  # may chage soonish
  
  my $self = bless {};    

  my $dsn = "DBI:".$::pipeConf{'DBI.driver'}.":database=".$::pipeConf{'ControlDB.name'}.";host=".$::pipeConf{'ControlDB.host'};
  my $dbh =  DBI->connect( $dsn, $::pipeConf{'ControlDB.user'}, $::pipeConf{'ControlDB.pass'}, {RaiseError => 1});
  $self->db( $dbh );
  return $self;
}


# make transition from one state to the other 
# for a given object_state_id
# line must be in transit for that
sub changeState {
  my $self = shift;
  my $obj_state_id = shift;
  my $to_nick = shift; 
  my @result;
   
  my $sth = $self->db->prepare( q{
  
   	SELECT o.inTransition, t.state_description_id, o.state_description_id 
	FROM object_state o, state_description s, state_description t
	WHERE o.object_state_id = ?
	  AND o.state_description_id = s.state_description_id 
	  AND s.object_class = t.object_class
	  AND t.state_nickname = ?
	} );

  $sth->execute( $obj_state_id, $to_nick );

  if( @result = $sth->fetchrow_array ) {
    if( $result[0] eq 'true' ) {
      $sth = $self->db->prepare( q{
        UPDATE object_state 
	SET last_change = NOW()
	    ,inTransition = 'false'
	    ,state_description_id = ?
	WHERE
	   object_state_id = ?
      } );
      $sth->execute( $result[1], $obj_state_id );
    } else {
      $self->throw( "Object not in transit, cant change state" );
    }
  } else {
    $self->throw( "Object not available, cant change state" );
  }
  
  # Log the change 
  # the log entry must already be there, just put the finish_transit
  # time into the row
  # objects should not pass the same state twice !!
  $sth = $self->db->prepare( q{
    UPDATE transition_log
       SET finish_transit = NOW()
     WHERE object_state_id = ?
       AND state_description_id = ?
       AND ISNULL( finish_transit )
  } );
  $sth->execute( $obj_state_id, $result[2] );
}

# get_by_criteria 
# many functions to get collection of hashes with 
# information about object_states.

# return all object_status which are in transit
# and have to be monitored. Used from TransitionManager.
sub get_inTransit_toMonitor { 
  my $self = shift;
  my @result;
  my $queryResult;
  
  my $cols = "object_state_id, object_id, object_class, state_nickname";
  
  my $query = qq{
    SELECT $cols
      FROM object_state o, state_description s
     WHERE o.inTransition = 'true'
       AND o.state_description_id = s.state_description_id
       AND s.needsMonitoring = 'true' };
  my $sth = $self->db->prepare( $query );
  $sth->execute;
  while( $queryResult = $sth->fetchrow_hashref ) {
    push( @result, $queryResult );
  }
  return @result; 
}

# return where the transit_module has to be started
sub get_nonFinal_nonTransit {
  my $self = shift;
  my @result;
  my $queryResult;
  
  my $cols = "object_state_id, object_id, object_class, state_nickname, transition_module, needsMonitoring";
  
  my $query = qq{
    SELECT $cols
      FROM object_state o, state_description s
     WHERE o.inTransition = 'false'
       AND o.state_description_id = s.state_description_id
       AND s.endState = 'false' };
  my $sth = $self->db->prepare( $query );
  $sth->execute;
  while( $queryResult = $sth->fetchrow_hashref ) {
    push( @result, $queryResult );
  }
  return @result; 
}

# used by transition monitors and aware transition modules
# to get what they want from the ControlDB
sub get_byIds_age {
  my $self = shift;
  my $minutes = shift;
  my $idListRef = shift;
  my $columnsListRef = shift;
  my @result;
  my $queryString;
  
  my $cols = "object_state_id, object_id, object_class, state_nickname, transition_module, needsMonitoring";
  if( scalar @$columnsListRef > 0 ) {
    $cols .= ", ".join( ", ",@$columnsListRef );
  }
  
  if( scalar @$idListRef > 0 ) {
    $queryString = "object_state_id in (".join(",",@$idListRef ).")";
  } else {
    return ();
  }
  
  my $query = qq{
    SELECT $cols
      FROM object_state o, state_description s
     WHERE $queryString
       AND o.state_description_id = s.state_description_id
       AND o.last_change < DATE_SUB( NOW(), INTERVAL $minutes MINUTE )
  };
 
  my $sth = $self->db->prepare( $query );
  $sth->execute;
  while( my $queryResult = $sth->fetchrow_hashref ) {
    push( @result, $queryResult );
  }
  return @result; 
}
  

# put an object from object_state table into transit mode
# log the time if requested.
sub object_toTransit {
  my $self = shift;
  my $obj_state_id_list = shift;
  my $resetTransit = shift;
  
  my $whereClause;
  if( ! defined $resetTransit ) {
    $resetTransit = "'true'";
  } else {
    $resetTransit = "'false'";
  }
  
  # check if the object is in transit ?
  # probably not, extra unnecessary query
  if( scalar @$obj_state_id_list == 1 ) {
    $whereClause = "object_state_id = ".$obj_state_id_list->[0];
  } else {
    $whereClause = "object_state_id in ( ".
    	join( ",",@$obj_state_id_list )." )";
  }
  my $sth = $self->db->prepare( qq{
    UPDATE object_state
       SET inTransition = $resetTransit
         , last_change = NOW()
     WHERE $whereClause
  } );
  $sth->execute;
  
  $sth = $self->db->prepare( qq{
    INSERT into transition_log
       ( object_state_id, object_id, state_description_id, start_transit )
    SELECT object_state_id, object_id, state_description_id, last_change
      FROM object_state
     WHERE $whereClause
  } );
  $sth->execute;
}

# very dirty, should seperate out the transit handling function

# here you can reset transits you have set before. Do so, when executing
# the TransitionModule start fails. If it starts, then dont!
sub object_reset_transit {
  my $self = shift;
  my $obj_state_id_list = shift;
  $self->object_toTransit( $obj_state_id_list, 'false' );
}


# put an object into the pipeline
# missing, a check so that the same object goes not in again.

sub submit {
  my $self = shift;
  my $object_id = shift;
  my $class = shift;
  my $state = shift;
  my @queryResult;
  
  # get the state_description_id
  my $sth = $self->db->prepare( q{
    SELECT state_description_id 
      FROM state_description
     WHERE object_class = ?
       AND state_nickname = ?
  } );
  $sth->execute( $class, $state );
  if( @queryResult = $sth->fetchrow_array ) {
    $sth = $self->db->prepare( q{
      INSERT INTO object_state( object_id, state_description_id, last_change, inTransition )
      VALUES( ?, ?, now(), 'false' )
    } );
    $sth->execute( $object_id, $queryResult[0] );
  }
}

sub db {
  my $self = shift;
  my $db = shift;
  $db && 
    ( $self->{_db} = $db );
  $self->{_db};
}


# convenience function to get rid of %pipeConf once
# maybe have params in the database?
sub get_param {
  my $self = shift;
  my $paramname = shift;
  
  return $::pipeConf{$paramname};
}
 
sub log_error {
  my $self = shift;
  my $obj_state_id = shift;
  my $msg = shift;
  
  my $sth = $self->db->prepare( q{
    INSERT into object_error( object_state_id, error_mesg )
    VALUES ( ?, ? )
  } );
  $sth->execute( $obj_state_id, $msg );
}

  
  
1;
