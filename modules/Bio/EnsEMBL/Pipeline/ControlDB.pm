# Copyright EMBL-EBI 2000
# Author: Arne Stabenau
# Creation: 11.07.2000


=head1 NAME

Bio::EnsEMBL::Pipeline::ControlDB - object_state table entries and convience functions

=head1 SYNOPSIS


=head1 DESCRIPTION

Handle object_state table entries and convience functions

=head1 CONTACT

Arne Stabenau
ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# module to handle object_state table entries
# contains convenience functions
# not yet objects for handeling the states.

package Bio::EnsEMBL::Pipeline::ControlDB;

use strict;
use DBI;
use Bio::Root::RootI;

use vars qw( @ISA );

# maybe import the configuration
BEGIN {
  require "Bio/EnsEMBL/Pipeline/pipeConf.pl";
}

@ISA = qw( Bio::Root::RootI );

=head2 new

 Title   : new
 Usage   : my $controldb = new Bio::EnsEMBL::Pipeline::ControlDB
 Function:
 Example :
 Returns : 
 Args    :

=cut

sub new {
    # without the pipeConf file nothing goes, so all parameters are there
    # may chage soonish
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);

    my $dsn = "DBI:".$::pipeConf{'DBI.driver'}.":database=".$::pipeConf{'ControlDB.name'}.";host=".$::pipeConf{'ControlDB.host'};
    my $dbh =  DBI->connect( $dsn, $::pipeConf{'ControlDB.user'}, $::pipeConf{'ControlDB.pass'}, {RaiseError => 1});
    $self->db( $dbh );
    return $self;
}


# make transition from one state to the other
# for a given object_state_id
# line must be in transit for that

=head2 changeState

 Title   : changeState
 Usage   : 
 Function: Transition object from one state to another
 Returns : 
 Args    : object state id, 
           nickname to transition to

=cut

sub changeState {
  my ($self,$obj_state_id, $to_nick) = @_;

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

=head2 get_inTransit_toMonitor

 Title   : get_inTransit_toMonitor
 Usage   : my @ids = $controldb->get_inTransit_toMonitor();
 Function: Transition object from one state to another
 Returns : all object states which are in transition
 Args    : none

=cut

sub get_inTransit_toMonitor {
  my ($self) = @_;
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

=head2 get_noFinal_nonTransit

 Title   : get_noFinal_nonTransit
 Usage   : my @ids = $controldb->get_noFinal_nonTransit();
 Function: Get where the transit_module has to be started
 Returns : ids where the transit_module has to be started
 Args    : none


=cut

sub get_nonFinal_nonTransit {
  my ($self) = @_;
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
=head2 get_byIds_age

 Title   : get_byIds_age
 Usage   : my @cols = $controldb->get_byIds_age($minutes,\@ids, \@cols); 
 Function: used by transition monitors and aware transition modules
           to get what they want from the ControlDB
 Returns : list of cols that match criteria
 Args    : minutes        - minutes old
           idListRef      - array ref of ids
           columnsListRef - array ref of columns to obtain

=cut

sub get_byIds_age {
  my ($self,$minutes,$idListRef,$columnsListRef) = @_;
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

=head2 object_toTransit

 Title   : object_toTransit
 Usage   : $controldb->object_toTransit(\@objectidlist, 1);
 Function: Set Objects transition state 
 Returns : none
 Args    : object stateid list - array ref of object state ids
           resetTransit        - boolean state to set transition

=cut

sub object_toTransit {
  my ($self,$obj_state_id_list,$resetTransit) = @_;

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

=head2 object_reset_transit

 Title   : object_reset_transit
 Usage   :
 Function: very dirty, should seperate out the transit handling function
           here you can reset transits you have set before. Do so, 
           when executing the TransitionModule start fails. 
           If it starts, then dont!
 Returns : none
 Args    : object stateid list - array ref of objects to reset to 
           transition to false

=cut

sub object_reset_transit {
  my ($self,$obj_state_id_list) = @_;
  $self->object_toTransit( $obj_state_id_list, 'false' );
}


# put an object into the pipeline
# missing, a check so that the same object goes not in again.

=head2 submit

 Title   : submit
 Usage   : $controldb->submit($objectid, $class,$state);
 Function: inserts an object in the pipeline - checks for duplicates
 Returns : none
 Args    : objectid,
           object class id,
           state nickname

=cut

sub submit {
  my ($self,$object_id,$class,$state) = @_;
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

=head2 db

 Title   : db
 Usage   : $obj->db($newval)
 Function: 
 Example : 
 Returns : value of db
 Args    : newvalue (optional)


=cut

sub db {
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'_db'} = $value;
    }
    return $obj->{'_db'};
}

# convenience function to get rid of %pipeConf once
# maybe have params in the database?
=head2 get_param

 Title   : get_param
 Usage   :
 Function: convenience function to get rid of %pipeConf once
           maybe have params in the database?
 Returns : config options
 Args    : paramname


=cut

sub get_param {
  my ($self,$paramname) = @_;

  return $::pipeConf{$paramname};
}

=head2 log_error

 Title   : log_error
 Usage   : $controldb->log_error($objectstateid, $msg);
 Function: logs a message in the database
 Returns : none
 Args    : object state id, 
           message

=cut

sub log_error {
  my ($self,$obj_state_id,$msg) = @_;

  my $sth = $self->db->prepare( q{
    INSERT into object_error( object_state_id, error_mesg )
    VALUES ( ?, ? )
  } );
  $sth->execute( $obj_state_id, $msg );
}

1;
