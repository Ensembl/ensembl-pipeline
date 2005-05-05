package Bio::EnsEMBL::Pipeline::DBSQL::PairwiseTrack;

use strict;
use DBI;
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

# While storing homology objects in the compara database it is
# occasionally necessary to avoid redundancy.  This adaptor manages
# a single table in the database that tracks the insertion of
# id pairs.  To use it, simply call the attempt_write_claim method
# with the candidate pair of member ids.  If these havent been
# previously inserted a true return value will be obtained, hence
# it is possible to go ahead and safely store the homology object.  
# If the pair of ids has previously been stored a false value
# will be returned, from which you will know better to attempt to
# store the homology object.
#
# The table this adaptor uses is very basic - a single column.  To 
# create it in any existing database, use this sql:
# CREATE TABLE pairwise_track (unique_id varchar(40) UNIQUE NOT NULL);


sub new {
  my ($class, @args) = @_;

  my $self = bless {}, $class;

  my ($db,
      $dbname,
      $dbhost,
      $dbuser,
      $dbpass,
      $dbport,
      $dbdriver) = &rearrange([qw(DB
			      DBNAME
			      DBHOST
			      DBUSER
			      DBPASS
			      DBPORT
			      DBDRIVER)
			     ], @args);

  if (defined $db) {

    unless ($db->isa("Bio::EnsEMBL::DBSQL::DBAdaptor")) {
      throw("Database adaptor is not a " .
	    "Bio::EnsEMBL::DBSQL::DBAdaptor it is a [$db]");
    }

    $dbdriver = $db->dbc->driver;
    $dbport   = $db->dbc->port;
    $dbuser   = $db->dbc->username;
    $dbpass   = $db->dbc->password;
    $dbhost   = $db->dbc->host;
    $dbname   = $db->dbc->dbname;

  }

  unless (defined $dbname && defined $dbhost && defined $dbuser) {
    throw("Need a minimum of database name [$dbname], host [$dbhost] and " . 
	  "username [$dbuser] to make a database connection.")
  }

  $dbdriver ||= 'mysql';
  $dbport   ||= 3306;

  my $dsn = "DBI:$dbdriver:database=$dbname;host=$dbhost;port=$dbport";

  my $dbh;

  eval{
    $dbh = DBI->connect($dsn,
			$dbuser,
			$dbpass,
			{'RaiseError' => 1});
  };

  if(!$dbh || $@ || !$dbh->ping()) {
    throw("Could not connect to database " . $dbname .
	  " as user " . $dbuser .
	  " using [$dsn] as a locator:\n" . $DBI::errstr);
  }

  $self->_dbh($dbh);

  return $self;
}


sub DESTROY {
  my $self = shift;

  $self->_dbh->disconnect;
}


sub _dbh {
  my $self = shift;

  if (@_) {
    $self->{_dbh} = shift;
  }

  return $self->{_dbh}
}


sub table_exists {
  my $self = shift;

  my $sth = 
    $self->_dbh->prepare("DESCRIBE pairwise_track");

  eval {
    $sth->execute
  };

  return 0 if $@;

  return 1
}


sub attempt_write_claim {
  my ($self, $id1, $id2) = @_;

  unless (defined $id1 && defined $id2){
    throw("Missing ids [$id1][$id2]");
  }

  my $sth = 
    $self->_dbh->prepare("INSERT INTO pairwise_track(unique_id) VALUES (?)");

  my @sorted_ids = sort {$a cmp $b} ($id1, $id2);

  my $pair_id = join '', @sorted_ids;

  return $sth->execute($pair_id)
}


sub previously_inserted {
  my ($self, $id1, $id2) = @_;

  unless (defined $id1 && defined $id2){
    throw("Missing ids [$id1][$id2]");
  }

  my $sth = 
    $self->_dbh->prepare("SELECT count(unique_id) FROM pairwise_track " . 
			 "WHERE unique_id = ?");

  my @sorted_ids = sort {$a cmp $b} ($id1, $id2);

  my $pair_id = join '', @sorted_ids;

  $sth->execute($pair_id);

  my ($count) = $sth->fetchrow_array();

  if ($count > 0) {
    return 1
  } else {
    return 0
  }
}

return 1
