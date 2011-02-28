#!/usr/local/ensembl/bin/perl -w

use strict;
use warnings;

use Data::Dumper;
use DBI qw( :sql_types );
use Getopt::Long qw( :config no_ignore_case );

my %master_tables = ( 'attrib_type'     => 1,
                      'external_db'     => 1,
                      'misc_set'        => 1,
                      'unmapped_reason' => 1 );
my @tables;

# Master database location:
my ( $mhost, $mport ) = ( undef, '3306' );
my ( $muser, $mpass );
my $mdbname = 'ensembl_production';

# User database location (default values):
my ( $host, $port ) = ( undef, '3306' );
my ( $user, $pass );
my $dbpattern;

my $drop_bak = 0;
my $verbose  = 0;

# Do command line parsing.
if ( !GetOptions( 'mhost|mh=s'     => \$mhost,
                  'mport|mP=i'     => \$mport,
                  'muser|mu=s'     => \$muser,
                  'mpass|mp=s'     => \$mpass,
                  'mdatabase|md=s' => \$mdbname,
                  'host|h=s'       => \$host,
                  'port|P=i'       => \$port,
                  'user|u=s'       => \$user,
                  'pass|p=s'       => \$pass,
                  'database|d=s'   => \$dbpattern,
                  'table|t=s'      => \@tables,
                  'drop|D!'        => \$drop_bak,
                  'verbose|v!'     => \$verbose )
     || !(    defined($host)
           && defined($user)
           && defined($pass)
           && defined($dbpattern)
           && defined($mhost)
           && defined($muser) ) )
{
  my $indent = ' ' x length($0);
  print <<USAGE_END;
This script copies the tables 'attrib_type', 'external_db', 'misc_set'
and 'unmapped_reason' from the production database into a user-defined
database.

Usage:

  $0 -h host [-P port] \\
  $indent -u user [-p password] -d database \\
  $indent -mh host [-mP port] \\
  $indent -mu user [-mp password] [-md database] \\
  $indent [-t table] [-t table] [-t ...] \\
  $indent [-D] [-v]

  -h / --host       User database server host
  -P / --port       User database server port (optional, default is 3306)

  -u / --user       User username (must have write-access)
  -p / --pass       User password

  -d / --database   User database name or pattern (Perl regular expression)
                    e.g. --database "(rnaseq|vega)_62"

  -mh / --mhost     Production database server host
  -mP / --mport     Production database server port
                    (optional, default is 3306)

  -mu / --muser     Production database username (no write-access required)
  -mp / --mpass     Production database password
                    (optional, default is undefined)

  -md / --mdatabase Production database name
                    (optional, default is 'ensembl_production')

  -t / --table      A specific table to update, may occur several times

  -D / --drop       Drop backup tables if they exists
  -v / --verbose    Be verbose, display every SQL statement as they
                    are executed

USAGE_END

  die("Need the following options: -h -u -p -d -mh -mu\n");

} ## end if ( !GetOptions( 'mhost|mh=s'...))

if (@tables) {
  foreach my $table (@tables) {
    if ( !exists( $master_tables{$table} ) ) {
      die( sprintf( "Invalid table specified: '%s'\n", $table ) );
    }
  }
} else {
  @tables = keys(%master_tables);
}

# Fetch all data from the master database.
my %data;
{
  my $dsn = sprintf( 'DBI:mysql:host=%s;port=%d;database=%s',
                     $mhost, $mport, $mdbname );
  my $dbh = DBI->connect( $dsn, $muser, $mpass,
                          { 'PrintError' => 1, 'RaiseError' => 1 } );

  foreach my $table (@tables) {
    my $sth = $dbh->prepare(
                   sprintf( 'SELECT * FROM %s',
                     $dbh->quote_identifier( undef, $mdbname, $table ) )
    );

    $sth->execute();

    while ( my $row = $sth->fetchrow_arrayref() ) {
      push( @{ $data{$table} }, [ @{$row} ] );
    }
  }

  $dbh->disconnect();
}

# Put all data into the specified database.
{
  my $dsn = sprintf( 'DBI:mysql:host=%s;port=%d', $host, $port );
  my $dbh = DBI->connect( $dsn, $user, $pass,
                          { 'PrintError' => 1, 'RaiseError' => 1 } );

  my $sth = $dbh->prepare('SHOW DATABASES');

  $sth->execute();

  my $dbname;
  $sth->bind_col( 1, \$dbname );

  while ( $sth->fetch() ) {
    if ( $dbname !~ /$dbpattern/ ) { next }

    print( '=' x 80, "\n" );
    printf( "\t%s\n", $dbname );
    print( '=' x 80, "\n" );

    foreach my $table ( keys(%data) ) {
      printf( "==> Inserting into %s\n", $table );

      my $full_table_name =
        $dbh->quote_identifier( undef, $dbname, $table );
      my $full_table_name_bak =
        $dbh->quote_identifier( undef, $dbname, $table . '_bak' );
      my $key_name = $table . '_id';

      # Drop backup table if it exists and if asked to do so.
      if ($drop_bak) {
        $dbh->do(
           sprintf( 'DROP TABLE IF EXISTS %s', $full_table_name_bak ) );
      }

      # Make a backup of any existing data.
      $dbh->do( sprintf( 'CREATE TABLE %s LIKE %s',
                         $full_table_name_bak, $full_table_name ) );
      $dbh->do( sprintf( 'INSERT INTO %s SELECT * FROM %s',
                         $full_table_name_bak, $full_table_name ) );

      # Truncate (empty) the table before inserting new data into it.
      $dbh->do( sprintf( 'TRUNCATE TABLE %s', $full_table_name ) );

      # Get column information.
      my $colinfo_sth =
        $dbh->column_info( undef, $dbname, $table, '%' );
      my $colinfo =
        $colinfo_sth->fetchall_hashref( ['ORDINAL_POSITION'] );

      my $numcols = scalar( keys( %{$colinfo} ) );

      # For each row read from the master table,
      # issue an INSERT statement.
      foreach my $row ( @{ $data{$table} } ) {
        my $insert_statement = sprintf(
          'INSERT INTO %s (%s) VALUES (%s)',
          $full_table_name,
          join( ', ',
                map { $colinfo->{$_}{'COLUMN_NAME'} } 1 .. $numcols ),
          join(
            ', ',
            map {
              $dbh->quote( $row->[ $_ - 1 ],
                           $colinfo->{$_}{'DATA_TYPE'} )
              } ( 1 .. $numcols ) ) );

        if ($verbose) {
          printf( "EXECUTING: %s\n", $insert_statement );
        }
        $dbh->do($insert_statement);
      }

      {
        my $statement = sprintf( 'SELECT %s '
                                   . 'FROM %s '
                                   . 'LEFT JOIN %s t USING (%s) '
                                   . 'WHERE t.%s IS NULL '
                                   . 'ORDER BY %s',
                                 $key_name,
                                 $full_table_name,
                                 $full_table_name_bak,
                                 $key_name,
                                 $key_name,
                                 $key_name );

        my $sth2 = $dbh->prepare($statement);

        if ($verbose) {
          printf( "EXECUTING: %s\n", $statement );
        }
        $sth2->execute();

        my $key;
        $sth2->bind_col( 1, \$key );

        my @keys;
        while ( $sth2->fetch() ) {
          push( @keys, $key );
        }

        if (@keys) {
          print("New data inserted:\n");
          printf( "SELECT * FROM %s WHERE %s_id IN (%s);\n",
                  $table, $table, join( ',', @keys ) );
          print("\n");
        }
      }
      {
        my $statement = sprintf( 'SELECT %s '
                                   . 'FROM %s '
                                   . 'LEFT JOIN %s t USING (%s) '
                                   . 'WHERE t.%s IS NULL '
                                   . 'ORDER BY %s',
                                 $key_name,
                                 $full_table_name_bak,
                                 $full_table_name,
                                 $key_name,
                                 $key_name,
                                 $key_name );

        my $sth2 = $dbh->prepare($statement);

        if ($verbose) {
          printf( "EXECUTING: %s\n", $statement );
        }
        $sth2->execute();

        my $key;
        $sth2->bind_col( 1, \$key );

        my @keys;
        while ( $sth2->fetch() ) {
          push( @keys, $key );
        }

        if (@keys) {
          print( '-' x 40, "\n" );
          print("!! Old data deleted:\n");
          printf( "SELECT * FROM %s WHERE %s_id IN (%s);\n",
                  $table . '_bak',
                  $table, join( ',', @keys ) );
          print( '-' x 40, "\n" );
          print("\n");
        }
      }

    } ## end foreach my $table ( keys(%data...))
    print("\n");
    print("Remember to drop the backup tables!\n\n");
  } ## end while ( $sth->fetch() )

  $dbh->disconnect();
}
