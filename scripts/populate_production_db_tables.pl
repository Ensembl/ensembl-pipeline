#!/usr/local/ensembl/bin/perl -w

use strict;
use warnings;

use Data::Dumper;
use DBI;
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
my $dbname;

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
                  'database|d=s'   => \$dbname,
                  'table|t=s'      => \@tables )
     || !(    defined($host)
           && defined($user)
           && defined($dbname)
           && defined($mhost)
           && defined($muser) ) )
{
  my $indent = ' ' x length($0);
  print <<USAGE_END;
Usage:

  $0 -h host [-P port] \\
  $indent -u user [-p password] -d database \\
  $indent -mh host [-mP port] \\
  $indent -mu user [-mp password] [-md database] \\
  $indent [-t table] [-t table] [-t ...]

  -h / --host       User database server host
  -P / --port       User database server port (optional, default is 3306)
  -u / --user       User username (must have write-access)
  -p / --pass       User password (optional, default is undefined)
  -d / --database   User database name

  -mh / --mhost     Production database server host
  -mP / --mport     Production database server port
                    (optional, default is 3306)
  -mu / --muser     Production database username (no write-access required)
  -mp / --mpass     Production database password
                    (optional, default is undefined)
  -md / --mdatabase Production database name
                    (optional, default is 'ensembl_production')

  -t / --table      A specific table to update, may occur several times

USAGE_END

  die("Need the following options: -h -u -d -mh -mu\n");

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
  my $dsn = sprintf( 'DBI:mysql:host=%s;port=%d;database=%s',
                     $host, $port, $dbname );
  my $dbh = DBI->connect( $dsn, $user, $pass,
                          { 'PrintError' => 1, 'RaiseError' => 1 } );

  foreach my $table ( keys(%data) ) {
    printf( "Inserting into %s.%s\n", $dbname, $table );

    # Truncate (empty) the table before inserting new data into it.
    $dbh->do(sprintf( 'TRUNCATE TABLE %s',
                      $dbh->quote_identifier( undef, $dbname, $table ) )
    );

    # Get column information.
    my $colinfo_sth = $dbh->column_info( undef, $dbname, $table, '%' );
    my $colinfo =
      $colinfo_sth->fetchall_hashref( ['ORDINAL_POSITION'] );

    my $numcols = scalar( keys( %{$colinfo} ) );

    # For each row read from the master table,
    # issue an INSERT statement.
    foreach my $row ( @{ $data{$table} } ) {
      my $insert_statement = sprintf(
        'INSERT INTO %s (%s) VALUES (%s)',
        $dbh->quote_identifier( undef, $dbname, $table ),
        join( ', ',
              map { $colinfo->{$_}{'COLUMN_NAME'} } 1 .. $numcols ),
        join(
          ', ',
          map {
            $dbh->quote( $row->[ $_ - 1 ], $colinfo->{$_}{'DATA_TYPE'} )
            } ( 1 .. $numcols ) ) );

      $dbh->do($insert_statement);
    }

  } ## end foreach my $table ( keys(%data...))

  $dbh->disconnect();
}
