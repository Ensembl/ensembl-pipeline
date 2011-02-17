#!/usr/local/ensembl/bin/perl -w

use strict;
use warnings;

use Data::Dumper;
use DBI;
use Getopt::Long qw( :config no_ignore_case );

my @tables =
  ( 'attrib_type', 'external_db', 'misc_set', 'unmapped_reason' );
my $user_table = '';


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
                  'mdatabase|md=s' => \$mdbname,
                  'host|h=s'       => \$host,
                  'port|P=i'       => \$port,
                  'user|u=s'       => \$user,
                  'pass|p=s'       => \$pass,
                  'database|d=s'   => \$dbname,
                  'table|t=s'     => \$user_table)
     || !(    defined($host)
           && defined($user)
           && defined($dbname)
           && defined($mhost)
           && defined($muser)
           && defined($mdbname) ) )
{
  my $indent = ' ' x length($0);
  print <<USAGE_END;
Usage:

  $0 -h host [-P port] \\
  $indent -u user [-p password] \\
  $indent -d database

  -h / --host       User database server host
  -P / --port       User database server port (optional, default is 3306)
  -u / --user       User username (must have write-access)
  -p / --pass       User password (optional, default is undefined)
  -d / --database   User database name

  -mh / --mhost       ensembl_production database server host
  -mP / --mport       ensembl_production database server port (optional, default is 3306)
  -mu / --muser       ensembl_production username (no write-access required)
  -md / --mdatabase   ensembl_production database name

  -t / --table       specify a table to update

USAGE_END
  die;
} ## end if ( !GetOptions( 'mhost|mh=s'...))

if ($user_table) {
  @tables = $user_table;
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
