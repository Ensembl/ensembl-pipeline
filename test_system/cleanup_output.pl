#!/usr/bin/env perl

# This script helps with cleaning up after a test run.

# To clean up only the DB:
#
#   cleanup_output.pl -dbhost xxx -dbuser xxx \
#       -dbpass xxx -dbport xxx -dbname xxx

# To clean up only the test output directory (and the files it
# contains):
#
#   cleanup_output.pl -output_dir xxx

# To clean up only the directory where the unzipped reference data files
# are:
#
#   cleanup_output.pl -sql_data_dir /path/to/the/data/directory

# For cleaning up any combination of the above, just use the relevant
# flags and provide the required information.

use strict;
use Getopt::Long;
use DBI;
use File::Path;

my $host;
my $user;
my $pass;
my $port = 3306;
my $dbname;
my $data_dir;
my $output_dir;
my $driver = 'mysql';

GetOptions( 'dbhost:s'       => \$host,
            'dbport:n'       => \$port,
            'dbuser:s'       => \$user,
            'dbpass:s'       => \$pass,
            'dbname:s'       => \$dbname,
            'dbdriver:s'     => \$driver,
            'sql_data_dir:s' => \$data_dir,
            'output_dir:s'   => \$output_dir,
  ) or
  throw("Can't parse command line arguemnts");

if ( defined($dbname) ) {
  my $dsn = sprintf( "DBI:%s:host=%s;port=%s", $driver, $host, $port );
  my $db = DBI->connect( $dsn, $user, $pass, { RaiseError => 1 } );
  $db->do("DROP DATABASE $dbname");
  $db->disconnect();
}

if ( defined($data_dir) && -d $data_dir ) {
  rmtree($data_dir);
}

if ( defined($output_dir) && -d $output_dir ) {
  rmtree($output_dir);
}
