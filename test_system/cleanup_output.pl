#This script's functionality is to make cleaning up databases and 
#directories left behind after a test run easier. 

#To clean up only the DB:
#cleanup_output.pl  -dbhost xxx -dbuser xxx -dbpass xxx -dbport xxx -dbname xxx

#To clean up only the test output directory (and the files it contains):
#cleanup_output.pl  -output_dir xxx

#To clean up only the directory where the unzipped reference data files are:
#cleanup_output.pl -sql_data_dir /path/to/the/data/directory

#For cleaning up any combination of the above, just use the relevant flags and
#provide the required information.


#!/usr/local/ensembl/bin/perl -w

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

&GetOptions(
            'dbhost:s'     => \$host,
            'dbport:n'     => \$port,
            'dbuser:s'     => \$user,
            'dbpass:s'     => \$pass,
            'dbname:s'     => \$dbname,
            'dbdriver:s' => \$driver,
            'sql_data_dir:s' => \$data_dir,
            'output_dir:s' => \$output_dir,
           ) or throw("Can't get opts");


if($dbname){
  my $locator = "DBI:".$driver.":host=".$host.";port=".$port;
  print "locator ".$locator."\n";
  my $db = DBI->connect($locator, $user, $pass, {RaiseError => 1});
  $db->do("Drop database ".$dbname);
}

if($data_dir){
  rmtree($data_dir);
}
if($output_dir){
  rmtree($output_dir);
}
