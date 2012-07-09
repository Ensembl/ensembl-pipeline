#!/usr/bin/env perl
# $Id: analysis_setup.pl,v 1.11 2012-07-09 13:56:16 ak4 Exp $

# POD documentation - main docs before the code

=pod

=head1 NAME

Analysis setup

=head1 SYNOPSIS
 
Script for writing an analysis table and generate the config needed to
generate an analysis table.

=head1 DESCRIPTION

This script will both take a config file and write to an analysis table
and take an analysis table and write a config file, it will work with
both straight core dbs and core dbs with pipeline tables


=head1 OPTIONS

Database options  

  -dbhost   host name for database (gets put as host= in locator)
  -dbport   what port to connect to (port= in locator) (optional)
  -dbname   what name to connect to (dbname= in locator)
  -dbuser   what username to connect as (dbuser= in locator)
  -dbpass   what password to use (dbpass= in locator)

Other options

  -read     Read a config file and write to the database.

  -write    Read the analysis table and write a config file.

  -file     Path to the file you either want to read from or write too.

  -pipeline_db  Indicates your database contains pipeline tables and the
                pipeline dbadaptor should be used, this is on by default
                but can be turned off with the -nopipeline_db option.

  -update   Update any existing analyses with the data in the config
            file.

  -help     Display help.

=head1 EXAMPLES

This will generate a config file based on the analysis table of the 
database it is pointed at

  ./analysis setup -dbhost ecs1b -dbuser ensadmin -dbpass **** \
                   -dbname my_pipeline_db -dbport 3306 -write \
                   -file analysis.conf 

This will fill in an analysis table based on the config file

  ./analysis setup -dbhost ecs1b -dbuser ensadmin -dbpass **** \
                   -dbname my_pipeline_db -dbport 3306 -read \
                   -file analysis.conf 

This will fill in an analysis table but leave the associated pipeline
tables alone

  ./analysis setup -dbhost ecs1b -dbuser ensadmin -dbpass **** \
                   -dbname my_pipeline_db -dbport 3306 -read \
                   -file analysis.conf -nopipeline_db


This is what a conf entry should like like.  The header should be the
logic_name, the rest of the columns should lie in key value pairs.

  [RepeatMask]
  db=repbase
  db_version=020713
  db_file=repbase
  program=RepeatMasker
  program_version=1
  program_file=RepeatMasker
  parameters=-low, -lib, /path/to/file/briggsae.lib 
  module=RepeatMasker
  module_version=1
  gff_source=RepeatMasker
  gff_feature=Repeat
  input_id_type=CONTIG

There is an example file in this directory called example_analysis.conf.

=cut

use strict;
use warnings;

use AnalysisCreation;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 3306;
my $dbname;
my $read;
my $write;
my $file;
my $help;
my $pipeline = 1;
my $update;

GetOptions( 'dbhost=s'     => \$dbhost,
            'dbname=s'     => \$dbname,
            'dbuser=s'     => \$dbuser,
            'dbpass=s'     => \$dbpass,
            'dbport=s'     => \$dbport,
            'read!'        => \$read,
            'write!'       => \$write,
            'pipeline_db!' => \$pipeline,
            'update!'      => \$update,
            'file=s'       => \$file,
            'h|help!'      => \$help,
  ) or
  usage();



if ( !($dbhost) || !($dbuser) || !($dbname) ) {
  print STDERR
    "need to pass in database arguments for script to work\n";
  print STDERR
    "-dbhost $dbhost -dbuser $dbuser -dbpass $dbpass -dbname" .
    " $dbname -dbport $dbport\n";
  $help = 1;
}

if ( ( !($read) && !($write) ) || ( $read && $write ) ) {
  print STDERR "you need to define either read or write on the " .
    "commandline but you shouldn't define both\n";
  $help = 1;
}

if ( !$file ) {
  print STDERR
    "You need to pass a file name which either represents a " .
    "analysis config file to be read and stored in the database or " .
    "written based on the analysis objects in the database\n";
  $help = 1;
}

if ($help) { usage() }

my $db;
if($pipeline){
  $db  = new Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
    (-host => $dbhost,
     -user => $dbuser,
     -pass => $dbpass,
     -dbname => $dbname,
     -port => $dbport,
    );
}else{
  $db  = new Bio::EnsEMBL::DBSQL::DBAdaptor
    (-host => $dbhost,
     -user => $dbuser,
     -pass => $dbpass,
     -dbname => $dbname,
     -port => $dbport,
    );
}



if($read){
  my @analyses = @{&parse_files($file)};
  &write_into_db($db, \@analyses, $update);
}

if($write){
  my $analyses = &read_db($db);
  $analyses = [sort {$a->dbID <=> $b->dbID} @$analyses];
  &write_file($file, $analyses);
}


sub usage {
  exec( 'perldoc', $0 );
  exit;
}

