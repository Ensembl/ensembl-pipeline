#!/usr/bin/env perl


# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


=head1 NAME

analysis_setup.pl

=head1 SYNOPSIS

Script for writing an analysis table given a configuration file, or the
other way around.

=head1 DESCRIPTION

This script will take a configuration file and write to an analysis table
or take an analysis table and write a configuration file.

=head1 OPTIONS

Database options:

  -dbhost    Host name for database.
  -dbport    Port to connect to (optional).
  -dbname    Database to connect to.
  -dbuser    Username to connect as.
  -dbpass    Password to use.

Other options:

  -read     Read a configuration file and write to the database.
  -write    Read the rule table and write a configuration file.
  -file     File to read from or write too.


  -pipeline_db  Your database contains pipeline tables and the pipeline
                database adaptor should be used, this is on by default
                but can be turned off with the -nopipeline_db option.

  -update   Update any existing analyses with the data in the
            configuration file.

  -help     Displays help text.

=head1 EXAMPLES

To generate a configuration file based on the analysis table of the
database it is pointed at:

  ./analysis setup -dbhost ecs1b -dbuser ensadmin -dbpass **** \
                   -dbname my_pipeline_db -dbport 3306 -write \
                   -file analysis.conf

To fill in an analysis table based on the configuration file:

  ./analysis setup -dbhost ecs1b -dbuser ensadmin -dbpass **** \
                   -dbname my_pipeline_db -dbport 3306 -read \
                   -file analysis.conf

To fill in an analysis table but leave the associated pipeline tables
alone:

  ./analysis setup -dbhost ecs1b -dbuser ensadmin -dbpass **** \
                   -dbname my_pipeline_db -dbport 3306 -read \
                   -file analysis.conf -nopipeline_db


This is what a conf entry should like like.  The header for each section
should be the analysis logic name, the rest is key-value pairs:

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

There is an example file in this directory called "example_analysis.conf".

=cut

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case);

use Bio::EnsEMBL::Pipeline::AnalysisCreation;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

sub usage {
  exec( 'perldoc', $0 );
  exit;
}

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

if ( !GetOptions( 'host|dbhost|h:s'     => \$dbhost,
                  'dbname|db|D:s'     => \$dbname,
                  'user|dbuser|u:s'     => \$dbuser,
                  'pass|dbpass|p:s'     => \$dbpass,
                  'port|dbport|P:s'     => \$dbport,
                  'read!'        => \$read,
                  'write!'       => \$write,
                  'pipeline_db!' => \$pipeline,
                  'update!'      => \$update,
                  'file=s'       => \$file,
                  'help!'      => \$help, ) ||
     $help )
{
  usage();
}

if ( !($dbhost) || !($dbuser) || !($dbname) ) {
  print STDERR
    "need to pass in database arguments for script to work\n";
  print STDERR
    "-dbhost $dbhost -dbuser $dbuser -dbpass $dbpass -dbname" .
    " $dbname -dbport $dbport\n";
  usage();
}

if ( ( !($read) && !($write) ) || ( $read && $write ) ) {
  print STDERR "you need to define either read or write on the " .
    "commandline but you shouldn't define both\n";
  usage();
}

if ( !$file ) {
  print STDERR
    "You need to pass a file name which either represents a " .
    "analysis config file to be read and stored in the database or " .
    "written based on the analysis objects in the database\n";
  usage();
}

my $db;

if ($pipeline) {
  $db =
    new Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor( -host   => $dbhost,
                                                  -user   => $dbuser,
                                                  -pass   => $dbpass,
                                                  -dbname => $dbname,
                                                  -port   => $dbport, );
}
else {
  $db =
    new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $dbhost,
                                        -user   => $dbuser,
                                        -pass   => $dbpass,
                                        -dbname => $dbname,
                                        -port   => $dbport, );
}

if ($read) {
  my @analyses = @{ parse_files($file) };
  write_into_db( $db, \@analyses, $update );
}
else {
  my $analyses = read_db($db);
  $analyses = [ sort { $a->dbID() <=> $b->dbID() } @{$analyses} ];
  write_file( $file, $analyses );
}
