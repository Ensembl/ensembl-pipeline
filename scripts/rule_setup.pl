#!/usr/bin/env perl
# $Source: /tmp/ENSCOPY-ENSEMBL-PIPELINE/scripts/rule_setup.pl,v $
# $Revision: 1.7 $

=head1 NAME

rule_setup.pl

=head1 SYNOPSIS

Script for writing the rule tables given a configuration file, or the
other way around.

=head1 DESCRIPTION

This script will take a configuration file and write to the rule tables
or take the rule tables and write a configuration file.  The database
you point at must alraedy contain the pipeline tables as defined in
ensembl-pipeline/sql/table.sql

=head1 OPTIONS

Database options:

  -dbhost   Host name for database.
  -dbport   Port to connect to (optional).
  -dbname   Database to connect to.
  -dbuser   Username to connect as.
  -dbpass   Password to use.

Other options:

  -read     Read a configuration file and write to the database.
  -write    Read the rule table and write a configuration file.
  -file     File to read from or write too.

  -help     Displays help text.

=head1 EXAMPLES

To generate a configuration file based on the rule table of the database
it is pointed at:

  ./rule_setup.pl -dbhost ecs1b -dbuser ensadmin -dbpass **** \
                  -dbname my_pipeline_db -dbport 3306 -write \
                  -file rule.conf

To fill in an rule table based on the configuration file passed in:

  ./rule_setup.pl -dbhost ecs1b -dbuser ensadmin -dbpass **** \
                  -dbname my_pipeline_db -dbport 3306 -read \
                  -file rule.conf

This is what a configuration entry should look like.  The header for
each section should be the goal analysis logic name and the conditions
should be key-value pairs:

  [RepeatMask]
  condition=SubmitContig

  [Pmatch]
  condition=SubmitChromosome

  [Pmatch_Wait]
  condition=Pmatch

  [BestPmatch]
  condition=Pmatch_Wait
  condition=SubmitGenome

Any line starting with a hash ('#') is treated as a comment.

See example_rule.conf for an example of a configuration file.

=cut

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case);

use RuleCreation;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;

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

if ( !GetOptions( 'host|dbhost|h:s'     => \$dbhost,
                  'dbname|db|D:s'     => \$dbname,
                  'user|dbuser|u:s'     => \$dbuser,
                  'pass|dbpass|p:s'     => \$dbpass,
                  'port|dbport|P:s'     => \$dbport,
                  'read|insert!' => \$read,
                  'write!'       => \$write,
                  'file=s'       => \$file,
                  'h|help!'      => \$help, ) ||
     $help )
{
  usage();
}

if ( !($dbhost) || !($dbuser) || !($dbname) ) {
  print STDERR "need to pass in database arguments for script to " .
    "work\n";
  print STDERR "-dbhost $dbhost -dbuser $dbuser -dbpass $dbpass " .
    "-dbname  $dbname -dbport $dbport\n";
  usage();
}

if ( ( !($read) && !($write) ) || ( $read && $write ) ) {
  print STDERR "you need to define either read or write on the " .
    "commandline but you shouldn't define both\n";
  usage();
}

if ( !$file ) {
  print STDERR "You need to pass a file name which either represents " .
    "a rule config file to be read and stored in the database or " .
    "written based on the rule objects in the database\n";
  usage();
}

my $db =
  new Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor( -host   => $dbhost,
                                                -user   => $dbuser,
                                                -pass   => $dbpass,
                                                -dbname => $dbname,
                                                -port   => $dbport, );

if ($read) {
  my $rule_hash = parse_files($file);
  my @rules = @{ create_rules( $db, $rule_hash ) };
  write_into_db( $db, \@rules );
}
else {
  my $analyses = read_db($db);
  $analyses = [
    sort {
      $a->goalAnalysis->dbID() <=> $b->goalAnalysis->dbID()
    } @{$analyses} ];
  write_file( $file, $analyses );
}
