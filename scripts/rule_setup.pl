#!/usr/local/ensembl/bin/perl -w

# POD documentation - main docs before the code

=pod

=head1 NAME

  rule setup

=head1 SYNOPSIS
 
  a script for writing the rule tables and generate the config needed
  to generate the rule tables rule_goal and rule_condition

=head1 DESCRIPTION

  this script will both take a config file and write to the rule tables
  and take the rule tables and write a config file. The database
  you point at must alraedy contain the pipeline tables as defined
  in ensembl-pipeline/sql/table.sql

=head1 OPTIONS

     Database options  

    -dbhost      host name for database (gets put as host= in locator)
    -dbport      For RDBs, what port to connect to (port= in locator)
    -dbname    For RDBs, what name to connect to (dbname= in locator)
    -dbuser    For RDBs, what username to connect as (dbuser= in locator)
    -dbpass    For RDBs, what password to use (dbpass= in locator)
    
     Other options

     -read this indicates to the script you want to read a config file
           and write to the database
     -write this indicates to the script you want to read the rule 
            table and write a config file
     -help prints out the perl docs
  
=head1 EXAMPLES

this will generate a config file based on the rule table of the 
database its is pointed at

./rule setup -dbhost ecs1b -dbuser ensadmin -dbpass **** 
                 -dbname my_pipeline_db -dbport 3306 -write 
                 -file rule.conf 

this will fill in an rule table based on the config file passed in
  
./rule setup -dbhost ecs1b -dbuser ensadmin -dbpass **** 
                 -dbname my_pipeline_db -dbport 3306 -read 
                 -file rule.conf 

this is what a conf entry should like like. The header should be the 
goal analysis logic_name , the conditions should lie in key value
pairs

[RepeatMask]
condition=SubmitContig

[Pmatch]
condition=SubmitChromosome

[Pmatch_Wait]
condition=Pmatch

[BestPmatch]
condition=Pmatch_Wait
condition=SubmitGenome

comments can be left if the line is started with a # symbol

see example_rule.conf for an example of a config file
=cut

use strict;
use RuleCreation;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
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


&GetOptions( 
	    'dbhost=s'      => \$dbhost,
	    'dbname=s'      => \$dbname,
	    'dbuser=s'      => \$dbuser,
	    'dbpass=s'      => \$dbpass,
	    'dbport=s'      => \$dbport,
	    'read!'         => \$read,
	    'write!'        => \$write,
	    'file=s'        => \$file,
	    'h|help!'       => \$help,
	    ) or useage();


if(!($dbhost) || !($dbuser) || !($dbname)){
  print STDERR "need to pass in database arguments for script to ".
    "work\n";
  print STDERR "-dbhost $dbhost -dbuser $dbuser -dbpass $dbpass ".
    "-dbname  $dbname -dbport $dbport\n";
  $help = 1;
}

if((!($read) && !($write)) || ($read && $write)){
  print STDERR "you need to define either read or write on the ".
    "commandline but you shouldn't define both\n";
  $help = 1;
}

if(!$file){
  print STDERR "You need to pass a file name which either represents ".
    "a rule config file to be read and stored in the database or ".
      "written based on the rule objects in the database\n";
  $help = 1;
}

if($help){
  useage();
}

my $db  = new Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
  (-host => $dbhost,
   -user => $dbuser,
   -pass => $dbpass,
   -dbname => $dbname,
   -port => $dbport,
  );

if($read){
  my $rule_hash = parse_files($file);
  my @rules = @{create_rules($db, $rule_hash)};
  &write_into_db($db, \@rules);
}

if($write){
  my $analyses = &read_db($db);
  $analyses = [sort {$a->goalAnalysis->dbID <=> $b->goalAnalysis->dbID} @{$analyses}];
  &write_file($file, $analyses);
}


sub useage{
  exec('perldoc', $0);
  exit;
}

