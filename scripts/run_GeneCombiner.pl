#!/usr/bin/env perl
# $Source: /tmp/ENSCOPY-ENSEMBL-PIPELINE/scripts/run_GeneCombiner.pl,v $
# $Revision: 1.8 $


=head1 NAME

run_GeneCombiner.pl

=head1 SYNOPSIS
 
  test_RunnableDB

=head1 DESCRIPTION


=head1 OPTIONS

    -dbhost      host name for database (gets put as host= in locator)

    -dbport      For RDBs, what port to connect to (port= in locator)

    -dbname    For RDBs, what name to connect to (dbname= in locator)

    -dbuser    For RDBs, what username to connect as (dbuser= in locator)

    -dbpass    For RDBs, what password to use (dbpass= in locator)

    -input_id  The input id for the RunnableDB

    -runnable  The name of the runnable module we want to run

    -analysis  The number of the analysisprocess we want to run
=cut

use warnings ;
use strict;
use Getopt::Long qw(:config no_ignore_case);

# this script connects to the db it is going to write to
use Bio::EnsEMBL::Pipeline::Config::GeneBuild::GeneCombiner qw (
                                        	REF_DBHOST
						REF_DBNAME
						REF_DBUSER
						REF_DBPASS
                                       );


use Bio::EnsEMBL::DBSQL::DBAdaptor;
#use Bio::EnsEMBL::DBLoader;
my $dbtype = 'rdb';
my $port   = '3306';
my $dbname = $REF_DBNAME;
my $dbuser = $REF_DBUSER;
my $dbpass = $REF_DBPASS;
my $host   = $REF_DBHOST;


my $runnable;
my $input_id;
my $write  = 0;
my $check  = 0;
my $params;
my $pepfile;
my $acc;
my $analysis_logic_name;

# can override db options on command line
GetOptions( 
	     'input_id:s'    => \$input_id,
	     'runnable:s'    => \$runnable,
	     'analysis:s'    => \$analysis_logic_name,
             'write'         => \$write,
             'check'         => \$check,
             'parameters:s'  => \$params,
             'dbname|db|D:s'      => \$dbname,
             'host|dbhost|h:s'      => \$host,
             'user|dbuser|u:s'      => \$dbuser,
             'pass|dbpass|p:s'      => \$dbpass,
             'port|dbport|P:n'      => \$port,
	     );

$| = 1;

die "No runnable entered" unless defined ($runnable);
(my $file = $runnable) =~ s/::/\//g;
require "$file.pm";

if ($check) {
   exit(0);
}

print STDERR "args: $host : $dbuser : $dbpass : $dbname\n";

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host             => $host,
    -user             => $dbuser,
    -dbname           => $dbname,
    -pass             => $dbpass,
    -port             => $port,
 );

die "No input id entered" unless defined ($input_id);

my $analysis = $db->get_AnalysisAdaptor->fetch_by_logic_name($analysis_logic_name);

my %hparams;
# eg -parameters param1=value1,param2=value2
if (defined $params){
  foreach my $p(split /,/, $params){
    my @sp = split /=/, $p;
    $sp[0] = '-' . $sp[0];
    $hparams{$sp[0]} =  $sp[1];
  }
}


print STDERR "run_GeneCombiner.pl: input_id: $input_id\n";
my $runobj = "$runnable"->new(-db       => $db,
			      -input_id => $input_id,
			      -analysis => $analysis,
                              %hparams,
			     );

#"$runnable"->dbobj($db);

$runobj->fetch_input;
$runobj->run;
 
my @out = $runobj->output;

if ($write) {
  $runobj->write_output;
}

