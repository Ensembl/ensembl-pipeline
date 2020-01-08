#!/usr/bin/env perl


# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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

run_EST_RunnableDB

=head1 SYNOPSIS
 
  test_RunnableDB

=head1 DESCRIPTION


=head1 OPTIONS

    -dbhost    host name for database (gets put as host= in locator)

    -dbport    For RDBs, what port to connect to (port= in locator)

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
use Bio::EnsEMBL::Pipeline::Config::cDNAs_ESTs::EST_GeneBuilder_Conf qw (
						     EST_REFDBNAME
						     EST_REFDBUSER
						     EST_REFDBPASS
						     EST_REFDBHOST
						     EST_GENE_DBHOST
						     EST_GENE_DBUSER
						     EST_GENE_DBNAME
						     EST_GENE_DBPASS
                                                      );


use Bio::EnsEMBL::DBSQL::DBAdaptor;
#use Bio::EnsEMBL::DBLoader;
my $dbtype = 'rdb';
my $port   = '3306';
my $dbname = $EST_REFDBNAME;
my $dbuser = $EST_REFDBUSER;
my $dbpass = $EST_REFDBPASS;
my $host   = $EST_REFDBHOST;


my $runnable;
my $input_id;
my $write  = 0;
my $check  = 0;
my $pepfile;
my $acc;
my $analysis;
my $use_label;

# can override db options on command line
GetOptions( 
	     'input_id:s'    => \$input_id,
	     'runnable:s'    => \$runnable,
             'analysis:s'    => \$analysis,
	     'write'         => \$write,
             'check'         => \$check,
             'dbname|db|D:s'      => \$dbname,
             'host|dbhost|h:s'      => \$host,
             'user|dbuser|u:s'      => \$dbuser,
             'pass|dbpass|p:s'      => \$dbpass,
	     'use_label'     => \$use_label,
           );

$| = 1;

die "No runnable entered" unless defined ($runnable);
(my $file = $runnable) =~ s/::/\//g;
require "$file.pm";

if ($check) {
   exit(0);
}

print STDERR "args: $host : $dbuser : $dbpass : $dbname\n";

die "No input id entered" unless defined ($input_id);

my $analysis_obj;
my $runobj;
{
 my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host             => $host,
    -user             => $dbuser,
    -dbname           => $dbname,
    -pass             => $dbpass,
    );

 eval{
   $analysis_obj = $db->get_AnalysisAdaptor->fetch_by_logic_name($analysis);
 };
 unless( $analysis_obj ){
   my $output_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
						     -host             => $EST_GENE_DBHOST,
						     -user             => $EST_GENE_DBUSER,
						     -dbname           => $EST_GENE_DBNAME,
						     -pass             => $EST_GENE_DBPASS,
						    );

   $analysis_obj = $output_db->get_AnalysisAdaptor->fetch_by_logic_name($analysis);
 }

  $runobj = "$runnable"->new(-db       => $db,
			      -input_id => $input_id,
			      -analysis => $analysis_obj,
        			     );

} # the db should go out of scope here?

if ($use_label){
 $runobj->_label($input_id);
}

$runobj->fetch_input;
$runobj->run;

 my @out = $runobj->output;

if ($write) {
  $runobj->write_output;
}

