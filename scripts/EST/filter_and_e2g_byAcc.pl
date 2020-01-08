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

  filter_and_e2g_byAcc.pl

=head1 SYNOPSIS
 
  filter_and_e2g_byAcc.pl
  Runs RunnableDB::FilterESTs_and_E2G over a chromosomal region and writes output to
  est database specified in ESTConf.pm

=head1 DESCRIPTION


=head1 OPTIONS

  -input_id in form chrname.start-end
  -acc
  all other options set in ESTConf.pm

=cut

use warnings ;
use strict;
use Getopt::Long qw(:config no_ignore_case);
use Bio::EnsEMBL::Pipeline::ESTConf qw (
					EST_RUNNER
					EST_FILTER_RUNNABLE
					EST_REFDBNAME
					EST_REFDBUSER
					EST_REFDBHOST
					EST_DBNAME
					EST_DBUSER
					EST_DBHOST
				       );

my $runner;
my $runnable;
my $refdbname;
my $refdbuser;
my $refdbhost;
my $estdbname;
my $estdbuser;
my $estdbhost;
my $input_id;
my $acc;

&get_variables();

my $command = "$runner -runnable $runnable -input_id $input_id -acc $acc -write";
#my $command = "$runner -runnable $runnable -input_id $input_id -acc $acc";

#print STDERR "command is $command (no write)\n";
print STDERR "command is $command\n";

my $output = `$command`;

# $output will contain anything that was printed to STDOUT during the run
if ($output ne ''){
  print "\noutput from FilterESTs_and_E2G: $output\n";
}

=head2 get_variables

  Title   : get_variables
  Usage   : get_variables
  Function: initialiases global variables according to input parameters and contents of EST_conf.pl 
            If required parameters are not provided, prints usgae statement and exits script.
  Returns : none - uses globals
  Args    : none - uses globals

=cut

sub get_variables {
  GetOptions( 
	      'input_id:s'      => \$input_id,
	      'acc:s'           => \$acc,
	     );

  $runner     = $EST_RUNNER;
  $runnable   = $EST_FILTER_RUNNABLE;
  $refdbname  = $EST_REFDBNAME;
  $refdbuser  = $EST_REFDBUSER;
  $refdbhost  = $EST_REFDBHOST;
  $estdbname  = $EST_DBNAME;
  $estdbuser  = $EST_DBUSER;
  $estdbhost  = $EST_DBHOST;

  if(!(defined $refdbhost && defined $refdbname & defined $refdbuser &&
       defined $estdbhost && defined $estdbname && defined $estdbuser &&
       defined $runner    && defined $runnable  &&
       defined $input_id)){
    print "Usage: filter_and_e2g.pl -input_id -acc\n" .
      "Additional options to be set in EST_conf.pl: runner, est_genome_runnable, refdbname, refdbuser, refdbhost, estdbname, estdbuser, estdbhost\n";
    exit (1);
  }
}
