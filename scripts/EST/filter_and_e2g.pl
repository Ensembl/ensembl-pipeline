#!/usr/local/bin/perl -w

=head1 NAME

  filter_and_e2g.pl

=head1 SYNOPSIS
 
  filter_and_e2g.pl
  Runs RunnableDB::FilterESTs_and_E2G over a chromosomal region and writes output to
  est database specified in EST_conf.pl

=head1 DESCRIPTION


=head1 OPTIONS

  -input_id in form chrname.start-end
  all other options set in EST_conf.pl

=cut

use strict;
use Getopt::Long;
require "Bio/EnsEMBL/Pipeline/EST_conf.pl";

my $runner;
my $runnable;
my $refdbname;
my $refdbuser;
my $refdbhost;
my $estdbname;
my $estdbuser;
my $estdbhost;
my $input_id;

&get_variables();

my $command = "$runner -runnable $runnable -input_id $input_id -write";

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
  &GetOptions( 
	      'input_id:s'      => \$input_id,
	     );

  $runner     = $::scripts_conf{'runner'};
  $runnable   = $::est_genome_conf{'est_genome_runnable'};
  $refdbname  = $::db_conf{'refdbname'};
  $refdbuser  = $::db_conf{'refdbuser'};
  $refdbhost  = $::db_conf{'refdbhost'};
  $estdbname  = $::db_conf{'estdbname'};
  $estdbuser  = $::db_conf{'estdbuser'};
  $estdbhost  = $::db_conf{'estdbhost'};

  if(!(defined $refdbhost && defined $refdbname & defined $refdbuser &&
       defined $estdbhost && defined $estdbname && defined $estdbuser &&
       defined $runner    && defined $runnable  &&
       defined $input_id)){
    print "Usage: filter_and_e2g.pl -input_id\n" .
      "Additional options to be set in EST_conf.pl: runner, est_genome_runnable, refdbname, refdbuser, refdbhost, estdbname, estdbuser, estdbhost\n";
    exit (1);
  }

}
