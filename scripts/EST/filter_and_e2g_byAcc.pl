#!/usr/local/ensembl/bin/perl -w

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

use strict;
use Getopt::Long;
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
  &GetOptions( 
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
