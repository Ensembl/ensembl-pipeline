#!/usr/local/bin/perl -w

BEGIN {
  # oooh this is not nice
  my $script_dir = $0;
  $script_dir =~ s/(\S+\/)\S+/$1/;
  use lib $script_dir;
  require "EST_conf.pl";
}

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

my %conf = %::EST_conf;

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

my $command = "$runner -runnable $runnable -dbuser $refdbuser -dbname $refdbname -host $refdbhost -input_id $input_id -parameters estdbname=$estdbname,estdbuser=$estdbuser,estdbhost=$estdbhost,seq_index=" . $conf{'estindex'} . " -write";

my $seqinx = "/data/blastdb/Ensembl/humanest";

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

  $runner     = $conf{'runner'};
  $runnable   = $conf{'est_genome_runnable'};
  $refdbname  = $conf{'refdbname'};
  $refdbuser  = $conf{'refdbuser'};
  $refdbhost    = $conf{'refdbhost'};
  $estdbname  = $conf{'estdbname'};
  $estdbuser  = $conf{'estdbuser'};
  $estdbhost    = $conf{'estdbhost'};

  if(!(defined $refdbhost && defined $refdbname & defined $refdbuser &&
       defined $estdbhost && defined $estdbname && defined $estdbuser &&
       defined $runner    && defined $runnable  &&
       defined $input_id)){
    print "Usage: filter_and_e2g.pl -input_id\n" .
      "Additional options to be set in EST_conf.pl: runner, est_genome_runnable, refdbname, refdbuser, refdbhost, estdbname, estdbuser, estdbhost\n";
    exit (1);
  }

}
