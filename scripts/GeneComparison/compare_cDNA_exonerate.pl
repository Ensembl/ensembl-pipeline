#!/usr/local/ensembl/bin/perl

=head1 NAME

compare_cDNA_exonerate.pl

=head1 SYNOPSIS
 
  compare_cDNA_exonerate.pl
  Runs Runnable::ExonerateEST ( and this runs Runnable::Exonerate )
  over an input cdnafile and a given transcript file/input id
  Writes output to /tmp
  Moves the /tmp output and errorfiles to specified output directory

=head1 DESCRIPTION


=head1 OPTIONS
   -chunkname name of the cDNA chunk to run with
   
   everything else is got from ESTConf.pm and GeneCompConf.pm configuration files

=cut

use strict;
use Getopt::Long;
use Bio::EnsEMBL::Pipeline::ESTConf qw (
					EST_EXONERATE
					EST_RUNNER
				       );

use Bio::EnsEMBL::Pipeline::GeneComparison::GeneCompConf qw (
							     COMP_EXONERATE_RUNNABLE
							     COMP_cDNA_CHUNKDIR
							     COMP_TRANSCRIPTS 
							     COMP_TMPDIR
							    );
my $runner;
my $runnable;
my $exonerate;
my $dbuser;
my $host;
my $estfiledir;
my $chunkname;
my $input_id;
my $outdir;
my $errfile;
my $tmperrfile;
my $outfile;
my $tmpoutfile;

&get_variables();

my $estfile = $estfiledir . "/" . $chunkname;

my $command = "$runner -runnable $runnable -input_id $input_id -parameters cdnas=$estfile,transcripts=$COMP_TRANSCRIPTS,exonerate=$exonerate  2>$tmperrfile | gzip -9 >$tmpoutfile";

#my $command = "$runner -runnable $runnable -input_id $input_id -parameters cdnas=$estfile,transcripts=$COMP_TRANSCRIPTS,exonerate=$exonerate  2>$tmperrfile  >$tmpoutfile";

print STDERR "command is $command\n";

my $output = `$command`;
#my $output = "";

# $output should be empty if redirect has worked
if ($output eq ''){
  my $mv1 = `mv $tmperrfile $errfile`;
  if($mv1 ne ""){
    warn "\nmessage from moving errfile: $mv1\n";
  }

  my $mv2 = `mv $tmpoutfile $outfile`;
  if($mv2 ne ""){
    print STDERR "message from moving outfile: $mv2\n";
  }
}
else {
  warn "\nproblem with exonerate redirect: $output\n";
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
	      'chunkname:s'      => \$chunkname,
	     );

  $runner     = $EST_RUNNER;
  $runnable   = $COMP_EXONERATE_RUNNABLE;
  $exonerate  = $EST_EXONERATE;
  $estfiledir = $COMP_cDNA_CHUNKDIR;
  $input_id   = $COMP_TRANSCRIPTS;
  $outdir     = $COMP_TMPDIR;

  if(!(defined $runner     && defined $runnable  && defined $exonerate &&
       defined $estfiledir && defined $chunkname && 
       defined $input_id   && defined $outdir)){
    print "Usage: compare_cDNA_exonerate.pl -chunkname\n" .
      "Additional options to be set in ESTConf.pm: EST_RUNNER, EST_EXONERATE_RUNNABLE, EST_EXONERATE, EST_REFDBNAME, EST_REFDBUSER, EST_REFDBHOST, EST_CHUNKDIR, EST_GENOMIC and EST_TMPDIR\n";
    exit (1);
  }

  # output directories have been created by make_bsubs.pl
  $outdir   .= "/exonerate_compare_cDNA/results/";
  my $errdir = $outdir . "stderr/";
  $outdir   .= "stdout/";

  die("can't open directory $errdir\n") unless opendir(DIR, $errdir);
  closedir DIR;
  
  die("can't open directory $outdir\n") unless opendir(DIR, $outdir);
  closedir DIR;

  my $err     =  "exest_"      . $chunkname . "_" . $$ . ".stderr";
  $errfile    = $errdir . $err;
  $tmperrfile = "/tmp/" . $err;
  #$tmperrfile = $err;

  my $out     = "exest_"      . $chunkname . "_" . $$ . ".stdout.gz";
  $outfile    = $outdir . $out;
  $tmpoutfile = "/tmp/" . $out;
  #$tmpoutfile = $out;

}
