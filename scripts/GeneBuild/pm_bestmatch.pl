#!/usr/local/bin/perl -w

BEGIN {
  # oooh this is not nice
  my $script_dir = $0;
  $script_dir =~ s/(\S+\/)\S+/$1/;
  unshift (@INC, $script_dir);
  require "GB_conf.pl";
}

=head1 NAME

  pm_bestmatch.pl

=head1 SYNOPSIS
 
  pm_bestmatch.pl

=head1 DESCRIPTION

  pm_bestmatch.pl processes data output from pm_filter.pl and returns the best hit(s) 
  in the genome for each pmatched protein.

  pm_bestmatch.pl takes input from STDIN and returns the best hit(s) for each of the 
  proteins whose results are passed in. If all the pm_filter.pl output files from 
  the various chromosomes are catted together, pm_bestmatch.pl will return the best 
  hit(s) in the genome.

  The top hit and any hits with coverage within 2% are returned.

  Typical usage:

  cat *.pm.out > allhits.pmatch
  ./pm_bestmatch.pl < allhits.pmatch
  

=head1 OPTIONS
  
  Options are to be set in GB_conf.pl
  The important ones for this script are:
     pm_output   directory to write filtered output files

     eg.
	    'pm_output'   => '/work2/vac/GeneBuild/',
  
  If pm_output directory is not provided, output file is written to current directory
    
=cut

use strict;
use pmatch_modules; 
use Bio::Seq;
use File::Find;
use Getopt::Long;

my %conf =  %::GB_conf; # configuration options

# global vars
my %plengths; # stores lengths of all proteins in $protfile 
my @hits;     # stores the hits from pmatch runs
my $outdir   = $conf{'pm_output'};
my $check    = 0;
my $outfile  = "pm_best.out";

&GetOptions( 
	    'check' => \$check
	   );

if ($check) {
  exit(0);
}

if (defined ($outdir) && $outdir ne '') {
  $outfile = $outdir . "/" . $outfile;
}

open (OUT, ">$outfile") or die "Can't open $outfile:$!\n";

# read pmatch results from STDIN
while(<>){

  #ctg12770:1140298,1146110:Q15486:18,140 75.8
  next unless /\S+/;
  next unless /^(\S+):(\d+),(\d+):(\S+):(\d+),(\d+)\s+(\S+)/;

  my $query   = $1;
  my $target  = $4;
  my $qstart  = $2;
  my $qend    = $3;
  my $tstart  = $5;
  my $tend    = $6;
  my $percent = $7;
  my $strand  = 1;

  if ($qend < $qstart) { $strand = -1; }

  my $cp = new CoordPair(
			 '-query'   => $query,
			 '-target'  => $target,
			 '-qstart'  => $qstart,
			 '-qend'    => $qend,
			 '-tstart'  => $tstart,
			 '-tend'    => $tend,
			 '-percent' => $percent,
			 '-strand'  => $strand,
			);

  my $mh = new MergedHit(
			 '-query'   => $query,
			 '-target'  => $target,
			 '-strand'  => $strand,
			 '-coverage'=> $percent,
			);

  $mh->add_CoordPair($cp);

  push (@hits, $mh);
}

 # find the best hit(s) in the genome for each protein in $protfile
 my $pmf2 = new Second_PMF( 
			   '-phits' => \@hits
			  );

 $pmf2->run;

 # output the results
 foreach my $hit($pmf2->output) {
   print OUT $hit->query  . ":" . 
         $hit->qstart . "," . 
         $hit->qend   . ":" . 
         $hit->target . ":" .
         $hit->tstart . "," .
         $hit->tend   . "\n";
 }

close (OUT) or die "Can't close $outfile:$!\n";

 ### END MAIN




