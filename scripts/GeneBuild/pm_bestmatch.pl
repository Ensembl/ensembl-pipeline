#!/usr/local/bin/perl

BEGIN {
  # oooh this is not nice
  my $script_dir = $0;
  $script_dir =~ s/(\S+\/)\S+/$1/;
  unshift (@INC, $script_dir);
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
  
  Options are to be set in GeneBuild config files
  The important ones for this script are:
     GeneBuild::Scripts:GB_PM_OUTPUT   directory to write filtered output files

  If GB_PM_OUTPUT is not provided, output file is written to current directory
    
  Database details:
     GeneBuild::Databases::GB_DBNAME
     GeneBuild::Databases::GB_DBHOST
     GeneBuild::Databases::GB_DBUSER
     GeneBuild::Databases::GB_DBPASS
=cut

use strict;
use pmatch_modules; 
use Bio::Seq;
use File::Find;
use Getopt::Long;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Databases qw (
							     GB_DBNAME
							     GB_DBHOST
							     GB_DBUSER
							     GB_DBPASS
							    );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Scripts qw (
							  GB_PM_OUTPUT
							 );


# global vars
my @hits;     # stores the hits from pmatch runs
my $outdir        = $GB_PM_OUTPUT;
my $check         = 0;
my $chromo_coords = 0;
my $outfile       = "pm_best.out";

&GetOptions( 
	    'check'         => \$check,
	    'chromo_coords' => \$chromo_coords,
	   );

if ($check) {
  exit(0);
}

if (defined ($outdir) && $outdir ne '') {
  $outfile = $outdir . "/" . $outfile;
}

open (OUT, ">$outfile") or die "Can't open $outfile:$!\n";

my $sgpa;
if($chromo_coords){
  $sgpa = &get_static_golden_path_adaptor;
}
else{
  print "Output will not be converted to chromosomal coords\n";
}


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

print STDERR "finished reading\n";

 # find the best hit(s) in the genome for each protein in $protfile
 my $pmf2 = new Second_PMF( 
			   '-phits' => \@hits
			  );

 $pmf2->run;

# need to get rid of duplicate lines - eg system sort -u on the output file?
 
 # output the results

 foreach my $hit($pmf2->output) {
  if($chromo_coords){
    my ($chr_name, $chrstart, $chrend) = $sgpa->convert_fpc_to_chromosome(
									 $hit->query, 
									 $hit->qstart,
									 $hit->qend
									);

    print OUT $chr_name       . ":" . $chrstart      . "," . $chrend        . ":" . 
	      $hit->target   . ":" . $hit->coverage . "\n";
  }
  
  else{
    print OUT $hit->query    . ":" . $hit->qstart   . "," . $hit->qend     . ":" . 
	      $hit->target   . ":" . $hit->coverage . "\n";
  }

}


close (OUT) or die "Can't close $outfile:$!\n";

 ### END MAIN

sub get_static_golden_path_adaptor {
  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host             => $GB_DBHOST,
    -user             => $GB_DBUSER,
    -dbname           => $GB_DBNAME,
    -pass             => $GB_DBPASS,
);


  my $sgpa = $db->get_SliceAdaptor();
  return $sgpa;
}
