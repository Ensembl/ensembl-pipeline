#!/usr/local/ensembl/bin/perl -w

=head1 NAME

  make_bsubs.pl

=head1 SYNOPSIS
 
  make_bsubs.pl
  Makes bsub entries for run_blat.pl, etc...
  bsubs can be submitted using submit.pl - they\'re not automatically 
  done from here as it\'s better to submit a few and check they come 
  back OK before sending a genome worth.

  Makes sure all the various scratch subdirectories needed are in place, 
  and makes them if necessary.

=head1 DESCRIPTION


=head1 OPTIONS

=cut

use strict;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Config::cDNAs_ESTs::Blat qw (
							 EST_TMPDIR
							 EST_REFDBHOST
							 EST_REFDBUSER
							 EST_REFDBNAME
							 LSF_OPTIONS
							 EST_SCRIPTDIR
							 EST_TMPDIR
							 EST_RUNNER
							 EST_BLAT_BSUBS
							 EST_BLAT_RUNNABLE
							 EST_BLAT_ANALYSIS
							 EST_CHUNKNUMBER
							 EST_CHUNKDIR
							 EST_FILE
							);

my %chrhash;

# declare these here so we can refer to them later
my $blat_bsubdir          = "blat_results/";

# make output directories
&make_directories();

# create jobs file for Exonerate
&make_blat_bsubs();

############################################################

sub make_directories {
  my $scratchdir =  $EST_TMPDIR ;

  # bsub output directories
  my $bsubdir = $scratchdir . "/" . $blat_bsubdir . "/";
  my $bsuberr = $bsubdir . "stderr/";
  my $bsubout = $bsubdir . "stdout/";
  makedir($bsubdir);
  makedir($bsuberr);
  makedir($bsubout);

}

############################################################

sub make_blat_bsubs {
  my $jobfile = $EST_BLAT_BSUBS;
  open (OUT, ">$jobfile") or die ("Can't open $jobfile for writing: $!");
  
  my $lsf_options   = $LSF_OPTIONS;
  my $scriptdir     = $EST_SCRIPTDIR;
  my $check         = $scriptdir . "/check_node.pl";
  my $blat          = $scriptdir . "/run_blat.pl";
  my $bsuberr       = $EST_TMPDIR . "/" . $blat_bsubdir . "/stderr/";
  my $bsubout       = $EST_TMPDIR . "/" . $blat_bsubdir . "/stdout/";
  my $runnable_db   = $EST_BLAT_RUNNABLE;
  my $analysis      = $EST_BLAT_ANALYSIS;
  
  my $estfile = $EST_FILE; # may be a full path
  my @path = split /\//, $estfile;
  $estfile = $path[$#path];
  $estfile .= "_chunk_";
  
  my $numchunks = $EST_CHUNKNUMBER;
  
  for(my $i = 0; $i < $numchunks; $i++){
    my $num = $i;
    while (length($num) < 7){
      $num = "0" . $num;
    }
    
    my $chunk     = $estfile . $num;
    my $outfile   = $bsubout . $chunk;
    my $errfile   = $bsuberr . $chunk;
    my $chunk_file = $EST_CHUNKDIR."/".$chunk;
    
    my $command = "bsub $lsf_options -o $outfile -e $errfile -E \"$check $chunk\" $blat -runnable $runnable_db -analysis $analysis -query_seq  $chunk_file  -write";
    
    print OUT "$command\n";
  }
  
  close (OUT) or die (" Error closing $jobfile: $!");
}

############################################################

sub makedir{
  my ($dir) = @_;
  if(opendir(DIR, $dir)){ closedir(DIR); }
  else{ system("mkdir $dir") == 0 or die "error creating $dir\n"; }
}


############################################################
