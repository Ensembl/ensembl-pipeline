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
use Bio::EnsEMBL::Pipeline::Config::cDNAs_ESTs::Exonerate qw (
							      EST_TMPDIR
							      EST_REFDBHOST
							      EST_REFDBUSER
							      EST_REFDBNAME
							      LSF_OPTIONS
							      EST_SCRIPTDIR
							      EST_TMPDIR
							      EST_RUNNER
							      EST_EXONERATE_BSUBS
							      EST_EXONERATE_RUNNABLE
							      EST_EXONERATE_ANALYSIS
							      EST_EXONERATE_OPTIONS
							      EST_CHUNKNUMBER
							      EST_CHUNKDIR
							      EST_FILE
							     );

my %chrhash;

# declare these here so we can refer to them later
my $exonerate_bsubdir          = "exonerate_results/";

# make output directories
&make_directories();

# create jobs file for Exonerate
&make_exonerate_bsubs();

############################################################

sub make_directories {
  my $scratchdir =  $EST_TMPDIR ;

  # bsub output directories
  my $bsubdir = $scratchdir . "/" . $exonerate_bsubdir . "/";
  my $bsuberr = $bsubdir . "stderr/";
  my $bsubout = $bsubdir . "stdout/";
  makedir($bsubdir);
  makedir($bsuberr);
  makedir($bsubout);

}

############################################################

sub make_exonerate_bsubs {
  my $jobfile = $EST_EXONERATE_BSUBS;
  open (OUT, ">$jobfile") or die ("Can't open $jobfile for writing: $!");
  
  my $lsf_options   = $LSF_OPTIONS;
  my $scriptdir     = $EST_SCRIPTDIR;
  my $check         = $scriptdir . "/check_node.pl";
  my $exonerate     = $scriptdir . "/run_exonerate.pl";
  my $bsuberr       = $EST_TMPDIR . "/" . $exonerate_bsubdir . "/stderr/";
  my $bsubout       = $EST_TMPDIR . "/" . $exonerate_bsubdir . "/stdout/";
  my $runnable_db   = $EST_EXONERATE_RUNNABLE;
  my $analysis      = $EST_EXONERTAE_ANALYSIS;
  
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
    
    my $command = "bsub $lsf_options -o $outfile -e $errfile -E \"$check $chunk\" $exonerate -runnable $runnable_db -analysis $analysis -query_seq  $chunk_file  -write";
    
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
