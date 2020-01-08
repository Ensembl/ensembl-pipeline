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

use warnings ;
use strict;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Config::cDNAs_ESTs::Exonerate qw (
							      EST_TMPDIR
							      EST_REFDBHOST
							      EST_REFDBUSER
							      EST_REFDBNAME
							      EST_DBHOST
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
my $exonerate_bsubdir          = "exonerate_results";

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
  
  my $refdb_host = $EST_REFDBHOST;
  my $est_host   = $EST_DBHOST;
  
  #$lsf_options .= " -R\"select[$refdb_host < 440 && $est_host < 440] rusage[$refdb_host=10]\" ";
  $lsf_options .= " -R\"select[$refdb_host < 440 && $est_host < 440] rusage[$refdb_host=10:duration=2:decay=1:$est_host=10:duration=2:decay=1]\"";
    
  my $scriptdir     = $EST_SCRIPTDIR;
  my $check         = $scriptdir . "check_node.pl";
  my $exonerate     = $scriptdir . "run_exonerate.pl";
  my $bsuberr       = $EST_TMPDIR . "/" . $exonerate_bsubdir . "/stderr/";
  my $bsubout       = $EST_TMPDIR . "/" . $exonerate_bsubdir . "/stdout/";
  my $runnable_db   = $EST_EXONERATE_RUNNABLE;
  my $analysis      = $EST_EXONERATE_ANALYSIS;
  
  my $estfile = $EST_FILE; # may be a full path
  my @path = split /\//, $estfile;
  $estfile = $path[$#path];
  $estfile =~ s/\.fa$//i; 
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
    my $query_chunk = "/tmp/".$chunk;

    #my $command = "bsub $lsf_options -f \" $chunk_file > $query_chunk \" -o $outfile -e $errfile -E \"$check $query_chunk\" $exonerate -runnable $runnable_db -analysis $analysis -query_seq  $query_chunk  -write";
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
