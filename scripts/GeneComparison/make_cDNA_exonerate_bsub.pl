#!/usr/bin/env perl


# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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

  make_cDNA_exonerate_bsub.pl

=head1 SYNOPSIS
 
  Makes bsub entries for compare_cDNA_exonerate.pl
  bsubs can be submitted using submit.pl or anything similar - they\'re not automatically 
  done from here as it\'s better to submit a few and check they come 
  back OK before sending a genome worth.

  Makes sure all the various output subdirectories needed are in place, 
  and makes them if necessary.

=head1 DESCRIPTION


=head1 OPTIONS

=cut

use warnings ;
use strict;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::ESTConf qw (
					EST_SCRIPTDIR
				       );

use Bio::EnsEMBL::Pipeline::GeneComparison::GeneCompConf qw (
							     COMP_EXONERATE_RUNNABLE
							     COMP_cDNA_CHUNKDIR
							     COMP_cDNA_CHUNKNUMBER
							     COMP_cDNA_FILE
							     COMP_TRANSCRIPTS 
							     COMP_TMPDIR
							     COMP_EXONERATE_BSUBS
							     COMP_QUEUE
							     COMP_SCRIPTDIR
							    );
my %chrhash;

# declare these here so we can refer to them later
my $ex_resultsdir       = "exonerate_compare_cDNA/results";
my $ex_bsubdir          = "exonerate_compare_cDNA/bsub";

# make output directories
&make_directories();

# create jobs file for Exonerate
&make_exonerate_bsubs();

=head2 make_directories

  Title   : make_directories
  Usage   : make_directories
  Function: makes sure needed output directories exist, and if not, makes them
  Returns : none - uses globals
  Args    : none - uses globals

=cut

sub make_directories {
  my $scratchdir =  $COMP_TMPDIR ;

  my @resdirs = split /\//, $ex_resultsdir;

  # exonerate_ests
  my $exoneratedir = $scratchdir . "/" . $resdirs[0] . "/";
  makedir($exoneratedir);

  # exonerate output directories
  my $exdir = $exoneratedir . $resdirs[1] . "/";  
  my $exerr = $exdir . "stderr/";
  my $exout = $exdir . "stdout/";
  makedir($exdir);
  makedir($exerr);
  makedir($exout);
  
  # bsub output directories
  my $bsubdir = $scratchdir . "/" . $ex_bsubdir . "/";
  my $bsuberr = $bsubdir . "stderr/";
  my $bsubout = $bsubdir . "stdout/";
  makedir($bsubdir);
  makedir($bsuberr);
  makedir($bsubout);
}

=head2 make_exonerate_bsubs

  Title   : make_exonerate_bsubs
  Usage   : make_exonerate_bsubs
  Function: makes bsubs to run exonerate_est.pl
  Returns : none - uses globals
  Args    : none - uses globals

=cut

sub make_exonerate_bsubs {
  my $jobfile = $COMP_EXONERATE_BSUBS;
  open (OUT, ">$jobfile") or die ("Can't open $jobfile for writing: $!");

  my $queue         = $COMP_QUEUE;
  my $scriptdir     = $COMP_SCRIPTDIR;
  my $check         = $COMP_SCRIPTDIR . "/check_node.pl";
  my $exonerate_est = $scriptdir . "/compare_cDNA_exonerate.pl";
  my $bsuberr       = $COMP_TMPDIR . "/" . $ex_bsubdir . "/stderr/";
  my $bsubout       = $COMP_TMPDIR . "/" . $ex_bsubdir . "/stdout/";
  
  my $estfile = $COMP_cDNA_FILE; # may be a full path
  my @path = split /\//, $estfile;
  $estfile = $path[$#path];
  $estfile .= "_chunk_";

  my $numchunks = $COMP_cDNA_CHUNKNUMBER;

  for(my $i = 0; $i < $numchunks; $i++){
    my $num = $i;
    while (length($num) < 7){
      $num = "0" . $num;
    }
    
    my $chunk     = $estfile . $num;
    my $outfile   = $bsubout . $chunk;
    my $errfile   = $bsuberr . $chunk;
    
    my $command = "bsub -q $queue -C0 -o $outfile -e $errfile -E \"$check $chunk\" $exonerate_est -chunkname $chunk";
    print OUT "$command\n";
  }
  
  close (OUT) or die (" Error closing $jobfile: $!");
}


=head2 makedir

  Title   : makedir
  Usage   : makedir
  Function: checks to see if the given directory exists, and makes it if it doesn\'t
  Returns : nothing, but will kill script if mkdir fails
  Args    : $dir - directory to be created

=cut

sub makedir{
  my ($dir) = @_;
  if(opendir(DIR, $dir)){ closedir(DIR); }
  else{ system("mkdir $dir") == 0 or die "error creating $dir\n"; }
}
