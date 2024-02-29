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


=head1 SYNOPSIS
 
  make_bsubs.pl
  Makes bsub entries for some runner script
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
use Bio::EnsEMBL::Pipeline::Config::PseudoGenes::PseudoGenes;

my %chrhash;

# declare these here so we can refer to them later
my $pseudogene_bsubdir          = "results/";

# make output directories
&make_directories();

# create jobs file for Exonerate
&make_pseudogene_bsubs();

############################################################

sub make_directories {
  my $scratchdir =  $TMPDIR ;

  makedir($scratchdir);
  # bsub output directories
  my $bsubdir = $scratchdir . "/" . $pseudogene_bsubdir . "/";
  my $bsuberr       = $TMPDIR . "/" . $pseudogene_bsubdir . "/stderr/";
  my $bsubout       = $TMPDIR . "/" . $pseudogene_bsubdir . "/stdout/";
  makedir($bsubdir);
  makedir($bsuberr);
  makedir($bsubout);
}

############################################################

sub make_pseudogene_bsubs {
  my $jobfile = $BSUBS_FILE;

  open (OUT, ">$jobfile") or die ("Can't open $jobfile for writing: $!");
  
  my $lsf_options   = $LSF_OPTIONS;
  my $pseudogene    = $PSEUDO_SCRIPT;
  my $pre_exec      = $PSEUDO_PRE_EXEC;
  my $bsuberr       = $TMPDIR . "/" . $pseudogene_bsubdir . "/stderr/";
  my $bsubout       = $TMPDIR . "/" . $pseudogene_bsubdir . "/stdout/";
  
  my $dir = $CHUNKS;
  open LS, "ls $dir |" or die "Can't open pipe from ls : $!";
  while(<LS>){
    chomp;
    my $chunk    = $_;
    print STDERR "chunk: $chunk\n";
    my $outfile   = $bsubout . $chunk;
    my $errfile   = $bsuberr . $chunk;
    
    my $command = "bsub $lsf_options -o $outfile -e $errfile -E \"$pre_exec\" $pseudogene -runnable Bio::EnsEMBL::Pipeline::RunnableDB::PseudoGeneFinder -analysis pseudogene -query_seq $chunk -write";
    
    print OUT "$command\n";
  }
  close(LS);
  close (OUT) or die (" Error closing $jobfile: $!");
}

############################################################

sub makedir{
  my ($dir) = @_;
  if(opendir(DIR, $dir)){ closedir(DIR); }
  else{ system("mkdir $dir") == 0 or die "error creating $dir\n"; }
}


############################################################
