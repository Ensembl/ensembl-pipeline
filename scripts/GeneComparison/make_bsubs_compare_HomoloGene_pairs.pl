#!/usr/bin/env perl


# Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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


use warnings ;
use strict;  
use Bio::EnsEMBL::Pipeline::GeneComparison::ObjectMap;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Bio::EnsEMBL::Pipeline::Runnable::Blast;
use Getopt::Long;

############################################################
# input is a file with the LL pairs
# LL.23 LL.435
# ...

my $input;
my $tmp_dir;
my $coding_exons;

# options
&GetOptions( 
	    'input:s'     => \$input,
	    'tmp_dir:s'   => \$tmp_dir,
	   );

unless( $input && $tmp_dir){
  print STDERR "Usage: $0 -tmp_dir /full/path/output_dir -input /full/path/locuslink_pairs\n";
  exit(0);
}


# make output directories
my ($bsuberr, $bsubout) = &make_directories();


my $jobfile = "Compare_isoforms_jobs";
open (OUT, ">$jobfile") or die ("Can't open $jobfile for writing: $!");

open ( IN, "<$input") or die ("cannot open $input");
while(<IN>){
  chomp;
  my @entries = split;
  &make_bsubs( $entries[0], $entries[3]);
}

close (OUT) or die (" Error closing $jobfile: $!");


############################################################

sub make_directories {
    
    # bsub output directories
    my $output  = $tmp_dir."/results_isoform_comparison";
    my $bsuberr = $output."/stderr/";
    my $bsubout = $output."/stdout/";
    makedir($output);
    makedir($bsuberr);
    makedir($bsubout);
    return ( $bsuberr, $bsubout );
}

############################################################

sub make_bsubs {
    my ( $human_id, $mouse_id ) = @_;
        
    #my $lsf_options   = "-q acari -C0 -m\"rlx_hosts ecs2_hosts ecs1_hosts\" ";
    my $lsf_options   = "-q acari -C0 -m\"rlx_hosts\"";
    #$lsf_options .= " -R\"select[myecs2f < 440]\" ";
    #$lsf_options .= " -R\"select[myecs2f < 440 && myecs1d < 440] rusage[myecs2f=10:myecs1d=10]\" ";
    
    my $check  = "/nfs/acari/eae/ensembl/ensembl-pipeline/scripts/GeneComparison/compare_HomoloGene_pair.pl";
    my $script = "/nfs/acari/eae/ensembl/ensembl-pipeline/scripts/GeneComparison/compare_HomoloGene_pair.pl";
    my $file_name = "comparison_".$human_id."_".$mouse_id;
    my $outfile   = $bsubout."/".$file_name;
    my $errfile   = $bsuberr."/".$file_name;
    my $command   = "bsub $lsf_options -o $outfile -e $errfile -E \"$check -check\" $script -refseq_mapping $input -gene_id1 $human_id -gene_id2 $mouse_id";
    
    if ( $coding_exons ){
      $command .= " -coding_exons ";
    }
    print OUT "$command\n";
    
  }


############################################################

sub makedir{
    my ($dir) = @_;
    if(opendir(DIR, $dir)){ closedir(DIR); }
    else{ system("mkdir $dir") == 0 or die "error creating $dir\n"; }
}


############################################################
