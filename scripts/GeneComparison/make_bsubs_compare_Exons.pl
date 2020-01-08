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


use warnings ;
use strict;  
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::GeneComparison::ObjectMap;
use Bio::EnsEMBL::Pipeline::GeneComparison::GeneCompConf;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Bio::EnsEMBL::Pipeline::Runnable::Blast;
use Getopt::Long;


my $tmp_dir;
my $coding_exons;

# options
&GetOptions( 
	    'tmp_dir:s'     => \$tmp_dir,
	    'coding_exons' => \$coding_exons,
	   );

unless ( $tmp_dir ){
    print STDERR "script to generate bsub lines for the compare_Exons.pl script\n";
    print STDERR "Usage: $0 -tmp_dir /full/path/output_dir/\n";
    print STDERR "\t\t\t[ -coding_exons ]\n";
    exit(0);
}

# make output directories
my ($bsuberr, $bsubout) = &make_directories();


############################################################
# if no file with homologs - read from the db

my $jobfile = "Compare_Exons";
open (OUT, ">$jobfile") or die ("Can't open $jobfile for writing: $!");
  

my @chr_names = &get_chrs;
foreach my $chr_name ( @chr_names ){
  &make_bsubs( $chr_name );
}    

close (OUT) or die (" Error closing $jobfile: $!");


############################################################

sub get_chrs{
  my $dbname = $DBNAME1;
  my $dbhost = $DBHOST1;
  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $dbhost,
					   -user  => 'ensro',
					   -dbname=> $dbname,
						   );
  my @chrs = @{$db->get_ChromosomeAdaptor->fetch_all};

  my @names = map { $_->chr_name } @chrs;
  return @names;
}

############################################################

sub make_directories {
    
  # bsub output directories
  my $output  = $tmp_dir."/results_exon_comparison";
  my $bsuberr = $output."/stderr/";
  my $bsubout = $output."/stdout/";
  makedir($output);
  makedir($bsuberr);
  makedir($bsubout);
  return ( $bsuberr, $bsubout );
}

############################################################

sub make_bsubs {
    my ( $chr_name ) = @_;
        
    #my $lsf_options   = "-q acari -C0 -m\"rlx_hosts ecs2_hosts ecs1_hosts\" ";
    my $lsf_options   = "-q acari -C0 ";
    #$lsf_options .= " -R\"select[myecs2f < 440]\" ";
    $lsf_options .= " -R\"select[myecs2b < 440 && myecs2a < 440] rusage[myecs2b=20:duration=15:myecs2a=20:duration=15]\" ";
    
    my $check  = "/nfs/acari/eae/ensembl/ensembl-pipeline/scripts/GeneComparison/compare_Exons.pl";
    my $script = "/nfs/acari/eae/ensembl/ensembl-pipeline/scripts/GeneComparison/compare_Exons.pl";
    my $file_name = "comparison_".$chr_name;
    my $outfile   = $bsubout."/".$file_name;
    my $errfile   = $bsuberr."/".$file_name;
    my $command   = "bsub $lsf_options -o $outfile -e $errfile -E \"$check -check\" $script -chr_name $chr_name";
    
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
