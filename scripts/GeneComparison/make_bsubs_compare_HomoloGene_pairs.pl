#!/usr/local/ensembl/bin/perl -w

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
