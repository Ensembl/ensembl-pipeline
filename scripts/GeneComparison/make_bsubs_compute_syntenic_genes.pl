#!/usr/local/ensembl/bin/perl -w

use strict;  
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::GeneComparison::ObjectMap;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Bio::EnsEMBL::Pipeline::Runnable::Blast;
use Getopt::Long;


# human_db
my $human_dbname = 'human_estgenes_eae';
my $human_dbhost = 'ecs2b';
my $human = 'Homo sapiens';

my $from_file;
my $gap_penalty;
my $tmp_dir;
my $coding_exons;

# options
&GetOptions( 
	    'from_file'     => \$from_file,
	    'gap_penalty:n' => \$gap_penalty,
	    'tmp_dir:s'     => \$tmp_dir,
	    'coding_exons' => \$coding_exons,
	   );

unless ( $tmp_dir ){
    print STDERR "script to generate bsub lines for the compute_syntenic_gene script\n";
    print STDERR "Usage: $0 -tmp_dir /full/path/output_dir/\n";
    print STDERR "\t\t\t[ -gap_penalty n -coding_exons -from_file < file_with_pairs]\n";
    exit(0);
}

unless( $gap_penalty ){
    $gap_penalty = -100;
}

# make output directories
my ($bsuberr, $bsubout) = &make_directories();


############################################################
# if no file with gene stable ids read it from db

my $jobfile = "Compute_syntenic_genes_jobs";
open (OUT, ">$jobfile") or die ("Can't open $jobfile for writing: $!");

unless( $from_file ){
    
    my $human_db;
    $human_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $human_dbhost,
						   -user  => 'ensro',
						   -dbname=> $human_dbname,
						   );
    
    my @human_ids     = @{$human_db->get_GeneAdaptor->list_stable_geneIds};
    
    foreach my $human_id ( @human_ids ){
      &make_bsubs( $human_id );
    }    
  }

if ( $from_file ){
  while(<>){
    chomp;
    my @entries = split;
    &make_bsubs( $entries[0]);
  }   
}

close (OUT) or die (" Error closing $jobfile: $!");


############################################################

sub make_directories {
    
    # bsub output directories
    my $output  = $tmp_dir."/results_syntenic_computation";
    my $bsuberr = $output."/stderr/";
    my $bsubout = $output."/stdout/";
    makedir($output);
    makedir($bsuberr);
    makedir($bsubout);
    return ( $bsuberr, $bsubout );
}

############################################################

sub make_bsubs {
    my ( $id) = @_;
        
    my $lsf_options   = "-q acari -C0";
    #$lsf_options .= " -R\"select[myecs2f < 440]\" ";
    $lsf_options .= " -R\"select[myecs2b < 440 && myecs2d < 440] rusage[myecs2b=20:duration=15:myecs2d=10:duration:15]\" ";
    
    my $check  = "/nfs/acari/eae/ensembl/ensembl-pipeline/scripts/GeneComparison/compute_syntenic_gene.pl";
    my $script = "/nfs/acari/eae/ensembl/ensembl-pipeline/scripts/GeneComparison/compute_syntenic_gene.pl";
    my $file_name = "syntenic_genes_".$id;
    my $outfile   = $bsubout."/".$file_name;
    my $errfile   = $bsuberr."/".$file_name;
    my $command   = "bsub $lsf_options -o $outfile -e $errfile -E \"$check -check\" $script -gene_id $id";
    
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
