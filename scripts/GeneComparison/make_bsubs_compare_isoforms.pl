#!/usr/local/ensembl/bin/perl

use strict;  
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::GeneComparison::ObjectMap;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Bio::EnsEMBL::Pipeline::Runnable::Blast;
use Getopt::Long;


# human_db
my $human_dbname = 'homo_sapiens_core_13_31';
my $human_dbhost = 'ecs2f';
my $human = 'Homo sapiens';

# human_dnadb
my $human_dnadbname = 'homo_sapiens_core_13_31';
my $human_dnadbhost = 'ecs2f';

# mouse_db
my $mouse_dbname = 'mus_musculus_core_13_30';
my $mouse_dbhost = 'ecs2f';
my $mouse = 'Mus musculus';

# mouse_dnadb
my $mouse_dnadbname = 'mus_musculus_core_13_30';
my $mouse_dnadbhost = 'ecs2f';

# compara_db
my $compara_dbname = 'ensembl_compara_13_1';
my $compara_dbhost = 'ecs2f';


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
    print STDERR "script to generate bsub lines for the compare_isoforms script\n";
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
# if no file with homologs - read from the db

my $jobfile = "Compare_isoforms_jobs";
open (OUT, ">$jobfile") or die ("Can't open $jobfile for writing: $!");

unless( $from_file ){
    
    my $human_dnadb;
    my $mouse_dnadb;
    my $mouse_db;
    my $human_db;
    my $compara_db;
    
    my %gene_pair;
    
    ############################################################
    # connect to the databases 
    $human_dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $human_dnadbhost,
						      -user  => 'ensro',
						      -dbname=> $human_dnadbname,
						      );
    
    $human_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $human_dbhost,
						   -user  => 'ensro',
						   -dbname=> $human_dbname,
						   -dnadb => $human_dnadb,
						   );
    
    
    $mouse_dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $mouse_dnadbhost,
						      -user  => 'ensro',
						      -dbname=> $mouse_dnadbname,
						      );
    
    $mouse_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $mouse_dbhost,
						   -user  => 'ensro',
						   -dbname=> $mouse_dbname,
						   -dnadb => $mouse_dnadb,
						   );
    
    
    #my $human_adaptor = $human_db->get_GeneAdaptor;
    #my $mouse_adaptor = $mouse_db->get_GeneAdaptor;
    
    print STDERR "Using compara_db\n";
    #my $compara_config =  '/nfs/acari/eae/ensembl/ensembl-compara/modules/Bio/EnsEMBL/Compara/Compara.conf';
    $compara_db = 
      Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(
						   -user      => 'ensro',
						   -dbname    => $compara_dbname,
						   -host      => $compara_dbhost,
						   #-conf_file => $compara_config,
						   );
    
    $compara_db->add_db_adaptor( $human_db);
    $compara_db->add_db_adaptor( $mouse_db);
    
    #my $mouse_db = $compara_db->get_db_adaptor($target_species,'NCBIM30');
    #my $human_db = $compara_db->get_db_adaptor($focus_species ,'NCBI31');
    
    my $homol_adaptor = $compara_db->get_HomologyAdaptor;
    my @human_ids     = $homol_adaptor->list_stable_ids_from_species($human);
    
    foreach my $human_id ( @human_ids ){
	
	my @mouse_homologs = 
	    $homol_adaptor->fetch_homologues_of_gene_in_species($human,$human_id,$mouse);
	
	foreach my $homology ( @mouse_homologs ){
	  #print STDERR "comparing $human_id and ".$homology->stable_id."\n";
	  
	  &make_bsubs( $human_id, $homology->stable_id );
	}    
      }
  }

if ( $from_file ){
  while(<>){
    chomp;
    my @entries = split;
    #print STDERR "making bsub for $entries[0] - $entries[1]\n";
    &make_bsubs( $entries[0], $entries[1]);
  }   
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
    my $lsf_options   = "-q acari -C0 ";
    #$lsf_options .= " -R\"select[myecs2f < 440]\" ";
    $lsf_options .= " -R\"select[myecs2f < 440] rusage[myecs2f=20:duration=15]\" ";
    
    my $check  = "/nfs/acari/eae/ensembl/ensembl-pipeline/scripts/GeneComparison/compare_isoforms.pl";
    my $script = "/nfs/acari/eae/ensembl/ensembl-pipeline/scripts/GeneComparison/compare_isoforms.pl";
    my $file_name = "comparison_".$human_id."_".$mouse_id;
    my $outfile   = $bsubout."/".$file_name;
    my $errfile   = $bsuberr."/".$file_name;
    my $command   = "bsub $lsf_options -o $outfile -e $errfile -E \"$check -check\" $script -gene_id1 $human_id -gene_id2 $mouse_id";
    
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
