#!/usr/local/ensembl/bin/perl -w


############################################################
# 
# script to check whether there are any orthologs
# for a given set of genes
#
# it uses heavily compara
#
############################################################


use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Tools::GeneUtils;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Bio::EnsEMBL::Pipeline::GeneComparison::ComparativeTools;
use Getopt::Long;

my $dbhost    = 'ecs2b';
my $dbuser    = 'ensro';
my $dbname    = 'ens_NCBI_31';
#my $port      = '19322';
my $dnadbhost;
my $dnadbname;

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    -host   => $dbhost,
					    -user   => $dbuser,
					    -dbname => $dbname,
					    #-port   => $port,
					    );

my $focus_species = 'Homo sapiens';
my $focus_dbname  = 'ens_NCBI_31';
my $focus_dbhost  = 'ecs2b';

my $target_species = 'Mus musculus';
my $target_dbname  = 'mus_musculus_core_12_3';
my $target_dbhost  = 'kaka';
my $target_dbuser  = 'anonymous';

my $target_species2 = 'Rattus norvegicus';
my $target_dbname2  = 'rattus_norvegicus_core_12_2';
my $target_dbhost2  = 'kaka';
my $target_dbuser2  = 'anonymous';



my $compara_dbname = 'ensembl_compara_12_1';
my $compara_dbhost = 'ecs2d';
my $compara_config =  '/nfs/acari/eae/ensembl/ensembl-compara/modules/Bio/EnsEMBL/Compara/Compara.conf';

my $compara_db = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(
							      -user      => 'ensro',
							      -dbname    => $compara_dbname,
							      -host      => $compara_dbhost,
							      -conf_file => $compara_config,
							     );

my $target_db  = $compara_db->get_db_adaptor($target_species,'MGSC3');
my $target_db2 = $compara_db->get_db_adaptor($target_species2,'RGSC2');
my $focus_db   = $compara_db->get_db_adaptor($focus_species ,'NCBI31');



unless( $focus_db){
  $focus_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
						 -host             => $focus_dbhost,
						 -user             => 'ensro',
						 -dbname           => $focus_dbname,
						);
}
unless( $target_db){
  $target_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
						     -host             => $target_dbhost,
						     -user             => $target_dbuser,
						     -dbname           => $target_dbname,
						    );
}

unless( $target_db2){
  $target_db2 = new Bio::EnsEMBL::DBSQL::DBAdaptor(
						   -host             => $target_dbhost2,
						   -user             => $target_dbuser2,
						   -dbname           => $target_dbname2,
						  );
}


my $id_file;

my $gene_adaptor = $db->get_GeneAdaptor;
my @ids = @{$gene_adaptor->list_geneIds};

my %gene_id;
my $threshold = 40;

foreach my $id (@ids){
  my $gene = $gene_adaptor->fetch_by_dbID( $id, 1);
  foreach my $transcript ( @{$gene->get_all_Transcripts} ){
    
    my $t_id = $transcript->stable_id || $transcript->dbID;
    $gene_id{$transcript} = $gene->stable_id || $gene->dbID;

    print STDERR "\n";
    print STDERR "testing transcript $t_id\n";
    #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Evidence( $transcript );
    Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript( $transcript );
    
    ############################################################
    # has it got homology in mouse?
    #print STDERR "Mouse orthology?\n";
    my $orthologues = Bio::EnsEMBL::Pipeline::GeneComparison::ComparativeTools
      ->test_for_orthology($transcript, $db, $focus_db, $focus_species, $compara_db, $target_db, $target_species, $threshold, \%gene_id );
    
    ############################################################
    # has it got homology in rat?
    #print STDERR "Rat orthology?\n";
    my $orthologues2 = Bio::EnsEMBL::Pipeline::GeneComparison::ComparativeTools
      ->test_for_orthology($transcript, $db, $focus_db, $focus_species, $compara_db, $target_db2, $target_species2, $threshold , \%gene_id);
    print STDERR "\n---------------------\n";
    
    

 #   ############################################################
#    # Does it break synteny?
#    my $synteny_breaking = $pseudogene_finder->test_for_synteny_breaking($transcript, $db, $focus_db, $focus_species, $compara_db, $target_db, $target_species );
#    if ( $synteny_breaking ){
#      print STDERR "synteny is broken\n";
#      $break_synteny = 1;
#    }
#    else{
#      print STDERR "synteny is preserved\n";
#      $break_synteny = 0;
#    }
    
#    ############################################################
#    # Does it overlap with repeats?
#    if ( $pseudogene_finder->overlap_repeats( $transcript, $focus_db ) ){
#      print STDERR "it overlaps with repeats\n";
#      $repeat = 1;
#    }
#    else{
#      print STDERR "it does not overlap with repeats\n";
#      $repeat = 0;
#    }
    
  }
}

############################################################

sub get_coverage{
  my ($transcript) = @_;
  my $coverage;
  foreach my $exon ( @{$transcript->get_all_Exons} ){
    foreach my $evi ( @{$exon->get_all_supporting_features} ){
      $coverage = $evi->score;
      if ($coverage > 0 ){
	return $coverage;
      }
    }
  }
  return 0;
}


############################################################

sub read_probabilities{
  my ($positivo, $negativo, $attributes ) = @_;
  
  my @attributes = @$attributes;

  ############################################################
  # target
  my %positivo;
  $positivo{frameshift}      = 0;
  $positivo{polyA}          = 0;
  $positivo{Met}            = 0;
  $positivo{mouse_homology} = 0;
  $positivo{break_synteny}  = 0;
  $positivo{repeat}         = 0;
  
  my $pos_count = 0;
  
  open( POSITIVO, "<$positivo") || die("cannot open $positivo");
  
  while (<POSITIVO>){
    chomp;
    my @att = split;
    next unless ( $att[1] eq '1' || $att[1] eq '0' );
    
    $positivo{frameshift}      += $att[0];
    $positivo{polyA}          += $att[1];
    $positivo{Met}            += $att[2];
    $positivo{mouse_homology} += $att[3];
    $positivo{break_synteny}  += $att[4];
    $positivo{repeat}         += $att[5];
    
    $pos_count++;
  }
  close(POSITIVO);
  
  ############################################################
  # target
  my %negativo;
  $negativo{frameshift}      = 0;
  $negativo{polyA}          = 0;
  $negativo{Met}            = 0;
  $negativo{mouse_homology} = 0;
  $negativo{break_synteny}  = 0;
  $negativo{repeat}         = 0;
  
  my $neg_count = 0;
  
  open( NEGATIVO, "<$negativo") || die("cannot open $negativo");
  
  while (<NEGATIVO>){
    chomp;
    my @att = split;
    next unless ( $att[1] eq '1' || $att[1] eq '0' );
    
    $negativo{frameshift}      += $att[0];
    $negativo{polyA}          += $att[1];
    $negativo{Met}            += $att[2];
    $negativo{mouse_homology} += $att[3];
    $negativo{break_synteny}  += $att[4];
    $negativo{repeat}         += $att[5];
    
    $neg_count++;
  }

  ############################################################
  # calculate the likelihoods
  
  my %Ppos;
  my %Pneg;
  
  foreach my $att ( @attributes ){
    $Ppos{$att}{0} = ( $pos_count - $positivo{$att} )/ $pos_count;
    $Ppos{$att}{1} = ( $positivo{$att} )/ $pos_count;
    $Pneg{$att}{0} = ( $neg_count - $negativo{$att} )/ $neg_count;
    $Pneg{$att}{1} = ( $negativo{$att} )/ $neg_count;
  }
  
  return ( $pos_count, \%Ppos, $neg_count, \%Pneg );
}

