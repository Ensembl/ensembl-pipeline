#!/usr/local/ensembl/bin/perl -w

use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Tools::GeneUtils;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Bio::EnsEMBL::Pipeline::Tools::PseudoGeneTests;
use Getopt::Long;

my $dbhost    = 'ecs2c'; 
my $dbuser    = 'ensro';
my $dbname    = 'human_finished_genewises';
my $port;
my $dnadbhost = 'ecs2e';
my $dnadbname = 'human_NCBI33_raw';

my $dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					       '-host'   => $dnadbhost,
					       '-user'   => $dbuser,
					       '-dbname' => $dnadbname,
					    );

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $dbhost,
					    '-user'   => $dbuser,
					    '-dbname' => $dbname,
					    '-dnadb'  => $dnadb,
					    );

my $focus_species = 'Homo sapiens';
my $focus_dbname  = 'ens_NCBI_31';
my $focus_dbhost  = 'ecs2b';

my $target_species = 'Mus musculus';
my $target_dbname  = 'mus_musculus_core_13_30';
#mus_musculus_core_12_3';
my $target_dbhost  = 'ecs2d';
#my $target_dbhost  = 'kaka';

my $target_species2 = 'Rattus norvegicus';
my $target_dbname2  = 'rattus_norvegicus_core_12_2';
my $target_dbhost2  = 'ecs2d';
my $target_dbuser2  = 'ensro';

my $compara_dbname = 'ensembl_compara_13_1';
my $compara_dbhost = 'ecs2d';
my $compara_config =  '/nfs/acari/eae/ensembl/ensembl-compara/modules/Bio/EnsEMBL/Compara/Compara.conf';

my $compara_db = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(
							      -user      => 'ensro',
							      -dbname    => $compara_dbname,
							      -host      => $compara_dbhost,
							      -conf_file => $compara_config,
							     );

my $target_db = $compara_db->get_db_adaptor($target_species,'NCBIM30');
my $focus_db  = $compara_db->get_db_adaptor($focus_species ,'NCBI31');
my $target_db2 = $compara_db->get_db_adaptor($target_species2,'RGSC2');

unless( $focus_db){
  $focus_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
						 -host             => $focus_dbhost,
						 -user             => 'ensro',
						 -dbname           => $focus_dbname,
						);
}


unless( $target_db){
  my $target_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
						     -host             => $target_dbhost,
						     -user             => 'ensro',
						     -dbname           => $target_dbname,
						    );
}

unless( $target_db2){
  my $target_db2 = new Bio::EnsEMBL::DBSQL::DBAdaptor(
						     -host             => $target_dbhost2,
						     -user             => 'ensro',
						     -dbname           => $target_dbname2,
						    );
}

print STDERR "databases target_db = $target_db  target_db2 = $target_db2\n";

my $transcript_id;
my $gene_id;
my $outfile;

&GetOptions(
	    'gene_id:s'          => \$gene_id,
	    'transcript_id:s'    => \$transcript_id,
	    #'outfile:s'          => \$outfile,
	   );


unless ( $gene_id && $transcript_id ){
  print STDERR "Usage: $0 -gene_id  -transcript_id\n";
  exit(0);
}

#my @attributes = ('frameshift', 'polyA', 'Met', 'spliced_elsewhere','mouse_homology', 'rat_homology','break_synteny_mouse','break_synteny_rat', 'repeat');
#my ( $pos_count, $Ppos, $neg_count, $Pneg ) = &read_probabilities(\@attributes );
#my %Ppos = %$Ppos;
#my %Pneg = %$Pneg;


my $pseudogene_tester = Bio::EnsEMBL::Pipeline::Tools::PseudoGeneTests->new();

my $threshold = 40;
  
############################################################
# get the gene in chr coordinates
my $gene = $db->get_GeneAdaptor->fetch_by_dbID( $gene_id,1);

#my $slice = $db->get_SliceAdaptor->fetch_by_chr_name( $gene->chr_name );
#$gene->transform($slice);
#Bio::EnsEMBL::Pipeline::Tools::GeneUtils->_print_Gene( $gene );

my %gene_ref;
foreach my $transcript ( @{$gene->get_all_Transcripts} ){
  #print STDERR "found transcript type: ".$transcript->type."\n";
  my $id = $transcript->dbID;
  next unless ( $id eq $transcript_id );

  $gene_ref{$transcript} = $gene_id;
  ############################################################
  #
  # PSEUDOGENE ATTRIBUTES 
  #
  ############################################################
  
  my ( $frameshift, $polyA, $Met, $spliced_elsewhere, $mouse_homology, $rat_homology, $break_synteny_mouse, $break_synteny_rat, $repeat ) = 
    $pseudogene_tester->pseudogene_test( $transcript, 
					 $db,
					 $compara_db, 
					 $focus_db, $focus_species, 
					 $target_db, $target_species,
					 $target_db2, $target_species2,
					 $threshold,\%gene_ref
				       );
  
  print STDERR "RESULT gene_id:$gene_id\ttranscript_id:$id\tframeshift:$frameshift\tpolyA:$polyA\tMet:$Met\tspliced_elsewhere:$spliced_elsewhere\tmouse_homology:$mouse_homology\trat_homology:$rat_homology\tmouse_synteny_break:$break_synteny_mouse\trat_synteny_break:$break_synteny_rat\trepeats:$repeat\n";

  
 # ############################################################
#  # to be a pseudogene:
#  my $P_pseudo = ( $pos_count / ( $pos_count  + $neg_count ) );
#  for ( my $i=0; $i<=$#attributes; $i++ ){
#    $P_pseudo *= $Ppos{$attributes[$i]}{$results[$i]};
#  }
  
#  ############################################################
#  # not to be a pseudogene:
#  my $P_not_pseudo = ( $neg_count / ( $pos_count  + $neg_count ) );
#  for ( my $i=0; $i<=$#attributes; $i++ ){
#    $P_not_pseudo *= $Pneg{$attributes[$i]}{$results[$i]};
#  }
  
#  ############################################################
#  # normalize it so that we can get the conditional probability
#  # that the target value is pseudogene/no-pseudogene given the
#  # observed attributes:
#  my $is_pseudo     = sprintf "%.3f", ($P_pseudo / ( $P_pseudo + $P_not_pseudo));
#  my $is_not_pseudo = sprintf "%.3f", ($P_not_pseudo / ( $P_pseudo + $P_not_pseudo));
  
#  print STDERR "PREDICTION Pseudogene:$is_pseudo\tNot-Pseudogene: $is_not_pseudo\n";
    
# 
  
}


############################################################

sub read_probabilities{
  my ($attributes ) = @_;
  
  my @attributes = @$attributes;


  ############################################################
  # target
  my %positivo;
  foreach my $att ( @attributes ){
    $positivo{$att} = 0;
  }
  
  my $pos_count = 0;
  
  open( POSITIVO, "</tmp/positivo") || die("cannot open /tmp/positivo");
  
  while (<POSITIVO>){
    chomp;
    my @att = split;
    next unless ( $att[1] eq '1' || $att[1] eq '0' );
    
    for( my $i=0; $i<=$#attributes; $i++ ){
      $positivo{$attributes[$i]} += $att[$i];
    }    
    $pos_count++;
  }
  close(POSITIVO);
  
  ############################################################
  # target
  my %negativo;
  foreach my $att ( @attributes ){
    $negativo{$att} = 0;
  }
    
  my $neg_count = 0;
  
  open( NEGATIVO, "</tmp/negative") || die("cannot open /tmp/negative");
  
  while (<NEGATIVO>){
    chomp;
    my @att = split;
    next unless ( $att[1] eq '1' || $att[1] eq '0' );
    
    for( my $i=0; $i<=$#attributes; $i++ ){
      $negativo{$attributes[$i]} += $att[$i];
    }   
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

