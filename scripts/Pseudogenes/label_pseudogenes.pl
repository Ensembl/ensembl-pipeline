#!/usr/local/bin/perl -w

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Bio::EnsEMBL::Pipeline::Tools::GeneUtils;
use Bio::EnsEMBL::Pipeline::Config::PseudoGenes::PseudoGenes;
use Bio::EnsEMBL::Pipeline::Tools::PseudoGeneTests;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Getopt::Long;
use strict;

my $gene_id;
my $check;

&GetOptions('gene_id:s'       => \$gene_id,
	    'check'     => \$check,
	   );

if ( $check ){
  exit(0);
}


unless( $gene_id ){
  print "script to test predicted genes for procesed-pseudogene properties\n";
  print "Usage: $0  -gene_id\n";
  exit(0);
}


# compara database
my $compara_db = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(
							      -user      => 'ensro',
							      -dbname    => $COMPARA_DBNAME,
							      -host      => $COMPARA_DBHOST,
							      );

# database where the dna and precomputes (repeats) are:
my $refdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					       '-host'   => $REF_DBHOST,
					       '-user'   => 'ensro',
					       '-dbname' => $REF_DBNAME,
					      );

# database where the predictions are
my $focus_species = $FOCUS_SPECIES;
my  $focus_db  = $compara_db->get_db_adaptor($focus_species, $LABEL_PATH);

unless( $focus_db ){
  $focus_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
						    '-host'   => $LABEL_DBHOST,
						    '-user'   => 'ensro',
						    '-dbname' => $LABEL_DBNAME,
						   );
  $compara_db->add_db_adaptor( $focus_db );
}
$focus_db->dnadb( $refdb );


my $gene = $focus_db->get_GeneAdaptor->fetch_by_dbID( $gene_id, 1);
my $db = $focus_db;


my ($target_db, $target_db2);
my ($target_species, $target_species2);

if ( @$COMPARATIVE_DBS ){
  if ( $COMPARATIVE_DBS->[0] ){
    $target_species = $COMPARATIVE_DBS->[0]->{SPECIES};
    
    $target_db = $compara_db->get_db_adaptor($target_species,$COMPARATIVE_DBS->[0]->{PATH});
    
    unless ( $target_db ){
      $target_db =  
	Bio::EnsEMBL::DBSQL::DBAdaptor
	  ->new(
		-user      => 'ensro',
		-dbname    => $COMPARATIVE_DBS->[0]->{DBNAME},
		-host      => $COMPARATIVE_DBS->[0]->{DBHOST},
	       );
      $target_db->assembly_type( $COMPARATIVE_DBS->[0]->{PATH} );
      $compara_db->add_db_adaptor( $target_db );
    }
  }
  if ( $COMPARATIVE_DBS->[1] ){
    $target_species2 = $COMPARATIVE_DBS->[1]->{SPECIES};
    $target_db2 = $compara_db->get_db_adaptor($target_species2,$COMPARATIVE_DBS->[1]->{PATH});
    
    unless( $target_db2){
      $target_db2 =  
	Bio::EnsEMBL::DBSQL::DBAdaptor
	  ->new(
		-user      => 'ensro',
		-dbname    => $COMPARATIVE_DBS->[1]->{DBNAME},
		-host      => $COMPARATIVE_DBS->[1]->{DBHOST},
	       );
      $target_db2->assembly_type( $COMPARATIVE_DBS->[1]->{PATH} );
      $compara_db->add_db_adaptor( $target_db2 );
    }
  }
}


print STDERR "databases target_db = $target_db  target_db2 = $target_db2\n";

my $pseudogene_tester = Bio::EnsEMBL::Pipeline::Tools::PseudoGeneTests->new();

my $threshold = 40;
my $evidence = 1;

print STDERR "\n ######## analysing gene $gene_id #######\n";
    
my $trans_count = scalar( @{ $gene->get_all_Transcripts} );

my %gene_ref;
foreach my $transcript ( @{$gene->get_all_Transcripts} ){

    #print STDERR "found transcript type: ".$transcript->type."\n";
  my $transcript_id = $transcript->stable_id || $transcript->dbID;
  
  $gene_ref{$transcript} = $gene_id;
  $transcript->type($gene->type);
  
  ############################################################
  #
  # PSEUDOGENE ATTRIBUTES 
  #
  ############################################################
  
  my ( $frameshift, $polyA, $Met, $spliced_elsewhere, $splice_sites_correct, $mouse_homology, $rat_homology, $break_synteny_mouse, $break_synteny_rat, $repeat ) = 
    $pseudogene_tester->pseudogene_test( $transcript, 
					 $db,
					 $compara_db, 
					 $focus_db, $focus_species, 
					 $target_db, $target_species,
					 $target_db2, $target_species2,
					 $threshold,\%gene_ref
				       );
  
  my $exon_number = scalar( @{$transcript->get_all_Exons} );
  
  if ($evidence){
    my $evidence_ids = &_get_evidence($transcript);
    print STDERR "RESULT gene_id:$gene_id\ttranscript_id:$transcript_id\texon_number:$exon_number\tframeshift:$frameshift\tpolyA:$polyA\tMet:$Met\tspliced_elsewhere:$spliced_elsewhere\tss_correct:$splice_sites_correct\thomology1:$mouse_homology\thomology2:$rat_homology\trepeats:$repeat\tevidence:$evidence_ids\ttrans_count:$trans_count\n";
  }
  else{
    $mouse_homology = 0 unless defined($mouse_homology);
    $rat_homology = 0 unless defined($rat_homology);
    
    print STDERR "RESULT gene_id:$gene_id\ttranscript_id:$transcript_id\texon_number:$exon_number\tframeshift:$frameshift\tpolyA:$polyA\tMet:$Met\tspliced_elsewhere:$spliced_elsewhere\tss_correct:$splice_sites_correct\tmouse_homology:$mouse_homology\trat_homology:$rat_homology\trepeats:$repeat\ttrans_count:$trans_count\n";
  }
}



sub _get_evidence{
    my ($transcript) = @_;
    my %evi;
    foreach my $exon (@{$transcript->get_all_Exons}){
	foreach my $evidence ( @{$exon->get_all_supporting_features} ){
	    $evi{$evidence->hseqname}++;
	}
    }
    my @evi_ids = keys %evi;
    my $ids = join ' ',@evi_ids;
    return $ids;
}
