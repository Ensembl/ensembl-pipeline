#!/usr/local/bin/perl -w
# script to take gene of a certain type to a database and
# transfer them to another database

# in this case we want to take all est_genomewise genes from
# 'manu_drosophila_et_complete' and put them together with
# similarity and targetted genewises in anew database

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;
use strict;


my $path    = 'fly';
my $type;
my $info;
my $help;

&GetOptions('type' => \$type, 
	    'info' => \$info,
	    'help' => \$help,
	   );

if ($help){
  print "script to transfer genes between two databases with the same assembly\n";
  print "Database details must be filled in in the script. Options:\n";
  print "-info (optional)\tprints extra info about the genes being transferred\n";
  print "-type (optional)\t restricts the transfer to this gene type\n";
  print "-help\tfor this help\n";
}

# dnadb adaptor
# dna database, we put this into the main database we read from, so if it needs any dna, 
# where to get it from
my $dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => 'ecs1d',
					    '-user'   => 'ensro',
					    '-dbname' => 'drosophila3',
					   );


# database where we read genes from
my $source_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => 'ecs1d',
					    '-user'   => 'ensro',
					    '-dbname' => 'manu_drosophila_est_complete',
					    '-dnadb'  => $dnadb,   
					   );


# big db adaptor  - already hold sim gw & targetted gw, we want to write
# est genes in here as well
my $target_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => 'ecs1e',
					    '-user'   => 'ensadmin',
					    '-dbname' => 'drosophila_genewise_genomewise_genes',
					    '-dnadb'  => $dnadb,   
					    '-pass'   => 'ensembl',
					    );


# fetch genes from est database
$source_db->static_golden_path_type($path);
my $source_sgp = $source_db->get_StaticGoldenPathAdaptor;

my @gene_ids = $source_db->get_GeneAdaptor->list_geneIds;

# get a gene-adaptors for the source and target dbs
my $source_gene_adaptor = $source_db->get_GeneAdaptor;
my $target_gene_adaptor = $target_db->get_GeneAdaptor;

foreach my $gene_id (@gene_ids){
  my $gene;
  eval {
    $gene = $source_gene_adaptor->fetch_by_dbID($gene_id);
  };
  if ( $@ ){
    print STDERR "Unable to read gene ".$gene->dbID."\n";
  }
  
  # check if we are using type
  if ($type){
    unless ( $gene->type eq $type ){
      next;
    }
  }

  
  ##############################
  # check few things first
  ##############################
  
  if ($info){
    print STDERR "gene: ".$gene->dbID."\n";
    my @transcripts = $gene->each_Transcript;
    
    foreach my $tran (@transcripts){
      print STDERR "transcript: ".$tran->dbID."\n";
      my @exons = $tran->get_all_Exons;
      
      foreach my $exon (@exons){
	print STDERR "exon: ".$exon->dbID." start: ".$exon->start." end: ".$exon->end."\n";
	my @evidence = $exon->each_Supporting_Feature();
	unless (@evidence){
	  print STDERR "no evidence retrieved!\n";
	}
      }
      my $translation = $tran->translation;
      if ( $translation ){
	print STDERR "translation found\n";
	print STDERR "start_exon_phase: ".$translation->start_exon->phase."\n";
	print STDERR "end_exon_phase  : ".$translation->end_exon->phase  ."\n";
      }
      else { 
	print STDERR "translation not found!\n";
      }
    }
  }
  ##############################
  # store gene, same assembly so no need to massage it
  ##############################
  eval{	
    print STDERR "storing gene ...\n";
    $target_gene_adaptor->store($gene);
  };
  if ( $@ ){
    print STDERR "Unable to store gene ".$gene->dbID."\n";
  }
  
  ##############################
  # info about the stored gene
  ##############################
  
  if ($info){
    print STDERR "gene: ".$gene->dbID."\n";
    my @newtranscripts = $gene->each_Transcript;
    
    foreach my $tran (@newtranscripts){
      print STDERR "transcript: ".$tran->dbID."\n";
      my @exons = $tran->get_all_Exons;
      print STDERR scalar( @exons )." exon(s) found\n";
      
      foreach my $exon (@exons){
	print STDERR "exon: ".$exon->dbID." start: ".$exon->start." end: ".$exon->end."\n";
	my @evidence = $exon->each_Supporting_Feature();
	unless (@evidence){
	  print STDERR "no evidence stored!\n";
	}
      }
      my $translation = $tran->translation;
      if ( $translation ){
	print STDERR "translation found\n";
	print STDERR "start_exon_phase: ".$translation->start_exon->phase."\n";
	print STDERR "end_exon_phase  : ".$translation->end_exon->phase  ."\n";
      }
    }
    print STDERR "\n";
  }
}


########################################################
