#!/usr/local/bin/perl -w
# script to take gene of a certain type to a database and
# transfer them to another database

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Getopt::Long;
use strict;


my $genetype;
my $info;
my $source_dbname;
my $source_dbhost;
my $target_dbname;
my $target_dbhost;

&GetOptions('info'             => \$info,
	    'genetype:s'       => \$genetype, 
	    'source_dbname:s'  => \$source_dbname,
            'source_dbhost:s'  => \$source_dbhost,
            'target_dbname:s'  => \$target_dbname,
            'target_dbhost:s'  => \$target_dbhost,
	    );

unless( $source_dbname && $source_dbhost && $target_dbhost && $target_dbname ){
    print "script to transfer genes between two databases with the same assembly\n";
    print "Usage: $0 -source_dbname -source_dbhost -target_dbhost -target_dbname\n";
    print "other options:\n";
    print "-info (optional)\tprints extra info about the genes being transferred\n";
    print "-genetype (optional)\t restricts the transfer to this gene type\n";
    exit(0);
}


# database where we read genes from
my $source_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
						   '-host'   => $source_dbhost,
						   '-user'   => 'ensro',
						   '-dbname' => $source_dbname,
						   );


# database where we put the genes
my  $target_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
						    '-host'   => $target_dbhost,
						    '-user'   => 'ensadmin',
						    '-dbname' => $target_dbname,
						    '-pass'   => 'ensembl',
						    );


# fetch genes from est database
my @gene_ids = @{$source_db->get_GeneAdaptor->list_geneIds};

# get a gene-adaptors for the source and target dbs
my $source_gene_adaptor = $source_db->get_GeneAdaptor;
my $target_gene_adaptor = $target_db->get_GeneAdaptor;

if ($genetype){
    print STDERR "transferring genes of type $genetype\n";
}
else{
    print STDERR "transferring all genes\n";
}

GENE:
foreach my $gene_id (@gene_ids){
  my $gene;
  eval {
    #print STDERR "fetching $gene_id\n";
    $gene = $source_gene_adaptor->fetch_by_dbID($gene_id);
  };
  if ( $@ ){
    my $id;
    if (defined $gene->dbID){
      $id = $gene->dbID;
    }
    else{
      $id = 'no id';
    }
    print STDERR "Unable to read gene ".$id." $@\n";
    next GENE;
  }
  
  # check if we are using type
  if ($genetype){
    unless ( $gene->type eq $genetype ){
      next GENE;
    }
  }
  
    
  ##############################
  # check few things first
  ##############################
  
  if ($info){
    my $id;
    if (defined $gene->dbID){
      $id = $gene->dbID;
    }
    else{
      $id = 'no id';
    }
    print STDERR "gene: ".$id."\n";
    my @transcripts = @{$gene->get_all_Transcripts};
    
    foreach my $tran (@transcripts){
      Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Evidence($tran);
      
      if ( defined $tran->translation){
	Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Translation($tran->translation);
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
    my $id;
    if (defined $gene->dbID){
      $id = $gene->dbID;
    }
    else{
      $id = 'no id';
    }
    print STDERR "Unable to store gene ".$id." $@\n";
    next GENE;
  }
  
  ##############################
  # info about the stored gene
  ##############################
  
  if ($info){
    my $id;
    if (defined $gene->dbID){
      $id = $gene->dbID;
    }
    else{
      $id = 'no id';
    }
    print STDERR "gene: ".$id."\n";
    my @transcripts = @{$gene->get_all_Transcripts};
    
    foreach my $tran (@transcripts){
      Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Evidence($tran);
      
      if ( defined $tran->translation){
	Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Translation($tran->translation);
      }
      else { 
	print STDERR "translation not found!\n";
      }
    }
  }
}



########################################################
