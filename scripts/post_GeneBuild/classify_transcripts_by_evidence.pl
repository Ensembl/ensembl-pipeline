#!/usr/local/ensembl/bin/perl -w

use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Getopt::Long;

my $dbhost;
my $dnadbhost;
my $dbuser    = 'ensro';
my $dbname;
my $dnadbname;
my $dbpass    = undef;

my $genetype;


$dbuser = "ensro";
GetOptions(
	   'dbname:s'    => \$dbname,
	   'dbhost:s'    => \$dbhost,
#	   'dnadbname:s'  => \$dnadbname,
#	   'dnadbhost:s'  => \$dnadbhost,
	   'genetype:s'  => \$genetype,
	  );

unless ( $dbname && $dbhost ){
  print STDERR "script to find a classification of transcripts by evidence\n";
  print STDERR "Usage: $0 -dbname -dbhost (-genetype)\n";
  exit(0);
}

#my $dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
#					       '-host'   => $dnadbhost,
#					       '-user'   => $dbuser,
#					       '-dbname' => $dnadbname,
#					       '-pass'   => $dbpass,
#					      );


my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $dbhost,
					    '-user'   => $dbuser,
					    '-dbname' => $dbname,
					    '-pass'   => $dbpass,
					    #'-dnadb'  => $dnadb,
					   );

print STDERR "connected to $dbname : $dbhost\n";

if ( $genetype ){
  print STDERR "checking genes of type $genetype\n";
}

my @gene_ids      = @{$db->get_GeneAdaptor->list_geneIds};
my $slice_adaptor = $db->get_SliceAdaptor;

#print "gene_id\ttran_id\texon_number\tframeshifts\ttest\n";
############################################################

GENE:
foreach my $gene_id ( @gene_ids){
  
  my $gene;
  eval{
    $gene = $db->get_GeneAdaptor->fetch_by_dbID($gene_id,1);
  };
  if ( $genetype ){
    next GENE unless ( $genetype eq $gene->type );
  }
  my $gene_id = $gene->stable_id() || $gene->dbID; 
  
  my $all_protein = 0;
  my $all_cdna    = 0;
  my $all_genscan = 0;
  my $all_other_protein = 0;
  
  my %transcripts;

 TRANS:
  foreach my $trans ( @{$gene->get_all_Transcripts} ) {
    
    my $tran_id   = $trans->stable_id || $trans->dbID;
    
    my ($protein,$other_protein,$cdna,$genscan) = &get_evidence($trans);
    if ($protein){
      push ( @{ $transcripts{protein} }, $trans );
    }
    if ($cdna){
      push ( @{ $transcripts{cdna} }, $trans );
    }
    if ( $other_protein ){
      push ( @{ $transcripts{other_protein} }, $trans );
    }
    if( $genscan ){
      push( @{ $transcripts{genscan} }, $trans );
    }
    
    print "TRAN:\t$gene_id\t$tran_id\tprotein:$protein\tother_protein:$other_protein\tcdna:$cdna\tgenscan:$genscan\n";
  }
  if ( $transcripts{protein} ){
    $all_protein = 1;
  }
  if ( $transcripts{cdna} ){
    $all_cdna = 1;
  }
  if ( $transcripts{other_protein} ){
    $all_other_protein = 1;
  }
  if ( $transcripts{genscan} ){
    $all_genscan = 1;
  }
  
  print "GENE:\t$gene_id\tprotein:$all_protein\tother_protein:$all_other_protein\tcdna:$all_cdna\tgenscan:$all_genscan\n";
}

############################################################


sub get_evidence{
  my ($tran) = @_;
  my @exons = @{$tran->get_all_Exons};
  
  my %cdna_evidence;
  my %protein_evidence;
  my %swall_evidence;
  my %other_evidence;

  foreach my $exon ( @exons ){
    foreach my $feature ( @{$exon->get_all_supporting_features}){
      #print STDERR "feature: ".$feature->analysis->logic_name."\n";
      if (  $feature->isa('Bio::EnsEMBL::DnaDnaAlignFeature') ){
	
	if ( $feature->analysis->logic_name =~/cdna/ ){
	  $cdna_evidence{$feature->hseqname}++;
	}
	elsif( $feature->analysis->logic_name =~/Vertrna/i 
	       || $feature->analysis->logic_name =~/Unigene/i ){
	  $other_evidence{$feature->hseqname}++;
	}
      }
      elsif (  $feature->isa('Bio::EnsEMBL::DnaPepAlignFeature') ){
	
	if ( $feature->analysis->logic_name =~/protein/i ){
	  $protein_evidence{$feature->hseqname}++;
	}
	elsif( $feature->analysis->logic_name =~/Swall/i ){
	  $swall_evidence{$feature->hseqname}++;
	}
      }
    }
  }
  my @cdnas          = keys %cdna_evidence;
  my @proteins       = keys %protein_evidence;
  my @swall_evidence = keys %swall_evidence;
  my @other_evidence = keys %other_evidence;
  
  my ($cdna,$protein,$other_protein,$genscan) = (0,0,0,0);
  if ( @cdnas ){
    $cdna = 1;
  }
  if ( @proteins ){
    $protein = 1;
  }
  if ( @swall_evidence && !@other_evidence ){
    $other_protein = 1;
  }
  if ( @other_evidence ){
    $genscan = 1;
  }
  return ($protein,$other_protein,$cdna,$genscan);
}
  
