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
  print STDERR "script to check splice sites\n";
  print STDERR "Usage: $0 -dbname -dbhost (-genetype)\n";
  exit(0);
}

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

my %intron_size;
my %trans_with_this_intron;

my @gene_ids      = @{$db->get_GeneAdaptor->list_geneIds};
my $slice_adaptor = $db->get_SliceAdaptor;

#print "gene_id\ttran_id\texon_number\tframeshifts\ttest\n";
############################################################

my $count = 0;
GENE:
foreach my $gene_id ( @gene_ids){
  my $gene = $db->get_GeneAdaptor->fetch_by_dbID($gene_id,1);
  #$count++;
  #last if $count >1000;
  if ( $genetype ){
    next GENE unless ( $genetype eq $gene->type );
  }
  my $gene__id = $gene->stable_id() || $gene->dbID; 
  
 TRANS:
  foreach my $trans ( @{$gene->get_all_Transcripts} ) {
    
    my %seen;
    my @exons = sort{ $a->start <=> $b->start } @{$trans->get_all_Exons};
    if ( scalar (@exons ) == 1 ){
      next TRANS;
    }
    for(my $i=0; $i<$#exons; $i++){
      my $intron = $exons[$i+1]->start - $exons[$i]->end - 1;
      $intron_size{$intron}++;
      unless ( $seen{$intron} ){
	$trans_with_this_intron{$intron}++;
      }
    }
  }
}


############################################################

my $introns = "intron_sizes";
my $transcripts = "transcripts_per_intron_size";

open (INTRONS,">$introns") or die;

print INTRONS "#intron sizes\n";
foreach my $size ( sort{ $a <=> $b } keys %intron_size ){
  print INTRONS $size."\t".$intron_size{$size}."\n";
}

close INTRONS;

open(TRANS,">$transcripts") or die;

print TRANS "#transcripts per intron size\n";
foreach my $size ( sort{ $a <=> $b } keys %trans_with_this_intron ){
  print TRANS $size."\t".$trans_with_this_intron{$size}."\n";
}
close TRANS;
