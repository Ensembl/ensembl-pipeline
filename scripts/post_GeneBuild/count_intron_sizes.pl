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
my $introns;
my $transcripts;

my $genetype;


$dbuser = "ensro";
GetOptions(
	   'dbname:s'    => \$dbname,
	   'dbhost:s'    => \$dbhost,
#	   'dnadbname:s'  => \$dnadbname,
#	   'dnadbhost:s'  => \$dnadbhost,
	   'genetype:s'  => \$genetype,
	   'intron_info_file:s' => \$introns,
	   'transcript_info_file:s' => \$transcripts,
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

my $total_intron_size = 0;
my $total_number_of_introns = 0;
my $one;
my $lesstwenty = 0;
my $lessfifty = 0;
my $maximum;
my $median;
my $mode = 0;
#my $introns = "/ecs2/work2/lec/code/rat_build/data/".$dbname.$dbhost."intron_size.txt";
#my $transcripts = "/ecs2/work2/lec/code/rat_build/data/".$dbname.$dbhost."transcripts_per_intron_size";

$introns = "intron_size.txt" if(!$introns);
$transcripts = "transcripts_per_intron_size.txt". if(!$transcripts);

my $total = scalar(keys(%intron_size));
my $middle = $total/2;
my $rounded = sprintf("%d", $middle);  
print "have opened ".$introns."\n";
open (INTRONS,">$introns") or die;
my $intron_count = 1;

print INTRONS "#intron sizes\n";
foreach my $size ( sort{ $a <=> $b } keys %intron_size ){
  print INTRONS $size."\t".$intron_size{$size}."\n";
  if($size == 1){
    $one = $intron_size{$size};
  }
  if($size <= 20){
    $lesstwenty += $intron_size{$size};
  }
  if($size <= 50){
    $lessfifty += $intron_size{$size};
  }
  if($intron_size{$size} > $mode){
    $mode = $intron_size{$size};
  }
  $total_intron_size += $size;
  $total_number_of_introns += $intron_size{$size};
  $maximum = $size;
  if($intron_count == $rounded){
    $median = $size;
  }
  $intron_count++;
}
close INTRONS;


print "1bp ".$one."\n<=20bp ".$lesstwenty."\n<=50bp ".$lessfifty."\n";

my $average = $total_intron_size/$total_number_of_introns;

print "average intron size = ".$average."\nmaximum intron size = ".
  $maximum."\nmedian intron size ".$median."\nmode intron size ".
  $mode."\n";

print "have opened ".$transcripts."\n";
open(TRANS,">$transcripts") or die;

print TRANS "#transcripts per intron size\n";
foreach my $size ( sort{ $a <=> $b } keys %trans_with_this_intron ){
  print TRANS $size."\t".$trans_with_this_intron{$size}."\n";
}
close TRANS;
