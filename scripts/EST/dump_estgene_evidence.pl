#!/usr/local/ensembl/bin/perl -w

# script that dumps a tab delimited file first column contains the
# transcript ids (defaulted to stable_id) and the second column contains
# the ids of the evidence (multiple in some cases) from the evidence

use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Bio::SeqIO;
use Getopt::Long;

my $file;

# database for the est genes:
my $dbhost = 'ecs2c';
my $dbuser = 'ensro';
my $dbname = 'ens_NCBI_31_estgene';
my $dbpass = undef;


my $genetype = 'genomewise';

&GetOptions(
	    'dbhost:s'        => \$dbhost,
	    'dbname:s'        => \$dbname,
	    'genetype:s'      => \$genetype,
	    'output:s'          => \$file,
	   );

unless ( $file){
  print STDERR "script to dump the evidence of any set of genes in the format \"transcript_id\\tevidence_id\"\n";
  print STDERR "Usage: $0 -output filename [ -dname -dbhost -genetype ]\n";
  exit(0);
}

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $dbhost,
					    '-user'   => $dbuser,
					    '-dbname' => $dbname,
					    '-pass'   => $dbpass,
					   );
print "connected to $dbname : $dbhost\n";

open( OUT, ">$file" ) || die("cannot open file $file");

my $gene_adaptor = $db->get_GeneAdaptor;
my @gene_ids = @{$gene_adaptor->list_geneIds()};
print STDERR "found ".scalar(@gene_ids)." genes\n";

foreach my $id ( @gene_ids ){
  my $gene;
  eval{
    $gene = $gene_adaptor->fetch_by_dbID($id);
  };
  unless ($gene ){
    print STDERR $@;
    next;
  }
  next unless ($gene->type eq $genetype);
  
  foreach my $tran ( @{$gene->get_all_Transcripts} ){
    
    my %evidence;
    
    foreach my $exon ( @{$tran->get_all_Exons} ){
      my @features =  @{$exon->get_all_supporting_features};
      #unless ( @features ){
      #	print STDERR "no evidence\n";
      #	print STDERR $exon->gffstring."\n";
      #      }
      foreach my $feature ( @{$exon->get_all_supporting_features} ){
	#print STDERR $feature->gffstring."\n";
	#print STDERR $feature->hseqname."\n";
	
	$evidence{ $feature->hseqname }++;
	
	
	
      }
    }
    
    my $id;
    if ( $tran->stable_id ){
      $id = $tran->stable_id;
    }
    else{
      $id = $tran->dbID;
    }

    my @names = keys %evidence;
    my $count = 0;
  OUT1:
    foreach my $name ( @names ){
      print OUT $id."\t".$name."\n";
      #$count++;
      #if ( $count > 9 ){
      #	last OUT1;
      #      }
    }
  }
}

close OUT;
