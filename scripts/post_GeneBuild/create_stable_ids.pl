#!/usr/local/bin/perl -w

use strict;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::SeqIO;
use Getopt::Long;
use Bio::EnsEMBL::Pipeline::ESTConf;

my $dbhost;
my $dbuser    = 'ensro';
my $dbname;
my $label;
my $dbport = '3306' ;
&GetOptions(
	    'dbname:s'    => \$dbname,
	    'dbhost:s'    => \$dbhost,
            'dbport:i'    => \$dbport,
	    'label:s'     => \$label,

	   );

unless ( $dbname && $dbhost && $label ){
  print STDERR "script to fake some stable ids for genes, transcripts, exons and translations\n";
  print STDERR "Usage: $0 -dbname -dbhost -label ( ENSEST, ENS, ENSRNO, ...)\n";
  exit(0);
}

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $dbhost,
					    '-user'   => $dbuser,
					    '-dbname' => $dbname,
                                            '-port'   => $dbport,
					   );

print STDERR "connected to databse  $dbname : $dbhost\n";

my $gene_file        = "gene_stable_id";
my $transcript_file  = "transcript_stable_id";
my $translation_file = "translation_stable_id";
my $exon_file        = "exon_stable_id";

open (GENE,">$gene_file")               or die("unable to open file $gene_file");
open (TRANSCRIPT,">$transcript_file")   or die("unable to open file $transcript_file");
open (TRANSLATION,">$translation_file") or die("unable to open file $translation_file");
open (EXON,">$exon_file")               or die("unable to open file $exon_file");

my $count_gene = 0;
my $count_transcript = 0;
my $count_exon = 0;

my %seen_gene;
my %seen_exon;
my %seen_transcript;

my $gene_adaptor = $db->get_GeneAdaptor;
foreach my $gene_id( @{$gene_adaptor->list_dbIDs} ) {
  
  #print STDERR "gene id $gene_id\n";
  eval {
    my $gene = $db->get_GeneAdaptor->fetch_by_dbID($gene_id);
    
    ####################
    # create gene id 
    ####################
    $count_gene++;
    my $gene_id = $gene->dbID;
    my $gene_stable_id = &make_stable_id($label."G",$gene);
    #  927491 | ENSestG00000000001 |       1 | 0000-00-00 00:00:00 | 0000-00-00 00:00:00 |
    print GENE "$gene_id\t$gene_stable_id\t1\t0\t0\n";
    
    #######################
    # create transcript id
    #######################
    foreach my $trans ( @{$gene->get_all_Transcripts} ) {
      $count_transcript++;
      my $transcript_id = $trans->dbID;
      my $transcript_stable_id = &make_stable_id($label."T",$trans);
      #  927491 | ENSestT00000000001 |       1 
      print TRANSCRIPT "$transcript_id\t$transcript_stable_id\t1\n";
      
      #######################
      # create translation id
      #######################
      my $translation           = $trans->translation;
      my $translation_id        = $translation->dbID;
      my $translation_stable_id = &make_stable_id($label."P",$translation);
      #  927491 | ENSestP00000000001 |       1 
      print TRANSLATION "$translation_id\t$translation_stable_id\t1\n";
    }

    ####################
    # create exon ids 
    ####################
    my @exons = @{$gene->get_all_Exons};
    foreach my $exon ( @exons ){
      $count_exon++;
      my $exon_id = $exon->dbID;
      my $exon_stable_id = &make_stable_id($label."E",$exon);
      #  927491 | ENSestE00000000001 |       1 
      print EXON "$exon_id\t$exon_stable_id\t1\t0\t0\n";
    }
  }
}

close EXON;
close TRANSCRIPT;
close TRANSLATION;
close GENE;


print STDERR "Files gene_stable_id, transcript_stable_id, translation_stable_id , exon_stable_id written\n"; 
 
sub make_stable_id{
  my ($label,$object) = @_;
  # $label is ENSestG, etc...

  my $number = $object->dbID;
  my $string = "$number";
  while ( length($string) < 11 ){
    $string = "0".$string;
  }
  $string = $label.$string;
  return $string;
}
