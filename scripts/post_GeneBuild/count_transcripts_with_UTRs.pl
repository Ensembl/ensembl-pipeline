#!/usr/local/ensembl/bin/perl -w

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::SeqFetcher::OBDAIndexSeqFetcher;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Getopt::Long;
use strict;

my $dbname;
my $dbhost;
my $dbuser = 'ensro';

my $genetype;

&GetOptions( 'dbhost:s'       => \$dbhost,
	     'dbname:s'       => \$dbname,
	     'genetype:s'     => \$genetype,
	   );

unless ( $dbhost && $dbname ){
  print STDERR "Usage: $0 -dbhost -dbname > & log_file\n";
  print STDERR "Optional: -genetype\n";
  exit(0);
}

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $dbhost,
					    -user   => $dbuser,
					    -dbname => $dbname,
					   );



print STDERR "connected to $dbname : $dbhost\n";
my $sa = $db->get_SliceAdaptor;

my  @ids = @{$db->get_GeneAdaptor->list_geneIds};

my $only_utr_three = 0;
my $only_utr_five  = 0;
my $both = 0;
my $no_utr = 0;
my $total;

GENE:
foreach my $gene_id( @ids) {
    
  # this gives the gene in chromosome coordinates:
  my $gene    = $db->get_GeneAdaptor->fetch_by_dbID( $gene_id, 1 );
  my $gene_id = $gene->stable_id || $gene->dbID;

  if ($genetype){
    next GENE unless( $genetype eq $gene->type );
  }

 TRANS:
  foreach my $trans ( @{$gene->get_all_Transcripts} ){
    
    $total++;

    my $five_seq  = $trans->five_prime_utr->seq;
    my $three_seq = $trans->three_prime_utr->seq;

    if ( $five_seq && $three_seq ){
      $both++;
      next TRANS;
    }
    
    if( $five_seq && !$three_seq ){
      $only_utr_five++;
      next TRANS;
    }

    if( !$five_seq && $three_seq ){
      $only_utr_three++;
      next TRANS;
    }

    $no_utr++;
  }
}

print STDERR "total transcripts: ".$total."\n";
print STDERR "no UTRs          : ".$no_utr."\n";
print STDERR "with 5'UTR only  : ".$only_utr_five."\n";
print STDERR "with 3'UTR only  : ".$only_utr_three."\n";
print STDERR "with both UTRs   : ".$both."\n";

############################################################
