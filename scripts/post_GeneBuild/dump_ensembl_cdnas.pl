#!/usr/local/ensembl/bin/perl

=head1 NAME

=head1 DESCRIPTION

dumps in fastaA format the cdnas of all the genes in a database specified

=head1 OPTIONS

=cut

use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::SeqIO;
use Getopt::Long;

my $file = 'ensembl_cdnas';

my $dbhost    = 'ecs1e';
my $dbuser    = 'ensro';
my $dbname    = 'mouse_whitehead_0401_denormalised';
my $dbpass    = undef;

my $dnadbhost = 'ecs1e';
my $dnadbuser = 'ensro';
my $dnadbname = 'mouse_whitehead_0401_denormalised';
my $dnadbpass = undef;

my $genetype;


&GetOptions(
	    'dbname:s'    => \$dbname,
	    'dbhost:s'    => \$dbhost,
	    'dnadbname:s' => \$dnadbname,
	    'dnadbhost:s' => \$dnadbhost,
	    'cdna_file:s'  => \$file,
	    'genetype:s'   => \$genetype,
	   );

unless ( $dbname && $dbhost && $dnadbname && $dnadbhost && $genetype){
  print STDERR "script to dump all the cdnas from the transcripts in a database\n";
 
  print STDERR "Usage: $0 -genetype -dbname -dbhost -dnadbname -dnadbhost -cdna_file \n";
  exit(0);
}


my $dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					       '-host'   => $dnadbhost,
					       '-user'   => $dnadbuser,
					       '-dbname' => $dnadbname,
					       '-pass'   => $dnadbpass,
					      );


my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $dbhost,
					    '-user'   => $dbuser,
					    '-dbname' => $dbname,
					    '-pass'   => $dbpass,
					    '-dnadb'  => $dnadb,
					   );


print STDERR "connected to $dbname : $dbhost\n";
my $sa = $db->get_StaticGoldenPathAdaptor();

open (OUT,">$file") or die("unable to open file $file");

my $seqio = Bio::SeqIO->new('-format' => 'Fasta' , -fh => \*OUT ) ;

my  @ids = $db->get_GeneAdaptor->list_geneIds;
foreach my $gene_id(@ids) {

  my $gene = $db->get_GeneAdaptor->fetch_by_dbID($gene_id);
  
  my $gene_id = $gene->dbID();
  
  foreach my $trans ( $gene->each_Transcript ) {
    my $gene_stable_id = $gene->stable_id;
    my $tran_stable_id = $trans->stable_id;
    
    eval {
      
      my $tran_seq = $trans->seq;
      $tran_seq->display_id("$tran_stable_id $gene_stable_id");
      
      my $result = $seqio->write_seq($tran_seq);
    };
    if( $@ ) {
      print STDERR "unable to process transcript $tran_stable_id, due to \n$@\n";
    }
  }
}

close (OUT);

