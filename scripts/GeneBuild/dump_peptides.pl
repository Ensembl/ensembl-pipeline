#!/usr/local/ensembl/bin/perl

=head1 NAME

  dump_translations.pl

=head1 DESCRIPTION

dump_peptides.pl dumps in fastaA format the translations of all the genes in a database specified in GeneConf.pm
It\'s a stripped down version of gene2flat.

Usage: dump_peptides.pl > ! peptide_file

it will output on STDERR the gene ids being dumped
It checks for '*' in the translation


=head1 OPTIONS

=cut

use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::SeqIO;
use Getopt::Long;

my $file;

my $dbhost;
my $dbuser    = 'ensro';
my $dbname;
my $dbpass    = undef;

my $dnadbhost;
my $dnadbuser = 'ensro';
my $dnadbname;
my $dnadbpass = undef;

my $genetype = "ensembl"; # default genetype


$dbuser = "ensro";
&GetOptions(
	    'dbname:s'    => \$dbname,
	    'dbhost:s'    => \$dbhost,
	    'dnadbname:s' => \$dnadbname,
	    'dnadbhost:s' => \$dnadbhost,
	    'peptide_file:s'  => \$file,
);

unless ( $dbname && $dbhost && $dnadbname && $dnadbhost && $genetype){
  print STDERR "script to dump all the translations from the transcripts in a database\n";
  print STDERR "Usage: $0 --dbname -dbhost -dnadbname -dnadbhost -peptide_file \n";
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
#golden path type is now obtained from the meta table
#$db->static_golden_path_type($path);
my $sa = $db->get_StaticGoldenPathAdaptor();


open (OUT,">$file") or die("unable to open file $file");

my $seqio = Bio::SeqIO->new('-format' => 'Fasta' , -fh => \*OUT ) ;

foreach my $gene_id($db->get_GeneAdaptor->list_geneIds) {
  print STDERR "gene id $gene_id\n";
  eval {
    my $gene = $db->get_GeneAdaptor->fetch_by_dbID($gene_id);
    my $gene_id = $gene->dbID();
    foreach my $trans ( $gene->each_Transcript ) {
      # get out first exon. Tag it to clone and gene on this basis
      my @exon = $trans->get_all_Exons;
      my $fe = $exon[0];
      
      my ($chr,$gene_start,$cdna_start) = find_trans_start($trans);
      
      my $tseq = $trans->translate();
      #print STDERR "translation->translate is a : $tseq\n";
      #print STDERR "translation: ".$tseq->seq()."\n";
      if ( $tseq->seq =~ /\*/ ) {
	print STDERR "translation of transcript: ".$trans->dbID." in chr $chr has stop codons. Skipping! (in clone". $fe->clone_id .")\n";
	next;
      }
      my $gene_version = $gene->version;
      
      $tseq->desc("Gene:$gene_id.$gene_version Clone:".$fe->clone_id . " Contig:" . $fe->contig_id . " Chr: " . $chr . " Pos: " . $cdna_start);
      my $result = $seqio->write_seq($tseq);
    }
  };
  
  if( $@ ) {
    print STDERR "unable to process $gene_id, due to \n$@\n";
    }
}

close (OUT);

sub  find_trans_start {
 my ($trans) = @_;

 my $start_pos;
 my $trans_pos;

 my $contig; 
 foreach my $exon ($trans->get_all_Exons) {
   if ($exon eq $trans->translation->start_exon) {
     $contig = $exon->contig_id;
     if ($exon->strand == 1) {
       $start_pos = $exon->start;
       $trans_pos = $exon->start + $trans->translation->start - 1;
     } else {
       $start_pos = $exon->end;
       $trans_pos = $exon->end - $trans->translation->start + 1;
     }
   }
 }
 if (!defined($start_pos)) {
   print STDERR "Couldn't find start exon for " . $trans->id . "\n";
   die;
 }

 my $query = "select chr_name,chr_start,chr_end,raw_start,raw_end,raw_ori from static_golden_path where raw_id = $contig";

 my $sth  = $db->prepare($query);
 my $res  = $sth->execute;
 
 my $row = $sth->fetchrow_hashref;
 
 my $chr = $row->{chr_name};
 my $chr_start = $row->{chr_start};
 my $chr_end   = $row->{chr_end};
 my $raw_start = $row->{raw_start};
 my $raw_end   = $row->{raw_end}; 
 my $raw_ori   = $row->{raw_ori};
 
 my $gene_start; 
 my $cdna_start;
 
 if ($raw_ori == 1) {
   $gene_start = $chr_start + ($trans_pos - $raw_start);
   $cdna_start = $chr_start + ($start_pos - $raw_start);
 } else {
   $cdna_start = $chr_end   - ($start_pos - $raw_start);
   $gene_start = $chr_end   - ($trans_pos - $raw_start);
 } 
 
 return ($chr,$gene_start,$cdna_start);
}

