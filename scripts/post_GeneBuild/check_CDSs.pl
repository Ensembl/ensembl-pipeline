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
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;

my $file = 'ensembl_CDSs';

my $dbhost;
my $dbuser    = 'ensro';
my $dbname;
my $dbpass    = undef;

my $dnadbhost;
my $dnadbuser = 'ensro';
my $dnadbname;
my $dnadbpass = undef;

my $genetype;


&GetOptions(
	    'dbname:s'    => \$dbname,
	    'dbhost:s'    => \$dbhost,
	    'dnadbname:s' => \$dnadbname,
	    'dnadbhost:s' => \$dnadbhost,
	    'file:s'  => \$file,
	    'genetype:s'   => \$genetype,
	   );

unless ( $dbname && $dbhost && $dnadbname && $dnadbhost ){
  print STDERR "script to check whether CDSs starts with ATG and end with stop (TAA|TGA|TAG)\n";
 
  print STDERR "Usage: $0 -dbname -dbhost -dnadbname -dnadbhost\n";
  print STDERR "Optional: -genetype -file (defaulted to ensembl_cdnas)\n";
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

open (OUT,">$file") or die("unable to open file $file");

my $seqio = Bio::SeqIO->new('-format' => 'Fasta' , -fh => \*OUT ) ;

my  @ids = @{$db->get_GeneAdaptor->list_geneIds};

my $start_correct = 0;
my $stop_correct  = 0;
my $both_correct  = 0;

GENE:
foreach my $gene_id(@ids) {
  
  my $gene = $db->get_GeneAdaptor->fetch_by_dbID($gene_id,1);
  if ($genetype){
    next GENE unless ( $gene->type eq $genetype );
  }


  my $gene_id = $gene->dbID();
  my $chr = $gene->chr_name;

 TRANS:
  foreach my $trans ( @{$gene->get_all_Transcripts} ) {
    my $gene_id = $gene->stable_id || $gene->dbID;
    my $tran_id = $trans->stable_id || $trans->dbID;
    my @evidence = &get_evidence($trans);
    
    my $strand = $trans->start_Exon->strand;
    my ($start,$end);
    my @exons;
    if ( $strand == 1 ){
      @exons = sort {$a->start <=> $b->end} @{$trans->get_all_translateable_Exons};
      $start = $exons[0]->start;
      $end   = $exons[$#exons]->end;
    }
    else{
      @exons = sort {$b->start <=> $a->end} @{$trans->get_all_translateable_Exons};
      $start = $exons[0]->end;
      $end   = $exons[$#exons]->start;
    }
    
    eval {      
      my $seq;
      foreach my $exon ( @exons ){
	$seq .= $exon->seq->seq;
      }
      my $first_codon = substr( $seq, 0, 3 );
      my $last_codon = substr( $seq, -3 );

      my $start = 0;
      my $end   = 0;
      if ( uc($first_codon) eq 'ATG' ){
	$start_correct++;
	$start =1 ;
      }
      if ( uc($last_codon) eq 'TAA' || uc($last_codon) eq 'TAG' || uc($last_codon) eq 'TGA' ){ 
	$stop_correct++;
	$end = 1;
      }
      if ( $start && $end ){
	$both_correct++;
      }
      
      print "Gene:$gene_id Transcript:$tran_id start:$start stop:$end\n";
      
      #my $tseq = $trans->translate();
      #if ( $tseq->seq =~ /\*/ ) {
      #	print STDERR "translation of ".$trans->dbID." has stop codons. Skipping!\n";
      #	Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Evidence($trans);
      #	next TRANS;
      #      }
      #my $tran_seq = Bio::Seq->new();
      #$tran_seq->seq($seq);
      #$tran_seq->display_id("Gene:$gene_id Transcript:$tran_id CODING SEQUENCE");
      #$tran_seq->desc("HMM:@evidence Chr:$chr Strand:$strand Start:$start End:$end");
      #my $result = $seqio->write_seq($tran_seq);
    };
    if( $@ ) {
      print STDERR "unable to process transcript $tran_id, due to \n$@\n";
    }
  }
}

print "start codons correct: $start_correct\n";
print "stop codons correct : $stop_correct\n";
print "both correct        : $both_correct\n";


close (OUT);

sub get_evidence{
  my ($trans) = @_;
  my %evi;
  foreach my $exon (@{$trans->get_all_Exons}){
    foreach my $evidence ( @{$exon->get_all_supporting_features} ){
      $evi{$evidence->hseqname}++;
   }
  }
  return keys %evi;
}
