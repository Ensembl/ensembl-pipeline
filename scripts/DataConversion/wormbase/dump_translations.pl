#!/usr/local/bin/perl

=head1 NAME

  dump_translations.pl

=head1 SYNOPSIS
 
  dump_translations.pl

=head1 DESCRIPTION

dump_translations.pl dumps out the translations of all the genes in a database specified in GeneConf.pm
It\'s a stripped down version of gene2flat.

=head1 OPTIONS

=cut

use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::SeqIO;


my $dbhost      = 'ecs2d';
my $dbuser    = 'ensro';
my $dbname    = 'caenorhabditis_elegans_core_10_93';
my $dbpass    = undef;







my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $dbhost,
					    '-user'   => $dbuser,
					    '-dbname' => $dbname,
					    '-pass'   => $dbpass,
					    
					   );


print STDERR "connected to $dbname : $dbhost\n";


my $seqio = Bio::SeqIO->new('-format' => 'Fasta' , -fh => \*STDOUT ) ;

foreach my $gene_id(@{$db->get_GeneAdaptor->list_geneIds}) {
  
  eval {
    my $gene = $db->get_GeneAdaptor->fetch_by_dbID($gene_id);
    my $gene_id = $gene->dbID();
    
    foreach my $trans ( @{$gene->get_all_Transcripts}) {
      # get out first exon. Tag it to clone and gene on this basis
      my @exon = @{$trans->get_all_Exons};
      my $fe = $exon[0];
      
      my ($chr,$gene_start,$cdna_start) = find_trans_start($trans);
      
      my $tseq = $trans->translate();
      if ( $tseq->seq =~ /\*/ ) {
	print STDERR "translation of ".$trans->id." has stop codons. Skipping! (in clone". $fe->contig->dbID .")\n";
	next;
      }
      my $gene_version = $gene->version;
      my $trans_id = $trans->translation->dbID;
      $tseq->display_id($trans_id);
      $tseq->desc("Translation id ".$trans_id." gene $gene_id Contig:" .$fe->contig->dbID. " Chr: " . $chr . " Pos: " . $cdna_start."\n");
      $seqio->write_seq($tseq);
    }
  };
  
  if( $@ ) {
    print STDERR "unable to process $gene_id, due to \n$@\n";
    }
}

sub  find_trans_start {
 my ($trans) = @_;
 #print STDERR "finding trans start\n";
 my $start_pos;
 my $trans_pos;

 my $contig; 
 foreach my $exon (@{$trans->get_all_Exons}) {
   if ($exon eq $trans->translation->start_Exon) {
     $contig = $exon->contig->dbID;
     if ($exon->strand == 1) {
       $start_pos = $exon->start;
       $trans_pos = $exon->start + $trans->translation->start - 1;
     } else {
       $start_pos = $exon->end;
       $trans_pos = $exon->end - $trans->translation->start + 1;
     }
    # print STDERR "have start ".$start_pos." trans ".$trans_pos."\n";
   }
 }
 if (!defined($start_pos)) {
   print STDERR "Couldn't find start exon for " . $trans->stable_id . "\n";
   die;
 }

 my $query = "select c.name, a.chr_start,a.chr_end,a.contig_start,a.contig_end,a.contig_ori from assembly a, chromosome c where contig_id = $contig and c.chromosome_id = a.chromosome_id";

 
 my $sth  = $db->prepare($query);
 my $res  = $sth->execute;
 
 my $row = $sth->fetchrow_hashref;
 
 my $chr = $row->{name};
 my $chr_start = $row->{chr_start};
 my $chr_end   = $row->{chr_end};
 my $raw_start = $row->{contig_start};
 my $raw_end   = $row->{contig_end}; 
 my $raw_ori   = $row->{contig_ori};
 
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

