#!/usr/local/ensembl/bin/perl


=head1 NAME

  dump_translations.pl

=head1 SYNOPSIS
 
  dump_translations.pl

=head1 DESCRIPTION

dump_translations.pl dumps out the translations of all the genes of a given gene type in a db.
This info is specified in the command line

=head1 OPTIONS

=cut

use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Bio::SeqIO;
use Getopt::Long;

my $dbhost;
my $dbuser    = 'ensro';
my $dbname;
my $dbpass    = undef;

my $dnadbhost;
my $dnadbuser = 'ensro';
my $dnadbname;
my $dnadbpass = undef;

my $genetype;
my $peptide_file = 'peptides.fa';

&GetOptions( 'dbhost:s'       => \$dbhost,
	     'dbname:s'       => \$dbname,
	     'dnadbhost:s'    => \$dnadbhost,
	     'dnadbname:s'    => \$dnadbname,
	     'peptide_file:s' => \$peptide_file,
	     'genetype:s'     => \$genetype,
	     );

unless ( $dbhost && $dbname && $dnadbname && $dnadbhost ){
    print STDERR "Usage: $0 -dbhost -dbname -dnadbhost -dnadbname > peptide_file\n";
    print STDERR "Optional: -peptide_file -genetype\n";
    print STDERR "peptide_file is defaulted to peptides.fa";
    print STDERR "if no genetype is specified all genetypes are retrieved\n";
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


print STDERR "connected to $dbname:$dbhost and $dnadbname:$dnadbhost\n";

open(OUT,">$peptide_file") or die("cannot open file $peptide_file");
my $seqio = Bio::SeqIO->new('-format' => 'Fasta' , -fh => \*OUT ) ;

my $gene_adaptor = $db->get_GeneAdaptor;

foreach my $gene_id(@{$db->get_GeneAdaptor->list_geneIds}) {
    eval {
      my $gene = $gene_adaptor->fetch_by_dbID($gene_id);
      if ( $genetype ){
	next unless ( $genetype eq $gene->type );
      }
      my $id = $gene->stable_id || $gene->dbID;
      print STDERR "fetching gene id $gene_id\n";
	
      TRANS:
	foreach my $trans ( @{$gene->get_all_Transcripts} ) {
	    my $t_id = $trans->stable_id || $trans->dbID;
	    unless ( defined $trans->translation ){
		print STDERR "transcript ".$t_id." has no translation. Skipping it.\n";
	      Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Evidence($trans);
	    }
	    
	    # get out first exon. Tag it to clone and gene on this basis
	    #my @exon = $trans->get_all_Exons;
	    #my $fe = $exon[0];
	    
	    #my ($chr,$gene_start,$cdna_start) = find_trans_start($trans);
	    
	    my $tseq = $trans->translate();
	    if ( $tseq->seq =~ /\*/ ) {
		print STDERR "translation of ".$trans->dbID." has stop codons. Skipping!\n";
	      Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Evidence($trans);
		next TRANS;
	    }
	    #my $gene_version = $gene->version;
	    
	    my $tr_id = $trans->translation->stable_id || $trans->translation->dbID;
	    $tseq->display_id("Translation:$tr_id");
	    $tseq->desc("Gene:$id Transcript:$t_id");
#$tseq->desc("Gene:$id.$gene_version Clone:".$fe->clone_id . " Contig:" . $fe->contig_id . " Chr: " . $chr . " Pos: " . $cdna_start);
	    $seqio->write_seq($tseq);
	}
    };
    
    if( $@ ) {
	print STDERR "unable to process $gene_id, due to \n$@\n";
    }
}

close(OUT);

sub  find_trans_start {
 my ($trans) = @_;
 
 my $start_pos;
 my $trans_pos;
 
 my $contig; 
 foreach my $exon (@{$trans->get_all_Exons}) {
     if ($exon == $trans->translation->start_Exon) {
	 $contig = $exon->contig;
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
     print STDERR "Couldn't find start exon for transcript dbID ". $trans->dbId . "\n";
     return;
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


