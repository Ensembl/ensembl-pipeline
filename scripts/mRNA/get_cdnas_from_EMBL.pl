#!/usr/local/bin/perl -w

=head1 NAME

  get_human_ests.pl

=head1 SYNOPSIS
 
  get_human_ests.pl

=head1 DESCRIPTION

  gets human ESTs from dbEST or cDNAs from embl/vertRNA and polyT/polyA clips them

=head1 OPTIONS

  -estfile
  -outfile

=cut


use strict; 
use Getopt::Long;
use Bio::Seq;
use Bio::SeqIO;
use Bio::EnsEMBL::Utils::PolyA;

#$| = 1; # disable buffering
#local $/ = '>';

my $estfile;
my $seqoutfile;
my $clip;
my $softmask;

&GetOptions( 
	    'cdnafile:s'     => \$estfile,
	    'clip'          => \$clip,
	    'softmask'      => \$softmask,
	   );

# usage
if(!defined $estfile  ){
  print  "USAGE: get_human_ests.pl -cdnafile cdnafile -outfile outfile -clip -softmask\n";
  exit(1);
}

my $seqin  = new Bio::SeqIO(-file   => "<$estfile",
			    -format => "Fasta",
			  );

############################################################
# check the type of cdna you want to get:
open ( PSEUDO, ">embl_pseudogenes.fa")       or die("cannot open embl_peudogenes.fa");
open ( CODING, ">embl_coding_rnas.fa")       or die("cannot open embl_coding_rnas.fa");
open ( NONCODING, ">embl_noncoding_rnas.fa") or die("cannot open embl_noncoding_rnas.fa");
open ( PARTIAL,   ">embl_partial_cds.fa")    or die("cannot open embl_partial_cds.fa");
open ( REST, ">embl_the_rest_of_rnas.fa")    or die("cannot open embl_the_rest_of_rnas.fa");

my $pseudo_out = new Bio::SeqIO(-fh   => \*PSEUDO,
				-format => "Fasta"
			       );
my $coding_out = new Bio::SeqIO(-fh   => \*CODING,
				-format => "Fasta"
			       );
my $noncoding_out = new Bio::SeqIO(-fh   => \*NONCODING,
				   -format => "Fasta"
			       );
my $partial_out = new Bio::SeqIO(-fh   => \*PARTIAL,
				-format => "Fasta"
			       );
my $rest_out = new Bio::SeqIO(-fh  => \*REST,
				-format => "Fasta"
			     );

SEQFETCH:
while( my $cdna = $seqin->next_seq ){
  
  my $display_id    = $cdna->display_id;
  my $description   = $cdna->desc;
  
  my $is_pseudogene = 0;
  my $is_noncoding  = 0;
  my $is_codinggene = 0;
  my $is_partial    = 0;

  # First select the species:
  next SEQFETCH unless (   $description =~ /Homo sapiens/ || $description =~ /DNA.*coding.*human/ );

  next SEQFETCH if ( $description =~ /RIKEN/ 
		     && ( $description =~/Mus musculus/
			||
			  $description =~/mouse/i 
			)
		   );
  
  if(  $description =~ /similar to/ || $description =~ /homolog/i ){
    
    next SEQFETCH unless ( $description =~ /Homo sapiens.*similar to/ 
			   || $description =~ /Homo sapiens.*homolog/i 
			 );
  }
  
  # pseudogenes
  if ( $description =~ /pseudogene/i ){
    $is_pseudogene = 1;
  }
  
  # structural genes
  if ( ( $description =~/snoRNA/i 
	 || $description =~/small nucleolar/i 
	 || $description =~/tRNA/i
	 #|| $description =~/small nuclear/i
       ) 
       && 
       ! (( $description =~/complete cds/i) 
	  ||
	  ( $description =~/cDNA/i )
	  ||
	  ( $description =~/synthetase/i )
	  ||
	  ( $description =~/mRNA/i )
	  ||
	  ( $description =~/tRNA-associated/i )
	  ||
	  ( $description =~/transferase/i )
	  ||
	  ( $description =~/partial cds/i )
	  
	 )
     ){
    $is_noncoding = 1;
  }
  if ( $description =~/non-coding/i || $description =~/non coding/i ){
    $is_noncoding = 1;
  }
  if ( $description =~/partial cds/i || $description =~/partial/i ){
    $is_partial = 1;
  }

  if ( $description =~/complete cds/i || $description =~ /protein/i ){
    $is_codinggene = 1;
  }
  
  #print STDERR "keeping $description\n";
  
  # GenBank
  if ( $display_id =~/gi\|\S+\|\S+\|(\S+\.\d+)\|/ || $description =~/gi\|\S+\|\S+\|(\S+\.\d+)\|/ ){
    $display_id = $1;
  }
  # EMBL vert-RNA
  else{
    my @labels = split /\s+/, $description;
    $display_id = $labels[0];
  }
  
  $cdna->display_id($display_id);
  $cdna->desc("");
  
  if($@){
    warn("can't parse sequence for [$description]:\n$@\n");
    next SEQFETCH;
  }

  ############################################################
  #################### clipping? 
  my $polyA_clipper = Bio::EnsEMBL::Utils::PolyA->new();
  my $new_cdna;
  if ($clip){
    #print STDERR "going to pass a $cdna\n";
    $new_cdna = $polyA_clipper->clip($cdna);
  }
  elsif( $softmask ){
    $new_cdna = $polyA_clipper->mask($cdna, 'soft');
  }
  else{
    $new_cdna = $cdna;
  }
  ############################################################
  
  if ( $is_partial ){
    print STDERR "PARTIAL: $description\n";
    $partial_out->write_seq($new_cdna);
  }
  elsif ($is_pseudogene ){
    print STDERR "PSEUDOGENE: $description\n";
    $pseudo_out->write_seq($new_cdna);
  }
  elsif( $is_noncoding ){
    print STDERR "NONCODING: $description\n";
    $noncoding_out->write_seq($new_cdna);
  }
  elsif( $is_codinggene ){
    print STDERR "CODING: $description\n";
    $coding_out->write_seq($new_cdna);
  }
  else{
    print STDERR "REST: $description\n";
    $rest_out->write_seq($new_cdna);
  }
}

close(REST);
close(PARTIAL);
close(CODING);
close(NONCODING);
close(PSEUDO);

