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
my $min_length = 60;

&GetOptions( 
	    'estfile:s'     => \$estfile,
	    'outfile:s'     => \$seqoutfile,
	    'clip'          => \$clip,
	    'softmask'      => \$softmask,
	    'min_length:s'    => \$min_length,
	   );

# usage
if(!defined $estfile    ||
   !defined $seqoutfile 
  ){
  print  "USAGE: get_human_ests.pl -estfile estfile -outfile outfile".
    " -min_length <min_est_length> -clip -softmask\n";
  exit(1);
}

my $seqin  = new Bio::SeqIO(-file   => "<$estfile",
			    -format => "Fasta",
			  );

my $seqout = new Bio::SeqIO(-file   => ">$seqoutfile", 
			    -format => "Fasta"
			   );

SEQFETCH:
while( my $cdna = $seqin->next_seq ){

  next unless $cdna->length > $min_length;
  
  my $display_id  = $cdna->display_id;
  my $description = $cdna->desc;

  # First select the species:
  next SEQFETCH unless (   $description =~ /Homo sapiens/
			   || $description =~ /DNA.*coding.*human/
		       );
  
  if(  $description =~ /similar to/ || $description =~ /homolog/i ){
    
    next SEQFETCH unless ( $description =~ /Homo sapiens.*similar to/ 
			   || $description =~ /Homo sapiens.*homolog/i 
			 );
  }
  #print STDERR "keeping $description\n";
  
  # GenBank
  if ( $display_id =~/gi\|\S+\|\S+\|(\S+\.\d+)\|/ || $description =~/gi\|\S+\|\S+\|(\S+\.\d+)\|/ ){
    $display_id = $1;
    $cdna->display_id($display_id);
  }
  
  $cdna->desc("");
  
  if($@){
    warn("can't parse sequence for [$description]:\n$@\n");
    next SEQFETCH;
  }

   #################### clipping? 
  my $polyA_clipper = Bio::EnsEMBL::Utils::PolyA->new();
  my $new_cdna;
  if ($clip){
#    print STDERR "going to pass a $cdna\n";
    $new_cdna = $polyA_clipper->clip($cdna);
  }
  elsif( $softmask ){
    $new_cdna = $polyA_clipper->mask($cdna, 'soft');
  }
  else{
    $new_cdna = $cdna;
  }

  # write sequence
  $seqout->write_seq($new_cdna);
}

