#!/usr/local/ensembl/bin/perl -w

# script to parse the header of refseq/embl file with cdnas
# it attempts to clip the polyA/polyT

use strict; 
use Getopt::Long;
use Bio::Seq;
use Bio::SeqIO;

$| = 1; # disable buffering
local $/ = '>';

my $rnafile;
my $seqoutfile;
my $clip;

&GetOptions( 
	    'mRNAfile:s'    => \$rnafile,
	    'outfile:s'     => \$seqoutfile,
	    'clip'          => \$clip,
	   );

# usage
if(!defined $rnafile    ||
   !defined $seqoutfile 
  ){
  print  "USAGE: get_human_ests.pl -mRNAfile rnafile -outfile outfile -clip\n";
  exit(1);
}

# fasta format
my $seqout = new Bio::SeqIO(-file => ">$seqoutfile", "-format" => "Fasta");

open(RNA, "<$rnafile") or die "Can't open rnafile [$rnafile]\n";

my %id_hash;

SEQFETCH:
while (<RNA>){
  my @lines = split /\n/;
  next SEQFETCH unless scalar(@lines) > 1;
  my $desc = shift (@lines);

  $desc      =~ s/^\>//;
  my @fields = split /\s+/, $desc;

  # is it a refseq?
  # >gi|8924222|ref|NM_018546.1| Homo sapiens hypothetical protein PRO2958 (PRO2958), mRNA
  
  if ( $fields[0] =~/\|ref\|(\S+\.\d+)\|/ ){
    $fields[0] = $1;
  }
  # else it might be embl:
  elsif ( $fields[0] =~/ri\|(\S+)\|(\S+)\|/ ){
    $fields[0] = $2;
  }
  else{
    $fields[0] = $fields[1];
  }
  
  # do not use the gene loci, i.e. the NG_ entries of RefSeq
  if ($fields[0] =~/^NG/ ){
    print STDERR "rejecting entry $fields[0]\n";
    next;
  }
  
  #print STDERR $fields[0]."\n";
  
  my $cdna = new Bio::Seq;

  my $seq;
  foreach my $line(@lines) {
    chomp $line;
    $seq .= $line;
  }

  #################### clipping? 
  my $clipped_seq;
  if ($clip){
    
    my $length = length($seq);
    print STDERR "seq length: $length\n";

    # is it a polyA or polyT?
    my $check_polyT = substr( $seq, 0, 10 );

    my $check_polyA = substr( $seq, -10 );
        
    my $t_count = $check_polyT =~ tr/Tt//;
    print STDERR "$t_count Ts in head\n";
    print STDERR "head: $check_polyT\n";
    
    my $a_count = $check_polyA =~ tr/Aa//;
    print STDERR "$a_count As in tail\n";
    print STDERR "tail: $check_polyA\n";


    #### polyA ####
    if ( $a_count >= 7 && $a_count > $t_count ){
           
      print STDERR "clipping PolyA tail:\n";
      # we calculate the number of bases we want to chop
      my $length_to_chop = 0;
      
      # we start with 3 bases
      my ($piece, $count ) = (3,0);

      # take 3 by 3 bases from the end
      while( $length_to_chop < $length ){
	my $chunk  = substr( $seq, ($length - ($length_to_chop + 3)), $piece);
	print STDERR "chunk: $chunk\n";
	$count = $chunk =~ tr/Aa//;
	if ( $count >= 2*( $piece )/3 ){
	  $length_to_chop += 3;
	  print STDERR "continue\n";
	}
	else{
	  print STDERR "last\n";
	  last;
	}
      }
      
      # do not chop the last base if it is not an A:
      my $last_base = substr( $seq, ($length - $length_to_chop), 1);
      unless ( $last_base eq 'A' || $last_base eq 'a' ){
	$length_to_chop--;
      }
      $clipped_seq = substr( $seq, 0, $length - $length_to_chop );
      
    }
    #### polyT ####
    elsif( $t_count >=7 && $t_count > $a_count ){
      
      # calculate the number of bases to chop
      my $length_to_chop = -3;

      # we start with 3 bases:
      my ($piece, $count) = (3,3);

      # take 3 by 3 bases from the beginning
      while ( $length_to_chop < $length ){
	my $chunk = substr( $seq, $length_to_chop + 3, $piece );
	print STDERR "chunk ".( $length_to_chop + 3)." : $chunk\n";
	$count = $chunk =~ tr/Tt//;
	if ( $count >= 2*( $piece )/3 ){
	  $length_to_chop +=3;
	  print STDERR "continue\n";
	}
	else{
	  print STDERR "last\n";
	  last;
	  
	}
      }
      if ( $length_to_chop > 0 ){
	# do not chop the last base if it is not a A:
	my $last_base = substr( $seq, ( $length_to_chop -1 ), 1 );
	unless ( $last_base eq 'T' || $last_base eq 't' ){
	  $length_to_chop--;
	}
	$clipped_seq = substr( $seq, $length_to_chop + 3);
      }
      else{
	$clipped_seq = $seq;
      }
    }
    else{
      # we cannot be sure of what it is
      # do not clip
      $clipped_seq = $seq;
    }
  }
  ### else, do not clip
  else{
    $clipped_seq = $seq;
  }

  ####################
  eval{
    $cdna->seq($clipped_seq);
    $cdna->display_id($fields[0]);
    $cdna->desc("");
  };
  
  if($@){
    warn("can't parse sequence for [$desc]:\n$@\n");
    next SEQFETCH;
  }
  
  # write sequence
  $seqout->write_seq($cdna);
  
}

close RNA or die "Can't close rnafile [$rnafile]\n";


