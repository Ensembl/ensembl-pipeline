#!/usr/local/bin/perl -w

=head1 NAME

  get_human_ests.pl

=head1 SYNOPSIS
 
  get_human_ests.pl

=head1 DESCRIPTION

  gets human ESTs from dbEST and polyT/polyA clips them ready for exonerate

=head1 OPTIONS

  -estfile
  -outfile

=cut


use strict; 
use Getopt::Long;
use Bio::Seq;
use Bio::SeqIO;

$| = 1; # disable buffering
local $/ = '>';

my $estfile;
my $seqoutfile;

&GetOptions( 
	    'estfile:s'     => \$estfile,
	    'outfile:s'     => \$seqoutfile,
	   );

# usage
if(!defined $estfile    ||
   !defined $seqoutfile 
  ){
  print  "USAGE: get_human_ests.pl -estfile estfile -outfile outfile\n";
  exit(1);
}

my $seqout = new Bio::SeqIO(-file => ">>$seqoutfile", "-format" => "Fasta");

open(EST, "<$estfile") or die "Can't open estfile [$estfile]\n";

SEQFETCH:
while (<EST>){
  my @lines = split /\n/;
  next SEQFETCH unless scalar(@lines) > 1;
  my $description = shift (@lines);
  next SEQFETCH unless $description =~ /Homo sapiens cDNA/;
  if($description =~ /similar to .* Homo sapiens cDNA/){
    next SEQFETCH unless $description =~ /Homo sapiens.*similar to.*Homo sapiens cDNA/;
    print STDERR "keeping $description\n";
  }

  my $cdna = new Bio::Seq;

  my $seq;
  foreach my $line(@lines) {
    chomp $line;
    $seq .= $line;
  }

  if (($description =~ /3\'/) && ($seq =~ /^T{3,}/)) {
    $seq =~ s/^T{3,}//;
    $description .= " polyTT"; 
  }
  else {
    $description .= " polyNO"; 
  }


  eval{
    $cdna->seq($seq);
  };

  if($@){
    warn("can't parse sequence for [$description]:\n$@\n");
    next SEQFETCH;
  }

  # modify description
  $cdna->display_id($description);

  # write sequence
  $seqout->write_seq($cdna);
  
}

close EST or die "Can't close estfile [$estfile]\n";


