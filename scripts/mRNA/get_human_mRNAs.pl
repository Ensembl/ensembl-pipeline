#!/usr/local/ensembl/bin/perl -w

=head1 NAME

  get_human_mRNAs.pl


=head1 SYNOPSIS
 
  get_human_mRNAs.pl -mrnafile rnafile -ouputfile outfile


=head1 DESCRIPTION

  gets human mRNAs from embl_vertrna and polyT/polyA clips them ready for exonerate

=head1 OPTIONS

  -mRNAfile
  -outfile

=cut


use strict; 
use Getopt::Long;
use Bio::Seq;
use Bio::SeqIO;

$| = 1; # disable buffering
local $/ = '>';

my $rnafile;
my $seqoutfile;

&GetOptions( 
	    'mRNAfile:s'     => \$rnafile,
	    'outfile:s'     => \$seqoutfile,
	   );

# usage
if(!defined $rnafile    ||
   !defined $seqoutfile 
  ){
  print  "USAGE: get_human_ests.pl -mRNAfile rnafile -outfile outfile\n";
  exit(1);
}

# append in fasta format
my $seqout = new Bio::SeqIO(-file => ">>$seqoutfile", "-format" => "Fasta");

open(RNA, "<$rnafile") or die "Can't open rnafile [$rnafile]\n";

my %id_hash;

SEQFETCH:
while (<RNA>){
  my @lines = split /\n/;
  next SEQFETCH unless scalar(@lines) > 1;
  my $desc = shift (@lines);
  next SEQFETCH unless $desc =~ /Homo sapiens/;
  
  ## if you want to avoid pseudogenes, repeats and retroviruses, uncomment the next line:
  #next SEQFETCH if ( $desc =~ /retrovirus/i || $desc =~ /repeat/i || $desc =~ /pseudogene/i );

  ## there are more mRNAs than those just labelled like RNA,cDNA or CDS, so this can be too restrictive:
  #next SEQFETCH unless ( $desc = ~/RNA/ || $desc =~ /cDNA/ || $desc =~ /cds/i );
  
  if($desc =~ /similar to .* Homo sapiens/i){
    next SEQFETCH unless $desc =~ /Homo sapiens.*similar to.*Homo sapiens/i;
    #print STDERR "keeping $desc\n";
  }

  $desc =~ s/^\>//;
  my @fields = split /\s+/, $desc;
  $fields[1] =~ s/\.\d+$//;
  #if ( !( $fields[0] eq  $fields[1]) ){
  #  print STDERR "$fields[0] not equal $fields[1]\n";
  #}
  if ( exists $id_hash{$fields[1]} ){
    print STDERR "Id: $fields[1] seen more than once\n";
  }

  $id_hash{$fields[1]}=1;

  my $cdna = new Bio::Seq;

  my $seq;
  foreach my $line(@lines) {
    chomp $line;
    $seq .= $line;
  }

## We don't cut the poly A tail if present
#
# if (($desc =~ /3\'/) && ($seq =~ /^T{3,}/)) {
#    $seq =~ s/^T{3,}//;
#    $desc .= " polyTT"; 
#  }
#  else {
#    $desc .= " polyNO"; 
#  }


  eval{
    $cdna->seq($seq);
  };

  if($@){
    warn("can't parse sequence for [$desc]:\n$@\n");
    next SEQFETCH;
  }

  # modify description
  $cdna->display_id($fields[1]);
  $cdna->desc("");
  

  # write sequence
  $seqout->write_seq($cdna);
  
}

close RNA or die "Can't close rnafile [$rnafile]\n";


