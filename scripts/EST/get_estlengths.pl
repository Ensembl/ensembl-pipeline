#!/usr/local/bin/perl -w
use strict;
$|=1;
use Getopt::Long;

my $estfile;
&GetOptions(
	    'estfile:s' => \$estfile,
	   );

unless ( $estfile ){
  print STDERR "script to dump to STDOUT all the lengths of ests/cdnas\n";
  
  print STDERR "Usage: $0 -estfile est-file/cdna-file\n";
  exit(0);
}

print STDERR "estfle: $estfile\n";
open(EST, "<$estfile") or die "Can't open [$estfile]\n";



$/ = '>'; 
SEQFETCH: while(<EST>){
  my @lines = split /\n/;
  next SEQFETCH unless scalar(@lines) > 1;
  my $est_id = shift (@lines);

  $est_id =~ /^(\S+)/; # find the first block of non-whitespace
  $est_id = $1;

  if($est_id =~ /\S+\|\S+\|\S+\|(\S+)\|\S+/){
    $est_id =~ s/\S+\|\S+\|\S+\|(\S+)\|\S+/$1/; # this is for dbEST format headers
  }

  my $seq;
  foreach my $line(@lines) {
    chomp $line;
    $seq .= $line;
  }
  print "\\N\t$est_id\t" . length($seq) . "\n";

}

$/ = "\n";
close EST or die "Can't close $estfile\n";




