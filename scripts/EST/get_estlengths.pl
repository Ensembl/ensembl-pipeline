#!/usr/local/bin/perl -w
use strict;
$|=1;
use Bio::EnsEMBL::Pipeline::ESTConf qw (
					EST_FILE
				       );

my $estfile = $EST_FILE;
#my $estfile = 'test.fa';
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




