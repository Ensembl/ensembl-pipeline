#!/usr/bin/env perl

use warnings ;
use strict;
use Getopt::Long qw(:config no_ignore_case);

my $file;
my @number;
my $backup = 1;
GetOptions(
            'file:s' => \$file,
	    'number_column:s@' => \@number,
	    'backup!' => \$backup,
	  ) or die("couldn't get opts");


if($backup){
  system("cp ".$file." ".$file.".bak");
}
open(FH, $file) or die "couldn't open $file";

my %replace;
foreach my $line(@number){
  my ($number, $column) = split /\:/, $line;
  $replace{$column} = $number;
}

my @keys = keys(%replace);
while(<FH>){
  chomp;
  my @values = split;
  foreach my $column(@keys){
    my $tmp = $values[$column];
    $values[$column] = $tmp + $replace{$column};
  }
  print join("\t", @values), "\n";
}

close(FH);
