#!/usr/local/ensembl/bin/perl -w

###############################
#
# Modification of the original id_increment.pl
# to be able to use ids with alfanumeric characters
# i.e. ENSHSG00000000012
#
# exon_id stable_id       version created_date    modified_date
# 1       ENSCFE00000000001       1       2005-11-14 15:57:05     2005-11-14 15:57:05
# 2       ENSCFE00000000002       1       2005-11-14 15:57:05     2005-11-14 15:57:05
# 3       ENSCFE00000000003       1       2005-11-14 15:57:05     2005-11-14 15:57:05
#
# perl alfanumeric_id_increment.pl -file exon_stable_id.txt -outfile exon_stable_id.table -number_column "100:0 100:1";
#  
# exon_id stable_id       version created_date    modified_date
# 101       ENSCFE00000000101       1       2005-11-14 15:57:05     2005-11-14 15:57:05
# 102       ENSCFE00000000102       1       2005-11-14 15:57:05     2005-11-14 15:57:05
# 103       ENSCFE00000000103       1       2005-11-14 15:57:05     2005-11-14 15:57:05 
#
###############################

use strict;
use Getopt::Long;

my $file;
my $numbers;
my $backup = 1;
my $outfile;

&GetOptions(
            'file=s'           => \$file,
	    'number_column:s' => \$numbers,
	    'backup!'          => \$backup,
	  ) or die("couldn't get opts");


if($backup){
  system("cp ".$file." ".$file.".bak");
}

$outfile = $file.".tmp";

open(FH, $file) or die "couldn't open $file";
open(OUT, ">>$outfile");

my %replace;
print $numbers,"\n";
my @number = split /\s/,$numbers;
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
    my ($prefix,$suffix) = $tmp =~ /([a-zA-Z]*)([0-9]+)/;
    if($suffix){
      my $new_id = $suffix + $replace{$column};
      if ($prefix){
        $values[$column] = sprintf "%s%011d",$prefix,$new_id;
      }else{
        $values[$column] = $new_id;
      }
    }
  }
  print join("\t", @values), "\n";
  print OUT join("\t", @values), "\n";

}

close(FH);
close(OUT);

system("mv ".$outfile." ".$file);

