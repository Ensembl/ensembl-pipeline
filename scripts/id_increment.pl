#!/usr/bin/env perl


# Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


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
