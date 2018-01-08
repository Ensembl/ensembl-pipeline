#!/usr/bin/env perl


# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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
$|=1;
use Getopt::Long qw(:config no_ignore_case);

my $estfile;
GetOptions(
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




