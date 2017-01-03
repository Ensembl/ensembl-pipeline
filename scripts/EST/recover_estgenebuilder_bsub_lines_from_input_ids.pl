#!/usr/bin/env perl


# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2017] EMBL-European Bioinformatics Institute
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

my $job_file;

GetOptions( 
	    'job_file:s'     => \$job_file,
	   );

unless ( $job_file ){
  print STDERR "$0 -job_file file_with_bsubs < input_ids\n";
  exit(0);
}


while(<>){
  chomp;
  my $input_id = $_;
  $input_id    =~ /(\S+\.\d+-\d+)/;
  #print $input_id."\n";
  
  #open (JOB,"<$job_file") or die("cannot open file $job_file");
  #while(<JOB>){
  #  chomp;
  #  my $line = $_;
  #  if ( $line =~/\s+$input_id/ ){
  #    print $line."\n";
  #  }
  #}
  #close (JOB);
  my $target = " $1";
    
  #my $command = "grep \"$target\" $job_file";
  #print $command."\n";;
  my $line = `grep "$target" $job_file`;
  print $line;

}
