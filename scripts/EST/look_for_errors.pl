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

open LS, 'ls |' or die "Can't open pipe from ls : $!";


FILES:
while(<LS>){
  chomp;
  my $file = $_;
  next FILES if ( $file =~/(\.)+\// );
  open(FILE, "<$file") or die "Can't open [$file]\n";
  
 THISFILE:
  while (<FILE>){
    chomp;
    my $line = $_;
    if ( $line =~ /glib/i     || 
	 $line =~ /ERROR/i    || 
	 $line =~ /exiting/i  ||
	 $line =~ /aborting/i || 
	 $line =~ /exit/i     || 
	 $line =~ /exception/i#|| 
	 #$line =~ /cannot/i 
       ){
      print $file."\n";
      print $line."\n\n";
      last THISFILE;
    }
  }
  close( FILE );
  #system("grep -i 'glib' $file > /work6a/eae.tmp/Mouse/jobs/error_file");
}
