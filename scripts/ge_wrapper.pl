#!/usr/bin/env perl


# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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


# A simple wrapper script for running jobs through GridEngine.
#

use warnings ;
use strict;  

print STDERR "Running on machine " . `hostname` . "\n"; 

my $ok = 1;
if (exists($ENV{PREEXEC})) {
  print "Pre exec = " . $ENV{PREEXEC} . "\n";
  if (system($ENV{PREEXEC})) {
    print STDERR "Failed pre exec " . $ENV{PREEXEC} . "\n";
    $ok = 0;
  }
}
if ($ok) {
  foreach my $arg (@ARGV) {
    print "Arg = $arg\n";
  }
  if (system(@ARGV)) {
    print STDERR "system @ARGV failed: $?";
  }
} else {
  print "Didn't run job\n";
}

print "Output file: " . $ENV{FINAL_STDOUT} . "\n";
print "Error  file: " . $ENV{FINAL_STDERR} . "\n";
print "Output path is " . $ENV{SGE_STDOUT_PATH} . "\n";
print "Error path is " . $ENV{SGE_STDERR_PATH} . "\n";
my $cmd = "rcp " . $ENV{SGE_STDOUT_PATH} ." " . $ENV{SGE_O_HOST} . ":" . $ENV{FINAL_STDOUT};
system ($cmd);
my $cmd = "rcp " . $ENV{SGE_STDERR_PATH} ." " . $ENV{SGE_O_HOST} . ":" .$ENV{FINAL_STDERR};
system ($cmd);

unlink $ENV{SGE_STDERR_PATH};
unlink $ENV{SGE_STDOUT_PATH};
