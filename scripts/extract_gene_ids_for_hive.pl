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


# Usage: ./extract_gene_ids_for_hive.pl -analysis_id <number - optional>
#                                       -status <status string - optional> 
#                                       -file <filename - fasta file from blast db>
#                                       > <output file>

use warnings ;
use strict;
use Getopt::Long qw(:config no_ignore_case);

my ($analysis_id, $status, $filename);

GetOptions("-analysis_id" => \$analysis_id,
	   "-status=s"    => \$status,
	   "-file=s"      => \$filename);

unless (defined $filename && $filename ne ''){
  print STDERR join ("\n", 
		     "Command-line options : ",
		     "-file          input file (required)",
		     "-status        hive status (optional - default is READY)",
		     "-analysis_id   hive analysis id (optional - default 1)",
		     "Output is written to STDOUT, so you might like to ",
		     "redirect it to a file.") . "\n";
  die "Command line options not specified."
}

$analysis_id = 1 unless $analysis_id;
$status = 'READY' unless $status;

$status = uc($status);

open(IN, $filename) or die "Unable to find or open [$filename]";

while (<IN>) {
  if (/>(\S+)/) {
    print STDOUT join("\t", ('','',$analysis_id,$1,'','',$status,'','','','','')) . "\n";
  }
}

close(IN);
