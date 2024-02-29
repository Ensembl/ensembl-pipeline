#!/usr/bin/env perl
use warnings ;
use strict;


# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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


#script to parse a fasta file and identify sequences with large strings of 'N's

#perl find_N.pl missing_fasta.out >many_n.out

my $percent = 2; #percentage of sequence which must be consecutive Ns
my $total_percent = 5 * $percent; #total Ns 

my $data = $ARGV[0];
my $a_count = 0;

local $/ = "\n>";

open(DATA, "<$data") or die ("Can't read $data $! \n");

while(<DATA>){ 
	#have a sequence:
	
	s/>//g;
	
	my $len = length $_;
	my $max_n = sprintf "%.0f", (($len / 100) * $percent); #threshold number of Ns which we want to flag 
	my $percent_n = 0;
	
	my ($name, $seq);
	if ($_=~/^([\w\.]+)\s+([\w\s]+)/){
		$name = $1;
		my @tmp = $2;
		
		for my $s (@tmp){
			$s =~s/\s//g;
		}
		$seq = join "", @tmp;
	}
	
	while ($_=~/(N+)/g){ #will match greedily
		$percent_n += length $1;
		if (length $1 >= $max_n && $percent_n >= $total_percent){
			
			print "$name\n"; #print the seq id
			last;
		}
	}
}

	
