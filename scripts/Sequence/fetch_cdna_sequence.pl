#!/usr/bin/env perl


# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016] EMBL-European Bioinformatics Institute
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


# dump_seq_into_fastA.pl
# it reads a bit of sequence and dump it into a fasA file, to eb able to view it in Apollo

use warnings ;
use strict;
use diagnostics;

use Bio::SeqIO;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long qw(:config no_ignore_case);



# get a contig with a piece-of/entire  chromosome

my $global_start;
my $global_end; 
my $chr_name;   
my $t_id;
my $dbhost;
my $dbname;
my $dnadbname;
my $dnadbhost;
my $reverse;


GetOptions(
	    't_id:s' => \$t_id,
	    'dbname|db|D:s'   => \$dbname,
	    'dbhost|host|h:s'   => \$dbhost,
	    'dnadbname:s'=> \$dnadbname,
	    'dnadbhost:s'=> \$dnadbhost,
	    'reverse'    => \$reverse,
	   );

unless ( $t_id  && $dbname && $dbhost && $dnadbname && $dnadbhost ){
  print STDERR "Usage: $0 -t_id (transcript dbID) -dbname -dbhost -dnadbname -dnadbhost\n";
  exit(0);
}


# connect to the database
my $dnadb= new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $dnadbhost,
					      -user  => 'ensro',
					      -dbname=> $dnadbname);

my $db= new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $dbhost,
					   -user  => 'ensro',
					   -dbname=> $dbname,
					   -dnadb => $dnadb,
					  );


my $tadaptor = $db->get_TranscriptAdaptor;

my $transcript = $tadaptor->fetch_by_dbID( $t_id );

my $seq = $transcript->seq;

my $outfile;
$outfile = "$t_id.fa";

open OUT, ">$outfile";
# get a Bio::SeqIO object
my $out = Bio::SeqIO->new(-format=>'Fasta',
			  -fh =>  \*OUT,
			 );


$seq->display_id($t_id);
$out->write_seq($seq);	       
close OUT;




