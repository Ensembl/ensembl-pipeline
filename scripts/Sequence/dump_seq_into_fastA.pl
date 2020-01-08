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
my $input_id;
my $masked;
my $softmasked;
my $dust;
my $dbhost;
my $dbname;
my $maskall;


GetOptions(
	    'input_id:s' => \$input_id,
	    'masked'     => \$masked,
	    'softmasked' => \$softmasked,
	    'maskall'    => \$maskall,
	    'dust'       => \$dust,
	    'dbname|db|D:s'   => \$dbname,
	    'dbhost|host|h:s'   => \$dbhost,
	   );

unless ( $input_id && $dbname && $dbhost ){
  print STDERR "Usage: $0 -input_id chr_name.chr_start-chr_end -dbname -dbhost ( -masked -dust -maskall -softmasked )\n";
  exit(0);
}


# connect to the database
my $db= new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $dbhost,
					   -user  => 'ensro',
					   -dbname=> $dbname);


my $sgp = $db->get_SliceAdaptor;

my $outfile = "$input_id.fa";
$input_id =~/(\S+)\.(\d+)-(\d+)/;
$chr_name     = $1;
$global_start = $2;
$global_end   = $3;

my $slice = $sgp->fetch_by_chr_start_end($chr_name,$global_start,$global_end);

# slice is a Bio::EnsEMBL::PrimarySeq


my @logic_names;
my $soft = 0;
my $string  = '';
if ($masked){
  push( @logic_names, 'RepeatMask' );
}
if ( $dust ){
  push (@logic_names, 'Dust' );
}
if ($softmasked){
  $soft = 1;
  $string = 'softmask';
}

open OUT, ">$outfile";
# get a Bio::SeqIO object
my $out = Bio::SeqIO->new(-format=>'Fasta',
			  -fh =>  \*OUT,
			 );


my $seq;
if ($maskall){
  $string .= " maskall ";
  $seq = $slice->get_repeatmasked_seq([''],$soft);
}
elsif( @logic_names ){
  print STDERR "getting sequence soft - repeatmaksed and dusted\n"; 
  $seq = $slice->get_repeatmasked_seq(\@logic_names,$soft);
}
else{
  $seq = $slice;
}

$seq->display_id($slice->display_id." @logic_names ".$string);
$out->write_seq($seq);	       
close OUT;




