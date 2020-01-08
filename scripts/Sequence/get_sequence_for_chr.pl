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

use Bio::SeqIO;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long qw(:config no_ignore_case);

# get a contig with a piece-of/entire  chromosome

my $mask;
my $softmask;
my $dust;
my $outfile;
my $name;
my $maskall;
my $coord_system;
my $dbhost;
my $dbuser = 'ensro';
my $dbname;
my $dbpass;
my $dbport = 3306;

GetOptions(
            'seq_region_name:s' => \$name,
            'dbname|db|D:s'       => \$dbname,
            'dbhost|host|h:s'       => \$dbhost,
            'dbpass|pass|p:s' => \$dbpass,
            'dbuser|user|u:s' => \$dbuser,
            'dbport|port|P:s' => \$dbport,
            'mask'           => \$mask,
            'dust'           => \$dust,
            'softmask'       => \$softmask,
            'maskall'        => \$maskall,
            'outfile:s'      => \$outfile,
            'coord_system|cs_name:s' => \$coord_system,
	    
	    );

unless( $dbhost && $dbname && $name && $outfile ){
  print STDERR ("Usage: -dbname $dbname -dbhost $dbhost -dbuser $dbuser ".
                "-seq_region_name $name -coord_system $coord_system ".
                "[ -mask -dust -softmask -maskall ] -outfile $outfile\n");
  exit(0);
}

# connect to the database
my $db= new Bio::EnsEMBL::DBSQL::DBAdaptor(
                                           -host  => $dbhost,
                                           -user  => $dbuser,
                                           -port => $dbport,
                                           -dbname=> $dbname,
                                           -dbpass => $dbpass,
                                          );


open OUT, ">$outfile";
my $out = Bio::SeqIO->new(-format=>'Fasta',
			  -fh =>  \*OUT,
			 );

my $slice = $db->get_SliceAdaptor->fetch_by_region($coord_system, $name );


my @logic_names;
my $soft = 0;
my $string  = '';
if ($mask){
    push( @logic_names, 'RepeatMask' );
}
if ( $dust ){
    push (@logic_names, 'Dust' );
}
if ($softmask){
    $soft = 1;
    $string = 'softmask';
}

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

############################################################
