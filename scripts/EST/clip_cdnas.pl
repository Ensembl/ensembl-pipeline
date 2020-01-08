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


# script to clip the ends of the given fasta file of cDNAs
# clipping can be hard (physically remove bases; default) or
# soft (replace clipped seq with Ns), and can be based
# on a given clip-range X i.e. clip X bases from each end)
# and/or on poly-A/T sequence occurring at the ends. If both
# range and polya clipping are asked for, the range clipping
# is performed first.  

use warnings ;
use strict; 
use Getopt::Long qw(:config no_ignore_case);
use Bio::Seq;
use Bio::SeqIO;

use Bio::EnsEMBL::Utils::PolyA;

my ($mask, $softmask, $poly_a_clip, $clip_len, $min_length, $outfile, $help);

GetOptions( 
	     'mask'         => \$mask,
	     'softmask'     => \$softmask,
	     'hardclip=s'   => \$clip_len,
	     'polyaclip'    => \$poly_a_clip,
	     'minlen=s'     => \$min_length,
             'outfile=s'    => \$outfile,
	     'help'         => \$help
	   );


&usage if $help;


die "You must give the name of an output filw with -outfile\n" 
    if not defined $outfile;

$min_length = 100 if not defined $min_length;

# fasta format

my $seqout = new Bio::SeqIO( '-format' => "fasta",
			     '-file' => ">$outfile" );

my $polyA_clipper = Bio::EnsEMBL::Utils::PolyA->new();

foreach my $file (@ARGV) {
    my $seqin = new Bio::SeqIO( '-format' => "fasta",
				'-file' => $file );

    while( my $cdna = $seqin->next_seq ){
	my $new_cdna;
	
	if ( $clip_len ){
	    my $seq = $cdna->seq;
	    my $seq_length = length( $seq );
	    
	    # skip it if you are going to clip more than the actual length of the EST
	    if ( 2*$clip_len >= $seq_length ){
		next SEQFETCH;
	    }

	    my $seq_left = substr( $seq, 0, $clip_len ); 
	    my $new_seq = substr( $seq, $clip_len, $seq_length - 2*$clip_len ); 
	    my $seq_right = substr( $seq, $clip_len + $seq_length - 2*$clip_len);

	    if ($mask) {
		$new_seq = "N" x $clip_len . $new_seq . "N" x $clip_len;
	    }
	    elsif ($softmask) {
		$new_seq = lc($seq_left) . $new_seq . lc($seq_right);
	    }
	    
	    # skip it if you are left with an EST of less than 100bp
	    if ( length( $new_seq ) < $min_length ){
		next SEQFETCH;
	    }
	    
	    $new_cdna = new Bio::Seq;
	    $new_cdna->display_id( $cdna->display_id );
	    $new_cdna->desc( $cdna->desc );
	    $new_cdna->seq($new_seq);
	}
	else{ 
	    $new_cdna = $cdna;
	}
	
	if ($poly_a_clip){
	    #print STDERR "going to pass ".$new_cdna->display_id."\n";
	    if ($mask or $softmask) {
		$new_cdna = $polyA_clipper->mask($new_cdna, $softmask);
	    }
	    else {
		$new_cdna = $polyA_clipper->clip($new_cdna);
	    }
	    
	    if (not $new_cdna or $new_cdna->length < $min_length) {
		next;
	    } 
	}
	
	# write sequence
	$seqout->write_seq($new_cdna);
    }
}



sub usage {
    print "Usage: clip_dnas.pl -outfile out.fa <-mask|-softmask> <-hardclip n> <-polyaclip> <-minlen n> file1.fa file2.fa ...\n\n";
    print "Recommended settings:\n   clip_cdnas.pl -outfile out.fa -polyaclip -minlen 100 file1.fa file.fa...\n";
    print "   (clips polyAs only, rejecting if result is less than 100 bp)\n";
    print "Other examples:\n";
    print "To softmask polyA:\n   clip_cdnas.pl -outfile out.fa -softmask -polyaclip file1.fa file2.fa ...\n";
    print "To hardmask polyA:\n   clip_cdnas.pl -outfile out.fa -mask -polyaclip file1.fa file2.fa ...\n";
    print "To hard clip 20bp from each end:\n   clip_cdnas.pl -outfile out.fa -hardclip 20 file1.fa file2.fa ...\n";
    print "Hard clip 20bp followed by polyA clip:\n   clip_cdnas.pl -outfile out.fa -hardclip 20 -polyA file1.fa file2.fa ...\n";
    print "To reject entries < 60bp after clipping:\n   clip_cdnas.pl -outfile out.fa -hardclip 20 -polyA -minlen 60 file1.fa file2.fa ...\n";
    print "To do nothing(!):\n   clip_cdnas.pl -outfile out.fa file1.fa file2.fa ...\n";
    exit(0);
}
