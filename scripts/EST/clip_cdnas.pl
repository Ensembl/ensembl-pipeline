#!/usr/local/ensembl/bin/perl

# script to clip the ends of the given fasta file of cDNAs
# clipping can be hard (physically remove bases; default) or
# soft (replace clipped seq with Ns), and can be based
# on a given clip-range X i.e. clip X bases from each end)
# and/or on poly-A/T sequence occurring at the ends. If both
# range and polya clipping are asked for, the range clipping
# is performed first.  

use strict; 
use Getopt::Long;
use Bio::Seq;
use Bio::SeqIO;

use Bio::EnsEMBL::Utils::PolyA;

my ($mask, $poly_a_clip, $clip_len, $min_length);

&GetOptions( 
	     'mask'         => \$mask,
	     'hardclip=s'   => \$clip_len,
	     'polya'        => \$poly_a_clip,
	     'minlen=s'     => \$min_length
	   );


$min_length = 100 if not defined $min_length;

# fasta format
my $seqin = new Bio::SeqIO( '-format' => "fasta",
			    '-fh' => \*STDIN );
my $seqout = new Bio::SeqIO( '-format' => "fasta",
			     '-fh' => \*STDOUT );

my $polyA_clipper = Bio::EnsEMBL::Utils::PolyA->new();

while( my $cdna = $seqin->next_seq ){
    my $new_cdna;

    if ( $clip_len ){
	my $seq = $cdna->seq;
	my $seq_length = length( $seq );
	
	# skip it if you are going to clip more than the actual length of the EST
	if ( 2*$clip_len >= $seq_length ){
	    next SEQFETCH;
	}

	my $new_seq = substr( $seq, $clip_len, $seq_length - 2*$clip_len ); 
	if ($mask) {
	    $new_seq = "N" x $clip_len . $new_seq . "N" x $clip_len;
	}

	# skip it if you are left with an EST of less than 100bp
	if ( length( $new_seq ) < $min_length ){
	    next SEQFETCH;
	}

	$new_cdna = new Bio::Seq;
	$new_cdna->display_id( $cdna->display_id );
	$new_cdna->description( $cdna->description );
	$new_cdna->seq($new_seq);
    }
    else{ 
	$new_cdna = $cdna;
    }
    
    if ($poly_a_clip){
	#print STDERR "going to pass ".$new_cdna->display_id."\n";
	if ($mask) {
	    $new_cdna = $polyA_clipper->mask($new_cdna);
	}
	else {
	    $new_cdna = $polyA_clipper->clip($new_cdna);
	}

	if (not $new_cdna or $new_cdna->length < $min_length) {
	  next SEQFETCH;  
	} 
    }

    # write sequence
    $seqout->write_seq($new_cdna);
}
