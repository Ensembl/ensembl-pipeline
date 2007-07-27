#!/usr/local/ensembl/bin/perl -w
use strict;

=pod

=head1 NAME

polA-clipping.pl

=head1 DESCRIPTION
         
script to parse a fasta file and identify sequences with polyA/T tails/heads
these are then clipped and stored in a file 

clipping:
  the non-A/T sequences at the ends must be <=10bp (set by $buffer)  
  the polyA/T string must be >4bp to be removed
  it only clips polyA tails or polyT heads using a sliding window of 3 bp
  the clipping is recursive but only clips one end of a sequence
  the head/tail is only clipped if the polyA/T string is longer than the non-polyA/T string at the end of the sequence

perl new_polyA_clipping.pl sequences.fasta polyat_clipped.out

=cut

use Bio::EnsEMBL::Analysis::Tools::PolyAClipping;
use Bio::SeqIO;

my $data = $ARGV[0]; # = '/path/to/unclipped/cdnas_unclipped.fa';
my $clipped_cdnas = $ARGV[1]; # = '/path/to/clipped/cdnas_clipped.fa';

my $seqin  = new Bio::SeqIO( -file   => "<$data",
                             -format => "Fasta",
                           );

my $seqout = new Bio::SeqIO( -file   => ">$clipped_cdnas",
                             -format => "Fasta"
                           );

SEQFETCH:
while ( my $unclipped = $seqin->next_seq ) {
  my ($clipped, $clip_end, $num_bases_removed) = clip_if_necessary($unclipped);
  $seqout->write_seq($clipped);
}  


