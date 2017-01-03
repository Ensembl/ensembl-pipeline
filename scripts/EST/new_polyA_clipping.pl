#!/usr/bin/env perl


# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2017] EMBL-European Bioinformatics Institute
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

=pod

=head1 NAME

  new_polyA_clipping.pl

=head1 DESCRIPTION

  Script to parse a fasta file and identify sequences with polyA/T
  tails/heads these are then clipped and stored in a file. Clipped,
  and potentially trimmed, sequences that are shorter than the
  specified minimum length are not stored.

=head2 Clipping

  - The non-A/T sequences at the ends must be <=10bp (set by $buffer)
  - The polyA/T string must be >4bp to be removed
  - It only clips polyA tails or polyT heads using a sliding window of 3 bp
  - The clipping is recursive but only clips one end of a sequence
  - The head/tail is only clipped if the polyA/T string is longer than
    the non-polyA/T string at the end of the sequence

=head1 OPTIONS

  -readfile     Input sequence file in fasta format
  -outfile      Output file in fasta format
  -min_length   Minimum sequence length (default 60)
  -trim         Flag for trimming
  -comments     sequence quality file

  Note that the -trim and -comments options are dependent and need
  to be specified together when trimming the sequences. The sequence
  quality file should contain the following tab columns:

    1./ accession number
    2./ high quality sequence start or stop
    3./ numeric value

  For human ESTs, the sequence quality file can be created by
  connecting to the Mole embl database and doing the following query,
  changing tax_division and data_class as appropriate:

    SELECT accession_version,comment_key,comment_value
    FROM entry e, comment c
    WHERE e.entry_id = c.entry_id
    AND e.tax_division = 'HUM'
    AND e.data_class = 'EST'

  The query takes about 5 mins for 1 million ESTs, depending on how
  busy the database is.

=head1 EXAMPLES

  For only polA clipping:

    perl new_polyA_clipping.pl -readfile sequences.fasta -outfile polyat_clipped.out

  For trimming and polyA clipping:

    perl new_polyA_clipping.pl -readfile sequences.fasta -outfile polyat_clipped.out -trim -comments comments.ls

  For trimming and polyA clipping with specified min length cutoff:

    perl new_polyA_clipping.pl -readfile sequences.fasta -outfile polyat_clipped.out -trim -comments comments.ls -min_length 50

=cut

use Bio::EnsEMBL::Analysis::Tools::PolyAClipping;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::Seq;
use Bio::SeqIO;
use Getopt::Long qw(:config no_ignore_case);

my (
    $data,
    $clipped_cdnas,
    $trim,
    $hqs_comment_file,
    $min_length,
   );

my $err_file = '';

#
# NOTE : This script is used in production code ( cdna_update ) - please do not
# change the calls / options, otherewise the cDNA update code breaks ...
#
GetOptions(
            'readfile:s'   => \$data,             # '/path/to/unclipped/cdnas_unclipped.fa';
            'outfile:s'    => \$clipped_cdnas,    # '/path/to/clipped/cdnas_clipped.fa';
            'trim'         => \$trim,             # trimming flag
            'comments:s'   => \$hqs_comment_file, # '/path/to/high-quality-sequence-comments-file';
            'min_length:s' => \$min_length,       # cutoff on sequence length
            'errfile:s'    => \$err_file,         # file to direct STDERR to
           );

# this is for backwards compatibility as options changed in r 1.10
$data = $ARGV[0] if $ARGV[0]; 
$clipped_cdnas = $ARGV[1] if $ARGV[1] ;

if($err_file){
  open STDERR, ">$err_file" or die "Can't redirect stderr\n";
}

if (defined $min_length) {
  print STDERR "Using minimum length of $min_length.\n";
} else {
  print STDERR "No minimum seq length specified. Using default of 60 bases.\n";
  $min_length = 60;
}

if ($trim && !defined $hqs_comment_file) {
  die "Please enter -comments file path if you'd like to trim sequences. "
    . "Please enter both or none of these options.\n";
}

if (!$trim && defined $hqs_comment_file) {
  die "You specified -comments file path but not -trim. "
    . "Please enter both or none of these options.\n";
}

my $seqin  = new Bio::SeqIO( -file   => "<$data",
                             -format => "Fasta",
                           );

my $seqout = new Bio::SeqIO( -file   => ">$clipped_cdnas",
                             -format => "Fasta"
                           );

my %comments;
if ($trim) {
  print STDERR "Reading in comments...\n";
  open(COMMENTS, "<$hqs_comment_file") or die "Cannot open $hqs_comment_file\n";
  while (my $line = <COMMENTS>) {
    chomp $line;
    my @fields = split('\t', $line);
    $comments{$fields[0]}{$fields[1]} = $fields[2];
  }
  close(COMMENTS);
}
print STDERR "Have ".scalar(keys %comments)." accessions with comments\n\nCleaning sequences...\n";

SEQFETCH:
while ( my $fullseq = $seqin->next_seq ) {

  # first trim, if necessary and if we have enough info
  my $unclipped;
  if ($trim) {
    ($unclipped) = trim_if_possible($fullseq, $min_length, \%comments);
  } else {
    $unclipped = $fullseq;
  }
  unless ($unclipped) {
   throw ("Sequence is undefined. Original seq was ". $fullseq->display_id."\n");
  }
  # now clip
  my ($clipped, $clip_end, $num_bases_removed) = clip_if_necessary($unclipped);
  if (defined $num_bases_removed) {
    # print STDERR "Clipped $num_bases_removed bases from seq. Writing to file.\n";
  }
  if (defined $clipped) {
    # Skip the sequence if length is shorter than $min_length after
    # it's been clipped.
    if ( $clipped->length < $min_length ) {
      print STDERR "After clipping, "
        . $unclipped->display_id
        . " is shorter than $min_length, it will be skipped.\n";
    } else {
      $seqout->write_seq($clipped);
    }
  } else {
    print STDERR "Sequence removed: ". $unclipped->display_id . "\n";
  }

}
if($err_file){
  close STDERR;
}


sub trim_if_possible {
  my ($seq, $min_length, $high_quality_seq) = @_;
  #my %high_quality_seq = %{$high_quality_seq};

  if (exists $high_quality_seq->{$seq->display_id}) {   
    my $start = 1;
    my $stop = $seq->length;

    foreach my $key (keys %{$high_quality_seq->{$seq->display_id}}) {
      if ($key =~ m/high quality sequence start/) {
        $start = $high_quality_seq->{$seq->display_id}{$key};
      } elsif ($key =~ m/high quality sequence stop/) {
        $stop = $high_quality_seq->{$seq->display_id}{$key};
      } else {
        die "EXCEPTION: Comment key '$key' not recognised. Exiting script.\n";
      }
    }
 
    # now let's see about trimming
    if ($start!=1 or $stop!=$seq->length){
      # we need to trim if it doesn't make seq too short
      if ($stop <= $seq->length) {
        if (($stop-$start+1 < $min_length) || ($stop == $start)) {
          # sometimes ESTs are all low-quality. The sequence after trimming
          # can become very short (shorter than $min_length), or in some
          # extreme cases, both 'high quality sequence start'
          # and 'high quality sequence stop' have the same value, and
          # trimming will result in a sequence with 1bp only. Therefore,
          # we should skip trimming and leave the sequence unchanged and
          # just print a note about it.
          #
          # The latter part of the if clause ($stop == $start) here is necessary 
          # to catch cases which are not worth trimming when $min_length 
          # is set to 0 by the user.
          print STDERR "Not trimming seq ".$seq->display_id.
                       " which would trim to length " . ($stop-$start+1) . "\n";
        } else {
          # do the trim
          $seq->seq($seq->subseq($start,$stop));
          print STDERR "Trim seq to $start - $stop\n";
          }
      } else {
        print STDERR "Not trimming seq ".$seq->display_id." which has invalid ".
                     "hq stop pos " . $start . " " . $stop . " " . $seq->length . "\n";
      }
    }

  } else {
    # no high_quality_seq so we just keep the original seq
    print STDERR "Seq ".$seq->display_id." has no high_quality_seq; cannot trim\n";
    }
  return $seq;
}
