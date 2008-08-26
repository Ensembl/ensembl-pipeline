#!/usr/local/ensembl/bin/perl -w
use strict;

=pod

=head1 NAME

polA-clipping.pl

=head1 DESCRIPTION
         
script to parse a fasta file and identify sequences with polyA/T tails/heads
these are then clipped and stored in a file 

Script will also trim the sequence before polyA clipping, if a path to the -comments
file is specified. This file should contain the following columns:
  1./ accession
  2./ high quality sequence start or stop
  3./ numeric value

clipping:
  the non-A/T sequences at the ends must be <=10bp (set by $buffer)  
  the polyA/T string must be >4bp to be removed
  it only clips polyA tails or polyT heads using a sliding window of 3 bp
  the clipping is recursive but only clips one end of a sequence
  the head/tail is only clipped if the polyA/T string is longer than the non-polyA/T string at the end of the sequence

For only polA clipping:
  perl new_polyA_clipping.pl -read sequences.fasta -out polyat_clipped.out
For trimming and polyA clipping 
  perl new_polyA_clipping.pl -read sequences.fasta -out polyat_clipped.out -trim -comments comments.ls

=cut

use Bio::EnsEMBL::Analysis::Tools::PolyAClipping;
use Bio::SeqIO;
use Getopt::Long;

my (
    $data,
    $clipped_cdnas,
    $trim,
    $hqs_comment_file,
    $min_length,
   );

&GetOptions(
            'readfile:s'   => \$data, # = '/path/to/unclipped/cdnas_unclipped.fa';
            'outfile:s'    => \$clipped_cdnas, # = '/path/to/clipped/cdnas_clipped.fa';
            'trim'         => \$trim, # flag
            'comments:s'   => \$hqs_comment_file, #  = '/path/to/high-quality-sequence-comments-file';
            'min_length:s' => \$min_length,
           );

if ($trim && !defined $hqs_comment_file) {
  die "Please enter -comment file path if you'd like to trim sequences\n";
}
if (!$trim && defined $hqs_comment_file) {
  die "You specified -comment file path but not -trim. Please enter both or none of these options.\n";
}
if (!defined $min_length) {
  $min_length = 60;
  print STDERR "No minimum seq length specified. Using default of 60 bases.\n";
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

  # now clip
  my ($clipped, $clip_end, $num_bases_removed) = clip_if_necessary($unclipped);
  if (defined $num_bases_removed) {
    print STDERR "Clipped $num_bases_removed bases from seq. Writing to file.\n";
  }
  if (defined $clipped) {
    # we currently do not specify a minimum seq length here, but it might be a good idea
    $seqout->write_seq($clipped);
  } else {
    print STDERR "Sequence removed: ".$unclipped->display_id."\n";
  }
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
        if ($stop-$start+1 < $min_length) {
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
