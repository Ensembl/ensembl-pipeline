#! /usr/local/ensembl/bin/perl -w

# A fairly hacky script to scan through a bunch of alignments
# and dig out the odd looking ones.

# This script contains a goodly amount of hard-coding.  Beware.

use strict;
use Bio::EnsEMBL::Pipeline::Alignment::EvidenceAlignment;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning info);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

use AlignmentSquizz;

# Handle command line options.

my ($chromosome_name,
    $dbname,
    $dbhost,
    $dbpass,
    $dbuser,
    $dbport,
    $dnadbname,
    $dnadbhost,
    $dnadbpass,
    $dnadbuser,
    $dnadbport,
    $transcript_stable_id,
    $transcript_dbid,
    $remove_introns,
    $type,
    $padding,
    $three_letter_codes,
    $fasta_line_length,
    $verbose,
    $output_file,
   ) = rearrange(['CHROMOSOME_NAME',
		  'DBNAME',
		  'DBHOST',
		  'DBPASS',
		  'DBUSER',
		  'DBPORT',
                  'DNADBNAME',
		  'DNADBHOST',
		  'DNADBPASS',
		  'DNADBUSER',
		  'DNADBPORT',
		  'TRANSCRIPT_STABLE_ID',
		  'TRANSCRIPT_DBID',
		  'REMOVE_INTRONS',
		  'TYPE',
		  'PADDING',
		  'THREE_LETTER_CODES',
		  'FASTA_LINE_LENGTH',
		  'VERBOSE',
		  'OUTPUT_FILE',
		 ], @ARGV);


# Check command line options for, at least, minimal sanity.

unless(defined $dbname &&
       defined $dbhost &&
       defined $dbuser){
  throw('Cannot connect to database without -dbname, -dbhost and ' .
	'-dbuser specified.')
}

#unless(defined $output_file){
#  throw('Please specify and output file where the alignment is ' .
#	'to be written.')
#}

if ($transcript_stable_id && $transcript_dbid){
  warning('Both a transcript_stable_id and transcript_dbid are ' .
	  'specified.  Only need one of these, but going with ' .
	  'the stable id.')
}

if (defined $type && ($type ne 'all' && 
		      $type ne 'nucleotide' && 
		      $type ne 'protein')){
  throw('Unrecognised alignment type.  The -type flag should be ' .
	'set to one of all, nucleotide or protein')

} elsif (! defined $type) {
  $type = 'all'
}

$padding = 10           unless defined $padding;
$three_letter_codes = 0 unless defined $three_letter_codes;
$fasta_line_length = 60 unless defined $fasta_line_length;
$verbose = 0            unless defined $verbose;


# Connect to database and go for it.

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
	    -dbname => $dbname,
	    -host   => $dbhost,
	    -port   => $dbport,
	    -user   => $dbuser,
	    -pass   => $dbpass);

if (defined($dnadbname)) {
  my $dnadb = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
  	    -dbname => $dnadbname,
  	    -host   => $dnadbhost,
  	    -port   => $dnadbport,
  	    -user   => $dnadbuser,
  	    -pass   => $dnadbpass);
  $db->dnadb($dnadb);
}

my $ta = $db->get_TranscriptAdaptor;
my $sa = $db->get_SliceAdaptor;
my $ga = $db->get_GeneAdaptor;

# Derive chromosome and work through a single transcript
# from each gene on this chromosome.

my $chr_slice = $sa->fetch_by_region('chromosome', $chromosome_name);

my $genes = $ga->fetch_all_by_Slice($chr_slice);


# Work through each gene.

foreach my $gene (@$genes) {

  next unless $gene->biotype eq 'protein_coding';

  my $transcript = $gene->get_all_Transcripts->[0];

#  if (defined $transcript_stable_id){
#    $transcript = $ta->fetch_by_stable_id($transcript_stable_id)
#  } else {
#    $transcript = $ta->fetch_by_dbID($transcript_dbid)
#  }

#  die "Unable to retrieve transcript.\n"
#    unless defined $transcript;

  my $seqfetcher = Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch->new();

  my $evidence_alignment = 
    Bio::EnsEMBL::Pipeline::Alignment::EvidenceAlignment->new(
      -transcript          => $transcript,
      -dbadaptor           => $db,
      -seqfetcher          => $seqfetcher,
      -padding             => $padding);

  my $alignment;

  eval {
    $alignment =
      $evidence_alignment->retrieve_alignment(
       '-type'            => $type,
       '-remove_introns'  => $remove_introns,
       '-three_letter_aa' => $three_letter_codes,);
  };


  if ($@ || ! $alignment) {
    print STDOUT "Alignment failed for transcript " . 
      $transcript->stable_id . "\n";
    next
  }

  my $align_squizz = AlignmentSquizz->new($alignment);
  $align_squizz->appraise;
  my $output = $align_squizz->output_array;

  foreach my $row (@$output) {
    if (($row->[6] > 0.05)||
	($row->[8] < 0.25)){
      print STDOUT join("\t", $transcript->stable_id, @$row) . "\n"
    }
  }
}
