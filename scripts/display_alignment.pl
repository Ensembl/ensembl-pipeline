#! /usr/local/ensembl/bin/perl -w

=head1 NAME

display_alignment.pl

=head1 DESCRIPTION

Simple script to drive the one-off reconsruction of the mutiple alignments 
for ensembl transcripts and their aligned evidence.  All options are 
stipulated at the command line.  This script basicallt just makes use of the 
Bio::EnsEMBL::Pipeline::Alignment::EvidenceAlignment module, which is where 
all the heavy lifting is conducted.  This module allows more complex alignments
to be attempted, specifically where an alignment is generated for a transcript
and the evidence for this transcript (or any other feature, really) resides
in another database.  See the docs in the EvidenceAlignment module.

=head1 OPTIONS

-dbname :
Name of core ensembl database for the species in question.

-dbhost :
Name of machine that hosts this database.

-dbuser :
Username for database access.

-dbpass :
(optional) Database password.

-dbport :
(optional) Database access port number.

-dnadbname :
(optional) Name of core ensembl database containing DNA for the species in question.

-dnadbhost :
(optional) Name of machine that hosts the DNA database.

-dnadbuser :
(optional) Username for database access to DNA database.

-dnadbpass :
(optional) Database password for DNA databse.

-dnadbport :
(optional) Database access port number for DNA databse.

-transcript_stable_id :
Stable identifier of transcript to display.

-transcript_dbid :
If you dont have a stable id, use the dbID.

-remove_introns :
(optional) Truncate intron sequences to display a
more compact alignment.

-padding :
(optional)  Number of bases padding upstream and
downstream of transcript.

-type :
(optional) Type of evidence sequences to display in
alignment, one of 'all', 'nucleotide' or 'protein'.

-three_letter_codes:
Display amino acids with three letter representation.

-fasta_line_length :
(optional) Fasta output file line length.

-verbose :
To get tabulated details of exon and feature coordinates.

-output_file :
File where alignment will be saved.

=cut

use strict;
use Bio::EnsEMBL::Pipeline::Alignment::EvidenceAlignment;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning info);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

# Handle command line options.

my ($dbname,
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
   ) = rearrange(['DBNAME',
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

unless (@ARGV) {
  print join("\n",
	     "Options are : ",
	     " -dbname :               Name of core ensembl database for the species",
	     "                           in question.",
	     " -dbhost :               Name of machine that hosts this database.",
	     " -dbuser :               Username for database access.",
	     " -dbpass :               (optional) Database password.",
	     " -dbport :               (optional) Database access port number.",
	     " -dnadbname :            (optional) Name of core ensembl database containing DNA.",
	     " -dnadbhost :            (optional) Name of machine that hosts DNA database.",
	     " -dnadbuser :            (optional) Username for DNA database access.",
	     " -dnadbpass :            (optional) DNA Database password.",
	     " -dnadbport :            (optional) DNA Database access port number.",
	     " -transcript_stable_id : Stable identifier of transcript to display.",
	     " -transcript_dbid :      If you dont have a stable id, use the dbID.",
	     " -remove_introns :       (optional) Truncate intron sequences to display a",
	     "                           more compact alignment.",
	     " -padding :              (optional) Number of bases padding upstream and",
	     "                           downstream of transcript.",
	     " -type :                 (optional) Type of evidence sequences to display",
	     "                           in alignment, one of 'all', 'nucleotide'",
	     "                           or 'protein'.",
	     " -three_letter_codes:    Display amino acids with three letter",
	     "                           representation.",
	     " -fasta_line_length :    (optional) Fasta output file line length.",
	     " -verbose :              To get tabulated details of exon and feature",
	     "                           coordinates.",
	     " -output_file :          File where alignment will be saved."
	    ) . "\n";
  die
}

unless(defined $dbname &&
       defined $dbhost &&
       defined $dbuser){
  throw('Cannot connect to database without -dbname, -dbhost and ' .
	'-dbuser specified.')
}

unless(defined $output_file){
  throw('Please specify and output file where the alignment is ' .
	'to be written.')
}

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

my $transcript;

if (defined $transcript_stable_id){
  $transcript = $ta->fetch_by_stable_id($transcript_stable_id)
} else {
  $transcript = $ta->fetch_by_dbID($transcript_dbid)
}

die "Unable to retrieve transcript.\n"
  unless defined $transcript;

my $seqfetcher = Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch->new();

my $evidence_alignment = 
  Bio::EnsEMBL::Pipeline::Alignment::EvidenceAlignment->new(
      -transcript          => $transcript,
      -dbadaptor           => $db,
      -seqfetcher          => $seqfetcher,
      -padding             => $padding);

my $alignment =
  $evidence_alignment->retrieve_alignment(
     '-type'            => $type,
     '-remove_introns'  => $remove_introns,
     '-three_letter_aa' => $three_letter_codes,);

$evidence_alignment->_print_tabulated_coordinates
  if $verbose;

my $align_seqs = $alignment->fetch_AlignmentSeqs;

open(OUT, ">$output_file") or die "Unable to write to output file [$output_file]";

foreach my $align_seq (@$align_seqs){
  my $seq = $align_seq->seq;
  $seq =~ s/(.{$fasta_line_length})/$1\n/g;

  print OUT ">" . $align_seq->name . "\n" . $seq . "\n";
}

close(OUT);
