#! /usr/local/ensembl/bin/perl -w

use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);


# Handle command line options.

my ($dbname,
    $dbhost,
    $dbpass,
    $dbuser,
    $dbport,
    $file) = rearrange(['DBNAME',
			  'DBHOST',
			  'DBPASS',
			  'DBUSER',
			  'DBPORT',
			  'PEPFILE',
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
             " -pepfile:               (optional) File to write to",
	    ) . "\n";
  die
}

unless(defined $dbname &&
       defined $dbhost &&
       defined $dbuser){
  die('Cannot connect to database without -dbname, -dbhost and ' .
	'-dbuser specified.')
}


$| = 1;
if (defined($file) && $file ne "stdout") {
  open FP,">$file";
} else {
  open FP,">-";
}

# Connect to database

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
                    -dbname => $dbname,
		    -host   => $dbhost,
		    -user   => $dbuser,
		    -pass   => $dbpass,
		    -port   => $dbport);

my $slice_adaptor     = $db->get_SliceAdaptor;
my $gene_adaptor      = $db->get_GeneAdaptor;

# Fetch all genes from all chromosomes

my @genes;

foreach my $chr (@{$slice_adaptor->fetch_all('chromosome')}){

  print STDERR "Chr " .  $chr->name . "\t";

  my $chr_slice = $slice_adaptor->fetch_by_name($chr->name);
  my $genes = $gene_adaptor->fetch_all_by_Slice($chr_slice);

  print STDERR scalar @$genes . " genes\n";

  push @genes, @$genes;
}


# Process each gene.

foreach my $gene (@genes) {

  next if (($gene->type =~ /pseudogene/i) ||
	   ($gene->type =~ /putative/i) ||
	   ($gene->type =~ /obsolete/i) ||
	   ($gene->type =~ /Novel_Transcript/i));

  my $transcript = &longest_transcript($gene);

  next unless $transcript;

  print STDERR join ("\t", ($gene->stable_id, 
			    $transcript->stable_id,
			    $transcript->translation->stable_id)
		    ) . "\n";

  print FP '>' . $gene->stable_id . "\t" . 
    $transcript->stable_id . "\t" . $transcript->translation->stable_id . "\n";

#  my $seq = $transcript->spliced_seq;
  my $seq = $transcript->translateable_seq;


  # Trim any stop codons TAA, TGA, TGG from sequence ends.
  $seq =~ s/TGA$//i;
  $seq =~ s/TAA$//i;
  $seq =~ s/TAG$//i;


  # Add line breaks
  $seq =~ s/(.{60})/$1\n/g;

  print FP uc($seq) . "\n";
}

sub longest_transcript ($) {
  my $gene = shift;

  my $longest = 0;
  my $longest_transcript;

  my $translatable = 0;

  foreach my $transcript (@{$gene->get_all_Transcripts}){

    eval {
      $transcript->translate;
    };

    next if $@;

    $translatable = 1;

#     my $length = $transcript->spliced_seq =~ tr/[ATGCatgcNn]//;
     my $length = $transcript->translateable_seq =~ tr/[ATGCatgcNn]//;

     if ($length > $longest){
       $longest = $length;
       $longest_transcript = $transcript
     }
  }


  unless ($longest_transcript){
    warn("Gene [" . $gene->stable_id . 
	 "] does not have a transcript with a translation.");
    $longest_transcript = 0;
  }


  return $longest_transcript
}
