#!/usr/local/ensembl/bin/perl


use strict;

use Getopt::Long;
use Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptPair;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;

############################################################
# refseq mapping is a file with the mappings of locus link
# ids to the refseq ids (both NM_ and NP_ entries)
my $refseq_mapping;

my $gene_id1;
my $gene_id2;
my $gap_penalty;
my $check;
my $coding_exons;

&GetOptions(
	    'refseq_mapping:s' => \$refseq_mapping,
	    'gene_id1:s'       => \$gene_id1,
	    'gene_id2:s'       => \$gene_id2,
	    'check'            => \$check,
	    'coding_exons'     => \$coding_exons,
	   );

if ( $check ){
    exit(0);
}

unless ( $refseq_mapping && $gene_id1 && $gene_id2 ){
  print STDERR "Usage: $0 -file -gene_id1 -gene_id2\n";
  exit(0);
}

open( MAP, "<$refseq_mapping") or die ("cannot open $refseq_mapping");

############################################################
# the map is of the format
# humanLL humanNM humanNP mouseLL mouseNM mouseNP
my %human_pair;
my %human_cdna2protein;
my %mouse_cdna2protein;
my %human_cdnas;
my %mouse_cdnas;
while(<MAP>){
  chomp;
  my @entries = split;
  my $humanLL = $entries[0];
  my $humanNM = $entries[1];
  my $humanNP = $entries[2];
  my $mouseLL = $entries[3];
  my $mouseNM = $entries[4];
  my $mouseNP = $entries[5];

  push( @{$human_pair{$humanLL}}, $mouseLL );
  push( @{$human_cdnas{$humanLL}}, $humanNM );
  push( @{$mouse_cdnas{$mouseLL}}, $mouseNM );
  $human_cdna2protein{$humanNM}=$humanNP;
  $mouse_cdna2protein{$mouseNP}=$mouseNP;
}
close(MAP);

my $this_humanLL = $gene_id1;
my $this_mouseLL = $gene_id2;

my %already_seen;
my @human_cdnas = @{$human_cdnas{$this_humanLL}};
my @mouse_cdnas = @{$mouse_cdnas{$this_mouseLL}};

print STDERR "comparing $this_humanLL (".scalar(@human_cdnas)." transcripts ) ".
  "and $this_mouseLL (".scalar(@mouse_cdnas)." transcripts )\n";

my $transcript_map = Bio::EnsEMBL::Pipeline::GeneComparison::ObjectMap->new();
my $protein_map = Bio::EnsEMBL::Pipeline::GeneComparison::ObjectMap->new();
my %seqobj;
my $seqfetcher = Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch->new();

my %target_coverage;
my %query_coverage;
my %perc_id;

HUMAN:
foreach my $humanNM ( @human_cdnas ){
 MOUSE:
  foreach my $mouseNM ( @mouse_cdnas ){ 
    
    next if $already_seen{$humanNM}{$mouseNM};

    ############################################################
    # compare transcripts
    my $humanNM_seq = $seqfetcher->get_Seq_by_acc($humanNM);
    my $mouseNM_seq = $seqfetcher->get_Seq_by_acc($mouseNM);
    unless ( $humanNM_seq ){
      print STDERR "Failed to find sequence $humanNM\n";
      next HUMAN;
    }
    unless ( $mouseNM_seq ){
      print STDERR "Failed to find sequence $mouseNM\n";
      next MOUSE;
    }

    my ($score,$best_features, $target_coverage, $query_coverage, $perc_id) = 
      Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptPair
	->blast_unmapped_transcripts( $humanNM_seq, $mouseNM_seq);
    if ( $score && $best_features ){
      $transcript_map->match($humanNM_seq, $mouseNM_seq, $score );
      $perc_id{$humanNM_seq}{$mouseNM_seq} = $perc_id;
      $target_coverage{$humanNM_seq}{$mouseNM_seq} = $target_coverage;
      $query_coverage{$humanNM_seq}{$mouseNM_seq} = $query_coverage;
    }
    unless ($score){
      $score = 0;
    }
    my $humanNP = $human_cdna2protein{$humanNM};
    my $mouseNP = $mouse_cdna2protein{$mouseNM};
    print STDERR "Pair [ $humanNM($humanNP)  , $mouseNM($mouseNP) ] score = $score\n";
    
    ############################################################
    # compare peptides
    my $humanNP_seq = $seqfetcher->get_Seq_by_acc($humanNP);
    my $mouseNP_seq = $seqfetcher->get_Seq_by_acc($mouseNP);
    unless ( $humanNP_seq ){
      print STDERR "Failed to find sequence $humanNP\n";
      next HUMAN;
    }
    unless ( $mouseNP_seq ){
      print STDERR "Failed to find sequence $mouseNP\n";
      next MOUSE;
    }
    
    my ($score,$best_features, $target_coverage, $query_coverage, $perc_id) = 
      Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptPair
	->blast_unmapped_proteins( $humanNP_seq, $mouseNP_seq);
    if ( $score && $best_features ){
      $protein_map->match($humanNP_seq, $mouseNP_seq, $score );
      $perc_id{$humanNP_seq}{$mouseNP_seq} = $perc_id;
      $target_coverage{$humanNP_seq}{$mouseNP_seq} = $target_coverage;
      $query_coverage{$humanNP_seq}{$mouseNP_seq} = $query_coverage;
    }
    unless ($score){
      $score = 0;
    }
    print STDERR "Pair [ $humanNP  , $mouseNP ] score = $score\n";
    
  }
}

############################################################
# pairs created:
my $best_transcript_pairs_object = $transcript_map->stable_marriage;

my $pair_count = scalar($best_transcript_pairs_object->list1);
print STDERR "Transcript pairs created: ".$pair_count."\n";
foreach my $element1 ( $best_transcript_pairs_object->list1 ){
  foreach my $partner ( $best_transcript_pairs_object->partners( $element1 ) ){
    # there should be only one
    my $id1 = $element1->display_id;
    my $id2 = $partner->display_id;
    my $protein_id1 =  $human_cdna2protein{$id1};
    my $protein_id2 =  $mouse_cdna2protein{$id2};
      
    print STDERR "Pair ( $id1 , $id2 ) with score: ".$best_transcript_pairs_object->score( $element1, $partner )."\n";
    
    ############################################################
    # summary line
    my $length1 = $element1->length;
    my $length2 = $partner->length;
    my $perc_id         = $perc_id{$element1}{$partner};
    my $target_coverage = $target_coverage{$element1}{$partner};
    my $query_coverage  = $query_coverage{$element1}{$partner};
    my $length_diff      = $length1 - $length2;
    print STDERR "TRANSCRIPT_PAIR\t".
      "$this_humanLL\t$id1\t$protein_id1\tlength1:$length1\t".
	"coverage1:$target_coverage\t".
	  "$this_mouseLL\t$id2\t$protein_id2\tlength2:$length2\t".
	    "coverage2:$query_coverage\t".
	      "perc_id:$perc_id\t".
		"length_diff:$length_diff\n";
  }
}

my $best_protein_pairs_object = $protein_map->stable_marriage;
my $pair_count = scalar($best_protein_pairs_object->list1);
print STDERR "Protein pairs created: ".$pair_count."\n";
foreach my $element1 ( $best_protein_pairs_object->list1 ){
  foreach my $partner ( $best_protein_pairs_object->partners( $element1 ) ){
    # there should be only one
    my $id1 = $element1->display_id;
    my $id2 = $partner->display_id;
    print STDERR "Pair ( $id1 , $id2 ) with score: ".$best_protein_pairs_object->score( $element1, $partner )."\n";
    
    ############################################################
    # summary line
    my $length1 = $element1->length;
    my $length2 = $partner->length;
    my $perc_id         = $perc_id{$element1}{$partner};
    my $target_coverage = $target_coverage{$element1}{$partner};
    my $query_coverage  = $query_coverage{$element1}{$partner};
    my $length_diff      = $length1 - $length2;
    print STDERR "PROTEIN_PAIR\t".
      "$this_humanLL\t$id1\tlength1:$length1\t".
	"coverage1:$target_coverage\t".
	  "$this_mouseLL\t$id2\tlength2:$length2\t".
	    "coverage2:$query_coverage\t".
	      "perc_id:$perc_id\t".
		"length_diff:$length_diff\n";
  }
}
