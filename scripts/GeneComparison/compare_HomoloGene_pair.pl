#!/usr/local/ensembl/bin/perl


use strict;

use Getopt::Long;
use Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptPair;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;
use Bio::EnsEMBL::Pipeline::GeneComparison::ObjectMap;

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
  print STDERR "Usage: $0 -refseq_mapping -gene_id1 -gene_id2\n";
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
my %seen_mouse;
my %seen_human;
my %seen_pair;

while(<MAP>){
  chomp;
  my @entries = split;
  my $humanLL = $entries[0];
  my $humanNM = $entries[1];
  my $humanNP = $entries[2];
  my $mouseLL = $entries[3];
  my $mouseNM = $entries[4];
  my $mouseNP = $entries[5];

  unless( $seen_pair{$humanLL}{$mouseLL} ){
    push( @{$human_pair{$humanLL}}, $mouseLL );
    $seen_pair{$humanLL}{$mouseLL} = 1;
  }
  unless( $seen_human{$mouseLL}{$humanNM} ){
    push( @{$human_cdnas{$humanLL}}, $humanNM );
    $human_cdna2protein{$humanNM}=$humanNP;
    $seen_human{$mouseLL}{$humanNM} = 1;
  }
  unless( $seen_mouse{$humanLL}{$mouseNM} ){
    push( @{$mouse_cdnas{$mouseLL}}, $mouseNM );
    $mouse_cdna2protein{$mouseNM}=$mouseNP;
    $seen_mouse{$humanLL}{$mouseNM} = 1;
  }
  
  }
close(MAP);

my $this_humanLL = $gene_id1;
my $this_mouseLL = $gene_id2;

print STDERR "check:\n";
foreach my $cdna ( @{$human_cdnas{$this_humanLL}} ){
  print STDERR "$this_humanLL\t$cdna\t".$human_cdna2protein{$cdna}."\n";
}
foreach my $cdna ( @{$mouse_cdnas{$this_mouseLL}} ){
  print STDERR "$this_mouseLL\t$cdna\t".$mouse_cdna2protein{$cdna}."\n";
}


my %already_seen;
my @human_cdnas = @{$human_cdnas{$this_humanLL}};
my @mouse_cdnas = @{$mouse_cdnas{$this_mouseLL}};

print STDERR "comparing $this_humanLL (".scalar(@human_cdnas)." transcripts ) ".
  "and $this_mouseLL (".scalar(@mouse_cdnas)." transcripts )\n";

my $transcript_map = Bio::EnsEMBL::Pipeline::GeneComparison::ObjectMap->new();
my $protein_map    = Bio::EnsEMBL::Pipeline::GeneComparison::ObjectMap->new();

my %seqobj;
my $seqfetcher = Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch->new();

my %target_coverage;
my %query_coverage;
my %target_gap;
my %query_gap;
my %perc_id;

my @human_cdna_seqs;
my %human_protein_seqs;
my %human_cdna_seqs;
my @mouse_cdna_seqs;
my %mouse_protein_seqs;
my %mouse_cdna_seqs;

HUMAN:
foreach my $humanNM ( @human_cdnas ){
  $humanNM =~/(\S+)\.\d*/;
  my $cdna_id = $1;
  my $humanNM_seq = $seqfetcher->get_Seq_by_acc($cdna_id);
  $humanNM_seq->display_id( $humanNM );
  unless ( $humanNM_seq ){
    print STDERR "Failed to find sequence $humanNM\n";
    next HUMAN;
  }
  push ( @human_cdna_seqs, $humanNM_seq );
  $human_cdna_seqs{ $humanNM } = $humanNM_seq;

  my $humanNP = $human_cdna2protein{$humanNM};
  $humanNP =~/(\S+)\.\d*/;
  my $prot_id = $1;
  my $humanNP_seq = $seqfetcher->get_Seq_by_acc($prot_id);
  $humanNP_seq->display_id( $humanNP );
  unless ( $humanNP_seq ){
    print STDERR "Failed to find sequence $humanNP\n";
    next HUMAN;
  }
  $human_protein_seqs{$humanNP} = $humanNP_seq;
}

MOUSE:
foreach my $mouseNM ( @mouse_cdnas ){
  $mouseNM =~/(\S+)\.\d*/;
  my $cdna_id = $1;
  my $mouseNM_seq = $seqfetcher->get_Seq_by_acc($cdna_id);
  $mouseNM_seq->display_id($mouseNM);
  unless ( $mouseNM_seq ){
    print STDERR "Failed to find sequence $mouseNM\n";
    next MOUSE;
  }
  push ( @mouse_cdna_seqs, $mouseNM_seq );
  $mouse_cdna_seqs{ $mouseNM } = $mouseNM_seq;
  my $mouseNP = $mouse_cdna2protein{$mouseNM};
  $mouseNP =~/(\S+)\.\d*/;
  my $prot_id = $1;
  my $mouseNP_seq = $seqfetcher->get_Seq_by_acc($prot_id);
  $mouseNP_seq->display_id( $mouseNP );
  unless ( $mouseNP_seq ){
    print STDERR "Failed to find sequence $mouseNP\n";
    next MOUSE;
  }
  $mouse_protein_seqs{$mouseNP} = $mouseNP_seq;
}

HUMAN:
foreach my $humanNM_seq ( @human_cdna_seqs ){
 MOUSE:
  foreach my $mouseNM_seq ( @mouse_cdna_seqs ){ 
    
    next if $already_seen{$humanNM_seq}{$mouseNM_seq};

    ############################################################
    # compare transcripts
    my ($score,$best_features, $target_coverage, $max_target_gap, $query_coverage, $max_query_gap, $perc_id) = 
      Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptPair
	->blast_unmapped_transcripts( $humanNM_seq, $mouseNM_seq);
    if ( $score && $best_features ){
      $transcript_map->match($humanNM_seq, $mouseNM_seq, $score );
      $perc_id{$humanNM_seq}{$mouseNM_seq} = $perc_id;
      $target_coverage{$humanNM_seq}{$mouseNM_seq} = $target_coverage;
      $target_gap{$humanNM_seq}{$mouseNM_seq} = $max_target_gap;
      $query_coverage{$humanNM_seq}{$mouseNM_seq} = $query_coverage;
      $query_gap{$humanNM_seq}{$mouseNM_seq} = $max_query_gap;
    }
    unless ($score){
      $score = 0;
    }
    my $humanNM = $humanNM_seq->display_id;
    my $mouseNM = $mouseNM_seq->display_id;
    my $humanNP = $human_cdna2protein{$humanNM};
    my $mouseNP = $mouse_cdna2protein{$mouseNM};
    print STDERR "Pair [ $humanNM($humanNP)  , $mouseNM($mouseNP) ] score = $score\n";
    
    ############################################################
    # compare peptides
    my $humanNP_seq = $human_protein_seqs{$humanNP};
    my $mouseNP_seq = $mouse_protein_seqs{$mouseNP};
    unless ( $humanNP_seq ){
      print STDERR "Failed to find sequence $humanNP\n";
      next HUMAN;
    }
    unless ( $mouseNP_seq ){
      print STDERR "Failed to find sequence $mouseNP\n";
      next MOUSE;
    }
    ($score,$best_features, $target_coverage, $max_target_gap, $query_coverage, $max_query_gap, $perc_id) = 
      Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptPair
	->blast_unmapped_proteins( $humanNP_seq, $mouseNP_seq);
    if ( $score && $best_features ){
      $protein_map->match($humanNP_seq, $mouseNP_seq, $score );
      $perc_id{$humanNP_seq}{$mouseNP_seq} = $perc_id;
      $target_coverage{$humanNP_seq}{$mouseNP_seq} = $target_coverage;
      $target_gap{$humanNP_seq}{$mouseNP_seq} = $max_target_gap;
      $query_coverage{$humanNP_seq}{$mouseNP_seq} = $query_coverage;
      $query_gap{$humanNP_seq}{$mouseNP_seq} = $max_query_gap;
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
my $transcript_pair_count = scalar($best_transcript_pairs_object->list1);
print STDERR "Transcript pairs created: ".$transcript_pair_count."\n";

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
    my $target_gap      = $target_gap{$element1}{$partner};
    my $query_coverage  = $query_coverage{$element1}{$partner};
    my $query_gap       = $query_gap{$element1}{$partner};
    my $length_diff     = $length1 - $length2;
    print STDERR "TRANSCRIPT_PAIR\t".
      "$this_humanLL\t$id1\t$protein_id1\tlength1:$length1\t".
	"coverage1:$target_coverage\tgap1:$target_gap\t".
	  "$this_mouseLL\t$id2\t$protein_id2\tlength2:$length2\t".
	    "coverage2:$query_coverage\tgap1:$query_gap\t".
	      "perc_id:$perc_id\t".
		"length_diff:$length_diff\n";
    
  }
}

############################################################
# protein pairs:
my %human_protein2cdna = reverse %human_cdna2protein;
my %mouse_protein2cdna = reverse %mouse_cdna2protein;


############################################################
# this requires a bit more of work - it can happen that
# 2 human proteins are exactly the same but have different transcripts
# so they have the same similarity to a mouse protein, leaving an ambiguity

# we need three graphs
# $stable  = with the stable (marriage) protein pairs (1-to-1)
# $correct = with the correct pairs (according to cdna pairs) (1-to-1)
# $all     = with all the comparisons (many-to-many)

# we bin the pairs in $all according to score
# from each bin we take the right pairs and if none
# we take the stable pairs

my $all     = $protein_map;
my $stable  = $protein_map->stable_marriage;
my $protein_pair_count = scalar($stable->list1);
print STDERR "Protein pairs created: ".$protein_pair_count."\n";
my $correct = Bio::EnsEMBL::Pipeline::GeneComparison::ObjectMap->new();
# build the correct pairs
foreach my $element1 ( $best_transcript_pairs_object->list1 ){
  foreach my $partner ( $best_transcript_pairs_object->partners( $element1 ) ){
    my $humanNP = $human_cdna2protein{ $element1->display_id};
    my $mouseNP = $mouse_cdna2protein{ $partner->display_id};
    my $humanNP_seq = $human_protein_seqs{$humanNP};
    my $mouseNP_seq = $mouse_protein_seqs{$mouseNP};
    $correct->match($humanNP_seq, $mouseNP_seq, $all->score($humanNP_seq,$mouseNP_seq) );
  }
}

my %bin_pairs;
foreach my $h ( $all->list1 ){
  foreach my $m ( $all->partners( $h ) ){
    my $score = $all->score($h,$m);
    push ( @{$bin_pairs{$score}}, [$h,$m] );
  }
}

my %seen;
my @selected_pairs;
foreach my $score ( sort { $b <=> $a } keys %bin_pairs ){
  my $has_correct = 0;
  foreach my $pair ( @{$bin_pairs{$score}} ){
    if ( $correct->score( $pair->[0], $pair->[1] ) ){
      $has_correct = 1;
    }
  }
  if ( $has_correct ){
    foreach my $pair ( @{$bin_pairs{$score}} ){
      if ( $correct->score( $pair->[0], $pair->[1] ) ){
	unless ( $seen{ $pair->[0] } || $seen{ $pair->[1] } ){
	  push ( @selected_pairs, $pair );
	  $seen{ $pair->[0] } = 1;
	  $seen{ $pair->[1] } = 1;
	}
      }
    }
  }
  else{
    foreach my $pair ( @{$bin_pairs{$score}} ){
      if ( $stable->score( $pair->[0], $pair->[1] ) ){
	unless ( $seen{ $pair->[0] } || $seen{ $pair->[1] } ){
	  push ( @selected_pairs, $pair );
	  $seen{ $pair->[0] } = 1;
	  $seen{ $pair->[1] } = 1;
	}
      }
    }
  }
}

print STDERR "pairs selected:\n";
foreach my $pair ( @selected_pairs ){
  my $humanNP = $pair->[0]->display_id;
  my $mouseNP = $pair->[1]->display_id;
  my $humanNM = $human_protein2cdna{ $humanNP };
  my $mouseNM = $mouse_protein2cdna{ $mouseNP };
  #print STDERR "$humanNP ($humanNM) <--> $mouseNP ($mouseNM)\n";
  
  ############################################################
  # summary line
  my $length1 = $pair->[0]->length;
  my $length2 = $pair->[1]->length;
  my $perc_id         = $perc_id{$pair->[0]}{$pair->[1]};
  my $target_coverage = $target_coverage{$pair->[0]}{$pair->[1]};
  my $target_gap      = $target_gap{$pair->[0]}{$pair->[1]};
  my $query_coverage  = $query_coverage{$pair->[0]}{$pair->[1]};
  my $query_gap       = $query_gap{$pair->[0]}{$pair->[1]};
  my $length_diff      = $length1 - $length2;
  print STDERR "PROTEIN_PAIR\t".
    "$this_humanLL\t$humanNM\t$humanNP\tlength1:$length1\t".
      "coverage1:$target_coverage\tgap1:$target_gap\t".
	"$this_mouseLL\t$mouseNM\t$mouseNP\tlength2:$length2\t".
	  "coverage2:$query_coverage\tgap2:$query_gap\t".
	    "perc_id:$perc_id\t".
	      "length_diff:$length_diff\n";

  ############################################################
  # do cdna alignments agree with the protein ones?
  my $humanNM_seq = $human_cdna_seqs{ $humanNM }; 
  my $mouseNM_seq = $mouse_cdna_seqs{ $mouseNM };
  if ( $best_transcript_pairs_object->score($humanNM_seq,$mouseNM_seq) ){
    print STDERR "PAIR_MATCH\n";
  }
  else{
    foreach my $actual_mouseNM_seq ( $best_transcript_pairs_object->partners($humanNM_seq) ){
      my $this_humanNM = $humanNM_seq->display_id;
      my $this_mouseNM = $actual_mouseNM_seq->display_id;
      my $this_humanNP = $human_cdna2protein{$this_humanNM};
      my $this_mouseNP = $mouse_cdna2protein{$this_mouseNM};
      print STDERR "PAIR_MISMATCH\t".
	"$this_humanLL\t$this_humanNM ($this_humanNP)\t$humanNP ($humanNM)\t".
	  "$this_mouseLL\t$this_mouseNM ($this_mouseNP)\tmouseNM\t$mouseNP ($mouseNM)\n";
    }
    foreach my $actual_humanNM_seq ( $best_transcript_pairs_object->partners($mouseNM_seq) ){
      my $this_humanNM = $actual_humanNM_seq->display_id;
      my $this_mouseNM = $mouseNM_seq->display_id;
      my $this_humanNP = $human_cdna2protein{$this_humanNM};
      my $this_mouseNP = $mouse_cdna2protein{$this_mouseNM};
      print STDERR "PAIR_MISMATCH\t".
	"$this_humanLL\t$this_humanNM ($this_humanNP)\t$humanNP ($humanNM)\t".
	  "$this_mouseLL\t$this_mouseNM ($this_mouseNP)\tmouseNM\t$mouseNP ($mouseNM)\n";
    }
  }

  ############################################################
  # check the missed transcript pair
  unless ( $best_transcript_pairs_object->partners($humanNM_seq) 
	   ||
	   $best_transcript_pairs_object->partners($mouseNM_seq)
	 ){
    my $humanNP = $human_cdna2protein{ $humanNM };
    my $mouseNP = $mouse_cdna2protein{ $mouseNM };
    print STDERR "TRANS_PAIR_MISSED\t".
      "$this_humanLL\t$humanNM ($humanNP)\t".
	"$this_mouseLL\t$mouseNM ($mouseNP)\n";
  }
}

foreach my $humanNM_seq ( $best_transcript_pairs_object->list1 ){
  foreach my $mouseNM_seq ( $best_transcript_pairs_object->partners( $humanNM_seq ) ){
    my $humanNM = $humanNM_seq->display_id;
    my $mouseNM = $mouseNM_seq->display_id;
    my $humanNP = $human_cdna2protein{ $humanNM };
    my $mouseNP = $mouse_cdna2protein{ $mouseNM };
    my $humanNP_seq = $human_protein_seqs{$humanNP};
    my $mouseNP_seq = $mouse_protein_seqs{$mouseNP};
    
    my $found = 0;
    foreach my $pair ( @selected_pairs ){
      if ( $pair->[0] == $humanNP_seq && $pair->[1] == $mouseNP_seq ){
	$found = 1;
      }
    }
    unless ( $found ){
      print STDERR "PROT_PAIR_MISSED\t".
	"$this_humanLL\t$humanNP($humanNM)\t".
	  "$this_mouseLL\t$mouseNP($mouseNM)\n";
    }
  }
}

if ( $protein_pair_count == 0 && $transcript_pair_count == 0 ){
  print STDERR "check:\n";
  foreach my $cdna ( @{$human_cdnas{$this_humanLL}} ){
    foreach my $mouse_cdna ( @{$mouse_cdnas{$this_mouseLL}} ){
      print STDERR "EMPTY_ORTHOLOGOUS_PAIR\t".
	"$this_humanLL\t$cdna\t".$human_cdna2protein{$cdna}."\t".
	  "$this_mouseLL\t$mouse_cdna\t".$mouse_cdna2protein{$mouse_cdna}."\n";
    }
  }
}
