#!/usr/local/ensembl/bin/perl

# load_align_features.pl
#
# Loads a 'extended' GFF file containing DNA/Protein align features 
# Alignments themselves will not be stored. 

use strict;
use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DnaDnaAlignFeature;


my (
    $dbname,
    $dbhost,
    $dbuser,
    $dbport,
    $dbpass,
    $gff_feature, 
    $gff_source
);


&GetOptions(
            'dbname=s' => \$dbname,
            'dbuser=s' => \$dbuser,
            'dbhost=s' => \$dbhost,
            'dbport=s' => \$dbport,
            'dbpass=s' => \$dbpass,
            'gff_feat=s' => \$gff_feature,
            'gff_source=s' => \$gff_source,
);

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
	'-dbname' => $dbname,
	'-host' => $dbhost,
	'-user' => $dbuser,
	'-port' => $dbport,
	'-pass' => $dbpass
);


my (%slices, %analyses, @all_features);

while (<>) {
  /^\#/ and next; 

  my @line = split /\t/;
  if ((not defined($gff_source) or $line[1] eq $gff_source) and
      (not defined($gff_feature) or $line[2] eq $gff_feature)) {

    my ($slice, $analysis, $f);

    my ($hitname)  = $line[8] =~ /hit_name\s+\s*\"([^\"]+)\"/;
    my ($hitstart) = $line[8] =~ /hit_start\s+\"([^\"]+)\"/;
    my ($hitend) = $line[8]   =~ /hit_end\s+\"([^\"]+)\"/;
    my ($hitstrand) = $line[8]   =~ /hit_strand\s+\"([^\"]+)\"/;
    if ($hitstrand eq "+") {
      $hitstrand = 1;
    } elsif ($hitstrand eq "-") {
      $hitstrand = -1;
    } else {
      # no strand, but we must set it to + because it seems you canot have strandless FeaturePairs
      $hitstrand = 1;
    }

    if (not $hitname or not $hitstart or not $hitend) {
      die "Line found without hit_name, hit_start and hit_end\n";
    }
    
    if (not ($slice = $slices{$line[0]})) {
      $slice = $db->get_SliceAdaptor->fetch_by_region('chromosome', $line[0]);
      $slices{$line[0]} = $slice;
    }

    if (not ($analysis = $analyses{$line[1]})) {
      if (not defined($analysis = $db->get_AnalysisAdaptor->fetch_by_logic_name($line[1]))) {
        $analysis  = Bio::EnsEMBL::Analysis->new(
                                                 -logic_name      => $line[1],
                                                 -gff_source      => $line[1],
                                                 -gff_feature     => 'similarity');
        
      }
      $analyses{$line[1]} = $analysis;
    }

    $f = Bio::EnsEMBL::FeaturePair->new();
    
    $f->score($line[5]);
    $f->start($line[3]);
    $f->end($line[4]);
    if ($line[6] eq "+") {
      $f->strand(1);
    } elsif ($line[6] eq "-") {
      $f->strand(-1);
    } else {
      # no strand, but we must set it to + because it seems you canot have strandless FeaturePairs
      $f->strand(1);
    }
    $f->hseqname($hitname);
    $f->hstart($hitstart);
    $f->hend($hitend);
    $f->hstrand($hitstrand);
    $f->slice($slice);
    $f->analysis($analysis);

    my $align_f = Bio::EnsEMBL::DnaDnaAlignFeature->new(-features => [$f]);

    push @all_features, $align_f;
  }
}



print STDERR "Created ", scalar(@all_features), " features; writing...\n";
$db->get_DnaAlignFeatureAdaptor->store(@all_features);

