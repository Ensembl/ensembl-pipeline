#!/usr/local/ensembl/bin/perl

# load_align_features.pl
#
# Loads a 'extended' GFF file containing DNA/Protein align features 
# Alignments themselves will not be stored. 

use strict;
use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::SimpleFeature;

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

  my @l = split /\t/;
  if ((not defined($gff_source) or $l[1] eq $gff_source) and
      (not defined($gff_feature) or $l[2] eq $gff_feature)) {

    my ($slice, $analysis, $f);

    my ($display_id)  = $l[8] =~ /display_id\s+\"([^\"]+)\"/;
    if (not $display_id) {
      die "Line found without display_id\n";
    }
    
    if (not ($analysis = $analyses{$l[1]})) {
      if (not defined($analysis = $db->get_AnalysisAdaptor->fetch_by_logic_name($l[1]))) {
        print STDERR "Warning; analysis not found in the database; creating\n";
        $analysis  = Bio::EnsEMBL::Analysis->new(
                                                 -logic_name      => $l[1],
                                                 -gff_source      => $l[1],
                                                 -gff_feature     => $l[2]);
        $analyses{$l[1]} = $analysis;
        
      }
    }

    if (not ($slice = $slices{$l[0]})) {
      $slice = $db->get_SliceAdaptor->fetch_by_region('chromosome', $l[0]);
      $slices{$l[0]} = $slice;
    }

    $f = Bio::EnsEMBL::SimpleFeature->new();    
    $f->score($l[5]);
    $f->start($l[3]);
    $f->end($l[4]);
    if ($l[6] eq "+") {
      $f->strand(1);
    } elsif ($l[6] eq "-") {
      $f->strand(-1);
    } else {
      $f->strand(0);
    }
    $f->display_label($display_id);
    $f->slice($slice);
    $f->analysis($analysis);
    push @all_features, $f;
  }
}

print STDERR "Created ", scalar(@all_features), " features; writing...\n";
$db->get_SimpleFeatureAdaptor->store(@all_features);


