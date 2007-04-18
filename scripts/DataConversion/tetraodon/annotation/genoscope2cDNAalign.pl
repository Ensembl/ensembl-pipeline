#!/usr/local/ensembl/bin/perl

# custom script to import the SQL dump of the Genoscope cDNA alignments
# and write them as DnaDnaAlignFeatures.

 
use strict;
 
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::DnaDnaAlignFeature;


my (
    $dbname,
    $dbhost,
    $dbuser,
    $dbport,
    $dbpass,
    $logic_name,
    $per_exon,
    $verbose,
    $debug,
    $log_file,
    $id_map_file,
    $id_map,
    %slices, 
    @all_features,
    );

&GetOptions(
            'dbname=s' => \$dbname,
            'dbuser=s' => \$dbuser,
            'dbhost=s' => \$dbhost,
            'dbport=s' => \$dbport,
            'dbpass=s' => \$dbpass,
            'perexon=s' => \$per_exon,
            'verbose'  => \$verbose,
            'debug'    => \$debug,
            'log=s'    => \$log_file,
            'idmapfile=s' => \$id_map_file,
);


my $logic_name = shift;
my $cdna_align_file = shift;

die "Usage: $0 cDNA_align_file.sql\n" 
    if not $logic_name or not $cdna_align_file;

if (defined $id_map_file) {
  $id_map = &process_exon_file($id_map_file);
}

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
	'-dbname' => $dbname,
	'-host' => $dbhost,
	'-user' => $dbuser,
	'-port' => $dbport,
	'-pass' => $dbpass
);


my $ana_obj;
if (not defined($ana_obj = $db->get_AnalysisAdaptor->fetch_by_logic_name($logic_name))) {
  $ana_obj  = Bio::EnsEMBL::Analysis->new(
                                         -logic_name      => $logic_name,
                                         -gff_source      => "Genoscope", 
                                         -gff_feature     => 'similarity');
}


$log_file and do {
  my $fname = $log_file;
  undef $log_file;
  open $log_file, ">$fname" or die "Could not open log file  $fname for writing\n";
};

$verbose and print STDERR "Processing alignment file for alignments...\n";

my $count = 0;
open(ALIGNS, $cdna_align_file) or die "Could not open 'cdna_align_file' for reading\n";
while(<ALIGNS>) {
  next unless /\S/;

  my @l = split;

  my $align_id = $l[0];

  if (defined $id_map) {
    next if not exists($id_map->{$align_id});
    $align_id = $id_map->{$align_id};
  }

  my $chr_name = $l[3]; $chr_name =~ s/^chr//;  
  my $cdna_name = $l[9];
  my ($chr_start, $chr_end, $chr_strand) = ($l[4], $l[5], $l[6]);
  my ($cdna_start, $cdna_end, $cdna_strand) = ($l[10], $l[11], $l[12]);
  my ($score, $pid) = ($l[15], $l[19]);
  my ($chr_gaps, $cdna_gaps, $intron_pos) = ($l[21], $l[22], $l[23]);

  $chr_strand = "+" if $chr_strand == 1;
  $chr_strand = "-" if $chr_strand == -1;
  
  
  if (not $slices{$chr_name}) {
    $slices{$chr_name} = $db->get_SliceAdaptor->fetch_by_region('chromosome', $chr_name);
  }

  $count++;

  # parse the intron information to extract the 
  my @align_blocks = &create_alignment($chr_gaps, $chr_end - $chr_start + 1, 
                                       $cdna_gaps, $cdna_end - $cdna_start + 1, 
                                       $intron_pos);

  
  if (not &transform_alignment( $chr_start, $chr_end, $chr_strand,
                                $cdna_start, $cdna_end, $cdna_strand, 
                                @align_blocks )) {
    print $log_file "Error parsing alignment $count\n";
    print $log_file join("\t", @l), "\n";
    &print_alignment($log_file, @align_blocks);
    next;
  }

  #$debug and do {
  #  print $log_file "GOOD ALIGNMENT\n";
  #  &print_alignment($log_file, @align_blocks);
  #};

  foreach my $exon (@align_blocks) {
    my @features;
    foreach my $block (@$exon) {
      my $fp = Bio::EnsEMBL::FeaturePair->new;
      $fp->seqname($chr_name);
      $fp->start($block->{gen_start});
      $fp->end($block->{gen_end});
      $fp->strand( $chr_strand eq "-" ? -1 : 1 );
      $fp->hseqname( $align_id );
      $fp->hstart($block->{rna_start});
      $fp->hend($block->{rna_end});
      $fp->hstrand($cdna_strand eq "-" ? -1 : 1 );
      $fp->score($score);
      $fp->percent_id($pid);

      push @features, $fp; 
    }
        
    my $align = Bio::EnsEMBL::DnaDnaAlignFeature->new(-features => \@features);

    $align->slice( $slices{$chr_name} );
    $align->analysis( $ana_obj );

    push @all_features, $align;
  }

  if ($count % 1000 == 0) {
    $verbose and print STDERR "Writing alignments up to and including $count...\n";
    $db->get_DnaAlignFeatureAdaptor->store(@all_features) if not $debug;
    @all_features = ();
  }

}
close(ALIGNS);

if (@all_features) {
  $verbose and print STDERR "Writing alignments up to and including $count...\n";
  $db->get_DnaAlignFeatureAdaptor->store(@all_features) if not $debug;
}



sub create_alignment {
  my ($gen_gaps,
      $gen_len,
      $rna_gaps,
      $rna_len,
      $intron_pos) = @_;

  my @introns;
  if ($intron_pos ne "NULL") {
    @introns = map { $_ =~ /(\d+)I(\d+)/;  [$1,$2]; } split /,/,$intron_pos;
  }
  my @inserts;
  if ($rna_gaps ne "NULL") {
    @inserts = map { $_ =~ /(\d+)\#(\d+)/; [$1,$2]; } split /,/,$rna_gaps;
  }
  my @deletes;
  if ($gen_gaps ne "NULL") {
    @deletes = map { $_ =~ /(\d+)\#(\d+)/; [$1,$2]; } split /,/,$gen_gaps;
  }

  # treat introns as inserts
  #push @inserts, @introns;
  #@introns = ();
  #@inserts = sort { $a->[0] <=> $b->[0] } @inserts;

  my @matches;

  my ($current_gen, $current_rna) = (1,1);
  # assume that we start with a match; we should do!
  push @matches, [ { gen_start => $current_gen, rna_start => $current_rna } ];

  while($current_gen <= $gen_len and $current_rna <= $rna_len) {
    $matches[-1]->[-1]->{gen_end} = $current_gen;
    $matches[-1]->[-1]->{rna_end} = $current_rna;

    # look for introns, inserts and deletes in this position and 
    # scan to the next match column; should only match one of the 
    # following 3 cases, if any
        
    if (@introns and ($introns[0]->[0] == $current_gen)) {

      # could be a delete at this position before the intron
      if (@deletes and $deletes[0]->[0] <= $current_gen) {
        # intron and genomic deletion at same position. 
        my $delete = shift @deletes;
        $current_rna += $delete->[1];
      }

      my $intron = shift @introns;

      if (not @inserts or $inserts[0]->[0] > $current_rna) {
        # finally, there could be a DELETE at the *end* of the intron!
        $current_gen += $intron->[1];

        if (@deletes and $deletes[0]->[0] <= $current_gen) {
          my $delete = shift @deletes; 
          $current_rna += $delete->[1];
        }
        $current_gen++;

      } else {
        my $insert = shift @inserts;
        $current_gen += $intron->[1] + $insert->[1] + 1;
      }

      $current_rna++;

      if ($per_exon) {
        push @matches, [{gen_start => $current_gen, rna_start => $current_rna}];
      } else {
        push @{$matches[-1]}, {gen_start => $current_gen, rna_start => $current_rna };
      }

    }
    elsif (@inserts and $inserts[0]->[0] == $current_rna) {
      my $insert = shift @inserts;
      $current_gen += $insert->[1] + 1;
      $current_rna++;

      push @{$matches[-1]}, { gen_start => $current_gen, rna_start=> $current_rna };

    } 
    elsif (@deletes and $deletes[0]->[0] == $current_gen) {
      my $delete = shift @deletes;
      $current_rna += $delete->[1] + 1;
      $current_gen++;

      push @{$matches[-1]}, { gen_start => $current_gen, rna_start => $current_rna };
    }
    else {      
      $current_gen++;
      $current_rna++;
    }
  }

  return @matches;
}




sub transform_alignment {
  my ($gen_start, $gen_end, $gen_strand,
      $rna_start, $rna_end, $rna_strand,
      @align_blocks) = @_;

  foreach my $exon (@align_blocks) {
    foreach my $block (@$exon) {
      my ($local_gen_start, $local_gen_end) = ($block->{gen_start}, $block->{gen_end});
      my ($local_rna_start, $local_rna_end) = ($block->{rna_start}, $block->{rna_end});

      if ($gen_strand eq "+") {
        $local_gen_start += $gen_start - 1;
        $local_gen_end   += $gen_start - 1;
      } else {
        my $save_local_gen_start = $local_gen_start;
        $local_gen_start = $gen_end - $local_gen_end   + 1;
        $local_gen_end   = $gen_end - $save_local_gen_start + 1;
      }

      if ($rna_strand eq "+") {
        $local_rna_start += $rna_start - 1;
        $local_rna_end   += $rna_start - 1;
      } else {
        my $save_local_rna_start = $local_rna_start;
        $local_rna_start = $rna_end - $local_rna_end   + 1;
        $local_rna_end   = $rna_end - $save_local_rna_start + 1;
      }

      $block->{gen_start} = $local_gen_start;
      $block->{gen_end} = $local_gen_end;
      $block->{rna_start} = $local_rna_start;
      $block->{rna_end} = $local_rna_end;
    }
  }

  # quick sanity check on the alignment

  my $test_gen_start = ($gen_strand eq "+" ) ? 
      $align_blocks[0]->[0]->{'gen_start'} :
      $align_blocks[-1]->[-1]->{'gen_start'};

  my $test_gen_end = ($gen_strand eq "+" ) ? 
      $align_blocks[-1]->[-1]->{'gen_end'} :
      $align_blocks[0]->[0]->{'gen_end'};

  my $test_rna_start = ($rna_strand eq "+" ) ? 
      $align_blocks[0]->[0]->{'rna_start'} :
      $align_blocks[-1]->[-1]->{'rna_start'};

  my $test_rna_end = ($rna_strand eq "+" ) ? 
      $align_blocks[-1]->[-1]->{'rna_end'} :
      $align_blocks[0]->[0]->{'rna_end'};

  if ($test_gen_start != $gen_start or
      $test_gen_end   != $gen_end   or
      $test_rna_start != $rna_start or
      $test_rna_end   != $rna_end) {

    return 0;
  } else {
    return 1;
  }

}



sub process_exon_file {
  my $exon_file = shift;

  my %id_map;

  $verbose and print STDERR "Processing exon file for id map...\n";

  open(EXONS, $exon_file) or die "Could not open '$exon_file' for reading;\n";
  while(<EXONS>) {
    next unless /\S/;
    my @l = split;
    
    next if $l[9] == 0;
    
    if (exists($id_map{$l[0]})) {
      # if strand is positive, take this one;
      if ($l[2] eq "+") {
        $id_map{$l[0]} = $l[12];
      }
      next;        
    }
    $id_map{$l[0]} = $l[12];
  }
  close(EXONS);
  
# should be 1-1
  %id_map = reverse %id_map;

  return \%id_map;
}



sub print_alignment {
  my $fh = shift;
  my @blocks = @_;

  foreach my $exon (@blocks) {
    print $fh "Exon:\n";
    foreach my $block (@$exon) {
      printf($fh "chr_start: %d chr_end: %d rna_start: %d rna_end: %d\n",
             $block->{gen_start}, $block->{gen_end},
             $block->{rna_start}, $block->{rna_end});
    }
  }

}
