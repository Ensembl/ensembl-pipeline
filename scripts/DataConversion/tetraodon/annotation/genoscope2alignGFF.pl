#!/usr/local/bin/perl
 
use strict;
 
use Getopt::Long;

my (@sources, %sources, %aligns);

&GetOptions("source=s@" => \@sources);
 
die "You must provide at least one source with -source" if not @sources;

map { $sources{$_} = 1 } @sources;

my (%aligns);

while(<>) {
  /^\#/ and next;

  my @l = split /\t/;

  $l[0] =~ s/^chr//;

  if (exists $sources{$l[1]}) {
    if ($l[2] eq "match" or $l[2] eq "transcript") {
      my ($align_id) = $l[8] =~ /$l[2]\s+(\S+)/;

      $aligns{$align_id}->{'chr'} = $l[0];
      $aligns{$align_id}->{'start'} = $l[3];
      $aligns{$align_id}->{'end'} = $l[4];
      $aligns{$align_id}->{'score'} = $l[5];
      $aligns{$align_id}->{'strand'} = $l[6];
      if ($l[1] =~ /^EG3/) {
        $aligns{$align_id}->{'strand'} = ".";
      }

      $aligns{$align_id}->{'source'} = $l[1];
      $aligns{$align_id}->{'feature'} = $l[2];
      if ($l[2] eq "transcript") {
        $aligns{$align_id}->{'is_protein'} = 1;
      }
    } elsif ($l[2] eq "HSP" or $l[2] = "exon" or $l[2]) {
      my ($align_id) = $l[8] =~ /$l[2]\s+(\S+)/;      

      my $exon = {
        chr      => $l[0],
        source   => $l[1],
        feature  => $l[2],            
        start    => $l[3],
        end      => $l[4],
        score    => $l[5],
        strand   => $l[6],
      };
      if ($l[1] =~ /^EG3/) {
        $exon->{'strand'} = ".";
      }

      push @{$aligns{$align_id}->{'exons'}}, $exon; 
    }    
  }

} 

foreach my $aid (keys %aligns) {
  my $strand = $aligns{$aid}->{'strand'};
  my @exons = sort {
    $a->{'start'} <=> $b->{'start'} 
  } @{$aligns{$aid}->{'exons'}};

  if ($strand eq "-") {
    @exons = reverse @exons;
  }
  
  my $hit_start = 1;
  foreach my $ex (@exons) {
    my $ex_len = $ex->{'end'} - $ex->{'start'} + 1;

    $ex->{'hit_strand'} = ".";
    $ex->{'hit_start'} = $hit_start;
    $ex->{'hit_end'} = $hit_start + $ex_len - 1;
    $ex->{'target_name'} = $aligns{$aid}->{'target_name'};
    $ex->{'hit_name'} = $aid;

    $hit_start += $ex_len;
  }
}

foreach my $aid (sort { 
  $aligns{$a}->{'source'} cmp $aligns{$b}->{'source'} or
  $aligns{$a}->{'chr'} cmp $aligns{$b}->{'chr'} or
  $aligns{$a}->{'start'} <=> $aligns{$b}->{'start'} } keys %aligns) {

  foreach my $ex (sort {$a->{'start'} <=> $b->{'start'}} @{$aligns{$aid}->{'exons'}}) {
    &print_gff($ex, "hit_name", "hit_start", "hit_end");
  }
}





sub print_gff {
  my ($gff_entry, @tag_list) = @_;

  printf("%s\t%s\t%s\t%d\t%d\t%.2f\t%s\t.",
         $gff_entry->{'chr'},
         $gff_entry->{'source'},
         $gff_entry->{'feature'},
         $gff_entry->{'start'},
         $gff_entry->{'end'},
         $gff_entry->{'score'},
         $gff_entry->{'strand'});

  if (@tag_list) {
    print "\t";
    my $first = 1;

    foreach my $tag (@tag_list) {
      printf("%s%s \"%s\";", 
             $first ? "" : " ",
             lc($tag), 
             $gff_entry->{$tag});
      $first = 0;
    }
  }
  print "\n";

}
