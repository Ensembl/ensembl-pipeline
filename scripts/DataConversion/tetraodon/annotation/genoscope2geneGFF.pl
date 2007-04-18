#!/usr/local/ensmbl/bin/perl

use strict;

use Getopt::Long;

my (%transcripts, @sources, %sources);

&GetOptions("sources=s@" => \@sources);

@sources = qw(GSTEN HOX CYT) if not @sources;
map { $sources{$_} = 1 } @sources;

while(<>) {
  /^\#/ and next;

  my @l = split /\t/;
  $l[6] = "+" if $l[6] == 1;
  $l[6] = "-" if $l[6] == -1;

  if ($sources{$l[1]}) {
    if ($l[2] eq "gene") {
      my ($gene_id) = $l[8] =~ /Gene\s+(\S+)/;

      my ($tran_id, $pep_id);
      if ($l[8] =~ /Note\s+(\S+)\s*;\s*Note\s+(\S+)/) {
        ($tran_id, $pep_id) = ($1, $2);
      }

      if (defined $tran_id) {
        # store global information if first time we have seen ref to this transcript. 
        if (not $transcripts{$tran_id}) {
          $l[0] =~ s/^chr//; 
          $transcripts{$tran_id}->{chromosome} = $l[0];
          $transcripts{$tran_id}->{strand} = $l[6];
          $transcripts{$tran_id}->{type} = $l[1];
          $transcripts{$tran_id}->{tran_id} = $tran_id;
        }
        $transcripts{$tran_id}->{gene_id} = $gene_id;
        $transcripts{$tran_id}->{pep_id} = $pep_id;
      }
    }
    elsif ($l[2] eq "mRNA") {
      my ($tran_id) = $l[8] =~ /mRNA\s+(\S+)/;

      # store global information if first time we have seen ref to this transcript. 
      if (not $transcripts{$tran_id}) {
        $l[0] =~ s/^chr//; 
        $transcripts{$tran_id}->{chromosome} = $l[0];
        $transcripts{$tran_id}->{strand} = $l[6];
        $transcripts{$tran_id}->{type} = $l[1];
        $transcripts{$tran_id}->{tran_id} = $tran_id;
      }

      my ($gene_id) = $l[8] =~ /Gene\s+(\S+)/;
      $transcripts{$tran_id}->{gene_id} = $gene_id;
      $transcripts{$tran_id}->{start} = $l[3];
      $transcripts{$tran_id}->{end} = $l[4];
      $transcripts{$tran_id}->{score} = $l[5];

    }
    elsif ($l[2] eq "transcript") {
      my ($tran_id) = $l[8] =~ /transcript\s+(\S+)/;

      # store global information if first time we have seen ref to this transcript. 
      if (not $transcripts{$tran_id}) {
        $l[0] =~ s/^chr//; 
        $transcripts{$tran_id}->{chromosome} = $l[0];
        $transcripts{$tran_id}->{strand} = $l[6];
        $transcripts{$tran_id}->{type} = $l[1];
        $transcripts{$tran_id}->{tran_id} = $tran_id;
      }

      $transcripts{$tran_id}->{start} = $l[3];
      $transcripts{$tran_id}->{end} = $l[4];
      $transcripts{$tran_id}->{score} = $l[5];
      if ($l[8] =~ /transcript\s+\S+\s*;\s*Note\s+\"(\S+)\"/) {
        $transcripts{$tran_id}->{source} = $1;
      }
    }
    elsif ($l[2] eq "UTR") {
      my ($tran_id) = $l[8] =~ /UTR\s+(\S+)/;
      my ($exon_id) = $l[8] =~ /UTR\s+\S+\s*;\s*Note\s+(\S+)/;
      
      # store global information if first time we have seen ref to this transcript. 
      if (not $transcripts{$tran_id}) {
        $l[0] =~ s/^chr//; 
        $transcripts{$tran_id}->{chromosome} = $l[0];
        $transcripts{$tran_id}->{strand} = $l[6];
        $transcripts{$tran_id}->{type} = $l[1];
        $transcripts{$tran_id}->{tran_id} = $tran_id;
      }
      
      my $utr_exon = { start  => $l[3],
                       end    => $l[4],
                       score  => $l[5],
                     };
      if (defined $exon_id) {
        $utr_exon->{stable_id} = $exon_id;
      }
      
      push @{$transcripts{$tran_id}->{utr_regions}}, $utr_exon;
    }
    elsif ($l[2] eq "CDS") {
      my ($tran_id) = $l[8] =~ /CDS\s+(\S+)/;
      my ($exon_id) = $l[8] =~ /CDS\s+\S+\s*;\s*Note\s+(\S+)/;

      # store global information if first time we have seen ref to this transcript. 
      if (not $transcripts{$tran_id}) {
        $l[0] =~ s/^chr//; 
        $transcripts{$tran_id}->{chromosome} = $l[0];
        $transcripts{$tran_id}->{strand} = $l[6];
        $transcripts{$tran_id}->{type} = $l[1];
        $transcripts{$tran_id}->{tran_id} = $tran_id;
      }

      my $cds_exon = { start  => $l[3],
                       end    => $l[4],
                       score  => $l[5],
                       phase  => $l[7],
                     };
      if (defined $exon_id) {
        $cds_exon->{stable_id} => $exon_id;
      }

      $transcripts{$tran_id}->{total_bps} += $cds_exon->{end} - $cds_exon->{start} + 1;
      push @{$transcripts{$tran_id}->{cds_regions}}, $cds_exon;

    }
    elsif ($l[2] eq "exon") {
      my ($tran_id) = $l[8] =~ /exon\s+(\S+)/;

      # store global information if first time we have seen ref to this transcript. 
      if (not $transcripts{$tran_id}) {
        $l[0] =~ s/^chr//; 
        $transcripts{$tran_id}->{chromosome} = $l[0];
        $transcripts{$tran_id}->{strand} = $l[6];
        $transcripts{$tran_id}->{type} = $l[1];
        $transcripts{$tran_id}->{tran_id} = $tran_id;
      }
            
      my $cds_exon = { start  => $l[3],
                       end    => $l[4],
                       score  => $l[5],
                       phase => $l[7],
                     };
      $transcripts{$tran_id}->{'total_bps'} += $cds_exon->{'end'} - $cds_exon->{'start'} + 1;
      push @{$transcripts{$tran_id}->{'cds_regions'}}, $cds_exon;
      
    }
  }
}

foreach my $tran_id (keys %transcripts) {
  my $tran = $transcripts{$tran_id};

  my @cds_regions = sort { $a->{start} <=> $b->{start} } @{$tran->{cds_regions}};
  my (@exons);
  
  if (exists($tran->{utr_regions})) {
    # meed to merge them into "exon"
    my @utr_regions = sort { $a->{start} <=> $b->{start} } @{$tran->{utr_regions}};
    my @tmp_exons = sort { $a->{start} <=> $b->{start} } (@utr_regions, @cds_regions);
    foreach my $ex (@tmp_exons) {
      if (not @exons or $exons[-1]->{end} + 1 < $ex->{start}) {

        push @exons, {
          start     => $ex->{start},
          end       => $ex->{end},
          score     => $ex->{score},
          stable_id => $ex->{stable_id}
        };
      } else {
        # merge; check the ids. 
        die "Inconsistent exon portions (" . $exons[-1]->{stable_id} . ", " . $ex->{stable_id} . ")\n" 
            if $exons[-1]->{stable_id} ne $ex->{stable_id}; 
        $exons[-1]->{end} = $ex->{end};
      }
    }

  } else {
    @exons = @cds_regions;
  }
  
  $tran->{'exons'} = \@exons;
  
  # finally, calculate phases for the CDS exons; GTF/GFF convention, NOT ensembl
  if ($tran->{strand} eq "-") {
    @cds_regions = reverse @cds_regions;
  }

  if ($cds_regions[0]->{phase} eq ".") {
    my $next_phase = 0;

    foreach my $exon (@cds_regions) {
      if ($exon->{phase} eq ".") {
        $exon->{phase} = $next_phase;
        $next_phase = (3 - ( ($exon->{end} - $exon->{start} + 1 - $next_phase) % 3)) % 3;
      }
    }
  }
}


foreach my $k (sort { 
  $transcripts{$a}->{type} cmp $transcripts{$b}->{type} or
      $transcripts{$a}->{chromosome} cmp $transcripts{$b}->{chromosome} or
      $transcripts{$a}->{start} <=> $transcripts{$b}->{start} } keys %transcripts) {
  
  &print_transcript($transcripts{$k});

}


sub print_transcript {
  my $tran = shift;

  my $chr = $tran->{chromosome};
  my $type = $tran->{type};
  my $strand = $tran->{strand};

  my $t_id = $tran->{tran_id};
  if (exists $tran->{source}) {
    $t_id .= ":" . $tran->{source};
  }
  
  my $g_id = $tran->{gene_id};
  $g_id = $t_id if not $g_id;

  my $p_id = $tran->{pep_id};
  $p_id = $t_id if not $p_id;

  my (@exons, @cds);

  if ($strand eq "-") {
    @exons = sort { $b->{start} <=> $a->{start} } @{$tran->{exons}};
    @cds =   sort { $b->{start} <=> $a->{start} } @{$tran->{cds_regions}};
  } else {    
    @exons = sort { $a->{start} <=> $b->{start} } @{$tran->{exons}};
    @cds =   sort { $a->{start} <=> $b->{start} } @{$tran->{cds_regions}};    
  }
  
  my $exon_count = 1;
  foreach my $exon (@exons) {
    my $e_id = $exon->{stable_id};
    $e_id = $t_id . "-" . sprintf("%03d", $exon_count++) if not $e_id;
    printf("%s\t%s\texon\t%d\t%d\t%s\t%s\t.\ttranscript_id \"%s\"; gene_id \"%s\"; exon_id \"%s\"\n",
           $chr, 
           $type,
           $exon->{start},
           $exon->{end},
           $exon->{score},
           $strand,
           $t_id,
           $g_id,
           $e_id);
  }

  foreach my $cds (@cds) {
    printf("%s\t%s\tCDS\t%d\t%d\t%s\t%s\t%d\ttranscript_id \"%s\"; gene_id \"%s\"; translation_id \"%s\"\;\n",
           $chr, 
           $type,
           $cds->{start},
           $cds->{end},
           $cds->{score},
           $strand,
           $cds->{phase},
           $t_id,
           $g_id,
           $p_id);
  }
}
