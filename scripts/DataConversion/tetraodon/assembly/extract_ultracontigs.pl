#!/usr/local/ensembl/bin/perl


while(<>) {
  my @line = split(/\s+/, $_);

  if ($line[4] =~ /N/ and $line[5] !~ /GAP_IS/) {
    push @{$ucs{$line[0]}}, [];
  }
  elsif ($line[4] =~ /F/) {
    if (not $ucs{$line[0]}) {
      push @{$ucs{$line[0]}}, [];
    }
    push @{$ucs{$line[0]}->[-1]}, \@line;
  }
}

foreach my $chr (sort keys %ucs) {
  my $uc_count = 1;
  #print "UCs on $chr = ", scalar(@{$ucs{$chr}}), "\n";
  #$count += scalar(@{$ucs{$chr}});
  foreach my $uc (sort {$a->[0]->[1] <=> $b->[0]->[1]} @{$ucs{$chr}}) {
    if (@$uc > 1) {
      print "#ULTRACONTIG\n";
      my $offset = $uc->[0]->[1] - 1;
      my $count = 1;
      foreach my $scaffold (sort { $a->[1] <=> $b->[1]} @$uc) {
        $scaffold->[1] -= $offset;
        $scaffold->[2] -= $offset;
        $scaffold->[3] = $count++;
        $scaffold->[0] .= "_UC_$uc_count";
        print join("\t", @$scaffold), "\n";
        
      }
      $uc_count++;
    }
  }
}
