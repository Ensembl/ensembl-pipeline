#!/usr/local/ensembl/bin/perl

my $repeat_source = "RepeatMasker";

while(<>) {
  /^\#/ and next;

  my @l = split /\t/;
  next if $l[1] ne $repeat_source;

  $l[0] =~ s/^chr//;
  $l[1] = "RepeatFeatures";
  $l[2] = "Repeat";

  if ($l[8] =~ /Note\s+\"([^\"]+)\"\s+\;\s+Note\s+\"([^\"]+)\"/) {
    my ($repeat_name, $repeat_class) = ($1, $2);
    $l[8] = "repeat_name \"$repeat_name\"; repeat_class \"$repeat_class\";";
  } elsif ($l[8] =~ /Note\s+\"([^\"]+)\"/) {
    my $repeat_name = $1;
    $l[8] = "repeat_name \"$repeat_name\"; repeat_class \"Tet_repeat\";";
  } 

  print join("\t", @l), "\n";
}
