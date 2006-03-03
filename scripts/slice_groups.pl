#!/usr/local/bin/perl

# This script generates groups of toplevel slices from a
# database, ensuring that each group is no bigger than
# maxslicesize in total DNA count and maxgroupsize in
# total group-member count.

# The resulting groups of slice names are written to files
# in the given directory


use strict;
use Getopt::Long;

use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Slice qw(split_Slices);

my (
    $dbname,
    $dbhost,
    $dbuser,
    $dbport,
    $dbpass,
    $slice_size,
    $group_size,
    $outdir,
);

$dbuser = 'ensro';

&GetOptions(
            'dbname=s' => \$dbname,
            'dbuser=s' => \$dbuser,
            'dbhost=s' => \$dbhost,
            'dbport=s' => \$dbport,
            'dbpass=s' => \$dbpass,
            'maxslicesize=s' => \$slice_size,
            'maxgroupsize=s' => \$group_size,
            'outdir=s'    => \$outdir,
);

my $db = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new(
	'-dbname' => $dbname,
	'-host' => $dbhost,
	'-user' => $dbuser,
	'-port' => $dbport,
	'-pass' => $dbpass
);


die "you must supply a positive slice size\n"
    if not defined $slice_size or $slice_size <= 0;

die "You must supply a valid output directory\n"
    if not defined $outdir or not -d $outdir;

my $csa = $db->get_CoordSystemAdaptor();
my $sa = $db->get_SliceAdaptor();

my @slices = @{$sa->fetch_all('toplevel')};

@slices = sort { $b->length <=> $a->length } @{split_Slices(\@slices,
                                                            $slice_size,
                                                            0)};

my @slice_groups;
foreach my $sl (@slices) {
  if (not @slice_groups 
      or $slice_groups[-1]->{total_len} + $sl->length > $slice_size 
      or (defined $group_size and scalar(@{$slice_groups[-1]->{slices}}) == $group_size)) {
    push @slice_groups, {
      total_len => $sl->length,
      slices => [$sl],
    };
  } else {
    push @{$slice_groups[-1]->{slices}}, $sl;
    $slice_groups[-1]->{total_len} += $sl->length;
  }
}

@slice_groups = sort { scalar(@{$a->{slices}}) <=> scalar(@{$b->{slices}}) } @slice_groups;

my $groupcount = 1;
foreach my $grp (@slice_groups) {
  open FILE, ">$outdir/slice_names_$groupcount"
      or die "Could not open $outdir/slice_names_$groupcount for writing\n";

  foreach my $sl (@{$grp->{slices}}) {
    print FILE $sl->name, "\n";
  }

  $groupcount++;
}
