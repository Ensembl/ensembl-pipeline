#!/usr/bin/env perl


# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2017] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


# This script generates groups of toplevel slices from a
# database, ensuring that each group is no bigger than
# maxslicesize in total DNA count and maxgroupsize in
# total group-member count.

# The resulting groups of slice names are written to files
# in the given directory


use warnings ;
use strict;
use Getopt::Long qw(:config no_ignore_case);

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
$dbport = '3306';
$dbuser = 'ensro';

GetOptions(
            'dbname|db|D:s' => \$dbname,
            'user|dbuser|u:s' => \$dbuser,
            'host|dbhost|h:s' => \$dbhost,
            'port|dbport|P:s' => \$dbport,
            'pass|dbpass|p:s' => \$dbpass,
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
