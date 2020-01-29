#!/usr/bin/env perl


# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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


=head1 NAME

set_toplevel.pl

=head1 SYNOPSIS

  set_toplevel.pl -dbhost host -dbport port -dbuser ensro -dbname homo_sapiens_core_20_34

=head1 DESCRIPTION

  This script flags the set of non-redundant sequence regions in a
  database as 'toplevel'. Doing this enables the use of 'toplevel'
  alias when retrieving all slices from the database and allows the
  generic construction of input ids for the pipeline regardless of
  what the actual assembly looks like. This script should be run
  before the pipeline is started.

  Before you can use this script the seq_region, assembly, and
  coord_system tables need to be fully populated. This is normally
  done via the B<load_seq_region.pl> and B<load_agp.pl> scripts.

=head1 OPTIONS

  -dbhost    host name for database (gets put as host= in locator)
  -dbname    what name to connect to (dbname= in locator)
  -dbuser    what username to connect as (dbuser= in locator)
  -dbpass    what password to use (dbpass= in locator)
  -help      displays this documentation with PERLDOC

=head1 EXAMPLES

  perl set_toplevel.pl -dbhost host -dbport port -dbuser user -dbname my_db -dbpass ****

=cut

use strict;
use warnings;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw);

use POSIX qw(ceil);

use Getopt::Long qw(:config no_ignore_case);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);


{ # block to avoid namespace pollution what are you talking about?
  my $dbhost = '';
  my $dbport = '3306';
  my $dbname = '';
  my $dbuser = '';
  my $dbpass = '';

  my $help;

  my @coord_system;

  GetOptions( 'host|dbhost|h:s'          => \$dbhost,
              'port|dbport|P:n'          => \$dbport,
              'user|dbuser|u:s'          => \$dbuser,
              'pass|dbpass|p:s'          => \$dbpass,
              'dbname|db|D:s'             => \$dbname,
              'ignore_coord_system:s@' => \@coord_system,
              'help'                 => \$help,
  ) or ( $help = 1 );

  if(!$dbhost || !$dbuser || !$dbname || !$dbpass){
    print STDERR "Can't store sequence without database details\n";
    print STDERR "-dbhost $dbhost -dbuser $dbuser -dbname $dbname ".
      " -dbpass $dbpass\n";
    $help = 1;
  }

  if ($help) {
    exec('perldoc', $0);
  }
  my %coord_systems_to_ignore;
  foreach my $cs(@coord_system){
    $coord_systems_to_ignore{$cs} = 1;
  }

  my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new
    (-dbname => $dbname,
     -host   => $dbhost,
     -user   => $dbuser,
     -port   => $dbport,
     -pass   => $dbpass);

  my $csa = $db->get_CoordSystemAdaptor();
  my $slice_adaptor = $db->get_SliceAdaptor();

  my $attrib_type_id = add_attrib_code($db);

  print STDERR "Deleting old toplevel attributes\n";
  my $sth = $db->dbc->prepare('DELETE FROM seq_region_attrib ' .
                         "WHERE attrib_type_id = ?");
  $sth->execute($attrib_type_id);
  $sth->finish();

  # get all of the regions that are components of other seq regions (assembly.cmp_seq_region_id)
  print STDERR "getting cmp_seq_region_id(s) from assembly table (and cs versions)\n";
  my $all_cmp_ids = all_component_ids($db);

  # prepare insert 'top_level' attributes
  $sth = $db->dbc->prepare('INSERT INTO seq_region_attrib ' .
                      'SET seq_region_id = ?,' .
                      '    attrib_type_id = ?,' .
                      '    value = ?');

  #get all of the seqlevel sequence regions and project them to 'toplevel'
  my $seqlevel_cs = $csa->fetch_by_name('seqlevel');

  if(!$seqlevel_cs) {
    throw("No 'seqlevel' CoordSystem has been defined." .
          "One must be specified using the load_seq_region.pl script.");
  }

  my $slices = $slice_adaptor->fetch_all($seqlevel_cs->name(),
                                         $seqlevel_cs->version);

  my %already_seen;

  my $total_number = @$slices;
  my $total_processed = 0;
  my $cs_added = {};
  my $five_percent = int($total_number/20);

  print STDERR "Adding new toplevel attributes to ".$total_number." slices.\n";

  foreach my $slice (@$slices) {
    my $projection = $slice->project('toplevel');
    foreach my $segment (@$projection) {
      my $proj_slice = $segment->[2];
      my $seq_region_id = $proj_slice->get_seq_region_id();

      if(!$already_seen{$seq_region_id}) {
        next if(keys(%coord_systems_to_ignore) &&
                $coord_systems_to_ignore{$proj_slice->coord_system->name});
        next if(exists $all_cmp_ids->{$proj_slice->get_seq_region_id});

        my $string = $proj_slice->coord_system->name();
        $string .= " " . $proj_slice->seq_region_name();
        print STDERR "Adding $string to toplevel\n";

        $sth->execute($seq_region_id,$attrib_type_id, 1);
        $already_seen{$seq_region_id} = 1;

        $cs_added->{$proj_slice->coord_system->name}++;
        $cs_added->{TOTAL}++;
      }
    }

    # print out progress periodically if there are a decent number of slices
    if($five_percent > 5){
      $total_processed++;
      if($total_processed % $five_percent == 0) {
	print STDERR ceil($total_processed/$total_number*100)."% complete\n";
      }
    }
  }

  if (%$cs_added) {
    print STDERR "toplevel regions ", join("\n\t", map {"$_:\t". $cs_added->{$_}} sort keys %$cs_added), "\n";
  }

  $sth->finish();
}

sub add_attrib_code {
  my $db = shift;
  # add a toplevel code to the attrib_type table if it is not there already

  my $sth = $db->dbc->prepare("SELECT attrib_type_id " .
                         "FROM attrib_type " .
                         "WHERE code = 'toplevel'");

  $sth->execute();

  if($sth->rows()) {
    my ($attrib_type_id) = $sth->fetchrow_array();
    $sth->finish();
    return $attrib_type_id;
  }
  $sth->finish();


  $sth = $db->dbc->prepare("INSERT INTO attrib_type " .
                    "SET code = 'toplevel', " .
                    "name = 'Top Level', " .
                    "description = 'Top Level Non-Redundant Sequence Region'");

  $sth->execute();
  my $attrib_type_id = $sth->{'mysql_insertid'};
  $sth->finish();

  return $attrib_type_id;
}

sub all_component_ids {
  # get all of the regions that are components of other seq regions (assembly.cmp_seq_region_id)
  my $db = shift;
  my $sth = $db->dbc->prepare('SELECT distinct cmp_seq_region_id FROM assembly');

  $sth->execute();

  my $result = {};
  if($sth->rows()) {
    while( my ($seq_region_id) = $sth->fetchrow_array) {
      $result->{$seq_region_id} = 0;
    }
  }

  $sth->finish();

  return $result;
}
