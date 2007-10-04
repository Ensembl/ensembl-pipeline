#!/usr/local/ensembl/bin/perl -w

=head1 NAME

set_toplevel.pl

=head1 SYNOPSIS

  set_toplevel.pl -dbhost host -dbuser ensro -dbname homo_sapiens_core_20_34

=head1 DESCRIPTION

This script flags the set of non-redundant sequence regions in a database as
'toplevel'.  Doing this enables the use of 'toplevel' alias when retrieving
all slices from the database and allows the generic construction of input
ids for the pipeline regardless of what the actual assembly looks like.
This script should be run before the pipeline is started.

Before you can use this script the seq_region, assembly, and coord_system
tables need to be fully populated.  This is normally done via the
B<load_seq_region.pl> and B<load_agp.pl> scripts.

here is an example commandline

  ./set_toplevel.pl -dbhost host -dbuser user -dbname my_db -dbpass ****

=head1 OPTIONS

    -dbhost    host name for database (gets put as host= in locator)
    -dbname    what name to connect to (dbname= in locator)
    -dbuser    what username to connect as (dbuser= in locator)
    -dbpass    what password to use (dbpass= in locator)
    -help      displays this documentation with PERLDOC

=cut

use strict;
use warnings;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw);

use POSIX qw(ceil);

use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);


{ # block to avoid namespace pollution
  my $host   = '';
  my $port   = '';
  my $dbname = '';
  my $dbuser = '';
  my $dbpass = '';
  my $help;
  my @coord_system;

  &GetOptions(
              'dbhost:s'   => \$host,
              'dbport:n'   => \$port,
              'dbname:s'   => \$dbname,
              'dbuser:s'   => \$dbuser,
              'dbpass:s'   => \$dbpass,
              'ignore_coord_system:s@' => \@coord_system,
              'h|help'     => \$help,
             ) or ($help = 1);

  if(!$host || !$dbuser || !$dbname || !$dbpass){
    print STDERR "Can't store sequence without database details\n";
    print STDERR "-dbhost $host -dbuser $dbuser -dbname $dbname ".
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
     -host   => $host,
     -user   => $dbuser,
     -port   => $port,
     -pass   => $dbpass);

  my $csa = $db->get_CoordSystemAdaptor();
  my $slice_adaptor = $db->get_SliceAdaptor();

  my $attrib_type_id = add_attrib_code($db);

  print STDERR "Deleting old toplevel attributes\n";
  my $sth = $db->dbc->prepare('DELETE FROM seq_region_attrib ' .
                         "WHERE attrib_type_id = ?");
  $sth->execute($attrib_type_id);
  $sth->finish();

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
        my $string = $proj_slice->coord_system->name();
        $string .= " " . $proj_slice->seq_region_name();
        print STDERR "Adding $string to toplevel\n";

        $sth->execute($seq_region_id,$attrib_type_id, 1);
        $already_seen{$seq_region_id} = 1;
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
