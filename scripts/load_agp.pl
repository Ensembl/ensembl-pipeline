#!/usr/local/ensembl/bin/perl

=head1 NAME

load_agp.pl


=head1 SYNOPSIS

  load_agp.pl <DB OPTIONS> \
    <OBJECT    COORDINATE SYSTEM OPTIONS> \
    <COMPONENT COORDINATE SYSTEM OPTIONS> \
    -agp_file <my.agp>


=head1 DESCRIPTION

Populates the assembly table from an AGP file

The AGP spec is described here:
http://www.ncbi.nlm.nih.gov/projects/genome/assembly/agp/AGP_Specification.shtml

Note, both the OBJECTS and COMPONENTS in the AGP should be loaded into
the seq_region table. This can be done with the load_seq_region.pl
script.

NB: What is called the 'object' in AGP terms is the assembled
sequence, composed of the individual gaps and components. We're not
talking about programming objects, but more like chromosomes or
scaffolds.


=head1 EXAMPLE USAGE

./load_agp.pl \
  -dbhost host -dbport 3601999 -dbuser user -dbname my_db -dbpass **** \
  -object_cs_name    chromosome -object_cs_version NCBI34 \
  -component_cs_name contig \
  -agp_file genome.agp


=head1 OPTIONS

    DB OPTIONS:
    -dbhost    Host name for the database
    -dbport    Port number for the database
    -dbuser    What username to connect as
    -dbpass    What password to use
    -dbname    What database to connect to

               For convenience, -host, -port, -user, -pass, and -D are
               also accepted as aliases for the above, respectively.


    COORDINATE SYSTEM OPTIONS:

    -object_cs_name
               The name of the object coordinate system.

    -object_cs_version [OPTIONAL]
               The version of the object coordinate system.

    -component_cs_name
               The name of the component coordinate system.

    -component_cs_version [OPTIONAL]
               The version of the component coordinate system.


    OTHER OPTIONS:

    -agp_file  The path to the the agp file.


    MISC OPTIONS:

    -help       Displays this documentation with PERLDOC.
    -verbose    Prints ... summin

=cut

use strict;
use warnings;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use Getopt::Long;

my ($dbhost, $dbport, $dbuser, $dbpass, $dbname);

my $assembled_cs_name;
my $assembled_cs_version;
my $component_cs_name;
my $component_cs_version;

my $agp_file;

my $help;
my $verbose = 0;

&GetOptions(
    'dbhost|host:s' => \$dbhost,
    'dbport|port:n' => \$dbport,
    'dbuser|user:s' => \$dbuser,
    'dbpass|pass:s' => \$dbpass,
    'dbname|D:s'    => \$dbname,

    'object_cs_name:s'       => \$assembled_cs_name,
    'object_cs_version:s'    => \$assembled_cs_version,
    'component_cs_name:s'    => \$component_cs_name,
    'component_cs_version:s' => \$component_cs_version,

    'agp_file:s' => \$agp_file,

    'verbose!' => \$verbose,
    'help'     => \$help,
) or ( $help = 1 );

if ( !$dbhost || !$dbuser || !$dbname || !$dbpass ) {
    print STDERR "Can't store sequence without database details\n";
    print STDERR "-dbhost $dbhost -dbuser $dbuser -dbname $dbname "
      . " -dbpass $dbpass\n";
    $help = 1;
}

if ( !$agp_file || !$assembled_cs_name || !$component_cs_name ) {
    print STDERR ( "Can't store assembly without an agp file or "
          . "coord system names for assembled and component pieces\n"
          . "-agp_file $agp_file -object_cs_name $assembled_cs_name "
          . "-component_cs_name $component_cs_name\n" );
    $help = 1;
}

if ($help) {
    exec( 'perldoc', $0 );
}



## Connect to the DB
my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -dbname => $dbname,
    -host   => $dbhost,
    -user   => $dbuser,
    -port   => $dbport,
    -pass   => $dbpass
);


## Get some adaptors
my $csa = $db->get_CoordSystemAdaptor();
my $sa  = $db->get_SliceAdaptor;



## Get the two existing coord systems
my $assembled_cs = $csa->
  fetch_by_name( $assembled_cs_name, $assembled_cs_version );
my $component_cs = $csa->
  fetch_by_name( $component_cs_name, $component_cs_version );


## Prepare the insert (isn't there an API way to do this?)
my $sql = <<EOS;
  INSERT IGNORE INTO assembly(
    asm_seq_region_id, asm_start, asm_end,
    cmp_seq_region_id, cmp_start, cmp_end, ori)
  VALUES(?, ?, ?, ?, ?, ?, ?)
EOS

my $sth = $db->dbc->prepare($sql);



## Process the AGP

open( FH, $agp_file )
  or throw("Can't open $agp_file");

my %component_reuse;
my $mapping_delimiter = '|';

while (<FH>) {
    next if /^\#/;
    chomp;

    my ($a_name, $a_start, $a_end, $ordinal, $type,
        $c_name, $c_start, $c_end, $ori ) = split;

    # Skip gap rows
    next if $type eq 'N';
    next if $type eq 'U';

    ## Convert AGP to Ensembl convention
    my $strand = $ori ne '-' ? +1 : -1;


    ## Look up the id for each object and component. Ditched the
    ## hasing code for clarity...

    my $a_slice = $sa->
      fetch_by_region( $assembled_cs->name, $a_name, undef, undef, undef,
                       $assembled_cs->version )
        or throw( "object    '$a_name' doesn't seem to exist in the database!\n" );
    my $c_slice = $sa->
      fetch_by_region( $component_cs->name, $c_name, undef, undef, undef,
                       $component_cs->version )
        or throw( "component '$c_name' doesn't seem to exist in the database!\n" );

    ## Sanity checks
    throw( "object    '$a_name' is longer than in the database! ($a_end vs. "
           . $a_slice->length. ")\n" )
      if $a_slice->length < $a_end;
    throw( "component '$c_name' is longer than in the database! ($c_end vs. "
           . $c_slice->length. ")\n" )
      if $c_slice->length < $c_end;


    ## I don't quite follow, but here you have it...
    if ( $component_reuse{$c_name}++ ) {
        $mapping_delimiter = '#';
        warning("You are already using component '$c_name' in another place "
                . "in your assembly. Are you sure you really want to?\n" );
    }

    ## INSERT
    $sth->execute( $sa->get_seq_region_id( $a_slice ), $a_start, $a_end,
                   $sa->get_seq_region_id( $c_slice ), $c_start, $c_end, $ori );
}

my $mapping_string = $assembled_cs->name;
$mapping_string .= ":" . $assembled_cs->version if ( $assembled_cs->version );
$mapping_string .= $mapping_delimiter . $component_cs->name;
$mapping_string .= ":" . $component_cs->version if ( $component_cs->version );

my $mc = $db->get_MetaContainer();
$mc->store_key_value( 'assembly.mapping', $mapping_string );

