#!/usr/local/ensembl/bin/perl -w


=head1 NAME

load_agp.pl

=head1 SYNOPSIS

  load_agp.pl

=head1 DESCRIPTION
 
This is a file for loading a standard agp file into the assembly table
The agp format is described here

http://www.sanger.ac.uk/Projects/C_elegans/DOCS/agp_files.shtml

Before you can use this script you need to have both the component
and assembled pieces loaded into the seq_region table. This can be
done with the load_seq_region.pl script

here is an example commandline

./load_agp.pl -dbhost host -dbuser user -dbname my_db -dbpass ****
-assembled_name chromosome -assembled_version NCBI34 
-component_name contig -agp_file genome.agp

=head1 OPTIONS

    -dbhost    host name for database (gets put as host= in locator)
    -dbname    For RDBs, what name to connect to (dbname= in locator)
    -dbuser    For RDBs, what username to connect as (dbuser= in locator)
    -dbpass    For RDBs, what password to use (dbpass= in locator)
    -assembled_name, the name of the coordinate system which represents
                   the assembled pieces
    -assembled_version, the version of the assembled coord system
    -component_name, the name of the coordinate system which represents
                     the component pieces
    -component_version, the version of the component coord system
    -agp_file path to the the agp file
    
    -help      displays this documentation with PERLDOC

=cut
use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

my $host   = '';
my $port   = '';
my $dbname = '';
my $dbuser = '';
my $dbpass = '';
my $assembled_name;
my $assembled_version;
my $component_name;
my $component_version;
my $agpfile;
my $help;

&GetOptions(
            'dbhost:s'   => \$host,
            'dbport:n'   => \$port,
            'dbname:s'   => \$dbname,
            'dbuser:s'   => \$dbuser,
            'dbpass:s'   => \$dbpass,
            'assembled_name:s' => \$assembled_name,
            'assembled_version:s' => \$assembled_version,
            'component_name:s' => \$component_name,
            'component_version:s' => \$component_version,
            'agp_file:s' => \$agpfile,
            'h|help'     => \$help,
            ) or ($help = 1);

if(!$host || !$dbuser || !$dbname || !$dbpass){
  print STDERR "Can't store sequence without database details\n";
  print STDERR "-dbhost $host -dbuser $dbuser -dbname $dbname ".
    " -dbpass $dbpass\n";
  $help = 1;
}

if(!$agpfile || !$assembled_name || !$component_name){
  print STDERR ("Can't store assembly without an agp file or ".
                "coord system names for assembled and component pieces ".
                "\n -agp_file $agpfile -assembled_name $assembled_name ".
                "-component_name $component_name\n");
  $help = 1;
}

if ($help) {
    exec('perldoc', $0);
}

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -dbname => $dbname,
    -host   => $host,
    -user   => $dbuser,
    -port   => $port,
    -pass   => $dbpass
);

my $csa = $db->get_CoordSystemAdaptor();
my $sa = $db->get_SliceAdaptor;

my $assembled_cs = $csa->fetch_by_name($assembled_name, 
                                       $assembled_version);
my $component_cs = $csa->fetch_by_name($component_name, 
                                       $component_version);

open(FH, $agpfile) or throw("Can't open $agpfile");

my %assembled_ids;
my %component_ids;
my $mapping_delimiter = '|';
LINE:while(<FH>){
  next if /^\#/;

  chomp;
  #cb25.fpc4250	119836	151061	13	W	c004100191.Contig2	1	31226	+
  #cb25.fpc4250	151062	152023	14	N	962	fragment	yes
  my ($a_name, $a_start, $a_end, $ordinal, $type, $c_name, $c_start, 
      $c_end, $ori) = split;
  if($type eq 'N'){
    next LINE; #skipping gaps as not stored in db
  }
  if(!$type){
    next LINE;
  }
  my $strand = 1;
  if($ori eq '-'){
    $strand = -1;
  }
  
  my ($a_id, $c_id);
  if(!$assembled_ids{$a_name}){
    if($a_name =~ /^chr(\S+)/){
      $a_name = $1;
    }
    my $a_piece = $sa->fetch_by_region($assembled_cs->name, $a_name,
                                       undef, undef, undef, 
                                       $assembled_cs->version);
    if(!$a_piece){
      throw($a_name." doesn't seem to exist in the database\n");
    }
    if($a_piece->length < $a_end){
      throw($a_name." apparent length ".$a_end. " is longer than ".
            "the length in the current database ".
            $a_piece->length."\n");
    }
    $a_id = $sa->get_seq_region_id($a_piece);
    $assembled_ids{$a_name} = $a_id;
  }else{
    $a_id = $assembled_ids{$a_name};
  }
  if($component_ids{$c_name}){
    warning("You are already using ".$c_name." in another place ".
         "in your assembly are you sure you want to\n");
    $mapping_delimiter = '#';
    $c_id = $component_ids{$c_name};
  }else{
    my $c_piece = $sa->fetch_by_region($component_cs->name, $c_name,
                                       undef, undef, undef, 
                                       $component_cs->version);
    if(!$c_piece){
      throw($c_name." doesn't seem to exist in the database\n");
    }
    if($c_piece->length < $c_end){
      throw($c_name." apparent length ".$c_end. " is longer than ".
            "the length in the current database ".
            $c_piece->length."\n");
    }
    $c_id = $sa->get_seq_region_id($c_piece);
    $component_ids{$c_name} = $c_id;
  }
  &insert_agp_line($a_id, $a_start, $a_end, $c_id, $c_start, $c_end, 
                   $strand, $db);
}

my $mapping_string = $assembled_cs->name;
$mapping_string .= ":".$assembled_cs->version if($assembled_cs->version);
$mapping_string .= $mapping_delimiter.$component_cs->name;
$mapping_string .= ":".$component_cs->version if($component_cs->version);

my $mc = $db->get_MetaContainer();
$mc->store_key_value('assembly.mapping', $mapping_string);

sub insert_agp_line{
  my ($chr_id, $chr_start, $chr_end, $contig, $contig_start, $contig_end, $contig_ori, $db) = @_;

  if(!$contig){
    #print STDERR "trying to insert into ".$chr_id." ".$chr_start." ".$chr_end."\n";
    die "contig id must be defined for this to work\n";
  }
  my $sql = "insert ignore into assembly(asm_seq_region_id, asm_start, asm_end, cmp_seq_region_id, cmp_start, cmp_end, ori) values(?, ?, ?, ?, ?, ?, ?)";
  
  my $sth = $db->dbc->prepare($sql);
  $sth->execute($chr_id, $chr_start, $chr_end, $contig, $contig_start, $contig_end, $contig_ori); 
}
