#!/usr/local/ensembl/bin/perl

=head1 NAME

load_seq_region.pl


=head1 SYNOPSIS

  load_seq_region.pl <DB OPTIONS> \
    <COORDINATE SYSTEM OPTIONS> \
    { -fasta file <my.fa> | -agp_file <my.agp> }


=head1 DESCRIPTION

This script can do three things:

1) Use entries in a FASTA file to load a set of seq_regions into the
   seq_region table.

2) The sequence from the FASTA file can be optionally added to the dna
   table.

3) It can load seq_regions that represent the objects in an AGP file.

In all cases, appropriate (configurable) entries will be added to the
coord_system_table.


Here are example usages:

This would load the *sequences* in the given FASTA file into the
database under a coord system called contig:

./load_seq_region.pl \
  -dbhost host -dbuser user -dbname my_db -dbpass **** \
  -coord_system_name contig -rank 4 -sequence_level \
  -fasta_file sequence.fa


This would just load seq_regions to represent the entries in this
FASTA file:

./load_seq_region.pl \
  -dbhost host -dbuser user -dbname my_db -dbpass **** \
  -coord_system_name clone -rank 3 \
  -fasta_file clone.fa


This will load the assembled pieces from the AGP file into the
seq_region table.

./load_seq_region \
  -dbhost host -dbuser user -dbname my_db -dbpass **** \
  -coord_system_name chromosome -rank 1 \
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

    -coord_system_name
               The name of the coordinate system being stored.

    -coord_system_version
               The version of the coordinate system being stored.

    -default_version
               Flag to denote that this version is the default version
               of the coordinate system.

    -rank      The rank of the coordinate system. The highest
               coordinate system should have a rank of 1 (e.g. the
               chromosome coordinate system). The nth highest should
               have a rank of n. There can only be one coordinate
               system for a given rank.

    -sequence_level
               Flag to denete that this coordinate system is a
               'sequence level'. This means that sequence will be
               stored from the FASTA file in the dna table. This
               option isn't valid for an agp_file.

    -components
               When loading an AGP, the default behaviour is to load
               seq_regions from the OBJECT level. Setting this flag
               will load seq_regions from the COMPONENT level. Note,
               it will not populate the assembly table, use
               load_agp.pl instead


    OTHER OPTIONS:

    -agp_file   The name of the agp file to be parsed.

    -fasta_file The name of the fasta file to be parsed. Without the
                presence of the -sequence_level option the sequence
                will not be stored.

    -regex      A regex to parse the desired sequence ID from the
                FASTA description line.


    MISC OPTIONS:

    -help       Displays this documentation with PERLDOC.
    -verbose    Prints ... summin

=cut

use strict;
use warnings;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use Bio::SeqIO;

use Getopt::Long;

my ($dbhost, $dbport, $dbuser, $dbpass, $dbname);
my $regex;

my $cs_name;
my $cs_version;
my $default = 0;
my $rank;
my $sequence_level = 0;
my $components = 0;

my $agp;
my $fasta;

my $help;
my $verbose = 0;

&GetOptions(
            'dbhost|host:s' => \$dbhost,
            'dbport|port:n' => \$dbport,
            'dbuser|user:s' => \$dbuser,
            'dbpass|pass:s' => \$dbpass,
            'dbname|D:s'    => \$dbname,

            'coord_system_name:s'    => \$cs_name,
            'coord_system_version:s' => \$cs_version,
            'default_version!'       => \$default,
            'rank:i'                 => \$rank,
            'sequence_level!'        => \$sequence_level,
            'components!'            => \$components,

            'agp_file:s'   => \$agp,
            'fasta_file:s' => \$fasta,
            'regex:s'      => \$regex,

            'verbose!' => \$verbose,
            'help'     => \$help,
           ) or ($help = 1);

if( !$dbhost || !$dbuser || !$dbname || !$dbpass ){
  print STDERR "Can't store sequence without database details\n";
  $help = 1;
}

if( !$cs_name || (!$fasta  && !$agp) ){
  print STDERR "Need a coord_system_name and a FASTA or an AGP file!\n";
  $help = 1;
}

if($agp && $sequence_level){
  print STDERR ("Can't use an AGP file $agp to store a ".
                "sequence level coordinate system!\n");
  $help = 1;
}

if(!$rank) {
  print STDERR "A rank for the coordinate system must be specified ".
    "with the -rank argument\n";
    $help = 1;
}

if ($help) {
    exec('perldoc', $0);
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
my $sa = $db->get_SliceAdaptor();


## Get an existing coord system...
my $cs;
eval{ $cs = $csa->fetch_by_name($cs_name, $cs_version) };

## or create a new one
if(!$cs){
  $cs = Bio::EnsEMBL::CoordSystem->new
    (
     -NAME            => $cs_name,
     -VERSION         => $cs_version,
     -DEFAULT         => $default,
     -SEQUENCE_LEVEL  => $sequence_level,
     -RANK            => $rank
    );
  $csa->store($cs);
}


## Process the FASTA file...
if($fasta){
  my $count_ambiguous_bases =
    &process_fasta( $fasta, $cs, $sa, $sequence_level, $regex );
  if ($count_ambiguous_bases) {
    throw("All sequences has loaded, but $count_ambiguous_bases slices have ambiguous bases - see warnings. Please change all ambiguous bases (RYKMSWBDHV) to N!");
  }
}


## or process the AGP file
if($agp){
  &process_agp($agp, $cs, $sa, $components);
}

warn "Done\n";





## SUBS

sub process_fasta{
  my ($filename, $cs, $sa, $store_seq, $regex) = @_;
  my $have_ambiguous_bases = 0;

  my $seqio = Bio::SeqIO->
    new( -format=>'Fasta',
         -file=>$filename,
       );
  
  while ( my $seq = $seqio->next_seq ) {
    
    my $name = $seq->id;
    
    # NOTE, the code used to generate the sequence name depends on the
    # format of your fasta headers and what id you want to use. In
    # this case, we use the first word of the description line, as
    # parsed by SeqIO. You may want the id, or a custom regular
    # experssion. You will need to check what this produces. If you
    # have checked your ids and you know what you are getting you may
    # want to comment out the warning about this.
    
    if ($regex) {
      ($name) = $name =~ /$regex/;
    }
    
    warning("You are going to store with name '$name'. Are you sure ".
            "this is what you wanted")
      if $verbose;
    
    my $slice =
      &make_slice($name, 1, $seq->length, $seq->length, 1, $cs);
    
    if($store_seq){
      # Check that we don't have ambiguous bases in the DNA sequence,
      # we are only allowed to load ATGCN
      if ($seq->seq =~ /[^ACGTN]+/i) {
        $have_ambiguous_bases++;
        warning("Slice '$name' has at least one non-ATGCN (RYKMSWBDHV) base.".
                "Please change to N.");
      }
      $sa->store($slice, \$seq->seq);
    }
    else{
      $sa->store($slice);
    }
  }
  return $have_ambiguous_bases;
}



sub process_agp{
  my ($agp_file, $cs, $sa, $components) = @_;
  
  my %sequence_length;
  
  open(FH, $agp_file)
    or throw("Can't open AGP file '$agp_file' : $!");
  
  while(<FH>){
    next if /^\#/;
    chomp;
    
    ## See:
    ## http://www.ncbi.nlm.nih.gov/projects/genome/assembly/agp/AGP_Specification.shtml
    my @value = split/\t/;
    
    #GL000001.1      1       615     1       F       AP006221.1      36117   36731   -
    #GL000001.1      616     167417  2       F       AL627309.15     103     166904  +
    #GL000001.1      167418  217417  3       N       50000   clone   yes
    
    #cb25.fpc4250       119836  151061  13      W       c004100191.Contig2      1       31226   +
    #cb25.fpc4250       151062  152023  14      N       962     telomere        yes
    
    
    ## Skip gap rows
    next if $value[4] eq 'U';
    next if $value[4] eq 'N';
    
    ## Collect values for the object level
    my ($name, $start, $end) = @value[0,1,2];
    
    ## Or use the component level
    if($components){
        ($name, $start, $end) = @value[5,6,7];
    }
    
    ## Get the length of the sequence
    $sequence_length{$name} = 0
      unless exists $sequence_length{$name};
    $sequence_length{$name} = $end
      if $end > $sequence_length{$name};
    
  }
  
  foreach my $name (keys %sequence_length){
    my $length = $sequence_length{$name};
    my $slice = &make_slice( $name, 1, $length, $length, 1, $cs );
    $sa->store($slice);
  }
  
  return scalar keys %sequence_length;
}



# Used by both process_fasta and process_agp

sub make_slice{
  my ( $name, $start, $end, $length, $strand, $coordinate_system ) = @_;

  my $slice = Bio::EnsEMBL::Slice->new
      (
       -seq_region_name   => $name,
       -start             => $start,
       -end               => $end,
       -seq_region_length => $length,
       -strand            => $strand,
       -coord_system      => $coordinate_system,
      );
  return $slice;
}
