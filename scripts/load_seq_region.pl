#!/usr/local/ensembl/bin/perl -w

=head1 NAME

load_seq_region.pl

=head1 SYNOPSIS

  load_seq_region.pl 

=head1 DESCRIPTION

This script can do three things. Use the entries in a fasta file to load
a set of seq_regions into the seq_region object. If desried the sequence
can also be stored. It can also load seq_regions which represent the
assembled pieces from an agp file. The assembled pieces are represented by 
the first 3 columns of the agp for further definition of the format see here

http://www.sanger.ac.uk/Projects/C_elegans/DOCS/agp_files.shtml

here are example commandlines

this would load the sequence in the given file into the database under
a coord system called contig. 

./load_seq_region.pl -dbhost host -dbuser user -dbname my_db -dbpass ****
-coord_system_name contig -rank 4 -sequence_level -fasta_file sequence.fa

this would just load seq_regions to represent the entries in this file

./load_seq_region.pl -dbhost host -dbuser user -dbname my_db -dbpass ****
-coord_system_name clone -rank 3 -fasta_file clone.fa

this will load the assembled pieces from the agp file into the seq_region
table. T
./load_seq_region -dbhost host -dbuser user -dbname my_db -dbpass ****
-coord_system_name chromosome -rank 1 -agp_file genome.agp





=head1 OPTIONS

    -dbhost    host name for database (gets put as host= in locator)
    -dbname    For RDBs, what name to connect to (dbname= in locator)
    -dbuser    For RDBs, what username to connect as (dbuser= in locator)
    -dbpass    For RDBs, what password to use (dbpass= in locator)

    -coord_system_name the name of the coordinate system being stored
    -coord_system_version the version of the coordinate system being stored
    -rank the rank of the coordinate system.  The highest coordinate system
          should have a rank of 1 (e.g. the chromosome coord system).  The nth
          highest should have a rank of n.  There can only be one coordinate 
          system for a given rank.
    -default_version shows this version is the default version of the 
                     coordinate system
    -sequence_level reflects this is a sequence level coordinate system and
    means sequence will be stored from the fasta file. This option isn't valid
    if an agp_file is being passed in'
    -agp_file the name of the agp file to be parsed
    -fasta_file, name of the fasta file to be parsed without the presence of
             the -sequence_level option the sequence will not be stored
    -verbose, prints the name which is going to be used can be switched 
              off with -noverbose
    -help      displays this documentation with PERLDOC

=cut

use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::SeqIO;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::CoordSystem;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
my $host   = '';
my $port   = '';
my $dbname = '';
my $dbuser = '';
my $dbpass = '';
my $help;
my $cs_name;
my $cs_version;
my $default = 0;
my $sequence_level = 0;
my $agp;
my $fasta;
my $rank;
my $verbose = 0;
my $regex;

&GetOptions(
            'dbhost:s'   => \$host,
            'dbport:n'   => \$port,
            'dbname:s'   => \$dbname,
            'dbuser:s'   => \$dbuser,
            'dbpass:s'   => \$dbpass,
            'coord_system_name:s' => \$cs_name,
            'coord_system_version:s' => \$cs_version,
            'rank:i' => \$rank,
            'sequence_level!' => \$sequence_level,
            'default_version!' => \$default,
            'agp_file:s' => \$agp,
            'fasta_file:s' => \$fasta,
            'verbose!' => \$verbose,
            'regex:s' => \$regex,
            'h|help'     => \$help,
           ) or ($help = 1);

if(!$host || !$dbuser || !$dbname || !$dbpass){
  print STDERR "Can't store sequence without database details\n";
  print STDERR "-dbhost $host -dbuser $dbuser -dbname $dbname ".
    " -dbpass $dbpass\n";
  $help = 1;
}
if(!$cs_name || (!$fasta  && !$agp)){
  print STDERR "Need coord_system_name and fasta/agp file to beable to run\n";
  print STDERR "-coord_system_name $cs_name -fasta_file $fasta -agp_file $agp\n";
  $help = 1;
}
if($agp && $sequence_level){
  print STDERR ("Can't use an agp file $agp to store a ".
                "sequence level coordinate system ".$cs_name."\n");
  $help = 1;
}
if(!$rank) {
  print STDERR "A rank for the coordinate system must be specified " .
    "with the -rank argument\n";
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

my $cs;
eval{
  $cs = $csa->fetch_by_name($cs_name, $cs_version);
};
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

my $sa  = $db->get_SliceAdaptor();


if($fasta){
  &parse_fasta($fasta, $cs, $sa, $sequence_level,$regex);
}

if($agp){
  &parse_agp($agp, $cs, $sa);
}

sub parse_fasta{
  my ($filename, $cs, $sa, $store_seq,$regex) = @_;

  my $seqio = new Bio::SeqIO(
                             -format=>'Fasta',
                             -file=>$filename
                            );
  
  while ( my $seq = $seqio->next_seq ) {
    
    #NOTE, the code used to generate the name very much depends on the 
    #format of yuor fasta headers and what id you want to use
    #In this case we use the first word of the sequences description as
    #parseed by SeqIO but you may want the id or you may want to use a
    #regular experssion to get the sequence you will need to check what 
    #this will produce, if you have checked your ids and you know what
    #you are getting you may want to comment out the warning about this
    #print STDERR "id ".$seq->id." ".$seq->desc."\n";
    #my @values = split /\s+/, $seq->desc;
    my $name = $seq->id;

    if ($regex) {
      ($name) = $name =~ /$regex/;
    }
    warning("You are going to store with name ".$name." are you sure ".
            "this is what you wanted") if($verbose);
    my $slice = &make_slice($name, 1, $seq->length, $seq->length, 1, $cs);
    if($store_seq){
      $sa->store($slice, \$seq->seq);
    }else{
      $sa->store($slice);
    }
  }
}


sub parse_agp{
  my ($agp_file, $cs, $sa) = @_;

  my %end_value;
  open(FH, $agp_file) or throw("Can't open ".$agp_file." ".$!);
 LINE:while(<FH>){   
    chomp;
    next if /^\#/;
    #cb25.fpc4250	119836	151061	13	W	c004100191.Contig2	1	31226	+
    #cb25.fpc4250	151062	152023	14	N	962	telomere	yes
    my @values = split;
    #if($values[4] eq 'N'){
    #  next LINE; 
    #}
    my $name = $values[0];
    my $end = $values[2];
    if(!$end_value{$name}){
      $end_value{$name} = $end;
    }else{
      if($end > $end_value{$name}){
        $end_value{$name} = $end;
      }
    }
  }
  foreach my $name(keys(%end_value)){
    my $end = $end_value{$name};
    my $slice = &make_slice($name, 1, $end, $end, 1, $cs);
    $sa->store($slice);
  }
  
  close(FH) or throw("Can't close ".$agp_file);
}

sub make_slice{
  my ($name, $start, $end, $length, $strand, $coordinate_system) = @_;

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


