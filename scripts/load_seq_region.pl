#!/usr/local/ensembl/bin/perl -w

=head1 NAME

load_scaffolds.pl

=head1 SYNOPSIS

  load_scaffolds.pl 

=head1 DESCRIPTION

this script will load a fasta file into a new schema database. It parses
the file using Bio::SeqIO and then turns each sequence into a Slice
before storing it

as standard this script takes what ever SeqIO parses as the id for the
slice name, you may need to edit this so it works properly and uses the 
names you actually want

here is an example commandline

./load_scaffolds.pl -dbhost host -dbuser user -dbname my_db -dbpass ****
-seqfile sequence.fa -coord_system_name contig 

=head1 OPTIONS

    -dbhost    host name for database (gets put as host= in locator)
    -dbname    For RDBs, what name to connect to (dbname= in locator)
    -dbuser    For RDBs, what username to connect as (dbuser= in locator)
    -dbpass    For RDBs, what password to use (dbpass= in locator)

    -seqfile   the fasta file which contains the sequence
    -coord_system_name the name of the coordinate system being stored
    -coord_system_version the version of the coordinate system being stored
    -default_version shows this version is the default version of the 
                     coordinate system
    -top_level  reflects this is the top level of sequence in the database
    -sequence_level reflects this is a sequence level coordinate system
    -agp_file shows the file you are passing in is an agp file and will 
                     be treated appropriately
    -fasta_file, shows the file you are passing in is a fasta file
                 the sequence will only be stored if the -sequence_level
                 option is also set the sequence will not be stored
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
my $seqfile;
my $help;
my $cs_name;
my $cs_version;
my $default = 0;
my $top_level = 0;
my $sequence_level = 0;
my $agp = 0;
my $fasta = 0;

&GetOptions(
            'dbhost:s'   => \$host,
            'dbport:n'   => \$port,
            'dbname:s'   => \$dbname,
            'dbuser:s'   => \$dbuser,
            'dbpass:s'   => \$dbpass,
            'seqfile:s'  => \$seqfile,
            'coord_system_name:s' => \$cs_name,
            'coord_system_version:s' => \$cs_version,
            'top_level!' => \$top_level,
            'sequence_level!' => \$sequence_level,
            'default_version!' => \$default,
            'agp_file!' => \$agp,
            'fasta_file!' => \$fasta,
            'h|help'     => \$help,
           ) or ($help = 1);

if(!$host || !$dbuser || !$dbname || !$dbpass){
  print STDERR "Can't store sequence without database details\n";
  print STDERR "-dbhost $host -dbuser $dbuser -dbname $dbname ".
    " -dbpass $dbpass\n";
  $help = 1;
}
if(!$cs_name || !$seqfile){
  print STDERR "Need coord_system_name and seqfile to beable to run\n";
  print STDERR "-coord_system_name $cs_name -seqfile $seqfile\n";
  $help = 1;
}
if($agp && $sequence_level){
  print STDERR ("Can't use an agp file ".$seq_file." to store a ".
                "sequence level coordinate system ".$cs_name."\n");
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
  $cs = $csa->fetch_by_name($cs_name);
};
if(!$cs){
  $cs = Bio::EnsEMBL::CoordSystem->new
    (
     -NAME            => $cs_name,
     -VERSION         => $cs_version,
     -DEFAULT         => $default,
     -SEQUENCE_LEVEL  => $sequence_level,
     -TOP_LEVEL       => $top_level,
    );
$csa->store($cs);
}

my $sa  = $db->get_SliceAdaptor();


if($fasta){
  &parse_fasta($seqfile, $cs, $sa, $sequence_level);
}

if($agp){
  &parse_agp($seqfile, $cs, $sa);
}

sub parse_fasta{
  my ($filename, $cs, $sa, $store_seq) = @_;

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
    
    my @values = split /\s+/, $seq->desc;
    my $name = $seq->id;
    warning("You are going to store with name ".$name." are you sure ".
            "this is what you wanted");
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
    #cb25.fpc4250	119836	151061	13	W	c004100191.Contig2	1	31226	+
    #cb25.fpc4250	151062	152023	14	N	962	fragment	yes
    my @values = split;
    if($values[4] eq 'N'){
      next LINE; 
    }
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


