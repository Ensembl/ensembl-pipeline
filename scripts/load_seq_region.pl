#!/usr/bin/env perl


# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016] EMBL-European Bioinformatics Institute
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

load_seq_region.pl

=head1 SYNOPSIS

  load_seq_region.pl 

=head1 DESCRIPTION

This script can do three things. Use the entries in a fasta file to load
a set of seq_regions into the seq_region object. If desired the sequence
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
table.
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

use warnings ;
use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::SeqIO;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::CoordSystem;
use Getopt::Long qw(:config no_ignore_case);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
my $host   = '';
my $port   = '3306';
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
my $name_file;
my $replace_ambiguous_bases = 0 ;

GetOptions(
            'host|dbhost|h:s'   => \$host,
            'port|dbport|P:n'   => \$port,
            'dbname|db|D:s'   => \$dbname,
            'user|dbuser|u:s'   => \$dbuser,
            'pass|dbpass|p:s'   => \$dbpass,
            'coord_system_name|cs_name:s' => \$cs_name,
            'coord_system_version:s' => \$cs_version,
            'rank:i' => \$rank,
            'sequence_level!' => \$sequence_level,
            'default_version!' => \$default,
            'agp_file:s' => \$agp,
            'fasta_file:s' => \$fasta,
            'verbose!' => \$verbose,
            'regex:s' => \$regex,
            'name_file:s' => \$name_file,
            'help'     => \$help,
            'replace_ambiguous_bases'  => \$replace_ambiguous_bases
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

my %acc_to_name;

if ($name_file){
  open(NF, $name_file) or throw("Can't open ".$name_file." ".$!);
  while(<NF>){   
    chomp;
    my @name_values = split(/\s+/,$_);
    $acc_to_name{$name_values[1]}=$name_values[0];

  }
}


if($fasta)
{
  my $count_ambiguous_bases = &parse_fasta($fasta, $cs, $sa, $sequence_level,$regex,$replace_ambiguous_bases);
  if ($count_ambiguous_bases) 
  {
    if( $replace_ambiguous_bases )
    {
        warn( "All sequences has loaded, but $count_ambiguous_bases slices had ambiguous bases - see warnings. Changed all ambiguous bases (RYKMSWBDHV) to N." ) ;
    }
    else
    {
        throw( "All sequences has loaded, but $count_ambiguous_bases slices had ambiguous bases - see warnings. Please change all ambiguous bases (RYKMSWBDHV) to N." ) ;
    }
  }
}

if($agp){
  &parse_agp($agp, $cs, $sa,%acc_to_name);
}

sub parse_fasta{
  my ($filename, $cs, $sa, $store_seq,$regex,$replace_ambiguous_bases) = @_;
  my $have_ambiguous_bases = 0;

  my $seqio = new Bio::SeqIO(
                             -format=>'Fasta',
                             -file=>$filename
                            );
  
  while ( my $seq = $seqio->next_seq ) 
  {
    
    #NOTE, the code used to generate the name very much depends on the 
    #format of your fasta headers and what id you want to use
    #In this case we use the first word of the sequences description as
    #parsed by SeqIO but you may want the id or you may want to use a
    #regular experssion to get the sequence you will need to check what 
    #this will produce, if you have checked your ids and you know what
    #you are getting you may want to comment out the warning about this
    #print STDERR "id ".$seq->id." ".$seq->desc."\n";
    #my @values = split /\s+/, $seq->desc;
    #my @name_vals = split /\|/, $seq->id;
    my $name = $seq->id;

    #my $name = $name_vals[3]; 
       
    if ($regex) 
    {
      ($name) = $name =~ /$regex/;
    }
    print( "Loading ".$name."\n" ) if($verbose) ;
    if( $name !~ /^\w+\.\d/ || length($name)>40 )
    {
        warning( "Name ".$name." does not look like a valid accession - are you sure ".
            "this is what you want?" )
    }

    my $slice = &make_slice($name, 1, $seq->length, $seq->length, 1, $cs);
    if($store_seq)
    {
      #ambiguous bases: detect, warn and substitute (if required) 
      if( $seq->seq =~ /[^ACGTN]+/i && $replace_ambiguous_bases )
      {
        warning( "Slice ".$name." had at least one non-ATGCN (RYKMSWBDHV) base. Changed all to N." ) ;
        $have_ambiguous_bases++ ;
        my $seq_clean = $seq->seq ;
        $seq_clean =~ tr/RYKMSWBDHV/N/ ;
        $seq->seq($seq_clean) ;
        $sa->store($slice, \$seq->seq);
      }
      elsif( $seq->seq =~ /[^ACGTN]+/i  )
      {
        warning( "Slice ".$name." had at least one non-ATGCN (RYKMSWBDHV) base." ) ;
        $have_ambiguous_bases++ ;        
        $sa->store( $slice, \$seq->seq ) ;
      }
      else
      {
          $sa->store($slice, \$seq->seq);
      }
    }
    else
    {
      $sa->store($slice);
    }
  }
  return $have_ambiguous_bases;
}


sub parse_agp{
  my ($agp_file, $cs, $sa,%acc_to_name) = @_;
  my %end_value;
  open(FH, $agp_file) or throw("Can't open ".$agp_file." ".$!);
 LINE:while(<FH>){   
    chomp;
    next if /^\#/;

    #GL000001.1      1       615     1       F       AP006221.1      36117   36731   -
    #GL000001.1      616     167417  2       F       AL627309.15     103     166904  +
    #GL000001.1      167418  217417  3       N       50000   clone   yes

    #cb25.fpc4250	119836	151061	13	W	c004100191.Contig2	1	31226	+
    #cb25.fpc4250	151062	152023	14	N	962	telomere	yes
    my @values = split;
    #if($values[4] eq 'N'){
    #  next LINE; 
    #}
    my $initial_name = $values[0];
   
    # remove the 'chr' string if it exists
    if ($initial_name =~ /^chr(\S+)/) {
      $initial_name = $1;
    }


    my $name;

    if ($acc_to_name{$initial_name}){
      $name = $acc_to_name{$initial_name};
    }else{
      $name =$initial_name;
    }

    print "Name: ",$name,"\n";
    
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


