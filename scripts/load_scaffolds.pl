#!/usr/local/bin/perl -w

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
    -default_version, reflects this is the default version of this 
                      coordinate system, this is true by default but can
                      be switched off qwith nodefault_version
    -top_level  reflects this is the top level of sequence in the database
                this is off by default
    -sequence_level reflects this is a sequence level coordinate system
                    this is on by default
    -help      displays this documentation with PERLDOC

=cut

use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::SeqIO;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::CoordSystem;
use Getopt::Long;

my $host   = '';
my $port   = '';
my $dbname = '';
my $dbuser = '';
my $dbpass = '';
my $seqfile;
my $help;
my $cs_name;
my $cs_version;
my $default = 1;
my $top_level = 0;
my $sequence_level = 1;

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


my $seqio = new Bio::SeqIO(
                           -format=>'Fasta',
                           -file=>$seqfile
                          );


while ( my $seq = $seqio->next_seq ) {

  my @values = split /\s+/, $seq->desc;
  my $name = $values[0];
  my $slice = Bio::EnsEMBL::Slice->new(
                                       -seq_region_name   => $name,
                                       -start             => 1,
                                       -end               => $seq->length,
                                       -seq_region_length => $seq->length,
                                       -strand            => 1,
                                       -coord_system      => $cs,
                                      );

  $sa->store($slice, \$seq->seq);

}
