#!/usr/local/ensembl/bin/perl -w

=head1 Synopsis

load_medaka.pl 

=head1 Description

Parses medaka genes out of the given GFF file and writes them to the database specified.

=head1 Config

All configuration is done through MedakaConf.pm

=cut

use strict;
use Carp;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::Analysis;
use Bio::SeqIO;
use MedakaConf;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(stack_trace_dump throw);
my $help;
my %opt; 

# options submitted with commandline override MedakaConf.pm 
GetOptions(
           \%opt ,
           '-h|help'    , 
           'dbhost=s' , 
           'dbuser=s' , 
           'dbpass=s' , 
           'dbport=i' , 
           'dbname=s' ,            
           'gene_type=s', 
           'logic_name=s' , 
           'gff_file=s' , 
           ) ; 

if ($opt{dbhost} && $opt{dbuser} && $opt{dbname} && $opt{dbpass} && $opt{dbport} ) {  
  $MED_DBHOST  = $opt{dbhost} ; 
  $MED_DBUSER = $opt{dbuser} ;  
  $MED_DBPASS = $opt{dbpass} ; 
  $MED_DBPORT = $opt{dbport} ; 
  $MED_DBNAME = $opt{dbname} ; 
}

$MED_GFF_FILE = $opt{gff_file} if $opt{gff_file} ; 
$MED_LOGIC_NAME =  $opt{logic_name} if $opt{logic_name} ; 
$MED_GENE_TYPE =  $opt{gene_type} if $opt{gene_type} ; 


unless ($MED_DBHOST && $MED_DBUSER && $MED_DBNAME && $MED_GFF_FILE && !$help){
  warn("Can't run without MedakaConf.pm values:
MED_DBHOST $MED_DBHOST 
MED_DBUSER $MED_DBUSER 
MED_DBNAME $MED_DBNAME
MED_DBPASS $MED_DBPASS
MED_GFF_FILE $MED_GFF_FILE
MED_LOGIC_NAME $MED_LOGIC_NAME
MED_GENE_TYPE $MED_GENE_TYPE
");
  $help = 1;
}

if ($help) {
    exec('perldoc', $0);
}
print STDERR "End of set-up. Time to parse using $MED_GFF_FILE.\n";
# end of set-up
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#PARSE GFF FILE FIRST

my (%gene)  = %{parse_gff()}; 
print STDERR "Done parsing. Now will print out.\n";

#######################################
# Parse GFF file. Use a hash of hashes.
# OuterKey = scaffold
# InnerKey = gene_stable_id

sub parse_gff {
  print "Parsing gff file\n";
  # example of GFF file.
  #	scaffold1       UTOLAPRE05100100001     initial-exon    129489  129606  .       +       0
  #	scaffold1       UTOLAPRE05100100001     internal-exon   129920  130027  .       +       2
  #	scaffold1       UTOLAPRE05100100001     internal-exon   130753  130839  .       +       2
  #	scaffold1       UTOLAPRE05100100001     final-exon      131859  132262  .       +       2
  #	scaffold6469    UTOLAPRE05100120178     single-exon     1604    2746    .       -       0

  # read in the file. 
  my %gff;
  open (GFF,$MED_GFF_FILE) or die "Cannot open gff file $MED_GFF_FILE\n";
  my $line;
  my $count;
  while (<GFF>){
    chomp;
    $line = $_;
    next if ($line =~ m/^\#/);
    my @fields = split/\t/, $line;
    if ($fields[2] eq 'initial-exon' || $fields[2] eq 'single-exon') {
      print STDERR "Resetting count to 0, ".$fields[2]."\n";
      $count = 0;      
    } else {
      $count++;
    }     
    #gene_id -> exon_number = (start, end, strand, frame, scaffold)
    @{$gff{$fields[1]}{$count}} = ($fields[3], $fields[4], $fields[6], $fields[7],$fields[0]);     
  }
  $line = '';
  close GFF; 
  return \%gff;
}



# check that it's all working so far
foreach my $g (sort keys %gene){
  foreach my $exon (sort {$a <=> $b} keys %{$gene{$g}}){
    print $g."\t".$exon."\t@{$gene{$g}{$exon}}\n";
  }
}


