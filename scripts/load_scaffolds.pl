#!/usr/local/bin/perl -w

=head1 NAME

load_scaffolds.pl

=head1 SYNOPSIS
 
  load_scaffolds.pl -pipe blabla.fa

=head1 DESCRIPTION

This script loads a fasta file into the ensembl clone/contig/dna schema. BEWARE! It is for gunshot sequencing projects like fugu, where one clone = one contig is the rule (and they are usually referred to as scaffolds. It does not deal with multiple contigs per clone.

If the option -pipe is set it also fills the ensembl pipeline InputIdAnalysis with each contig id and value 1 (i.e. ready to be analysed)

=head1 OPTIONS


    -host    host name for database (gets put as host= in locator)
    -dbname    For RDBs, what name to connect to (dbname= in locator)
    -dbuser    For RDBs, what username to connect as (dbuser= in locator)
    -dbpass    For RDBs, what password to use (dbpass= in locator)
    -module    Module name to load (Defaults to Bio::EnsEMBL::DBSQL::DBAdaptor)
    -help      displays this documentation with PERLDOC
=cut


use strict;
use vars qw($USER $PASS $DB $HOST $DSN);
use Bio::EnsEMBL::DBLoader;
use Bio::SeqIO;
use Bio::EnsEMBL::PerlDB::Contig;
use Bio::EnsEMBL::PerlDB::Clone;
use Getopt::Long;

my $dbtype = 'rdb';
my $host   = 'localhost';
my $port   = '';
my $dbname = 'main_trunk_pipeline';
my $dbuser = 'root';
my $dbpass = undef;
my $module = 'Bio::EnsEMBL::DBSQL::DBAdaptor';

my $help = 0;
my $pipe = 0;
my $verbose = 0;

&GetOptions( 
	     'dbtype:s'   => \$dbtype,
	     'host:s'     => \$host,
	     'port:n'     => \$port,
	     'dbname:s'   => \$dbname,
	     'dbuser:s'   => \$dbuser,
	     'dbpass:s'   => \$dbpass,
	     'module:s'   => \$module,
             'pipe'       => \$pipe,
	     'verbose'    => $verbose,
	     'h|help'     => \$help
	     );
if ($help) {
    exec('perldoc', $0);
}

$SIG{INT} = sub {my $sig=shift;die "exited after SIG$sig";};

my $locator = "$module/host=$host;port=;dbname=$dbname;user=$dbuser;pass=$dbpass";
my $db =  Bio::EnsEMBL::DBLoader->new($locator);

my ($seqfile) = @ARGV;

if( !defined $seqfile ) { die 'cannot load because sequence file to load sequences from was not passed in as argument';}

my $seqio = new Bio::SeqIO(-format=>'Fasta', 
			   -file=>$seqfile);
my $count = 0;
while ( my $seq = $seqio->next_seq ) {
    my $cloneid= $seq->id;
    my $contigid = $cloneid.".1";
    $verbose && print STDERR ("Parsed contig $contigid : contig length ".$seq->length."\n");
    my $clone     = new Bio::EnsEMBL::PerlDB::Clone;
    my $contig    = new Bio::EnsEMBL::PerlDB::Contig;
    $clone->htg_phase(3);
    $clone->id($cloneid);
    $contig->id($contigid);
    $contig->internal_id($count++);
    $contig->seq($seq);    
    $clone->add_Contig($contig);
    eval { 
       $db->write_Clone($clone);
       $verbose && print STDERR "Written $clone scaffold into db\n";
    };
    if( $@ ) {
      print STDERR "Could not write clone into database, error was $@\n";
    }
    elsif ($pipe) {
      $db->prepare("insert into InputIdAnalysis (inputId,class,analysisId,created) values(".$contigid.",contig,1,now())");
    }
}
