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
    -help      displays this documentation with PERLDOC

=cut

use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::SeqIO;
use Bio::EnsEMBL::RawContig;
use Bio::EnsEMBL::Clone;
use Getopt::Long;

my $dbtype = 'rdb';
my $host   = '';
my $port   = '';
my $dbname = '';
my $dbuser = '';
my $dbpass = '';
my $write  = 0;

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
             'pipe'       => \$pipe,
	     'verbose'    => \$verbose,
	     'write'      => \$write,
	     'h|help'     => \$help
	     );
if ($help) {
    exec('perldoc', $0);
}

$SIG{INT} = sub {my $sig=shift;die "exited after SIG$sig";};

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -dbname => $dbname,
    -host   => $host,
    -user   => $dbuser,
    -port   => $port,
    -pass   => $dbpass
);

my ($seqfile) = @ARGV;

if( !defined $seqfile ) { die 'cannot load because sequence file to load sequences from was not passed in as argument';}

my $seqio = new Bio::SeqIO(-format=>'Fasta',
			   -file=>$seqfile);
my $count = 0;
while ( my $seq = $seqio->next_seq ) {
    my $cloneid= $seq->id;
    $verbose && print STDERR ("Parsed contig $contigid : contig length ".$seq->length."\n");
    my $clone     = new Bio::EnsEMBL::Clone;
    my $contig    = new Bio::EnsEMBL::RawContig;
    $clone->htg_phase(-1);
    $clone->id($cloneid);
    $clone->version(1);
    $clone->embl_id("");
    $clone->embl_version(0);
    my $now = time;
    $clone->created($now);
    $clone->modified($now);
    $contig->seq($seq->seq);
    my $bp = $seq->length;   #
    $contig->length($bp);  
    my $contigid = $cloneid.".1.".$bp;
    $contig->name($contigid);
    $contig->embl_offset(1);
    $clone->add_Contig($contig);
    if ($write) {
       eval {
          $db->get_CloneAdaptor->store($clone);
          $verbose && print STDERR "Written ".$clone->id." scaffold into db\n";
       };
       if( $@ ) {
         print STDERR "Could not write clone into database, error was $@\n";
       }
       if ($pipe) {
          my $sth = $db->prepare("insert into input_id_analysis (input_id,class,analysis_id,created) values('".$contigid."','contig',1,now())");
         $sth->execute;
         $verbose && print STDERR "Written input_id_analysis entry for ".$clone->id."\n";
       }
    }
}
