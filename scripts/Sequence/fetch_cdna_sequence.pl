#!/usr/local/bin/perl -w

# dump_seq_into_fastA.pl
# it reads a bit of sequence and dump it into a fasA file, to eb able to view it in Apollo

use strict;
use diagnostics;

use Bio::SeqIO;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;



# get a contig with a piece-of/entire  chromosome

my $global_start;
my $global_end; 
my $chr_name;   
my $t_id;
my $dbhost;
my $dbname;
my $dnadbname;
my $dnadbhost;
my $reverse;


&GetOptions(
	    't_id:s' => \$t_id,
	    'dbname:s'   => \$dbname,
	    'dbhost:s'   => \$dbhost,
	    'dnadbname:s'=> \$dnadbname,
	    'dnadbhost:s'=> \$dnadbhost,
	    'reverse'    => \$reverse,
	   );

unless ( $t_id  && $dbname && $dbhost && $dnadbname && $dnadbhost ){
  print STDERR "Usage: $0 -t_id (transcript dbID) -dbname -dbhost -dnadbname -dnadbhost\n";
  exit(0);
}


# connect to the database
my $dnadb= new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $dnadbhost,
					      -user  => 'ensro',
					      -dbname=> $dnadbname);

my $db= new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $dbhost,
					   -user  => 'ensro',
					   -dbname=> $dbname,
					   -dnadb => $dnadb,
					  );


my $tadaptor = $db->get_TranscriptAdaptor;

my $transcript = $tadaptor->fetch_by_dbID( $t_id );

my $seq = $transcript->seq;

my $outfile;
$outfile = "$t_id.fa";

open OUT, ">$outfile";
# get a Bio::SeqIO object
my $out = Bio::SeqIO->new(-format=>'Fasta',
			  -fh =>  \*OUT,
			 );


$seq->display_id($t_id);
$out->write_seq($seq);	       
close OUT;




