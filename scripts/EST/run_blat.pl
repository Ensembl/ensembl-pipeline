#!/usr/local/ensembl/bin/perl -w

=head1 NAME

run_blat

=head1 SYNOPSIS
 

=head1 DESCRIPTION


=head1 OPTIONS

    -input_id  The input id for the RunnableDB
    -runnable  The name of the runnable module we want to run
    -analysis  The number of the analysisprocess we want to run

=cut

use strict;
use Getopt::Long;
use Bio::SeqIO;

use Bio::EnsEMBL::Pipeline::Config::cDNAs_ESTs::Blat qw (
							 EST_DBNAME
							 EST_DBUSER
							 EST_DBPASS
							 EST_DBHOST
							 EST_BLAT_GENOMIC
							 EST_BLAT_OPTIONS
							 EST_REFDBNAME  
							 EST_REFDBHOST
							 EST_REFDBUSER
							 EST_REFDBPASS
							);


use Bio::EnsEMBL::DBSQL::DBAdaptor;

# this is the db where the resulting genes will be written
my $dbname = $EST_DBNAME;
my $dbuser = $EST_DBUSER;
my $dbpass = $EST_DBPASS;
my $host   = $EST_DBHOST;


my $runnable;
my $input_id;
my $write  = 0;
my $check  = 0;
my $params;
my $pepfile;
my $analysis;
my $query_seq;

# can override db options on command line
&GetOptions( 
	    'runnable:s'    => \$runnable,
	    'analysis:s'    => \$analysis,
	    'write'         => \$write,
	    'check'         => \$check,
	    'query_seq:s'   => \$query_seq,
	   );

$| = 1;

die "No runnable entered" unless defined ($runnable);
(my $file = $runnable) =~ s/::/\//g;
require "$file.pm";

if ($check) {
   exit(0);
}

print STDERR "args: $host : $dbuser : $dbpass : $dbname\n";

my $dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					       -host             => $EST_REFDBHOST,
					       -user             => $EST_REFDBUSER,
					       -dbname           => $EST_REFDBNAME,
					       );

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    -host             => $host,
					    -user             => $dbuser,
					    -dbname           => $dbname,
					    -pass             => $dbpass,
					    -dnadb            =>  $dnadb,
					   );



my $analysis_obj = $db->get_AnalysisAdaptor->fetch_by_logic_name($analysis);

die "No input entered" unless( defined ($query_seq));
			       

# convert the multiple fasta file with query est/rna sequences into an array of seq features
my @sequences;
my $seqio= Bio::SeqIO->new(
			   -format => "Fasta",
			   -file   => "$query_seq",
			  );
while( my $seq = $seqio->next_seq() ){
    if (defined($seq) && !($seq->seq eq '') && !($seq->display_id eq '') ){
	push( @sequences, $seq );
    }
    else{
	print STDERR "problems getting sequence ".$seq->display_id."\n".$seq->seq."\n";
    }
}

print STDERR "got ".scalar(@sequences)." sequence objects\n";


my $runobj = "$runnable"->new(-db         => $db,
			      -input_id   => \@sequences,
			      -rna_seqs   => \@sequences,
			      -analysis   => $analysis_obj,
			      -database   => $EST_BLAT_GENOMIC,
			      -query_type => 'dna',
			      -target_type=> 'dna',
			      -options    => "$EST_BLAT_OPTIONS",
			     );

$runobj->fetch_input;
$runobj->run;

if ($write) {
  print STDERR "run_blat.pl : wrtting!\n";
  $runobj->write_output;
}
