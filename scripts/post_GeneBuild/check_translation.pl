#!/usr/local/bin/perl

=head1 NAME

  dump_translations.pl

=head1 SYNOPSIS
 
  dump_translations.pl

=head1 DESCRIPTION

dump_translations.pl dumps out the translations of all the genes in a database specified in GeneConf.pm
It\'s a stripped down version of gene2flat.

=head1 OPTIONS

=cut

use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::SeqIO;
use Getopt::Long;

my $dbhost;
my $dbuser = 'ensro';
my $dbname;
my $dbpass = undef;

my $dnadbhost;
my $dnadbuser = 'ensro';
my $dnadbname;
my $dnadbpass = undef;

my $genetype;


my $tstable_id;
my $t_id;
my $gstable_id;
my $g_id;

&GetOptions(
	    'tstable_id:s' => \$tstable_id,
	    't_id:s' => \$t_id,
	    'g_id:s' => \$g_id,
	    'gstable_id:s' => \$gstable_id,
	    'dbhost:s'        => \$dbhost,
	    'dbname:s'        => \$dbname,
	    'dnadbhost:s'     => \$dnadbhost,
	    'dnadbname:s'     => \$dnadbname,
	    'genetype:s'      => \$genetype,
	    );

unless ( $dbhost && $dbname && $dnadbhost && $dnadbname && ( $t_id || $tstable_id || $g_id || $gstable_id) ){
  print STDERR "script to check the translation from the transcripts or genes in a database\n";
  
  print STDERR "Usage: $0 [-dbname -dbhost -dnadbname -dnadbhost -genetype] -t_id -tstable_id -g_id -gstable_id\n";
  exit(0);
}




my $dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					       '-host'   => $dnadbhost,
					       '-user'   => $dnadbuser,
					       '-dbname' => $dnadbname,
					       '-pass'   => $dnadbpass,
					      );


my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $dbhost,
					    '-user'   => $dbuser,
					    '-dbname' => $dbname,
					    '-pass'   => $dbpass,
					    '-dnadb'  => $dnadb,
					   );



print STDERR "connected to $dbname : $dbhost\n";

my $seqio = Bio::SeqIO->new('-format' => 'Fasta' , -fh => \*STDOUT ) ;

if ( $t_id){

  # first method
  my $tadaptor = $db->get_TranscriptAdaptor;
  my $trans    = $tadaptor->fetch_by_dbID($t_id);
  my $tseq     = $trans->translate();
  $tseq->desc("Transcript dbID: $t_id, transcript from TranscriptAdaptor");
  $seqio->write_seq($tseq);
  
  
}



