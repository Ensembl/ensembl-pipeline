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

my $dbhost    = '';
my $dbuser    = '';
my $dbname    = '';
my $dbpass    = undef;
my $dnadbhost    = '';
my $dnadbuser    = '';
my $dnadbname    = '';
my $dnadbpass    = undef;
my $dbport    = 3306;
my $stable_id = 0;
my $db_id = 0;
my $file;

GetOptions(
	   'dbhost=s'    => \$dbhost,
	   'dbname=s'    => \$dbname,
	   'dbuser=s'    => \$dbuser,
	   'dbpass=s'    => \$dbpass,
	   'dbport=s'    => \$dbport,
	   'dnadbhost=s'    => \$dnadbhost,
	   'dnadbname=s'    => \$dnadbname,
	   'dnadbuser=s'    => \$dnadbuser,
	   'dnadbpass=s'    => \$dnadbpass,
	   'stable_id!' => \$stable_id,
	   'db_id!' => \$db_id,
	   'file=s' => \$file,
)
or die ("Couldn't get options");

if(!$dbhost || !$dbuser || !$dbname){
  die ("need to pass database settings in on the commandline -dbhost -dbuser -dbname -dbpass");
}
if(!$dnadbhost || !$dnadbname || !$dnadbuser){
  print STDERR "you have passed in no dnadb settings you db must have dna in it otherwise this won't work\n";
}
if(!$stable_id && !$db_id){
  die "need to specify to use either stable_id or dbId for the header line";
}elsif($stable_id && $db_id){
  print STDERR "you have defined both stable_id and db_id your identifier will have the format db_id.stable_id\n";
}

my $dnadb;
 
$dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $dnadbhost,
					    '-user'   => $dnadbuser,
					    '-dbname' => $dnadbname,
					    '-pass'   => $dnadbpass,
					    '-port'   => $dbport,
					   ) unless(!$dnadbhost || !$dnadbname || !$dnadbuser);


#print STDERR "have dnadb ".$dnadb."\n";
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $dbhost,
					    '-user'   => $dbuser,
					    '-dbname' => $dbname,
					    '-pass'   => $dbpass,
					    '-port'   => $dbport,
					    '-dnadb' => $dnadb,
					   );


print STDERR "connected to $dbname : $dbhost going to write to file ".
  " $file\n";

my $fh;
if($file){
  open (FH, '>'.$file) or die "couldn't open file ".$file." $!";
  $fh = \*FH;
} else{
  $fh = \*STDOUT;
}



my $seqio = Bio::SeqIO->new('-format' => 'Fasta' , -fh => $fh ) ;

my $gene_adaptor = $db->get_GeneAdaptor();
my $gene_ids = $gene_adaptor->list_dbIDs();

foreach my $gene (@{$gene_adaptor->fetch_all_by_dbID_list($gene_ids)}) {
  my $gene_id = $gene->dbID();

  foreach my $trans ( @{$gene->get_all_Transcripts}) {
    next if (!$trans->translation);

    my $identifier;
    if($db_id){
      $identifier = $trans->translation->dbID;
    }
    if($stable_id){
      if(!$db_id){
        $identifier = $trans->stable_id;
      } else {
        $identifier .= ".".$trans->stable_id;
      }
    }
    my $tseq = $trans->translate();
    if ( $tseq->seq =~ /\*/ ) {
      print STDERR "Translation of $identifier has stop codons ",
        "- Skipping! (in ",$trans->slice->name(),")\n";
      next;
    }

    $tseq->display_id($identifier);
    $tseq->desc("Translation id $identifier gene $gene_id");
    $seqio->write_seq($tseq);
  }
}
close($fh);
