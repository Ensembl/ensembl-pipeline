#!/usr/local/ensembl/bin/perl -w

use strict;
use IO::File;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;

my $fh = new IO::File;

# connect to the database
my $species = 'mouse';
#my $species = 'rat';
my $dbhost = 'ecs2f';
my $dbuser = 'ensro';
my $dbname = "genewisedb_".$species;
my $outdir;
my $threshold = 50;


&GetOptions(
	    'dbname:s'       => \$dbname,
	    'dbhost:s'       => \$dbhost,
	    'outdir:s'       => \$outdir,
	    );

my $db= new Bio::EnsEMBL::DBSQL::DBAdaptor(
					   -host  => $dbhost,
					   -user  => $dbuser,
					   -dbname=> $dbname
					  );


my $input_ids = &get_input_ids($db,$threshold);
my %ids = %$input_ids;

foreach my $hmm ( keys %ids ){
  foreach my $dna_id ( keys %{$ids{$hmm}} ){
    
    my $command = "bsub -q acari -C0 -f  \"/ecs2/work1/eae/GeneWiseHMM/HMMs/$hmm.hmm > /tmp/$hmm.hmm\" -o /ecs2/work1/eae/GeneWiseHMM/genewiseHMM_".$species."_jobs_dir/$hmm-$dna_id.out -e /ecs2/work1/eae/GeneWiseHMM/genewiseHMM_".$species."_jobs_dir/$hmm-$dna_id.err -E \"/nfs/acari/eae/ensembl/ensembl-pipeline/scripts/HMMs/run_genewiseHMM.pl -check \" /nfs/acari/eae/ensembl/ensembl-pipeline/scripts/HMMs/run_genewiseHMM.pl -input_id $dna_id  -hmm /tmp/$hmm.hmm";
    print $command."\n";
  }    
}

sub get_input_ids{
  my $db = shift;
  my $threshold = shift;
  my $q = "SELECT protein_id, dna_id 
           FROM   genewisedb
           WHERE  bit_score > $threshold
             AND  bit_score < 10000";
  
  my $sth = $db->prepare($q) || $db->throw("can't prepare: $q");
  my $res = $sth->execute || $db->throw("can't execute: $q");
  
  my %ids;
  while( my ($protein_id,$dna_id) = $sth->fetchrow_array) {
    $ids{$protein_id}{$dna_id} = 1;
  }
  return \%ids;
}
