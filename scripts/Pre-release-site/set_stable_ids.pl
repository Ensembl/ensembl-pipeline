#!/usr/local/bin/perl -w
=head1 NAME

  set_stable_ids.pl

=head1 SYNOPSIS
 
  Script to create stable ids for the pre-release targetted gene set.
  Sets stable ids to be the name of the evidence used to buikld the 
  structure, plus a .version in case more than one gene is built from the 
  same protein.

=head1 DESCRIPTION


=head1 OPTIONS

    -dbhost      host name for database (gets put as host= in locator)

    -dbport      For RDBs, what port to connect to (port= in locator)

    -dbname      For RDBs, what name to connect to (dbname= in locator)

    -dbuser      For RDBs, what username to connect as (dbuser= in locator)

    -dbpass      For RDBs, what password to use (dbpass= in locator)

=cut
use strict;
use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;

my $host;
my $port;
my $name;
my $user;
my $pass;

&GetOptions( 
	     'dbhost:s'      => \$host,
	     'dbport:n'      => \$port,
	     'dbname:s'      => \$name,
	     'dbuser:s'      => \$user,
	     'dbpass:s'      => \$pass,
	     );

my %proteins;

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $host,
					    '-user'   => $user,
					    '-dbname' => $name,
					    '-pass'   => $pass,
					    '-port'   => $port,
					   );

my $gene_query = "insert into gene_stable_id values (?,?,?,".
                 "FROM_UNIXTIME(?), FROM_UNIXTIME(?))";
my $gene_sth = $db->dbc->prepare($gene_query);
my $trans_query = "insert into transcript_stable_id values (?,?,?,".
                 "FROM_UNIXTIME(?), FROM_UNIXTIME(?))";
my $trans_sth = $db->prepare($trans_query);
my $tln_query = "insert into translation_stable_id values (?,?,?,".
                 "FROM_UNIXTIME(?), FROM_UNIXTIME(?))";
my $tln_sth = $db->prepare($tln_query);
my $exon_query = "insert into exon_stable_id values (?,?,?,".
                 "FROM_UNIXTIME(?), FROM_UNIXTIME(?))";
my $exon_sth = $db->prepare($exon_query);

my $time = time;

GENE:
foreach my $gene_id(@{$db->get_GeneAdaptor->list_dbIDs}) {
  print "gene id $gene_id\n";

  my $gene = $db->get_GeneAdaptor->fetch_by_dbID($gene_id);
  my @exons = @{$gene->get_all_Exons};
  my $protein_id;
 EXON:foreach my $e(@exons){
    my @sf = @{$e->get_all_supporting_features};
    foreach my $f(@sf){
      last EXON if defined $protein_id;
      eval{
        $protein_id = $f->hseqname;
      }
    }
  }
  if(!$protein_id){
    print "Gene ".$gene->dbID." seems to have no supporting features\n";
    next GENE;
  }
  if(!$proteins{$protein_id}){
    $proteins{$protein_id} = 1;
  }
  my $stable_id = "Built_from_".$protein_id;
  print "Stable_id ".$stable_id."\n";
  $gene_sth->execute($gene_id,$stable_id, $proteins{$protein_id}, 
                $time, $time);
  my $exon_count = 1;
  foreach my $exon(@exons){
    $exon_sth->execute($exon->dbID, $stable_id."_".$exon_count, 
                  $proteins{$protein_id}, $time, $time);
    $exon_count++;
  }
  foreach my $transcript(@{$gene->get_all_Transcripts}){
    $trans_sth->execute($transcript->dbID, $stable_id, 
                        $proteins{$protein_id}, $time, $time);
    if($transcript->translation){
      $tln_sth->execute($transcript->translation->dbID, $stable_id,
                        $proteins{$protein_id}, $time, $time);
    }
  }
  $proteins{$protein_id}++;
}

