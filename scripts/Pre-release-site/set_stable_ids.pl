#!/usr/local/bin/perl -w
=head1 NAME

  set_stable_ids.pl

=head1 SYNOPSIS
 
  Script to create stable ids for the pre-release targetted gene set.
  Sets stable ids to be the name of the evidence used to buikld the structure, 
  plus a .version in case more than one gene is built from the same protein.

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
use Bio::EnsEMBL::Utils::Exception qw( deprecate throw warning );

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
my %transcript_proteins;
my %exons;
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $host,
					    '-user'   => $user,
					    '-dbname' => $name,
					    '-pass'   => $pass,
					    '-port'   => $port,
					   );

my $query = "insert into gene_stable_id values (?,?,?,now(),now())";
my $sth = $db->dbc->prepare($query);
my $trans_query = "insert into transcript_stable_id ".
  "values (?,?,?,now(),now())";
my $trans_sth = $db->dbc->prepare($trans_query);
my $translation_query = "insert into translation_stable_id ".
  "values (?,?,?,now(),now())";
my $translation_sth = $db->dbc->prepare($translation_query);
my $exon_query = "insert into exon_stable_id ".
  "values (?,?,?,now(),now())";
my $exon_sth = $db->dbc->prepare($exon_query);
GENE:
foreach my $gene_id(@{$db->get_GeneAdaptor->list_dbIDs}) {
  #print STDERR "gene id $gene_id\n";
  my $gene = $db->get_GeneAdaptor->fetch_by_dbID($gene_id);
  next GENE if $gene->stable_id;
  my $gene_protein_id = get_gene_protein_id($gene);
  #print "Have protein id ".$gene_protein_id."\n";
  throw("Gene ".$gene_id." has no protein id") unless defined $gene_protein_id;
  $proteins{$gene_protein_id}++;
  eval{
    $sth->execute($gene_id,$gene_protein_id."_".$proteins{$gene_protein_id}, $proteins{$gene_protein_id});
    foreach my $transcript(@{$gene->get_all_Transcripts}){
      my $transcript_protein_id = get_transcript_protein_id($transcript);
      $transcript_proteins{$transcript_protein_id}++;
      $trans_sth->execute($transcript->dbID, $transcript_protein_id."_".$transcript_proteins{$transcript_protein_id}, 
                          $transcript_proteins{$transcript_protein_id});
      my $translation = $transcript->translation;
      $translation_sth->execute($translation->dbID, $transcript_protein_id."_".$transcript_proteins{$transcript_protein_id},
                                $transcript_proteins{$transcript_protein_id});
      my $exon_count = 1;
      EXON:foreach my $exon(@{$transcript->get_all_Exons}){
        next EXON if($exons{$exon->dbID});
        my $stable_id = $transcript_protein_id."_".$transcript_proteins{$transcript_protein_id}.".".$exon_count;
        $exon_count++;
        #print "Storing ".$exon->dbID." with ".$stable_id." ".$transcript_proteins{$transcript_protein_id}."\n";
        $exon_sth->execute($exon->dbID, $stable_id, 
                           $transcript_proteins{$transcript_protein_id});
        $exons{$exon->dbID} = $stable_id;
      }
    }
  };
  if($@){
    throw("Stable id insertion for gene ".$gene_id." failed $@");
  }
}


sub get_gene_protein_id{
  my ($gene) = @_;
 
  my @exons = @{$gene->get_all_Exons};
  my $protein_id;
  
  foreach my $exon(@exons){
    foreach my $sf(@{$exon->get_all_supporting_features}){
      $protein_id = $sf->hseqname;
      return $protein_id if($protein_id);
    }
  }
  if(!$protein_id){
    throw("Found no protein id for ".$gene->dbID);
  }
}

sub get_transcript_protein_id{
  my ($transcript) = @_;
  my $protein_id;
  foreach my $sf(@{$transcript->get_all_supporting_features}){
    $protein_id = $sf->hseqname;
      return $protein_id if($protein_id);
  }
  if(!$protein_id){
    throw("Found no protein id for ".$transcript->dbID);
  }
}


