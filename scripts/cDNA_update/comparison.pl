#!/usr/local/ensembl/bin/perl

=pod

=head1 NAME

cDNA_setup - comparion

=head1 SYNOPSIS / OPTIONS

perl comparison.pl
         -chrom chomosome to ananlyse
         -oldname oldFeatureName
         -newname newFeatureName
         -dir dataDIR for output files
         -olddbhost db with older alignments (host)
         -olddbport db with older alignments (port)
         -olddbname db with older alignments (dbname)
         -newdbhost db with new alignments (host)
         -newdbport db with new alignments (port)
         -newdbname db with new alignments (dbname)
         -olddnadbhost db with older dna (host)
         -olddnadbport db with older dna (port)
         -olddnadbname db with older dna (dbname)
         -newdnadbhost db with new dna (host)
         -newdnadbport db with new dna (port)
         -newdnadbname db with new dna (dbname)

All arguments are required.

=head1 DESCRIPTION

This is a script to compare (cDNA-alignment) features in two databases.
It is used a a health-check / result-analysis after running the cDNA_update script,
which adds the latest cDNA information on an existing gene-build.
This script writes out a file for every chromosome containing the location of
features that were missed and that were added respectively in the new run.
It bases the comparison on genes and tries to trace the hits using the supporting evidence.

=head1 CONTACT

ensembl-dev@ebi.ac.uk

=cut


use warnings;
use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;

my ($chromosome, $oldFeatureName, $newFeatureName, $resultfilebase);
my ($olddbhost, $olddbport, $olddbname, $newdbhost, $newdbport, $newdbname);
my ($olddnadbhost, $olddnadbport, $olddnadbname, $newdnadbhost, $newdnadbport, $newdnadbname);

&GetOptions(
            'chrom=s'        => \$chromosome,
            'oldname=s'      => \$oldFeatureName,
            'newname=s'      => \$newFeatureName,
            'dir=s'          => \$resultfilebase,
            'olddbhost=s'    => \$olddbhost,
            'olddbport=s'    => \$olddbport,
            'olddbname=s'    => \$olddbname,
            'newdbhost=s'    => \$newdbhost,
            'newdbport=s'    => \$newdbport,
            'newdbname=s'    => \$newdbname,
            'olddnadbhost=s' => \$olddnadbhost,
            'olddnadbport=s' => \$olddnadbport,
            'olddnadbname=s' => \$olddnadbname,
            'newdnadbhost=s' => \$newdnadbhost,
            'newdnadbport=s' => \$newdnadbport,
            'newdnadbname=s' => \$newdnadbname
	   );

#connect dbs
#new dna_db
my $db1     = connect_db($newdnadbhost, $newdnadbport, $newdnadbname, "ensro");
my $sa1     = $db1->get_SliceAdaptor();
#new alignments
my $db2     = connect_db($newdbhost, $newdbport, $newdbname, "ensro", $db1);
my $sa2     = $db2->get_SliceAdaptor();

#old dna_db
my $db3      = connect_db($olddnadbhost, $olddnadbport, $olddnadbname, "ensro");
#old alignments
my $db4      = connect_db($olddbhost, $olddbport, $olddbname, "ensro", $db3);
my $sa4      = $db4->get_SliceAdaptor();

my( $gene, $chr_slice);
my (%results1, %results2);
my $num_alis_1 = 0;
my $num_alis_2 = 0;
my @missed_locations = ();
my @new_locations    = ();

my $resultfile  = $resultfilebase."/comparison_result_chr_".$chromosome.".txt";
open WP,">$resultfile" or die "cant open $resultfile.";

print WP "Type\tHit-name\tGene-ID\tChromosome\tGene-Start\tGene-End\tAligned Features & (new/old) Locations\n";

#new cDNA-update-db
$chr_slice = $sa2->fetch_by_region('chromosome', $chromosome);
#get genes
foreach $gene ( @{$chr_slice->get_all_Genes($newFeatureName)} ) {
  #check if theres a gene in older cDNA db here
  my $mini_slice = $sa4->fetch_by_region('chromosome', $chromosome, $gene->start, $gene->end);
  if (!scalar @{$mini_slice->get_all_Genes($oldFeatureName)}) {
    push(@new_locations, $gene->dbID);
  }
}
#get feature stats
$num_alis_1 = scalar @{ $chr_slice->get_all_DnaAlignFeatures($newFeatureName) };

#get infos about features at new locations (write to file)
_pursue_evidence(\@new_locations, "new", $db2, $db4, $oldFeatureName);


#old cdna-db
$chr_slice = $sa4->fetch_by_region('chromosome', $chromosome);
# get genes
foreach $gene ( @{$chr_slice->get_all_Genes($oldFeatureName)} ) {
  #check if theres a gene in the new cDNA db here
  my $mini_slice = $sa2->fetch_by_region('chromosome', $chromosome, $gene->start, $gene->end);
  if (!scalar @{$mini_slice->get_all_Genes($newFeatureName)}) {
    push(@missed_locations, $gene->dbID);
  }
}
#get feature stats
$num_alis_2 = scalar @{ $chr_slice->get_all_DnaAlignFeatures($oldFeatureName) };

#get infos about features at previous locations (write to file)
_pursue_evidence(\@missed_locations, "missed", $db4, $db2, $newFeatureName);

close WP;

$resultfile  = $resultfilebase."/comparison_result.txt";
open WP,">>$resultfile" or die "cant open global $resultfile.";

#print global results to shared file
print WP "\nChromosome $chromosome\n".
         "\tNumber of alignments (old): ".$num_alis_2.
         "\n\tNumber of alignments (new): ".$num_alis_1.
         "\n\tnew locations covered: ".scalar @new_locations.
         "\n\told locations missed: ".scalar @missed_locations."\n";

close WP;

########################################


sub _pursue_evidence{
  my $gene_ids      = shift;
  my $option        = shift;
  my $db1           = shift;
  my $db2           = shift;
  my $genetype      = shift;
  my $sql;
  my (%gene_hit_hash, %prev_gene_hit_hash);
  my $ga1 = $db1->get_GeneAdaptor;
  my $ga2 = $db2->get_GeneAdaptor;
  my $found_previous;

  $sql = ("select distinct(hit_name), gene_id ".
	  "from dna_align_feature, supporting_feature, ".
	  "exon_transcript, transcript ".
	  "where dna_align_feature_id = feature_id ".
	  "and feature_type = 'dna_align_feature' ".
	  "and supporting_feature.exon_id = exon_transcript.exon_id ".
	  "and exon_transcript.transcript_id = transcript.transcript_id ".
	  "and gene_id = ?");
  my $sth1 = $db1->dbc->prepare($sql) or die "sql error 1";

  $sql = ("SELECT distinct(gene.gene_id), hit_name ".
	  "FROM gene, transcript, exon_transcript, ".
	  "supporting_feature, dna_align_feature ".
	  "WHERE hit_name = ? ".
	  "AND dna_align_feature_id = feature_id ".
	  "AND feature_type = 'dna_align_feature' ".
	  "AND supporting_feature.exon_id = exon_transcript.exon_id ".
	  "AND exon_transcript.transcript_id = transcript.transcript_id ".
	  "AND transcript.gene_id = gene.gene_id ".
	  "AND gene.biotype = '".$genetype."'");
  my $sth2 = $db2->dbc->prepare($sql) or die "sql error 2";

  #print previous and new location
  foreach my $gene_id (@$gene_ids){
    $found_previous = 0;
    $sth1->execute($gene_id);
    my ($hit_name, $gene_id) = $sth1->fetchrow_array;
    $gene_hit_hash{$gene_id} = $hit_name;
    my $gene1 = $ga1->fetch_by_dbID($gene_id);
    print WP "\n".$option."\t".$hit_name."\t".$gene_id."\t".$gene1->slice->seq_region_name."\t".$gene1->start."\t".$gene1->end."\t";

    $sth2->execute($hit_name);
    while(my ($prev_gene_id, $hit_name) = $sth2->fetchrow_array){
      $found_previous++;
      $prev_gene_hit_hash{$hit_name} = $prev_gene_id;
      my $gene2 = $ga2->fetch_by_dbID($prev_gene_id);
      print WP $gene2->dbID."\t".$gene2->slice->seq_region_name."\t".$gene2->start."\t".$gene2->end;
    }

    if(!$found_previous){
      print WP "new feature";
    }
  }

  return 1;
}


sub connect_db{
  my $host       = shift;
  my $port       = shift;
  my $dbname     = shift;
  my $user       = shift;
  my $dnadb      = shift;
  my $dbObj;

  if($dnadb){
    $dbObj      = new Bio::EnsEMBL::DBSQL::DBAdaptor(
						     -host   => $host,
						     -port   => $port,
						     -user   => $user,
						     -dbname => $dbname,
						     -dnadb  => $dnadb
						    );
  }
  else{
    $dbObj      = new Bio::EnsEMBL::DBSQL::DBAdaptor(
						     -host    => $host,
						     -port    => $port,
						     -user    => $user,
						     -dbname  => $dbname
						    );
  }
  if(!$dbObj){
    die "\ncould not connect to \"$dbname\".\n";
  }
  return $dbObj;
}


1;
