#!/usr/local/ensembl/bin/perl

=head1 NAME

cluster_Genes.pl
 
=head1 DESCRIPTION

script to cluster transcripts into genes. It reads genes and takes all the transcripts and recluster
them into genes.


=head1 OPTIONS


=cut

use strict;  
use diagnostics;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::GeneComparison::GeneComparison;
use Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster;
use Getopt::Long;

## load all the parameters
use Bio::EnsEMBL::Pipeline::GeneComparison::GeneCompConf;

my $dbhost = 'ecs1d';
my $dbname = 'homo_sapiens_estgene_9_30';
my $dbuser = 'ensro'; 

my $dnadbhost	= 'ecs1d';
my $dnadbname   = 'homo_sapiens_core_9_30';
my $dnadbuser   = 'ensro';

my $input_id;
my $write  = 0;
my $check  = 0;
my $pepfile;
my $gff_file;

my $genome;
my @genetypes;

# options
&GetOptions( 
	    'input_id:s'  => \$input_id,
	    'dbname:s'    => \$dbname,
	    'dbhost:s'    => \$dbhost,
	    'genome'      => \$genome,
	    'genetypes:s'  => \@genetypes,
	   );

unless( $input_id || $genome ){     
  print STDERR "Usage: $0 -dbname -dbhost\n";
  print STDERR "          -input_id < optional: chrname.chrstart-chrend >\n";
  print STDERR "          -genome ( for all the genes )\n";
  print STDERR "          -genetypes ( optional: omit for all genetypes)\n";
  exit(0);
}
    
# connect to the databases 
# only use one of them
#my $dna_db= new Bio::EnsEMBL::DBSQL::DBAdaptor(
#					       -host  => $dnadbhost,
#					       -user  => $dnadbuser,
#					       -dbname=> $dnadbname,
#					      );

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    -host  => $dbhost,
					    -user  => $dbuser,
					    -dbname=> $dbname,
#					    -dnadb => $dna_db,
					   );


my $sa = $db->get_SliceAdaptor;

my @slices;
my @all_clusters;
my $transcript_count = 0;

if ( $input_id ){
  # get genomic region 
  my $chr      = $input_id;
  $chr         =~ s/\.(.*)-(.*)//;
  my $chrstart = $1;
  my $chrend   = $2;
  
  unless ( $chr && $chrstart && $chrend ){
    print STDERR "bad input_id option, try something like 20.1-5000000\n";
  }
  
  print STDERR "Fetching region $chr, $chrstart - $chrend\n";
  my $slice = $sa->fetch_by_chr_start_end($chr,$chrstart,$chrend);
  push (@slices, $slice );
}
elsif( $genome ){
  my %chr_lengths = %{ &get_chrlengths( $db) };
    
  foreach my $chr_name ( keys %chr_lengths ){
    my $start = 1;
    my
 $end   = 1000000;
    
    while ( $start < $chr_lengths{$chr_name} ){
      print STDERR "calling $chr_name . $start - $end\n";
      my $slice = $sa->fetch_by_chr_start_end($chr_name,$start,$end);
      $start += 1000000;
      $end   += 1000000;
      if ( $end > $chr_lengths{$chr_name} ){
	$end = $chr_lengths{$chr_name};
      }
      push ( @slices, $slice );
    }
  }
}

# get the genes of type @$type2

foreach my $slice ( @slices ){

  my @genes;
  my @transcripts;
  
  my @all_the_genes = @{$slice->get_all_Genes};
  
  if ( @genetypes ){
    foreach my $gene ( @all_the_genes ){
      my $type = $gene->type;
      if ( grep /^$type$/, @genetypes ){
	push ( @genes, $gene );
      }
    }
  }
  else{
    @genes = @all_the_genes;
  }
  foreach my $gene ( @genes ){
    push ( @transcripts, @{$gene->get_all_Transcripts} );
  }
  print STDERR scalar(@genes)." genes found\n";
  print STDERR "with ".scalar(@transcripts)." transcripts\n";
  my $clusters = &cluster_Transcripts( \@transcripts );
  print STDERR scalar(@$clusters)." clusters found\n";
  push (@all_clusters, @$clusters );
}

my %order;
print STDERR "Total number of clusters: ".scalar(@all_clusters)."\n";
foreach my $cluster ( @all_clusters ){
  # store for each key=number of transcripts, how many clusters have this number of transcripts:
  push( @{ $order{ scalar( @{$cluster->get_Transcripts} ) } }, $cluster );
}

foreach my $key ( sort{ $a <=> $b } keys %order ){
  print STDERR "number of clusters with $key transcript(): ".scalar( @{ $order{$key} } )."\n";
}

############################################################

sub get_transcript_start_end_strand {
  my ($transcript) = @_;
  my $start;
  my $end;
  
  my $start_exon = $transcript->start_Exon;
  my $end_exon = $transcript->end_Exon;
  
  if ($start_exon->strand == 1) {
    $start = $start_exon->start;
    $end   = $end_exon->end;
  } else {
    $end   = $start_exon->end;
    $start = $end_exon->start;
  }
  return ($start, $end, $start_exon->strand);
}

sub by_transcript_high {
  my $alow;
  my $blow;
  my $ahigh;
  my $bhigh;
  
  if ($a->start_Exon->strand == 1) {
    $alow = $a->start_Exon->start;
    $ahigh = $a->end_Exon->end;
  } else {
    $alow = $a->end_Exon->start;
    $ahigh = $a->start_Exon->end;
  }

  if ($b->start_Exon->strand == 1) {
    $blow = $b->start_Exon->start;
    $bhigh = $b->end_Exon->end;
  } else {
    $blow = $b->end_Exon->start;
    $bhigh = $b->start_Exon->end;
  }

  if ($ahigh != $bhigh) {
    return $ahigh <=> $bhigh;
  } else {
    return $alow <=> $blow;
  }
}
############################################################

sub cluster_Transcripts {

  my $transcripts_unsorted = shift;
  my @transcripts_unsorted = @$transcripts_unsorted;

  my @transcripts = sort by_transcript_high @transcripts_unsorted;

  my @clusters;

  # clusters transcripts by whether or not any exon overlaps with an exon in
  # another transcript (came from prune in GeneBuilder)
  foreach my $tran (@transcripts) {
    my @matching_clusters;
    my ($trans_start, $trans_end, $trans_strand) = get_transcript_start_end_strand($tran);
    
    #print "transcript limits: $trans_start $trans_end \n";

  CLUSTER: 
    foreach my $cluster (@clusters) {
      
      #print "Testing against cluster with limits " . $cluster->start ." to " . $cluster->end . "\n";
      
      if (!($trans_start > $cluster->end || $trans_end < $cluster->start) &&
           $trans_strand == $cluster->strand) {
	#print "In range\n";
        foreach my $cluster_transcript (@{$cluster->get_Transcripts()}) {
          foreach my $exon1 (@{$tran->get_all_Exons}) {
	    
            foreach my $cluster_exon (@{$cluster_transcript->get_all_Exons}) {
              if ($exon1->overlaps($cluster_exon) &&
                  $exon1->strand == $cluster_exon->strand) {
                push (@matching_clusters, $cluster);
                next CLUSTER;
              }
            }
          }
        }
      }
    }
    
    if (scalar(@matching_clusters) == 0) {
      #print STDERR "Found new cluster for " . $tran->dbID . "\n";
      my $newcluster = new Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster;
      $newcluster->put_Transcripts($tran);
      push(@clusters,$newcluster);
      
    } 
    elsif (scalar(@matching_clusters) == 1) {
      #print STDERR "Adding to cluster for " . $tran->dbID . "\n";
      $matching_clusters[0]->put_Transcripts($tran);

    } 
    else {
      # Merge the matching clusters into a single cluster
      #print STDERR "Merging clusters for " . $tran->dbID . "\n";
      my @new_clusters;
      my $merged_cluster = new Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster;
      foreach my $clust (@matching_clusters) {
        $merged_cluster->put_Transcripts(@{$clust->get_Transcripts});
      }
      $merged_cluster->put_Transcripts($tran);
      push @new_clusters,$merged_cluster;
      # Add back non matching clusters
      foreach my $clust (@clusters) {
        my $found = 0;
        MATCHING: 
	foreach my $m_clust (@matching_clusters) {
          if ($clust == $m_clust) {
            $found = 1;
            last MATCHING;
          }
        }
        if (!$found) {
          push @new_clusters,$clust;
        }
      }
      @clusters = @new_clusters;
    }
  }
  
  #Safety checks
  my $ntrans = 0;
  my %trans_check_hash;
  foreach my $cluster (@clusters) {
    $ntrans += scalar(@{$cluster->get_Transcripts});
    foreach my $trans (@{$cluster->get_Transcripts}) {
      if (defined($trans_check_hash{"$trans"})) {
        print STDERR ("Transcript " . $trans->dbID . " added twice to clusters\n");
      }
      $trans_check_hash{"$trans"} = 1;
    }
    if (!scalar(@{$cluster->get_Transcripts})) {
      print STDERR ("Empty cluster");
    }
  }
  
  if ($ntrans != scalar(@transcripts)) {
    print STDERR ("Not all transcripts have been added into clusters $ntrans and " . scalar(@transcripts). " \n");
  }
  #end safety checks

  print "Started with " . scalar(@transcripts) . " transcripts and ended with " . scalar(@clusters) . " clusters\n";

  return \@clusters;
}


############################################################

sub get_chrlengths{
  my $db   = shift;
  my $type = shift;

  my %chrhash;

  my $q = qq( SELECT chr.name, max(ass.chr_end)  
	      FROM   chromosome chr, assembly ass
              WHERE  chr.chromosome_id = ass.chromosome_id
	    GROUP BY ass.chromosome_id
            );

  my $sth = $db->prepare($q) || $db->throw("can't prepare: $q");
  my $res = $sth->execute || $db->throw("can't execute: $q");

  while( my ($chr, $length) = $sth->fetchrow_array) {
    $chrhash{$chr} = $length;
  }
  return \%chrhash;
}

