#!/usr/local/bin/perl -w

=head1 NAME

  exonerate2mysqltab.pl

=head1 SYNOPSIS
 
  exonerate2mysqltab
  filters exonerate results based on percent_id, and makes a tab
  delimited file suitable for loading into feature table.
  Reads input from STDIN.

=head1 DESCRIPTION


=head1 OPTIONS

=cut

use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::ESTConf qw (
					EST_REFDBHOST
					EST_REFDBUSER
					EST_REFDBNAME
				       );


$| = 1;

my %cids;
my $analysis = 1; # from ens_UCSC_0801_est
my $name ='exonerate';

&get_golden_contigs;

while(<>){
  #print;
  chomp;
  #c006201126.1.55275      41440   41695   673.00  74      1       AB006208        142     400     -1      176M1I7M2D1M1D1M1D2M1D4M1I62M
  my ($contig_id, $contig_start, $contig_end, $score, $percent_id, $contig_strand, 
      $est_id, $est_start, $est_end, $est_strand, $cigar) = split();
  
  next unless $percent_id >= 90;

  $est_id =~ s/\S+\|\S+\|\S+\|(\S+)\|\S+/$1/;

  my $contig = $cids{$contig_id};
  #print STDERR "have feature from ".$contig."\n";
  next unless defined $contig;
  
  # if both contig and est strands are the same, convention is to set both to be 1
  # if they differ, convention is to set contig strand to -1, est strand to 1
  if($contig_strand == $est_strand){
    $contig_strand = 1;
    $est_strand    = 1;
  }
  else{
    $contig_strand = -1;
    $est_strand = 1;
  }
  
  # need to print out tab delimited line suitable for mysql. Autogenerate feature internal_id
  print "\\N\t$contig\t$contig_start\t$contig_end\t$contig_strand\t$analysis_id\t$est_start\t$est_end\t$est_strand\t$est_strand\t$cigar\t0\t$percent_id\t$score\n";

}

=head2 get_golden_contigs

  Title   : get_golden_contigs
  Usage   : get_golden_contigs
  Function: fetches internal ids and ids of golden contigs from reference database
  Returns : none - uses globals
  Args    : none - uses globals

=cut

sub get_golden_contigs{
  # get information about golden contigs from reference database
  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					      -host   => $EST_REFDBHOST,
					      -user   => $EST_REFDBUSER,
					      -dbname => $EST_REFDBNAME,
					     );
    
    my $query = "select contig.name, contig.contig_id from contig, assembly, meta where contig.contig_id=assembly.contig_id and assembly.type= meta.meta_value and meta.meta_key='assembly.default'";
  #print STDERR $query."\n";
  my $sth = $db->prepare($query) || $db->throw("can't prepare: $query");
  my $res = $sth->execute        || $db->throw("can't execute: $query");
  
  while( my ($internal_id, $id) = $sth->fetchrow_array) {
    #print STDERR "have ".$id." ".$internal_id."\n";
    $cids{$internal_id} = $id;
  }
  
}
