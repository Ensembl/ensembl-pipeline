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
use Bio::EnsEMBL::Pipeline::ESTConf qw (
					EST_REFDBHOST
					EST_REFDBUSER
					EST_REFDBNAME
					EST_GOLDEN_PATH
				       );


$| = 1;

my %cids;
my $analysis = 1; # from ens_UCSC_0801_est
my $name ='exonerate';

&fetch_golden_contigs;

while(<>){
  # AC004073.1.1.79612      75445   76005   2293.00 89      1       gi|12779912|emb|AL516419.1|AL516419     296     861     1
  my ($contig_id, $contig_start, $contig_end, $score, $percent_id, $contig_strand, 
      $est_id, $est_start, $est_end, $est_strand) = split();
  
  next unless $percent_id >= 90;

  $est_id =~ s/\S+\|\S+\|\S+\|(\S+)\|\S+/$1/;

  my $contig = $cids{$contig_id};
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
  print "\\N\t$contig\t$contig_start\t$contig_end\t$score\t$contig_strand\t$analysis\t$name\t$est_start\t$est_end\t$est_id\tNULL\t$percent_id\t0\t0\n";

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
					     )
    
    my $query = "select contig.id, contig.internal_id from contig, static_golden_path where contig.internal_id=static_golden_path.raw_id and static_golden_path.type=$EST_GOLDEN_PATH";
  
  my $sth = $db->prepare($query) || $db->throw("can't prepare: $query");
  my $res = $sth->execute        || $db->throw("can't execute: $query");
  
  while( my ($internal_id, $id) = $sth->fetchrow_array) {
    $cids{$internal_id} = $id;
  }
  
}
