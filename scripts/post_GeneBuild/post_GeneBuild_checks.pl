#!/usr/local/bin/perl

=head1 NAME

   post_GeneBuild_checks.pl

=head1 DESCRIPTION

prototype for a suite of checks after the gene build. The checks are based on virtual
contigs of 5Mb (this size should be the same one as the one used during the genebuild). It checks so far:

1.- all exons in a gene are in the same strand
2.- also checks for folded transcripts
3.- it flags also single exon genes and from these, the ones which are longer than 50000 bases
4.- it also checks some function calls. Not all may be relevant for the genome at hand

For mouse denormalised contigs: in order to know the size of the chromosomes, 
you must have a file with the internal_ids and raw_contig ids.
For standard contigs this should be changed. 

It ought to read the database parameters from GeneConf but it all depends on whether your final genes
are in the same database as the one you put GeneConf.

=head1 OPTIONS

=cut

use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::SeqIO;
use Bio::EnsEMBL::Pipeline::GeneConf qw (
					 GB_FINALDBHOST
					 GB_FINALDBNAME
					 GB_DBHOST
					 GB_DBNAME
					 GB_FINAL_GENETYPE
					);




use Bio::EnsEMBL::Utils::Eprof('eprof_start','eprof_end','eprof_dump');

my $dbhost      = $GB_FINALDBHOST;
my $dbuser    = 'ensro';
my $dbname    = $GB_FINALDBNAME;
my $dbpass    = undef;

my $path      = 'CHR';

my $dnadbhost = $GB_DBHOST;
my $dnadbuser = 'ensro';
my $dnadbname = $GB_DBNAME;
my $dnadbpass = undef;

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
#$db->static_golden_path_type($path);
my $sa = $db->get_StaticGoldenPathAdaptor();

my $genetype= $GB_FINAL_GENETYPE;

# get genomic regions
my $chr_lengths = get_chrlengths($db,$db->static_golden_path_type);
my %chr_lengths = %$chr_lengths;

# take 1MB vc's and check genes for each vc:
foreach my $chr ( keys( %chr_lengths ) ){
  
  my $start   = 1;
  my $end     = 5000000;
  my $stop  = 0;

  while ( $start < $chr_lengths{$chr} ){
    my $vc =  $sa->fetch_VirtualContig_by_chr_start_end($chr,$start,$end);
    print STDERR "vc: $chr.$start-$end\n";
    my @genes;
    eval{
      @genes = $vc->get_Genes_by_Type($genetype);
    };
    if ($@){
      print STDERR "trouble trying to get genes by type $genetype in this vc\n";
      print STDERR $@;
    }
  GENE: 
    foreach my $gene ( @genes ){
      
      # check all exons are on the same strand
      my $strand;
      EXON:
	foreach my $exon($gene->get_all_Exons){
	  if(!defined $exon->strand || ($exon->strand != 1 && $exon->strand != -1)){
	    print STDERR "Exon " . $exon->dbID . " has no strand!!\n";
	    last EXON;
	  }
	  if(!defined $strand){ $strand = $exon->strand; }
	  if($exon->strand != $strand){
	    print STDERR "strand problem with gene " . $gene->dbID . "\n";
	    last EXON;
	  }
	}
	
      # check there are no folded transcripts
      TRANSCRIPT:
	foreach my $transcript($gene->each_Transcript){
	  my @exons = $transcript->get_all_Exons;
	  
	  # check also for long single-exons genes
	  if ( scalar( @exons ) == 1 ){
	    print STDERR "single exon transcript : ".$transcript->dbID."\n";
	    my $start  = $exons[0]->start;
	    my $end    = $exons[0]->end;
	    my $length = $end-$start+1;
	    if ( $length > 50000 ){
	      print STDERR "long exon : ".$exons[0]->dbID." length: $length\n";
	    }
	  }

	  my $i;
	  for ($i = 1; $i < $#exons; $i++) {
	    if ($exons[0]->strand == 1) {
	      if ($exons[$i]->start < $exons[$i-1]->end) {
		print STDERR "ERROR:  Transcript folds back on itself. Transcript : " . $transcript->dbID . "\n";
		foreach my $exon ( @exons ){
		  print STDERR $exon->start.":".$exon->end."\t";
		}
		print STDERR "\n";
		next TRANSCRIPT;
	      } 
	    } elsif ($exons[0]->strand == -1) {
	      if ($exons[$i]->end > $exons[$i-1]->start) {
		print STDERR "ERROR:  Transcript folds back on itself. Transcript : " . $transcript->dbID . "\n";
		foreach my $exon ( @exons ){
		  print STDERR $exon->start.":".$exon->end."\t";
		}
		print STDERR "\n";
		next TRANSCRIPT;
	      } 
	    } else {
	      print STDERR "EEEK:In transcript  " . $transcript->dbID . " No strand for exon - can't check for folded transcript\n";
	    }
	  }
	}
      }

    # check some function calls
    eval{
      my @genes     = $vc->get_all_Genes_exononly();
      my @genes     = $vc->get_all_Genes('evidence');
      
      my @tmpgenes    = $vc->get_Genes_by_Type('ensembl','evidence');
      my @preds     = $vc->get_all_PredictionFeatures;
      
      my @features  = $vc->get_all_SimilarityFeatures_above_score('cpg',25);
      @features     = $vc->get_all_SimilarityFeatures_above_score('trna',80);
      @features     = $vc->get_all_SimilarityFeatures_above_score('unigene.seq',80);
      @features     = $vc->get_all_SimilarityFeatures_above_score('embl_vertrna',80);
      @features     = $vc->get_all_SimilarityFeatures_above_score('dbEST',1,0);
      @features     = $vc->get_all_SimilarityFeatures_above_score('swall',1);
      @features     = $vc->get_all_SimilarityFeatures_above_score('human_mrna',1);
    };
    if ($@){
      print STDERR "Error testing a functionality:\n$@\n";
    }
    
    $start += 5000000;
    $end   += 5000000;
    if ( $end > $chr_lengths{$chr} ){
      $end = $chr_lengths{$chr};
    }
  }  
}


########################################################


sub get_chrlengths{
  my $db   = shift;
  my $type = shift;

  if (!$db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor')) {
    die "get_chrlengths should be passed a Bio::EnsEMBL::DBSQL::DBAdaptor\n";
  }

  my %chrhash;

  my $q = qq( SELECT chr_name,max(chr_end) FROM static_golden_path as sgp
               WHERE sgp.type = '$type' GROUP BY chr_name
            );

  my $sth = $db->prepare($q) || $db->throw("can't prepare: $q");
  my $res = $sth->execute || $db->throw("can't execute: $q");

  while( my ($chr, $length) = $sth->fetchrow_array) {
    $chrhash{$chr} = $length;
  }
  return \%chrhash;
}
