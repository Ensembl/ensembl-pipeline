#!/usr/local/ensembl/bin/perl -w

use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;

my $dbhost;
my $dnadbhost;
my $dbuser    = 'ensro';
my $dbname;
my $dnadbname;
my $dbpass    = undef;

my $genetype;


$dbuser = "ensro";
GetOptions(
	   'dbname:s'    => \$dbname,
	   'dbhost:s'    => \$dbhost,
	   'dnadbname:s'  => \$dnadbname,
	   'dnadbhost:s'  => \$dnadbhost,
	   'genetype:s'  => \$genetype,
	  );

unless ( $dbname && $dbhost && $dnadbhost && $dnadbname ){
  print STDERR "script to check splice sites\n";
  print STDERR "Usage: $0 -dbname -dbhost -dnadbname -dnadbhost (-genetype)\n";
  exit(0);
}

my $dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					       '-host'   => $dnadbhost,
					       '-user'   => $dbuser,
					       '-dbname' => $dnadbname,
					       '-pass'   => $dbpass,
					      );


my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $dbhost,
					    '-user'   => $dbuser,
					    '-dbname' => $dbname,
					    '-pass'   => $dbpass,
					    '-dnadb'  => $dnadb,
					   );


print STDERR "connected to $dbname : $dbhost\n";

if ( $genetype ){
  print STDERR "checking genes of type $genetype\n";
}

my @gene_ids      = @{$db->get_GeneAdaptor->list_geneIds};
my $slice_adaptor = $db->get_SliceAdaptor;

print "gene_id\ttran_id\tcanonical\tnon_canonical\twrong\tevidence\n";
############################################################

GENE:
foreach my $gene_id ( @gene_ids){
  
  my $gene = $db->get_GeneAdaptor->fetch_by_dbID($gene_id,1);

  if ( $genetype ){
    next GENE unless ( $genetype eq $gene->type );
  }

 TRANS:
  foreach my $trans ( @{$gene->get_all_Transcripts} ) {
    #my @mouse_evidence       = &get_mouse_only_evidence( $trans );
    #unless ( @mouse_evidence ){
    #  next TRANS;
    #}
    my @evidence = &get_evidence($trans);

    my $tranid   = $trans->stable_id || $trans->dbID;
    my $geneid   = $gene->stable_id  || $gene->dbID;
    my ($canonical,$non_canonical,$wrong) = &check_splice_sites($trans);
    print $geneid."\t".$tranid."\t".$canonical."\t".$non_canonical."\t".$wrong."\t@evidence\n";
  }
}

############################################################

sub get_mouse_only_evidence{
  my ($t) = @_;
  my %evidence;
  my $mouse = 0;
  my $other = 0;
  foreach my $exon ( @{$t->get_all_Exons} ){
    foreach my $feature ( @{$exon->get_all_supporting_features}){
      if (  $feature->isa('Bio::EnsEMBL::DnaDnaAlignFeature') 
	    && 
	    $feature->analysis->logic_name eq 'mouse_cdna' ){
	
	#print STDERR "logic_name:".$feature->analysis->logic_name. " evidence: ".$feature->hseqname."\n";
	$mouse = 1;
	$evidence{ $feature->hseqname } = 1;
      }
      else{
	$other = 1;
      }
    }
  }
  unless ( $other ==  1){
    return keys %evidence;
  }
  return ();
}

############################################################

sub get_evidence{
  my ($t) = @_;
  my %evidence;
  my $mouse = 0;
  my $other = 0;
  foreach my $exon ( @{$t->get_all_Exons} ){
    foreach my $feature ( @{$exon->get_all_supporting_features}){
      $evidence{ $feature->hseqname } = 1;
    }
  }
  return keys %evidence;
}

############################################################

sub check_splice_sites{
  my ($transcript) = @_;
  $transcript->sort;
  
  my $strand = $transcript->start_Exon->strand;
  my @exons  = @{$transcript->get_all_Exons};
  
  my $introns  = scalar(@exons) - 1 ; 
  if ( $introns <= 0 ){
    return (0,0,0);
  }
  
  my $correct  = 0;
  my $wrong    = 0;
  my $other    = 0;
  
  # all exons in the transcripts are in the same seqname coordinate system:
  my $slice = $transcript->start_Exon->contig;
  
  if ($strand == 1 ){
    
  INTRON:
    for (my $i=0; $i<$#exons; $i++ ){
      my $upstream_exon   = $exons[$i];
      my $downstream_exon = $exons[$i+1];
      my $upstream_site;
      my $downstream_site;
      eval{
	$upstream_site = 
	  $slice->subseq( ($upstream_exon->end     + 1), ($upstream_exon->end     + 2 ) );
	$downstream_site = 
	  $slice->subseq( ($downstream_exon->start - 2), ($downstream_exon->start - 1 ) );
      };
      unless ( $upstream_site && $downstream_site ){
	print STDERR "problems retrieving sequence for splice sites\n$@";
	next INTRON;
      }
      
      #print STDERR "upstream $upstream_site, downstream: $downstream_site\n";
      ## good pairs of upstream-downstream intron sites:
      ## ..###GT...AG###...   ...###AT...AC###...   ...###GC...AG###.
      
      ## bad  pairs of upstream-downstream intron sites (they imply wrong strand)
      ##...###CT...AC###...   ...###GT...AT###...   ...###CT...GC###...
      
      if (  ($upstream_site eq 'GT' && $downstream_site eq 'AG') ||
	    ($upstream_site eq 'AT' && $downstream_site eq 'AC') ||
	    ($upstream_site eq 'GC' && $downstream_site eq 'AG') ){
	$correct++;
      }
      elsif (  ($upstream_site eq 'CT' && $downstream_site eq 'AC') ||
	       ($upstream_site eq 'GT' && $downstream_site eq 'AT') ||
	       ($upstream_site eq 'CT' && $downstream_site eq 'GC') ){
	$wrong++;
      }
      else{
	$other++;
      }
    } # end of INTRON
  }
  elsif ( $strand == -1 ){
    
    #  example:
    #                                  ------CT...AC---... 
    #  transcript in reverse strand -> ######GA...TG###... 
    # we calculate AC in the slice and the revcomp to get GT == good site
    
  INTRON:
    for (my $i=0; $i<$#exons; $i++ ){
      my $upstream_exon   = $exons[$i];
      my $downstream_exon = $exons[$i+1];
      my $upstream_site;
      my $downstream_site;
      my $up_site;
      my $down_site;
      eval{
	$up_site = 
	  $slice->subseq( ($upstream_exon->start - 2), ($upstream_exon->start - 1) );
	$down_site = 
	  $slice->subseq( ($downstream_exon->end + 1), ($downstream_exon->end + 2 ) );
      };
      unless ( $up_site && $down_site ){
	print STDERR "problems retrieving sequence for splice sites\n$@";
	next INTRON;
      }
      ( $upstream_site   = reverse(  $up_site  ) ) =~ tr/ACGTacgt/TGCAtgca/;
      ( $downstream_site = reverse( $down_site ) ) =~ tr/ACGTacgt/TGCAtgca/;
      
      #print STDERR "upstream $upstream_site, downstream: $downstream_site\n";
      if (  ($upstream_site eq 'GT' && $downstream_site eq 'AG') ||
	    ($upstream_site eq 'AT' && $downstream_site eq 'AC') ||
	    ($upstream_site eq 'GC' && $downstream_site eq 'AG') ){
	$correct++;
      }
      elsif (  ($upstream_site eq 'CT' && $downstream_site eq 'AC') ||
	       ($upstream_site eq 'GT' && $downstream_site eq 'AT') ||
	       ($upstream_site eq 'CT' && $downstream_site eq 'GC') ){
	$wrong++;
      }
      else{
	$other++;
      }
      
    } # end of INTRON
  }
  return ( $correct, $other, $wrong );
}
