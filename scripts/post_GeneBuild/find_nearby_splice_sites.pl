#!/usr/local/ensembl/bin/perl -w

use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;


my $dbhost;
my $dbuser    = 'ensro';
my $dbname;
my $dbpass    = undef;
my $dnadbhost;
my $dnadbname;

my $genetype; # default genetype


$dbuser = "ensro";
GetOptions(
	   'dbname:s'    => \$dbname,
	   'dbhost:s'    => \$dbhost,
	   'dnadbhost:s' => \$dnadbhost,
	   'dnadbname:s' => \$dnadbname,
	   'genetype:s'  => \$genetype,
	  );

unless ( $dbname && $dbhost && $dnadbname && $dnadbhost ){
  print STDERR "script to check splice sites. If one site is not canonical it will try to\n";
  print STDERR "find a canonical one by shifting 1 and then 2 bases on either side\n";
  print STDERR "Usage: $0 -dbname -dbhost (-genetype)\n";
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

print STDERR "checking genes of type $genetype\n";


my @gene_ids      = @{$db->get_GeneAdaptor->list_geneIds};
my $slice_adaptor = $db->get_SliceAdaptor;

print "gene_stable_id\ttran_stable_id\ttotal_introns\tcorrect\tone_corrected\tboth_corrected\n";

GENE:
foreach my $gene_id ( @gene_ids){
  my $gene = $db->get_GeneAdaptor->fetch_by_dbID($gene_id);
  
  if ($genetype){
    next GENE unless ( $gene->type eq $genetype );
  }

  #unless ($gene->is_known){
  #  next GENE;
  #}

  #my $gene_stable_id = $gene->stable_id(); 
  
  
 TRANS:
  foreach my $trans ( @{$gene->get_all_Transcripts} ) {
    
    #print STDERR "trans_id: ".$trans->dbID."\n";
    my $slice = $slice_adaptor->fetch_by_transcript_id( $trans->dbID );
    #print STDERR "slice: ".$slice." ".$slice->name."\n";

    
    my $newgene = Bio::EnsEMBL::Gene->new();
    $newgene->add_Transcript($trans);
    $newgene->transform($slice);
    
    #my @mouse_evidence       = &get_mouse_only_evidence( $trans );
    #unless ( @mouse_evidence ){
    #  next TRANS;
    #}
    #my $tran_stable_id   = $trans->stable_id;
    
    my @transcripts = @{$newgene->get_all_Transcripts};
    my $transcript = $transcripts[0];

    my ($intron_count,$correct,$one_corrected,$both_corrected) = &check_splice_sites($transcript);
    unless( $intron_count ){
      next TRANS;
    }
    print $gene_id."\t".$trans->dbID."\t".$intron_count."\t".$correct."\t".$one_corrected."\t".$both_corrected."\n";
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



sub check_splice_sites{
  my ($transcript) = @_;
  $transcript->sort;
  
  my $strand = $transcript->start_Exon->strand;
  my @exons  = @{$transcript->get_all_Exons};
  
  my $introns  = scalar(@exons) - 1 ; 
  if ( $introns <= 0 ){
    return (0,0,0,0);
  }
  
  my $correct        = 0;
  my $one_corrected  = 0;
  my $both_corrected = 0;
  my $wrong          = 0;
  my $other          = 0;
  
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
      elsif (  ($upstream_site eq 'GT' && &find_downstream_in_forward_strand('AG',$downstream_exon,$slice) ) ||
	       (&find_upstream_in_forward_strand('GT',$upstream_exon,$slice) && $downstream_site eq 'AG') ||

	       ($upstream_site eq 'AT' && &find_downstream_in_forward_strand('AC',$downstream_exon,$slice) ) ||
	       (&find_upstream_in_forward_strand('AT',$upstream_exon,$slice) && $downstream_site eq 'AC') ||

	       ($upstream_site eq 'GC' && &find_downstream_in_forward_strand('AG',$downstream_exon,$slice) ) ||
	       (&find_upstream_in_forward_strand('GC',$upstream_exon,$slice) && $downstream_site eq 'AG') 
	    ){
	$one_corrected++;
      }
      elsif (  ( &find_upstream_in_forward_strand('GT',$upstream_exon,$slice) 
		 && 
		 &find_downstream_in_forward_strand('AG',$downstream_exon,$slice) ) 
	       ||
	       ( &find_upstream_in_forward_strand('GT',$upstream_exon,$slice) 
		 && 
		 &find_downstream_in_forward_strand('AC',$downstream_exon,$slice) ) 
	       ||
	       ( &find_upstream_in_forward_strand('GC',$upstream_exon,$slice)
		 && 
		 &find_downstream_in_forward_strand('AG',$downstream_exon,$slice) ) 
	    ){
	$both_corrected++;
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
      elsif (  ($upstream_site eq 'GT' && &find_downstream_in_reverse_strand('AG',$downstream_exon,$slice) ) ||
	       (&find_upstream_in_reverse_strand('GT',$upstream_exon,$slice) && $downstream_site eq 'AG') ||
	       
	       ($upstream_site eq 'AT' && &find_downstream_in_reverse_strand('AC',$downstream_exon,$slice) ) ||
	       (&find_upstream_in_reverse_strand('AT',$upstream_exon,$slice) && $downstream_site eq 'AC') ||
	       
	       ($upstream_site eq 'GC' && &find_downstream_in_reverse_strand('AG',$downstream_exon,$slice) ) ||
	       (&find_upstream_in_reverse_strand('GC',$upstream_exon,$slice) && $downstream_site eq 'AG') 
	    ){
	$one_corrected++;
      }
      elsif (  ( &find_upstream_in_reverse_strand('GT',$upstream_exon,$slice) 
		 && 
		 &find_downstream_in_reverse_strand('AG',$downstream_exon,$slice) ) 
	       ||
	       ( &find_upstream_in_reverse_strand('GT',$upstream_exon,$slice) 
		 && 
		 &find_downstream_in_reverse_strand('AC',$downstream_exon,$slice) ) 
	       ||
	       ( &find_upstream_in_reverse_strand('GC',$upstream_exon,$slice)
		 && 
		 &find_downstream_in_reverse_strand('AG',$downstream_exon,$slice) ) 
	    ){
	$both_corrected++;
      }
    } # end of INTRON
  }
  return ( $introns, $correct, $one_corrected, $both_corrected );
}


############################################################

sub find_downstream_in_forward_strand{
  my ($site,$downstream_exon,$slice) = @_;
  my @shifts = (-1,+1,-2,+2);
  foreach my $shift (@shifts){
    my $downstream_site = 
      $slice->subseq( ($downstream_exon->start - 2 + $shift), ($downstream_exon->start - 1 + $shift) );
    if ( $downstream_site eq $site ){
      return 1;
    }
  }
  return 0;
}

############################################################

sub find_upstream_in_forward_strand{
  my ($site,$upstream_exon,$slice) = @_;
  my @shifts = (-1,+1,-2,+2);
  foreach my $shift (@shifts){
    my $upstream_site = 
      $slice->subseq( ($upstream_exon->end     + 1 + $shift), ($upstream_exon->end     + 2 + $shift) );
    if ( $upstream_site eq $site ){
      return 1;
    }
  }
  return 0;
}

############################################################

sub find_downstream_in_reverse_strand{
  my ($site,$downstream_exon,$slice) = @_;
  my @shifts = (-1,+1,-2,+2);
  foreach my $shift (@shifts){
    my $downstream_site;
    my $down_site = 
      $slice->subseq( ($downstream_exon->end + 1 + $shift), ($downstream_exon->end + 2 + $shift) );
    ( $downstream_site = reverse( $down_site ) ) =~ tr/ACGTacgt/TGCAtgca/;
    if ( $downstream_site eq $site ){
      return 1;
    }
  }
  return 0;
}

############################################################

sub find_upstream_in_reverse_strand{
  my ($site,$upstream_exon,$slice) = @_;
  my @shifts = (-1,+1,-2,+2);
  foreach my $shift (@shifts){
    my $upstream_site;
    my $up_site = 
      $slice->subseq( ($upstream_exon->start - 2 + $shift), ($upstream_exon->start - 1 + $shift) );
    ( $upstream_site   = reverse(  $up_site  ) ) =~ tr/ACGTacgt/TGCAtgca/;
    if ( $upstream_site eq $site ){
      return 1;
    }
  }
  return 0;
}

############################################################
