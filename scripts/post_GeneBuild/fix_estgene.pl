#!/usr/local/bin/perl -w

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Bio::EnsEMBL::Pipeline::Tools::GeneUtils;

use Getopt::Long;
use strict;


my $genetype;
my $info;
my $source_dbname = 'human_estgenes_eae';
my $source_dbhost = 'ecs2b';
my $source_dnadbname = 'homo_sapiens_core_15_33';
my $source_dnadbhost = 'ecs2f';
my $target_dbname;
my $target_dbhost;
my $gene_id;
my $check;

&GetOptions('gene_id:s'          => \$gene_id,
	    'target_dbname:s'  => \$target_dbname,
            'target_dbhost:s'  => \$target_dbhost,
	    'check'            => \$check,
	    );

# database where we read genes from
my $source_dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
						      '-host'   => $source_dnadbhost,
						      '-user'   => 'ensro',
						      '-dbname' => $source_dnadbname,
						      );

my $source_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
						   '-host'   => $source_dbhost,
						   '-user'   => 'ensro',
						   '-dbname' => $source_dbname,
						   '-dnadb'  => $source_dnadb,
						   );

if ( $check ){
    exit(0);
}

unless( $gene_id ){
    print "script to fix the splice sites and translations of estgenes\n";
    print "Usage: $0 -gene_id -target_dbname -target_dbhost -check\n";
    exit(0);
}


# database where we put the genes
my  $target_db;
if ( $target_dbhost && $target_dbname ){
    $target_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
						    '-host'   => $target_dbhost,
						    '-user'   => 'ensadmin',
						    '-dbname' => $target_dbname,
						    '-pass'   => 'ensembl',
						    );
}

if ( $check ){
    exit(0);
}

# fetch gene
my $source_gene_adaptor = $source_db->get_GeneAdaptor;

my $target_gene_adaptor;
if ( $target_db ){
    $target_gene_adaptor = $target_db->get_GeneAdaptor;
}

my $gene;
eval {
    print STDERR "fetching $gene_id\n";
    $gene = $source_gene_adaptor->fetch_by_stable_id($gene_id,1);
};
if ( $@ ){
    print STDERR "Unable to read gene $gene_id $@\n";
}

############################################################
# clone the gene 
############################################################
print STDERR "gene : ".scalar(@{$gene->get_all_Transcripts})." transcripts, ".
	    scalar(@{$gene->get_all_Exons})." exons\n";

my $cloned_gene = Bio::EnsEMBL::Pipeline::Tools::GeneUtils->_clone_Gene($gene);
print STDERR "cloned gene : ".scalar(@{$cloned_gene->get_all_Transcripts})." transcripts, ".
	    scalar(@{$cloned_gene->get_all_Exons})." exons\n";

my $newgene = &_clone_Gene_only($gene);

############################################################
# fix the transcripts
############################################################

 TRANS:
    foreach my $transcript (@{$cloned_gene->get_all_Transcripts}){
	my $trans_id = $transcript->stable_id || $transcript->dbID;
	
	print STDERR "\n --- checking transcript $trans_id ---\n";
	
	############################################################
	my @evidence = &get_evidence( $transcript );
	#print STDERR "going to check transcript $trans_id:\n";
	#Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript($transcript);
	
	############################################################
	# modify splice sites (not more than 2bp) if necessary
	############################################################
	print STDERR "before fix_splice_sites:\n";
      Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript( $transcript );

	my ($newtranscript1,$intron_count,$correct,$one_corrected,$both_corrected) = &fix_splice_sites($transcript);
	print STDERR "after fix_splice_sites: end_exon ".$newtranscript1->translation->end_Exon."length: ".$newtranscript1->translation->end_Exon->length."\n";
	if ( $transcript->translation ){
	    $newtranscript1->translation( $transcript->translation );
	}
	
	if ( $one_corrected || $both_corrected ){
	    print STDERR "SS_CORRECTED ".$gene_id."\t".$trans_id."\t".$intron_count."\t".$correct."\t".$one_corrected."\t".$both_corrected."\t@evidence\n";
		
	    if ( $intron_count == $correct + $one_corrected + $both_corrected ){
		print STDERR "all introns are correct now\n";
	    }
	}
	if ( Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->check_splice_sites($newtranscript1)  ){
	    print STDERR "check: introns correct\n";
	}
	else{
	    print STDERR "not all splice sites are correct\n";
	}
	print STDERR "after fix_splice_sites:\n";
      Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript( $newtranscript1 );

	############################################################
	# modify translation if necessary
	############################################################
	my $newtranscript2;
	if ( Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils
	     ->_check_Translation($newtranscript1) ){
	    print STDERR "good-translation for $trans_id - not modifying\n";
	    $newtranscript2 = $newtranscript1;
	}
	else{
	    print STDERR "stop-codons for $trans_id - modifying\n";
	    $newtranscript2 = &fix_translation( $newtranscript1 );
	    print STDERR "after fix_translaton: end_exon ".$newtranscript2->translation->end_Exon." length: ".$newtranscript2->translation->end_Exon->length."\n";
	    
	}

	############################################################
	# extend the translation object to include a stop if there is one:
	############################################################
	my $newtranscript3 = Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils
	    ->set_stop_codon($newtranscript2);
	print STDERR "after set_stop_codon: end_exon length: ".$newtranscript3->translation->end_Exon->length."\n";
	
	############################################################
	# check translation - just in case
	############################################################
	unless (Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_check_Translation($newtranscript3) ){
	    print STDERR "FAILED: translation has been corrupted\n";
	}
	############################################################
	# add the transcript to the cloned gene
	############################################################
	if ( $transcript->stable_id ){
	    $newtranscript3->stable_id( $transcript->stable_id );
	    #$newtranscript2->stable_id( $transcript->stable_id );
	}
	
	$newgene->add_Transcript($newtranscript3);

    }

&prune_Exons($newgene);

print STDERR "newgene : ".scalar(@{$newgene->get_all_Transcripts})." transcripts, ".
    scalar(@{$newgene->get_all_Exons})." exons\n";

unless ( scalar(@{$newgene->get_all_Transcripts}) > 0 ){
    print STDERR "gene $gene_id has been left without good transcripts, skipping\n";
    exit(0);
}

############################################################

if ( $target_db ){
    
    ##############################
    # store gene, need to transform
    ##############################
    eval{
	$newgene->transform;
	print STDERR "storing gene $gene_id...\n";
	$target_gene_adaptor->store($newgene);
    };
    if ( $@ ){
	print STDERR "Unable to store newgene ".$newgene->dbID." $@\n";
    }
    
    ##############################
    # info about the stored gene
    ##############################
    
    if ($info){
	print STDERR "stored gene: ".$newgene->dbID."\n";
	my @transcripts = @{$newgene->get_all_Transcripts};
	
	foreach my $tran (@transcripts){
	  Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Evidence($tran);
	    
	    if ( defined $tran->translation){
	      Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Translation($tran->translation);
	    }
	    else { 
		print STDERR "translation not found!\n";
	    }
	}
	    }
}
else{
    print STDERR "Not storing gene\n";
}






########################################################


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


sub fix_splice_sites{
  my ($transcript) = @_;

  my $verbose = 1;

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
  
  #print STDERR "=== strand = $strand ===\n";
  
  if ($strand == 1 ){
    
    INTRON:
      for (my $i=0; $i<$#exons; $i++ ){
	  # upstream  : $exons[$i];
	  # downstream: $exons[$i+1];
	  my $upstream_site;
	  my $downstream_site;
	  eval{
	      $upstream_site = 
		  $slice->subseq( ($exons[$i]->end     + 1), ($exons[$i]->end     + 2 ) );
	      $downstream_site = 
		  $slice->subseq( ($exons[$i+1]->start - 2), ($exons[$i+1]->start - 1 ) );
	  };
	  unless ( $upstream_site && $downstream_site ){
	      print STDERR "problems retrieving sequence for splice sites\n$@";
	      next INTRON;
	  }
	  
	  print STDERR "fix_splice_sites: upstream $upstream_site, downstream: $downstream_site\n";
	  ## good pairs of upstream-downstream intron sites:
	  ## ..###GT...AG###...   ...###AT...AC###...   ...###GC...AG###.
	  
	  ## bad  pairs of upstream-downstream intron sites (they imply wrong strand)
	  ##...###CT...AC###...   ...###GT...AT###...   ...###CT...GC###...
	  
	  if (  ($upstream_site eq 'GT' && $downstream_site eq 'AG') ||
		($upstream_site eq 'AT' && $downstream_site eq 'AC') ||
		($upstream_site eq 'GC' && $downstream_site eq 'AG') ){
	      print STDERR "correct\n" if $verbose;
	      $correct++;
	  }
	  elsif (  $upstream_site eq 'GT' ){
	      my $shift = &find_downstream_in_forward_strand('AG',$exons[$i+1],$slice); 
	      if ( $shift != 0 ){
		  my $start = $exons[$i+1]->start;
		  $exons[$i+1]->start( $start + $shift );
		  $one_corrected++;
		  print STDERR "$exons[$i+1] corrected to AG, shift = $shift\n" if $verbose;
	      }
	  }
	  elsif( $downstream_site eq 'AG' ){
	      my $shift = &find_upstream_in_forward_strand('GT',$exons[$i],$slice);
	      if ( $shift != 0 ){
		  my $end = $exons[$i]->end;
		  $exons[$i]->end( $end + $shift );
		  $one_corrected++;
		  print STDERR "$exons[$i] corrected to GT, shift = $shift\n" if $verbose;
	      }
	  }
	  elsif( $upstream_site eq 'AT' ){
	      my $shift = &find_downstream_in_forward_strand('AC',$exons[$i+1],$slice);
	      if ( $shift != 0 ){
		  my $start = $exons[$i+1]->start;
		  $exons[$i+1]->start( $start + $shift );
		  $one_corrected++;
		  print STDERR "$exons[$i+1] corrected to AC, shift = $shift\n" if $verbose;
	      }
	  }
	  elsif(   $downstream_site eq 'AC' ){
	      my $shift = &find_upstream_in_forward_strand('AT',$exons[$i],$slice);
	      if ( $shift != 0 ){ 
		  my $end = $exons[$i]->end;
		  $exons[$i]->end( $end + $shift ); 
		  $one_corrected++;
		  print STDERR "$exons[$i] corrected to AT, shift = $shift\n" if $verbose;
	      }
	  }
	  elsif( $upstream_site eq 'GC' ){
	      my $shift = &find_downstream_in_forward_strand('AG',$exons[$i+1],$slice);
	      if ( $shift != 0 ){
		  my $start = $exons[$i+1]->start;
		  $exons[$i+1]->start( $start + $shift );
		  $one_corrected++;
		  print STDERR "$exons[$i+1] corrected to AG, shift = $shift\n" if $verbose;
	      }
	  }
	  elsif( $downstream_site eq 'AG' ){
	      my $shift = &find_upstream_in_forward_strand('GC',$exons[$i],$slice);
	      if ( $shift != 0 ){
		  my $end = $exons[$i]->end;
		  $exons[$i]->end( $end + $shift ); 
		  $one_corrected++;
		  print STDERR "$exons[$i] corrected to GC, shift = $shift\n" if $verbose;
	      }
	  }
	  else {
	      my $shift_GT =  &find_upstream_in_forward_strand('GT',$exons[$i],$slice);
	      my $shift_AG =  &find_downstream_in_forward_strand('AG',$exons[$i+1],$slice);
	      my $shift_AT =  &find_upstream_in_forward_strand('AT',$exons[$i],$slice);
	      my $shift_AC =  &find_downstream_in_forward_strand('AC',$exons[$i+1],$slice);
	      my $shift_GC =  &find_upstream_in_forward_strand('GC',$exons[$i],$slice);
	      if ( $shift_GT != 0 && $shift_AG != 0 ){
		  my $start = $exons[$i+1]->start;
		  $exons[$i+1]->start( $start + $shift_AG );
		  my $end = $exons[$i]->end;
		  $exons[$i]->end( $end + $shift_GT ); 
		  $both_corrected++;
	      }
	      elsif( $shift_AT != 0 && $shift_AC != 0 ){
		  my $start = $exons[$i+1]->start;
		  $exons[$i+1]->start( $start + $shift_AC );
		  my $end = $exons[$i]->end;
		  $exons[$i]->end( $end + $shift_AT ); 
		  $both_corrected++;
	      }
	      elsif( $shift_GC != 0 && $shift_AG != 0 ){
		  my $start = $exons[$i+1]->start;
		  $exons[$i+1]->start( $start + $shift_AG );
		  my $end = $exons[$i]->end;
		  $exons[$i]->end( $end + $shift_GC ); 
		  $both_corrected++;
	      }
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
	  # $upstream_exon   = $exons[$i];
	  # $downstream_exon = $exons[$i+1];
	  my $upstream_site;
	  my $downstream_site;
	  my $up_site;
	  my $down_site;
	  eval{
	      $up_site = 
		  $slice->subseq( ($exons[$i]->start - 2), ($exons[$i]->start - 1) );
	      $down_site = 
		  $slice->subseq( ($exons[$i+1]->end + 1), ($exons[$i+1]->end + 2 ) );
	  };
	  unless ( $up_site && $down_site ){
	      print STDERR "problems retrieving sequence for splice sites\n$@";
	      next INTRON;
	  }
	  ( $upstream_site   = reverse(  $up_site  ) ) =~ tr/ACGTacgt/TGCAtgca/;
	  ( $downstream_site = reverse( $down_site ) ) =~ tr/ACGTacgt/TGCAtgca/;
	  
	  print STDERR "fix_splice_sites: upstream $upstream_site, downstream: $downstream_site\n";
	  if (  ($upstream_site eq 'GT' && $downstream_site eq 'AG') ||
		($upstream_site eq 'AT' && $downstream_site eq 'AC') ||
		($upstream_site eq 'GC' && $downstream_site eq 'AG') ){
	      $correct++;
	      print STDERR "correct\n" if $verbose;
	  }
	  elsif (  $upstream_site eq 'GT' ){
	      my $shift = &find_downstream_in_reverse_strand('AG',$exons[$i+1],$slice); 
	      print STDERR "shift: $shift\n";
	      if ( $shift != 0 ){
		  print STDERR "shift: $shift\n";
		  my $end = $exons[$i+1]->end;
		  $exons[$i+1]->end( $end + $shift );
		  $one_corrected++;
		  print STDERR "$exons[$i+1] corrected to AG, shift = $shift\n" if $verbose;
	      }
	   }
	  elsif( $downstream_site eq 'AG' ){
	      my $shift = &find_upstream_in_reverse_strand('GT',$exons[$i],$slice);
	      print STDERR "shift: $shift\n";
	      if ( $shift != 0 ){
		  print STDERR "shift: $shift\n";
		    my $start = $exons[$i]->start;
		  $exons[$i]->start( $start + $shift );
		  $one_corrected++;
		  print STDERR "$exons[$i] corrected to GT, shift = $shift\n" if $verbose;
	      }
	  }
	  elsif( $upstream_site eq 'AT' ){
	      my $shift = &find_downstream_in_reverse_strand('AC',$exons[$i+1],$slice);
	      if ( $shift != 0 ){
		  my $end = $exons[$i+1]->end;
		  $exons[$i+1]->end( $end + $shift );
		  $one_corrected++;
		  print STDERR "$exons[$i+1] corrected to AC, shift = $shift\n" if $verbose;
	      }
	  }
	  elsif(   $downstream_site eq 'AC' ){
	      my $shift = &find_upstream_in_reverse_strand('AT',$exons[$i],$slice);
	      if ( $shift != 0 ){ 
		  my $start = $exons[$i]->start;
		  $exons[$i]->start( $start + $shift ); 
		  $one_corrected++;
		  print STDERR "$exons[$i] corrected to AT, shift = $shift\n" if $verbose;
	      }
	  }
	  elsif( $upstream_site eq 'GC' ){
	      my $shift = &find_downstream_in_reverse_strand('AG',$exons[$i+1],$slice);
	      if ( $shift != 0 ){
		  my $end = $exons[$i+1]->end;
		  $exons[$i+1]->end( $end + $shift );
		  $one_corrected++;
		  print STDERR "$exons[$i+1] corrected to AG, shift = $shift\n" if $verbose;
	      }
	  }
	  elsif( $downstream_site eq 'AG' ){
	      my $shift = &find_upstream_in_reverse_strand('GC',$exons[$i],$slice);
	      if ( $shift != 0 ){
		  my $start = $exons[$i]->start;
		  $exons[$i]->start( $start + $shift ); 
		  $one_corrected++;
		  print STDERR "$exons[$i] corrected to GC, shift = $shift\n" if $verbose;
	      }
	  }
	  else {
	      my $shift_GT =  &find_upstream_in_reverse_strand('GT',$exons[$i],$slice);
	      my $shift_AG =  &find_downstream_in_reverse_strand('AG',$exons[$i+1],$slice);
	      my $shift_AT =  &find_upstream_in_reverse_strand('AT',$exons[$i],$slice);
	      my $shift_AC =  &find_downstream_in_reverse_strand('AC',$exons[$i+1],$slice);
	      my $shift_GC =  &find_upstream_in_reverse_strand('GC',$exons[$i],$slice);
	      if ( $shift_GT != 0 && $shift_AG != 0 ){
		  my $end = $exons[$i+1]->end;
		  $exons[$i+1]->end( $end + $shift_AG );
		  my $start = $exons[$i]->start;
		  $exons[$i]->start( $start + $shift_GT ); 
		  $both_corrected++;
		  print STDERR "modified [$i] start: $start -> ". $exons[$i]->start.
		      " [$i+1] end: $end -> ".$exons[$i+1]->end."\n";
	      }
	      elsif( $shift_AT != 0 && $shift_AC != 0 ){
		  my $end = $exons[$i+1]->end;
		  $exons[$i+1]->end( $end + $shift_AC );
		  my $start = $exons[$i]->start;
		  $exons[$i]->start( $start + $shift_AT ); 
		  $both_corrected++;
	      }
	      elsif( $shift_GC != 0 && $shift_AG != 0 ){
		  my $end = $exons[$i+1]->end;
		  $exons[$i+1]->end( $end + $shift_AG );
		  my $start = $exons[$i]->start;
		  $exons[$i]->start( $start + $shift_GC ); 
		  $both_corrected++;
	      }
	  }
      } # end of INTRON
  }
  
  $transcript->flush_Exons;
  foreach my $exon ( @exons ){
      if ( $exon->stable_id ){
	  $exon->version(1);
	  $exon->created(0);
	  $exon->modified(0);
      }
      ############################################################
      # re-set the sequence or delete the current on
      # so that it gets re-calculated next time it is called
      delete $exon->{'_seq_cache'};
      #my $seq_string = $exon->contig->subseq( $exon->start, $exon->end, $exon->strand );
      #my $exon_seq = Bio::Seq->new(
      #-DISPLAY_ID => $exon->stable_id || $exon->dbID,
      #-MOLTYPE    => 'dna',
      #-SEQ        => $seq_string,
      #);
      #$exon->seq($exon_seq);
      $transcript->add_Exon($exon);
  }
  
  
  return ( $transcript,$introns, $correct, $one_corrected, $both_corrected );
}


############################################################

sub find_downstream_in_forward_strand{
    my ($site,$downstream_exon,$slice) = @_;
    my @shifts = (-1,+1,-2,+2,-3,+3);
    foreach my $shift (@shifts){
	my $downstream_site = 
	    $slice
		->subseq( ($downstream_exon->start - 2 + $shift), ($downstream_exon->start - 1 + $shift) );
	print STDERR "downstream:$downstream_site site:$site shift:$shift\n";
	if ( $downstream_site eq $site ){
	    print STDERR "Bingo!\n";
	    return $shift;
	}
    }
    return 0;
}

############################################################

sub shift_downstream_translation{
    my $value = shift;
    return 2 if $value == 1 || $value == -2;
    return 1 if $value == 2 || $value == -1;
    return 0;
}

sub shift_upstream_translation{
    my $value = shift;
    return 1 if $value == 1 || $value == -2;
    return 2 if $value == 2 || $value == -1;
    return 0;
}


############################################################

sub find_upstream_in_forward_strand{
    my ($site,$upstream_exon,$slice) = @_;
    my @shifts = (-1,+1,-2,+2);
    foreach my $shift (@shifts){
	my $upstream_site = 
	    $slice
		->subseq( ($upstream_exon->end     + 1 + $shift), ($upstream_exon->end     + 2 + $shift) );
	print STDERR "upstream:$upstream_site site:$site shift:$shift\n";
	if ( $upstream_site eq $site ){
	    print STDERR "Bingo!\n";
	    return $shift;
	}
    }
    return 0;
}

############################################################

sub find_downstream_in_reverse_strand{
  my ($site,$downstream_exon,$slice) = @_;
  my @shifts = (-1,+1,-2,+2,-3,+3);
  foreach my $shift (@shifts){
    my $downstream_site;
    my $down_site = 
	$slice->subseq( ($downstream_exon->end + 1 + $shift), ($downstream_exon->end + 2 + $shift) );
    ( $downstream_site = reverse( $down_site ) ) =~ tr/ACGTacgt/TGCAtgca/;
    print STDERR "downstream:$downstream_site site:$site shift:$shift\n";
    if ( $downstream_site eq $site ){
	print STDERR "Bingo\n";
	return $shift;
    }
  }
  return 0;
}

############################################################

sub find_upstream_in_reverse_strand{
  my ($site,$upstream_exon,$slice) = @_;
  my @shifts = (-1,+1,-2,+2,-3,+3);
  foreach my $shift (@shifts){
      my $upstream_site;
      my $up_site = 
	  $slice->subseq( ($upstream_exon->start - 2 + $shift), ($upstream_exon->start - 1 + $shift) );
      ( $upstream_site   = reverse(  $up_site  ) ) =~ tr/ACGTacgt/TGCAtgca/;
      print STDERR "upstream:$upstream_site site:$site shift:$shift\n";
      if ( $upstream_site eq $site ){
	  print STDERR "Bingo\n";
	  return $shift;
      }
  }
  return 0;
}

############################################################


############################################################


sub fix_translation{
    my $trans = shift;

    my $verbose = 0;

    my @met_predictions   = &run_translate( $trans,1);
    my @nomet_predictions = &run_translate( $trans );
    
    my $count = 0;
    while ( $count < 2 && $met_predictions[$count] ){
	my @entry = @{$met_predictions[$count]};
	print STDERR "MET length:$entry[0] start:$entry[1] end:$entry[2]\n";
	$count++;
    }
    $count = 0;
    while ( $count < 2 && $nomet_predictions[$count] ){
	my @entry = @{$nomet_predictions[$count]};
	print STDERR "NO_MET length:$entry[0] start:$entry[1] end:$entry[2]\n";
	$count++;
    }
    my $translation = $trans->translation;
    my $length = $trans->seq->length;
    my $best;
    if ( @met_predictions && @nomet_predictions ){
	my $met_best   = $met_predictions[0];
	my $nomet_best = $nomet_predictions[0];
	if ( $nomet_best->[0] > 2*$met_best->[0] ){
	    $best = $nomet_best;
	}
	else{
	    $best = $met_best;
	}
    }
    elsif( @met_predictions ){
	$best = $met_predictions[0];
    }
    elsif( @nomet_predictions ){
	$best = $nomet_predictions[0];
    }
    my @entry = @{$best};
    my $orf_start = $entry[1];
    my $orf_end   = $entry[2];
    print STDERR "BEST length:$entry[0] start:$entry[1] end:$entry[2]\n";
    my @exons;
    my $strand = $trans->start_Exon->strand;
    if ( $strand == 1 ){
	@exons = sort { $a->start <=> $b->start } @{$trans->get_all_Exons};
    }
    else{
	@exons = sort { $b->start <=> $a->start } @{$trans->get_all_Exons};
    }
    my $transl_start;
    my $transl_end;
    my $transl_start_Exon;
    my $transl_end_Exon;
    my $exon_count = 0;
    print STDERR "transcript length: $length\n";
    my $pos = 1;
    foreach my $exon ( @exons ){
	$exon_count++;
	print STDERR "exon:$exon_count exon_length:".$exon->length." pos:$pos orf_start:$orf_start orf_end:$orf_end pos+:".($pos + $exon->length - 1)."\n" if $verbose;
	if ( $orf_start >= $pos && $orf_start <= $pos + $exon->length - 1 ){
	    $transl_start_Exon = $exon;
	    $transl_start      = $orf_start - $pos + 1;
	    print STDERR "start found\n" if $verbose;
	    $exon->phase(0);
	}
	if ( $orf_end >= $pos && $orf_end <= $pos + $exon->length - 1 ){
	    $transl_end_Exon   = $exon;
	    $transl_end        = $orf_end - $pos + 1;
	    print STDERR "end found\n" if $verbose;
	}
	$pos += $exon->length;
	#last if ( $pos > $length );
    }
    print STDERR "original translation:\n";
  Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Translation($trans);

    my $newtranslation;
    if ( $transl_start && $transl_end &&  $transl_start_Exon && $transl_end_Exon ){
	$newtranslation = Bio::EnsEMBL::Translation->new();
	$newtranslation->start( $transl_start );
	$newtranslation->end( $transl_end );
	$newtranslation->start_Exon( $transl_start_Exon );
	$newtranslation->end_Exon( $transl_end_Exon );
	$trans->translation($newtranslation);
	print STDERR "modified translation:\n";
      Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Translation($trans);
    }
    else{
	print STDERR "problem modifying the translation\n";
    }
    $trans->flush_Exons;
    foreach my $exon ( @exons ){
	$trans->add_Exon($exon);
    }
    return $trans;
}

############################################################

sub run_translate{
    my ($trans,$met) = @_;
    
    my $verbose = 0;

    my $trans_id = $trans->stable_id || $trans->dbID;
    my $seq = $trans->seq;
    unless ( $seq->display_id ){
	$seq->display_id( $trans_id );
    }
    my $length = $seq->length;
    
    ############################################################
    # create file
    my $file = "/tmp/"."cdna_".$$.".fa";
    open ( SEQ, ">$file" ) || die("could not open file $!");
    my $seqout = Bio::SeqIO->new('-format' => 'Fasta',
				 '-fh'     => \*SEQ);
    $seqout->write_seq($seq);
    close(SEQ);
    
    my $command;
    if ( $met){
	$command = "/usr/local/ensembl/bin/translate -m $file |";
    }
    else{
	$command = "/usr/local/ensembl/bin/translate $file |";
    } 
    open ( ORF, $command ) || die( "Error running translate" );
    ############################################################
    # out put is of the form:
    #> gi|20070124|ref|NM_000918.2|.44    length 62, nt 2236..2051
    #AHDRRRSPGLREGEGPGLCRAPGLAATSSSSRHGGHPDRIRKSPFTQKCKSHDQSWRHCRRY
    #> gi|20070124|ref|NM_000918.2|.45    length 34, nt 2047..1946
    #VTMSSPAPSLPHGGQASPRRPGQGGTNTLMSKNV
    #> gi|20070124|ref|NM_000918.2|.46    length 34, nt 1942..1841
    #KSHRRNFQKEEKPPAGGRQRDSEHGSKHSGQTHV
    
    my @orf_predictions;
  ORF:
    while ( <ORF> ){
	chomp;
	next ORF unless /\>/;
	print STDERR "$_\n" if $verbose;
	my @entries = split;
	my $id = $entries[1];
	my $orf_length = $entries[3];
	$orf_length =~s/\,//;
	$entries[5] =~/(\d+)\.\.(\d+)/;
	my $orf_start = $1;
	my $orf_end   = $2;
	next ORF if $orf_start>=$orf_end;
	print STDERR "id:$id\torf_length:$orf_length\tstart:$orf_start\tend:$orf_end\n" if $verbose;
	my @prediction = ($orf_length,$orf_start,$orf_end);
	push( @orf_predictions, \@prediction );
    }
    my @sorted_predictions = 
	map { $_->[1] } sort { $b->[0] <=> $a->[0] } map { [$_->[0], $_] } @orf_predictions;
    return @sorted_predictions;
}


############################################################



sub prune_Exons {
  my ($gene) = @_;
  my @unique_Exons; 
  
  # keep track of all unique exons found so far to avoid making duplicates
  # need to be very careful about translation->start_exon and translation->end_exon
  foreach my $tran (@{$gene->get_all_Transcripts}) {
    my @newexons;
    foreach my $exon (@{$tran->get_all_Exons}) {
      my $found;
    UNI:
      foreach my $uni (@unique_Exons) {
	if ($uni->start     == $exon->start  &&
	    $uni->end       == $exon->end    &&
	    $uni->strand    == $exon->strand &&
	    $uni->phase     == $exon->phase  &&
	    $uni->end_phase == $exon->end_phase
	   ) {
	  $found = $uni;
	  last UNI;
	}
      }
      
      if (defined($found)) {
	  push(@newexons,$found);
	  if ($exon == $tran->translation->start_Exon){
	      $tran->translation->start_Exon($found);
	  }
	  if ($exon == $tran->translation->end_Exon){
	      $tran->translation->end_Exon($found);
	  }
      } 
      else {
	  push(@newexons,$exon);
	  push(@unique_Exons, $exon);
      }
    }          
    $tran->flush_Exons;
    foreach my $exon (@newexons) {
      $tran->add_Exon($exon);
    }
  } # end of TRANSCRIPT
  return;
}


############################################################

sub _clone_Gene_only{
  my ($gene) = @_;
  
  my $newgene = new Bio::EnsEMBL::Gene;
  if ($gene->type){
      $newgene->type( $gene->type);
  }
  #if ( defined $gene->dbID ){
  #    $newgene->dbID($gene->dbID);
  #}
  if ( defined $gene->stable_id ){
      $newgene->stable_id( $gene->stable_id );
      $newgene->created($gene->created);
      $newgene->modified($gene->modified);
      $newgene->version($gene->version);
  }
  if ( defined $gene->analysis ){
      $newgene->analysis($gene->analysis);
  }
  return $newgene;
}
