package WormBase;
require Exporter;


our @ISA = qw(Exporter);
our @EXPORT = qw(get_seq_ids get_sequences_pfetch agp_parse parse_gff write_genes translation_check make_Clone make_Contig insert_agp_line display_exons non_translate);

use strict;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;

sub get_seq_ids{
  my ($fh) = @_;

  my @seq_ids;
  my @non_ids;
  while(<$fh>){
   chomp;
   #I	47490	107680	3	F	AC024796.1	1	60191	+
   #print;
   #print "\n";
   my ($status, $contig) =
    (split)[4, 5];
   if($status eq 'N'){
     next;
   }
   if($contig eq '.'){
     next;
   }
   if(!$contig =~ /\S+\.\d+/){
     push(@non_ids, $contig);
     next;
   }
   push(@seq_ids, $contig)
  }
  return \@seq_ids, \@non_ids;
}


sub get_sequences_pfetch{
  my ($seq_ids, $seqfetcher) = @_;
  unless($seqfetcher->isa("Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch")){
    die("seqfetcher ".$seqfetcher." needs to be a pfetch for this too work");
  }
  my %seqs;
  foreach my $id(@$seq_ids){
    my $seq;
    eval{
      $seq = $seqfetcher->get_Seq_by_acc($id);
    };
    if($@){
      warn "$id isn't most recent sequence trying archive\n";
      $seqfetcher->options('-a');
      $seq = $seqfetcher->get_Seq_by_acc($id);
    }
    if($seq){
      $seqs{$id} = $seq;
    }else{
      warn "sequence ".$id." wasn't found\n";
    }
  }
  return(\%seqs);
}


sub agp_parse{
  my ($fh, $chr_id, $agp_type) = @_;
  my $chr_hash = {};
  #print STDERR "parsing agp\n";
  while(<$fh>){
    chomp;
    #I	47490	107680	3	F	AC024796.1	1	60191	+
    #print;
    #print "\n";
    my ($chr, $chr_start, $chr_end, $gap,  $contig, $raw_start, $raw_end, $raw_ori) =
      (split)[0, 1, 2, 4, 5, 6, 7, 8];
    if($gap eq 'N'){
      next;
    }
    if($contig eq '.'){
      next;
    }
    if(!$contig =~ /\S+\.\d+/){
      next;
    }
    if ($raw_ori eq '+') {
      $raw_ori = 1;
    }
    elsif ($raw_ori eq '-') {
      $raw_ori = -1;
    }
    else {
      $raw_ori = 1;
      #print "$chr Contig  $contig  $chr_start \n";
      #print "Warning assumed default orientation for $contig\n";
    }
    #print "have raw coords ".$raw_start." ".$raw_end." ".$raw_ori."\n";
    if($chr_hash->{$contig}){
      die "contig ".$contig." has been found twice something odd is going on\n";
    }else{
      $chr_hash->{$contig} = { 'chromosome_id' => $chr_id,
			       'chr_start' => $chr_start,
			       'chr_end' => $chr_end,
			       'superctg_name' => $chr,
			       'superctg_start' => $chr_start,
			       'superctg_end' => $chr_end,
			       'superctg_ori' => 1,
			       'contig_start' => $raw_start,
			       'contig_end' => $raw_end,
			       'contig_ori' => $raw_ori,
			       'type' => $agp_type,
			     }

    }
  }
  return $chr_hash;
}




sub parse_gff{
  my ($file, $seq, $analysis) = @_;

  print STDERR "opening ".$file."\n";
  open(FH, $file) or die"couldn't open ".$file." $!";

  die " seq ".$seq." is not a Bio::Seq " unless($seq->isa("Bio::SeqI") || 
						$seq->isa("Bio::Seq")  || 
						$seq->isa("Bio::PrimarySeqI"));
  my @genes;
  my $transcripts = undef;
  $transcripts = &process_file(\*FH);
  #print "there are ".keys(%$transcripts)." distinct transcripts\n";
  my $processed_transcripts = undef;
  $processed_transcripts = &process_transcripts($transcripts, $seq, $analysis);
  #print "there are ".keys(%$processed_transcripts)." transcript\n";
  my $genes = undef;
  $genes = &create_transcripts($processed_transcripts);
  #print "PARSE GFF there are ".keys(%$genes)." genes\n";
  foreach my $gene_id(keys(%$genes)){
    my $transcripts = $genes->{$gene_id};
    my $unpruned = &create_gene($transcripts, $gene_id);
    #print STDERR "gene ".$unpruned."\n";
    my $gene = &prune_Exons($unpruned);
    push(@genes, $gene);
  }
  close(FH);
  #print "PARSE_GFF ".@genes." genes\n";
  return \@genes;
}





sub process_file{
  my ($fh) = @_;
  
  my %transcripts;
 LOOP: while(<$fh>){
#    print;
    chomp;
    my($chr, $status, $type, $start, $end, $score, $strand, $frame, $sequence, $gene) = split;
    my $element = $_;
    if($chr =~ /sequence-region/){
      #print STDERR $_;
      next LOOP;
    }
    if(!$status && !$type){
      #print "status and type no defined skipping\n";
      next LOOP;
    }
    my $line = $status." ".$type;
#    print "line ".$line."\n";
    if($line ne 'curated CDS'){
      next LOOP;
    }
    $gene =~ s/\"//g;
    if(!$transcripts{$gene}){
      $transcripts{$gene} = [];
      push(@{$transcripts{$gene}}, $element);
    }else{
      push(@{$transcripts{$gene}}, $element);
    }
    
  }
  return \%transcripts;
}


sub process_transcripts{
  my ($transcripts, $slice, $analysis) = @_;
  
  my %genes;
  my %transcripts = %$transcripts;
  my @names = keys(%transcripts);
  #print STDERR "PROCESSING TRANSCRIPTS \n";
  foreach my $name(@names){
    my @lines = @{$transcripts{$name}};
    $transcripts{$name} = [];
    my @exons;
    foreach my $line(@lines){
      my($chr, $status, $type, $start, $end, $score, $strand, $frame, $sequence, $gene) = split /\s+/, $line;
      $chr =~ s/CHROMOSOME_//;
      if($start == $end){
	next;
      }
     
      my $exon = new Bio::EnsEMBL::Exon;
      my $phase = (3 - $frame)%3;
      #print "have ".$phase." phase\n";
      $exon->start($start);
      $exon->end($end);
      $exon->analysis($analysis);
      $exon->contig($slice);
      $exon->phase($phase);
      my $end_phase = ($phase + ($exon->end-$exon->start) + 1)%3;
      $exon->end_phase($end_phase);
      if($strand eq '+'){
	$exon->strand(1);
      }else{
	$exon->strand(-1);
      }
      $exon->score(100);
      push(@exons, $exon);
    }
    if($exons[0]->strand == -1){
      @exons = sort{$b->start <=> $a->start} @exons;
    }else{
      @exons = sort{$a->start <=> $b->start} @exons;
    }
    #print STDERR $name." transcript\n";
    my $phase = 0;
    foreach my $e(@exons){
      #print "phase ".$phase."\n";
      #$e->phase($phase);
      #my $end_phase = ($phase + ($e->end-$e->start) + 1)%3;
      #$e->end_phase($end_phase);
      #print "end phase ".$end_phase."\n";
      #$phase = $end_phase;
      push(@{$transcripts{$name}}, $e);
    }
  }
  #print STDERR "FINISHED PROCESSING TRANSCRIPTS\n";
  return \%transcripts;

}

sub create_transcripts{
  my ($transcripts) = @_;
  #print STDERR "CREATING TRANSCRIPTS\n";
  my %transcripts = %$transcripts;
  my @non_translate;
  my %genes;
  my $gene_name;
  my $transcript_id;
  foreach my $transcript(keys(%transcripts)){
    my $time = time;
    #print STDERR "transcript ".$transcript." \n" if($transcript =~ /C07F11\.1/);
    my @exons = @{$transcripts{$transcript}};
    if($transcript =~ /\w+\.\d+[a-z A-Z]/){
      #print STDERR "parsing gene name\n";
     ($gene_name) = $transcript =~ /(\w+\.\d+)[a-z A-Z]/;
     $transcript_id = $transcript;
    }else{
      #print STDERR "taking gene name as is\n";
      $gene_name = $transcript;
      $transcript_id = $transcript;
    }
    #print STDERR "have gene name ".$gene_name." and transcript id ".$transcript_id."\n" if($transcript =~ /C07F11\.1/);
    my $transcript = new Bio::EnsEMBL::Transcript;
    my $translation = new Bio::EnsEMBL::Translation;
    my @sorted_exons;
    if($exons[0]->strand == 1){
      @sorted_exons = sort{$a->start <=> $b->start} @exons
    }else{
      @sorted_exons = sort{$b->start <=> $a->start} @exons  
    }
    my $exon_count = 1;
    my $phase = 0;
    foreach my $exon(@sorted_exons){
      $exon->created($time);
      $exon->modified($time);
      $exon->version(1);
      $exon->stable_id($transcript_id.".".$exon_count);
      #$exon->phase($phase);
      #my $end_phase = ($phase + ($end-$start) + 1)%3; 
      #$exon->end_phase($end_phase);
      #$phase = $end_phase;
      $exon_count++;
      $transcript->add_Exon($exon);
    }
    $translation->start_Exon($sorted_exons[0]);
    $translation->end_Exon  ($sorted_exons[$#sorted_exons]);
 
    if ($sorted_exons[0]->phase == 0) {
      $translation->start(1);
    } elsif ($sorted_exons[0]->phase == 1) {
      $translation->start(3);
    } elsif ($sorted_exons[0]->phase == 2) {
      $translation->start(2);
    }
    $translation->end  ($sorted_exons[$#sorted_exons]->end - $sorted_exons[$#sorted_exons]->start + 1);
#    $translation->version(1);
    $translation->stable_id($transcript_id);
    $transcript->translation($translation);
    $transcript->version(1);
    $transcript->stable_id($transcript_id);
   # my $mrna = $transcript->translateable_seq;
#    my $last_codon = substr($mrna, -3);
#    if($last_codon eq 'TAG' || $last_codon eq 'TGA' || $last_codon eq 'TAA'){
#      my $translation_end =   ($sorted_exons[$#sorted_exons]->end - $sorted_exons[$#sorted_exons]->start + 1)-3;
#      if($translation_end <= 0){
#	my $exon_number = $#sorted_exons
#      }
#      $translation->end  (($sorted_exons[$#sorted_exons]->end - $sorted_exons[$#sorted_exons]->start + 1)-3); 
#    }
    my $peptide = $transcript->translate->seq;
    
    if(!$genes{$gene_name}){
      $genes{$gene_name} = [];
      push(@{$genes{$gene_name}}, $transcript);
    }else{
      push(@{$genes{$gene_name}}, $transcript);
    }
  }

  
  #print STDERR "FINISHED CREATEING TRANSCRIPTS\n";
  return \%genes;

}


sub create_gene{
  my ($transcripts, $name) = @_;
  #print STDERR "CREATING GENE\n";
  my $time = time;
  my $gene = new Bio::EnsEMBL::Gene; 
  my $exons = $transcripts->[0]->get_all_Exons;
  my $analysis = $exons->[0]->analysis;
  $gene->analysis($analysis);
  $gene->type($analysis->logic_name);
  $gene->created($time);
  $gene->modified($time);
  $gene->version(1);
  $gene->stable_id($name);
  foreach my $transcript(@$transcripts){
    $gene->add_Transcript($transcript);
  }
  
  return $gene;
}


sub prune_Exons {
  my ($gene) = @_;
  
  my @unique_Exons; 
  
  # keep track of all unique exons found so far to avoid making duplicates
  # need to be very careful about translation->start_Exon and translation->end_Exon
  
  foreach my $tran (@{$gene->get_all_Transcripts}) {
    my @newexons;
    foreach my $exon (@{$tran->get_all_Exons}) {
      my $found;
      #always empty
    UNI:foreach my $uni (@unique_Exons) {
	if ($uni->start  == $exon->start  &&
	    $uni->end    == $exon->end    &&
	    $uni->strand == $exon->strand &&
	    $uni->phase  == $exon->phase  &&
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
      } else {
	push(@newexons,$exon);
	push(@unique_Exons, $exon);
      }
    }          
    $tran->flush_Exons;
    foreach my $exon (@newexons) {
      $tran->add_Exon($exon);
    }
  }
  return $gene;
}



sub write_genes{
  my ($genes, $db) = @_;

  #print "WRITE GENES there are ".@$genes." genes\n";
  my %non_translating;
  my %non_transforming;
  
 GENE: foreach my $gene(@$genes){
    eval{
      $gene->transform;
    };
    if($@){
      warn("gene ".$gene->stable_id." wouldn't transform ".$@);
      my ($clone_name) = $gene->stable_id =~ /(\S+)\.\S+/;
      if(!$non_transforming{$clone_name}){
	$non_transforming{$clone_name} = [];
	push(@{$non_transforming{$clone_name}}, $gene);
	next GENE;
      }else{
	push(@{$non_transforming{$clone_name}}, $gene);
	next GENE;
      }
    }
 
  
    my $gene_adaptor = $db->get_GeneAdaptor;
    eval{
      $gene_adaptor->store($gene);
    };
    if($@){
      die "couldn't store ".$gene->stable_id." problems ".$@;
    }
  }

  return \%non_transforming;
}

sub translation_check{
  my ($gene) = @_;
  
  #print STDERR "checking translation of ".$gene->stable_id."\n";
  my @transcripts = @{$gene->get_all_Transcripts};
  foreach my $t(@transcripts){
    my $pep = $t->translate->seq;
    if($pep =~ /\*/){
      print STDERR "transcript ".$t->stable_id." doesn't translate\n";
      return undef;
    }
  }
  return $gene;
  
}

sub make_Contig{
  my ($name, $seq, $length) = @_;

  my $contig = Bio::EnsEMBL::RawContig->new();

  $contig->name($name);
  $contig->seq($seq);
  $contig->length($length);

  return $contig;
  
}

sub make_Clone{
  my ($name, $version, $embl_acc, $embl_version, $htg_phase, $contig, $created, $modified) = @_;

  my $clone = Bio::EnsEMBL::Clone->new();
  $clone->id($name);
  $clone->version($version);
  $clone->embl_id($embl_acc);
  $clone->embl_version($embl_version);
  $clone->htg_phase($htg_phase);
  $clone->add_Contig($contig);
  $clone->created($created);
  $clone->modified($modified);

  return $clone;
  
}

sub insert_agp_line{
  my ($chr_id, $chr_start, $chr_end, $superctg_name, $superctg_start, $superctg_end, $superctg_ori, $contig, $contig_start, $contig_end, $contig_ori, $type, $db) = @_;

  my $sql = "insert into assembly(chromosome_id, chr_start, chr_end, superctg_name, superctg_start, superctg_end, superctg_ori, contig_id, contig_start, contig_end, contig_ori, type) values(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)";
  my $sth = $db->prepare($sql);
  $sth->execute($chr_id, $chr_start, $chr_end, $superctg_name, $superctg_start, $superctg_end, $superctg_ori, $contig, $contig_start, $contig_end, $contig_ori, $type); 
}


sub display_exons{
  my (@exons) = @_;

  @exons = sort{$a->start <=> $b->start || $a->end <=> $b->end} @exons if($exons[0]->strand == 1);

  @exons = sort{$b->start <=> $a->start || $b->end <=> $a->end} @exons if($exons[0]->strand == -1);
  
  foreach my $e(@exons){
       print $e->stable_id."\t ".$e->start."\t ".$e->end."\t ".$e->strand."\t ".$e->phase."\t ".$e->end_phase."\n";
    }
  
}

sub non_translate{
  my (@transcripts) = @_;
  
  foreach my $t(@transcripts){
    
    my @exons = @{$t->get_all_Exons};
#    print "transcript sequence :\n".$t->seq."\n";
    foreach my $e(@exons){
      my $seq = $e->seq;
      my $pep0 = $seq->translate('*', 'X', 0);
      my $pep1 = $seq->translate('*', 'X', 1);
      my $pep2 = $seq->translate('*', 'X', 2);
      print "exon sequence :\n".$e->seq->seq."\n\n";
      print $e->seqname." ".$e->start." : ".$e->end." translation in 0 frame\n ".$pep0->seq."\n\n";
      print $e->seqname." ".$e->start." : ".$e->end." translation in 1 phase\n ".$pep2->seq."\n\n";
      print $e->seqname." ".$e->start." : ".$e->end." translation in 2 phase\n ".$pep1->seq."\n\n";
      print "\n\n";
      
    }
    
  }
}



1;
