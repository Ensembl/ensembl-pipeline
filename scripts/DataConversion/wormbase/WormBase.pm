#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code
# mod fsk

=pod 

=head1 NAME

WormBase

=head1 SYNOPSIS

the methods are all used by wormbase_to_ensembl.pl script which should be in the same directory

=head1 DESCRIPTION

parse gff and agp files from the wormbase database for caenorhabditis elegans to create ensembl database.

=head1 CONTACT

ensembl-dev@ebi.ac.uk about code issues
wormbase-hel@wormbase.org about data issues

=head1 APPENDIX

=cut

use lib '~/PerlCode/ensembl-pipeline/scripts/DataConversion/wormbase';

package WormBase;
require Exporter;


our @ISA = qw(Exporter);
our @EXPORT = qw(get_seq_ids get_sequences_pfetch agp_parse parse_gff write_genes translation_check insert_agp_line display_exons non_translate process_file parse_operons write_simple_features parse_rnai parse_expr parse_SL1 parse_SL2 parse_pseudo_gff store_coord_system store_slice parse_tRNA parse_rRNA_genes parse_tRNA_genes parse_pseudo_files);

use strict;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::SimpleFeature;
use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::Slice;

$|=1; #turns off buffering on STDOUT

=head2 get_seq_ids

  Arg [1]   : filehandle to an agp file
  Function  : retrives sequence ids from an agp file
  Returntype: array ref of seq ids and array ref of ides which don't fit the format'
  Exceptions: non
  Caller    : 
  Example   : &get_seq_ids($fh);

=cut

sub get_seq_ids{
  my ($fh) = @_;

  my @seq_ids;
  my @non_ids;
  while(<$fh>){
   chomp;
   #I	47490	107680	3	F	AC024796.1	1	60191	+
   ##print;
   #print "\n";
   my ($status, $contig) =
    (split)[4, 5];
   if(!$contig =~ /\S+\.\d+/){
     #print STDERR "contig doesn't match ".$contig." accepted format\n";
     push(@non_ids, $contig);
     next;
   }
   if($contig =~ /\S+\.$/){
     #print STDERR "contig doesn't match ".$contig." accepted format\n";
     push(@non_ids, $contig);
     next;
   }
   #print STDERR "contig being added ".$contig."\n";
   push(@seq_ids, $contig)
  }
  return \@seq_ids, \@non_ids;
}


=head2 get_sequences_pfetch

  Arg [1]   : array ref of sequence ids which will be recognised by pfetch
  Arg [2]   : a Bio::EnsEMBL::Pipeline::Seqfetcher::Pfetch object
  Function  : gets sequences for the ids passed to it using pfetch
  Returntype: hash keyed on seq id containing Bio::Seq object
  Exceptions: throws if seqfetcher passed to it isn't pfetch'
  Caller    : 
  Example   : %get_sequences_pfetch($seq_ids, $seqfetcher);

=cut

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


=head2 agp_parse

  Arg [1]   : filehandle to an agp file
  Arg [2]   : chromosome id from chromsome table
  Arg [3]   : agp type for assembly table
  Function  : parses an agp file into a suitable format for inserting into the assembly table
  Returntype: hash ref, keyed on contig id, each hash element is itself a hash which contains all the info needed for the assembly table
  Exceptions: dies if the same contig id is found twice
  Caller    : 
  Example   : 

=cut

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
    if(!$contig =~ /\S+\.\d+/){
      next;
    }
    if($contig =~ /\S+\.$/){
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



=head2 parse_gff

  Arg [1]   : filename of gff file
  Arg [2]   : Bio::Seq object
  Arg [3]   : Bio::EnsEMBL::Analysis object
  Function  : parses gff file given into genes
  Returntype: array ref of Bio::EnEMBL::Genes
  Exceptions: dies if can't open file or seq isn't a Bio::Seq 
  Caller    : 
  Example   : 

=cut

sub parse_gff{
  my ($file, $seq, $analysis) = @_;

  use Storable qw(store retrieve freeze thaw dclone);
  
  #print STDERR "opening ".$file."\n";
  open(FH, $file) or die "couldn't open ".$file." $!";
  
  die " seq ".$seq." is not a Bio::Seq " unless($seq->isa("Bio::SeqI") || 
						$seq->isa("Bio::Seq")  || 
						$seq->isa("Bio::PrimarySeqI"));
  my @genes;

  my ($transcripts, $five_prime, $three_prime) = &process_file(\*FH);
  print "there are ".keys(%$transcripts)." distinct transcripts\n";
  my $genes = undef;

  my ($processed_transcripts, $five_start, $three_end, $trans_start_exon, $trans_end_exon) = &generate_transcripts($transcripts, $seq, $analysis, $five_prime, $three_prime);

  print "\nthere are ".keys(%$processed_transcripts)." transcripts\n";
  #print keys(%$five_start)." transcripts have 5' UTRs and ".keys(%$three_end)." have 3' UTRs\n";

  $genes = undef;
  $genes = &create_transcripts($processed_transcripts, $five_start, $three_end, $trans_start_exon, $trans_end_exon);

  print "\nPARSE GFF has ".keys(%$genes)." genes\n";
  foreach my $gene_id(keys(%$genes)){
    my $transcripts = $genes->{$gene_id};
    my $unpruned = &create_gene($transcripts, $gene_id);
    #print STDERR "gene ".$unpruned."\n";
    my $gene = &prune_Exons($unpruned);
    push(@genes, $gene);
  }

  close(FH);
  print "\nPARSE_GFF got ".@genes." genes\n";
  return \@genes;
}


=head2 process_file

  Arg [1]   : filehandle pointing to a gff file
  Function  : parses out lines for exons
  Returntype: hash keyed on transcript id each containig array of lines for that transcript
              and two hashes of hashes with arrays containing 5 prime and 3 prime UTR lines
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub process_file{
  my ($fh) = @_;

  my %genes;
  my $transcript;
  my %five_prime;
  my %three_prime;
 LOOP: while(<$fh>){
    #CHROMOSOME_I    curated three_prime_UTR 11696828        11697110        .       -       .       Transcript "T22H2.5a"
    #CHROMOSOME_I    curated three_prime_UTR 11697157        11697230        .       -       .       Transcript "T22H2.5a"
    #CHROMOSOME_I    curated five_prime_UTR  11698944        11698954        .       -       .       Transcript "T22H2.5a"

    chomp;
    my($chr, $status, $type, $start, $end, $score, $strand, $frame, $sequence, $gene) = split;
    my $element = $_;
    #print STDERR $element."\n" if($type eq 'exon');
    if($chr =~ /sequence-region/){
      #print STDERR $_;
      next LOOP;
    }
    if(/^##/){
      next LOOP;
    }
    if(!$status && !$type){
      #print "status and type no defined skipping\n";
      next LOOP;
    }
    my $line = $status." ".$type;
    if( ($line eq 'Coding_transcript five_prime_UTR') or ($line eq 'Coding_transcript three_prime_UTR') ){
      $gene =~ s/\"//g;
      $transcript = $gene;
      #remove transcript-specific part: Y105E8B.1a.2
      $gene =~ s/(\.\w+)\.\d+$/$1/;
      my ($position) = $type;
      if($position =~/^five/){
	#print STDERR "have 5 prime utr ".$gene." (".$transcript.")\n";
	if(!$five_prime{$gene}){
	  $five_prime{$gene} = {};
	  if(!$five_prime{$gene}{$transcript}){
	    $five_prime{$gene}{$transcript} = [];
	  }
	  push(@{$five_prime{$gene}{$transcript}}, $element);
	}
	else{
	  if(!$five_prime{$gene}{$transcript}){
	    $five_prime{$gene}{$transcript} = [];
	  }
	  push(@{$five_prime{$gene}{$transcript}}, $element);
	}
      }elsif($position =~/^three/){
	#print STDERR "have 3 prime utr ".$gene." (".$transcript.")\n";
	if(!$three_prime{$gene}){
	  $three_prime{$gene} = {};
	  if(!$three_prime{$gene}{$transcript}){
	    $three_prime{$gene}{$transcript} = [];
	  }
	  push(@{$three_prime{$gene}{$transcript}}, $element);
	}
	else{
	  if(!$three_prime{$gene}{$transcript}){
	    $three_prime{$gene}{$transcript} = [];
	  }
	  push(@{$three_prime{$gene}{$transcript}}, $element);
	}
      }
      next LOOP;
    }elsif($line ne 'curated coding_exon'){
      next LOOP;
    }
    $gene =~ s/\"//g;
    if(!$genes{$gene}){
      $genes{$gene} = [];
      push(@{$genes{$gene}}, $element);
    }else{
      push(@{$genes{$gene}}, $element);
    }
  }
  print STDERR "Have ".keys(%genes). " genes (CDS), ".
    keys(%five_prime)." have 5' UTR and ".keys(%three_prime)." have 3' UTR information\n";
  return \%genes, \%five_prime, \%three_prime;
}


=head2 generate_transcripts
      
  Arg [1]   : hash ref (as returned by process_file, containing
              information about the CDS)
  Arg [2]   : Bio::EnsEMBL::Slice
  Arg [3]   : Bio::EnsEMBL::Analysis
  Arg [4]   : ref to hash of hashes (as returned by process_file, containing
              information about the 5' UTR regions
  Arg [5]   : ref to hash of hashes (as returned by process_file, containing
              information about the 3' UTR regions
  Function  : takes line representing a transcript and creates an exon for each one
  Returntype: hash ref hash keyed on transcript id containing an array of exons
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub generate_transcripts{
  my ($genesRef, $slice, $analysis, $five_prime, $three_prime) = @_;
  my %genes;
  my %transcripts;
  my %temp_transcripts;
  my %five_trans_start;
  my %three_trans_end;
  my %trans_start_exon;
  my %trans_end_exon;
  my $translation_end;
  my $genecount = 0;
  my @global_exons;
  my %overlapcheck;

  use Bio::EnsEMBL::Pipeline::Tools::ExonUtils;

  #go through all genes
  GENE: foreach my $gene_name(keys(%$genesRef)){
    #print "\nGENE $gene_name";
    #create gene-hash entry
    $genes{$gene_name} = [];
    my $transcriptcount = 0;
    %temp_transcripts = ();


    ## collect all "curated_coding_exons" for this gene ##
    my @lines = @{$$genesRef{$gene_name}}; #is this right?
    my @global_exons = ();
    my %three_prime_exons;
    my %five_prime_exons;

    foreach my $line(@lines){
      my($chr, $status, $type, $start, $end, $score, $strand, $frame, $sequence, $gene) = split /\s+/, $line;
      $chr =~ s/CHROMOSOME_//;
      # we re currently not using singe nucleotide exons
      if($start == $end){
	next;
      }
      my $exon = new Bio::EnsEMBL::Exon;
      my $phase = (3 - $frame)%3; # wormbase gff cotains frame which is effectively the opposite of phase 
                                  # for a good explaination of phase see the Bio::EnsEMBL::Exon documentation
      #print STDERR "phase calculated to be ".$phase." \n";
      $exon->start($start);
      $exon->end($end);
      $exon->analysis($analysis);
      $exon->slice($slice);
      $exon->phase($phase);
      my $end_phase = ($phase + ($exon->end-$exon->start) + 1)%3;
      #print STDERR "end phase calculated to be ".$end_phase."\n";
      $exon->end_phase($end_phase);
      if($strand eq '+'){
	$exon->strand(1);
      }else{
	$exon->strand(-1);
      }
      #$exon->score(100);
      push(@global_exons, $exon);
    }

    #sort exons for this gene
    if($global_exons[0]->strand == -1){
      @global_exons = sort{$b->start <=> $a->start} @global_exons;
    }else{
      @global_exons = sort{$a->start <=> $b->start} @global_exons;
    }

    #save information if there is not further UTR info
    if(!defined($$five_prime{$gene_name}) and !defined($$three_prime{$gene_name})){
      #print "\nno alternative transcripts...";
      $transcripts{$gene_name} = \@global_exons;
      #print "\nCOUNT: ".keys %transcripts;
      next GENE;
    }


    ## check different transcripts using UTR information ##
    #collect 5' UTRs
    foreach my $transcript_name ( keys(%{$$five_prime{$gene_name}}) ){
      #print "\nchecking 5' transcript $transcript_name. ";
      my @five_prime_exons = ();
      %overlapcheck = ();
      #more than one transcript at 5 prime level
      $temp_transcripts{$transcript_name} = 1;
      #get UTR lines
      foreach my $line(@{$$five_prime{$gene_name}{$transcript_name}}){
	my($chr, $status, $type, $start, $end, $score, $strand, $frame, $sequence, $gene) = split(/\s+/, $line);
	#avoid saving two identical exons
	if(defined $overlapcheck{$start}){
	  #print "\nexon already defined";
	  next;
	}
	$overlapcheck{$start} = 1;
	my $exon = new Bio::EnsEMBL::Exon;
	my $phase = -1;
	$exon->start($start);
	$exon->end($end);
	$exon->analysis($analysis);
	$exon->slice($slice);
	$exon->phase($phase);
	my $end_phase = -1;
	$exon->end_phase($end_phase);
	if($strand eq '+'){
	  $exon->strand(1);
	}else{
	  $exon->strand(-1);
	}
	push(@five_prime_exons, $exon);
      }
      #sort exons for this transcript
      if($five_prime_exons[0]->strand == -1){
	@five_prime_exons = sort{$b->start <=> $a->start} @five_prime_exons;
      }else{
	@five_prime_exons = sort{$a->start <=> $b->start} @five_prime_exons;
      }
      #save them to transcript
      #print "\nhave ".scalar @five_prime_exons." UTR lines. ";
      $five_prime_exons{$transcript_name} = \@five_prime_exons;
    }


    ## collect 3' UTRs ##
    foreach my $transcript_name ( keys(%{$$three_prime{$gene_name}}) ){
      #print "\nchecking 3' transcript $transcript_name. ";
      my @three_prime_exons = ();
      %overlapcheck = ();
      #more than one transcript at the 3 prime level, save the name
      $temp_transcripts{$transcript_name} = 1;
      #get UTR lines
      foreach my $line(@{$$three_prime{$gene_name}{$transcript_name}}){
	my($chr, $status, $type, $start, $end, $score, $strand, $frame, $sequence, $gene) = split /\s+/, $line;
	#avoid saving two identical exons
	if(defined $overlapcheck{$start}){
	  #print "\nexon already defined";
	  next;
	}
	$overlapcheck{$start} = 1;
	my $exon = new Bio::EnsEMBL::Exon;
	my $phase = -1;
	$exon->start($start);
	$exon->end($end);
	$exon->analysis($analysis);
	$exon->slice($slice);
	$exon->phase($phase);
	my $end_phase = -1;
	$exon->end_phase($end_phase);
	if($strand eq '+'){
	  $exon->strand(1);
	}else{
	  $exon->strand(-1);
	}
	push(@three_prime_exons, $exon);
      }
      #sort exons for this transcript
      if($three_prime_exons[0]->strand == -1){
	@three_prime_exons = sort{$b->start <=> $a->start} @three_prime_exons;
      }else{
	@three_prime_exons = sort{$a->start <=> $b->start} @three_prime_exons;
      }
      #save them to transcript
      #print "\nhave ".scalar @three_prime_exons." UTR lines. ";
      $three_prime_exons{$transcript_name} = \@three_prime_exons;
    }

    ## combine exons, 5' and 3' for every transcript ##
    foreach my $transcript_name (keys %temp_transcripts){
      print "transcript $transcript_name\n";
      $transcriptcount++;
      my @exons = ();
      foreach my $temp_exon (@global_exons){
	push(@exons, Bio::EnsEMBL::Pipeline::Tools::ExonUtils->_clone_Exon($temp_exon));
      }
      my $translation_start = 1;
      my $first = 1;
      #set default translation range
      $trans_start_exon{$transcript_name} = 0;
      $trans_end_exon{$transcript_name} = $#exons;
      print "\ntrans-exons: ".$trans_start_exon{$transcript_name}." - ".$trans_end_exon{$transcript_name}." (".scalar @exons.")";

      #check 5' exons
      if(defined($five_prime_exons{$transcript_name})){
	my @five_prime_exons = @{$five_prime_exons{$transcript_name}};
	print "\nworking on 5' of $transcript_name (".scalar @five_prime_exons.") ";
	#is last 5' UTR exon part of first coding exon?

	FIVEUTR: while(my $five_prime_exon = shift(@five_prime_exons)){
	  
	  my $start = $five_prime_exon->start;
	  my $end = $five_prime_exon->end;
	  my $strand = $five_prime_exon->strand;
	  
	  #print "\n- 5'exon: $start - $end.";
	  if($exons[$trans_start_exon{$transcript_name}]->strand == 1 and $strand == 1){
	    #forward strand
	    if($start > $end){
	      #print "\n>>strange 5' UTR (+) exon: ".$end." - ".$start;
	      next FIVEUTR;
	    }
	    if($end == ($exons[$trans_start_exon{$transcript_name}]->start)-1){
	      #combine exons, adjust translation start
	      #print "\ncombine exons, adjust translation start...";
	      $translation_start = $exons[$trans_start_exon{$transcript_name}]->start - $start + 1;
	      #print $translation_start."(+) ";
	      $five_trans_start{$transcript_name} = $translation_start;
	      $exons[$trans_start_exon{$transcript_name}]->start($start);
	    }
	    elsif($end < $exons[$trans_start_exon{$transcript_name}]->start -1){
	      #additional non-coding exon
	      #add to exon array, keep translation start on last coding exon
	      $trans_start_exon{$transcript_name}++;
	      $trans_end_exon{$transcript_name}++;
	      unshift(@exons, $five_prime_exon);
	      print "\nadditional non-coding exon (+) ".$start." - ".$end;
	    }
	    else{
	      print STDERR "\n>>$transcript_name strange 5' UTR exon (+): $start - $end with 1.exons at ".
		$exons[$trans_start_exon{$transcript_name}]->start.
		" - ".$exons[$trans_start_exon{$transcript_name}]->end;
	      next FIVEUTR;
	    }
	  }
	  elsif($exons[$trans_start_exon{$transcript_name}]->strand == -1 and $strand == -1){
	    #reverse strand
	    if($start > $end){
	      #print "\n>>strange 5' UTR (-) exon: ".$end." - ".$start;
	      next FIVEUTR;
	    }
	    if($start == ($exons[$trans_start_exon{$transcript_name}]->end)+1){
	      #combine exons, adjust translation start
	      #print "\ncombine exons, adjust translation start...";
	      $translation_start = ($end - $exons[$trans_start_exon{$transcript_name}]->end + 1);
	      #print $translation_start."(-) ";
	      $five_trans_start{$transcript_name} = $translation_start;
	      $exons[$trans_start_exon{$transcript_name}]->end($end);
	    }
	    elsif($start > ($exons[$trans_start_exon{$transcript_name}]->end)+1){
	      #additional non-coding exon
	      #add to exon array, keep translation start on last coding exon
	      $trans_start_exon{$transcript_name}++;
	      $trans_end_exon{$transcript_name}++;
	      unshift(@exons, $five_prime_exon);
	      #print "\nadditional 5' non-coding exon (-)".$start." - ".$end;
	    }
	    else{
	      print "\n>>$transcript_name strange 5' UTR exon (-): $start - $end with 1.exons at ".
		$exons[$trans_start_exon{$transcript_name}]->start.
		" - ".$exons[$trans_start_exon{$transcript_name}]->end;
	      next FIVEUTR;
	    }
	  }
	  else{
	    print STDERR "\n>>strand switch in UTR / coding!";
	  }
	}
      }

      #check 3' exons
      if(defined($three_prime_exons{$transcript_name})){
	my @three_prime_exons = @{$three_prime_exons{$transcript_name}};
	#print "\nworking on 3' of $transcript_name (".scalar @three_prime_exons.") ";
	#is first 3' UTR exon part of last coding exon?

	THREEUTR: while(my $three_prime_exon = shift(@three_prime_exons)){
	  
	  my $start = $three_prime_exon->start;
	  my $end = $three_prime_exon->end;
	  my $strand = $three_prime_exon->strand;

	  #print "\n- 3'exon: $start - $end.";
	  if($exons[$trans_end_exon{$transcript_name}]->strand == 1 and $strand == 1){
	    #forward strand
	    if($start > $end){
	      print STDERR "\n>>$transcript_name strange 3' UTR (+) exon: ".$start." - ".$end;
	      next THREEUTR;
	    }
	    if($start == (($exons[$trans_end_exon{$transcript_name}]->end)+1)){
	      #combine exons, adjust translation start
	      #print "\ncombine exons, keep current translation end...";
	      $translation_end = (($exons[$trans_end_exon{$transcript_name}]->end - $exons[$trans_end_exon{$transcript_name}]->start) + 1);
	      #print $translation_end."(+) ";
	      $three_trans_end{$transcript_name} = $translation_end;
	      $exons[$trans_end_exon{$transcript_name}]->end($end);
	    }
	    elsif($start > (($exons[$trans_end_exon{$transcript_name}]->end)+1)){
	      #additional non-coding exon
	      #add to exon array
	      push(@exons, $three_prime_exon);
	      #print "\nadditional 3'  non-coding exon (+)";
	    }
	    else{
	      print STDERR "\n$transcript_name strange 3' UTR exon (+): $start - $end with 1.exons at ".$exons[$trans_end_exon{$transcript_name}]->start;
	      next THREEUTR;
	    }
	  }
	  elsif($exons[$trans_end_exon{$transcript_name}]->strand == -1 and $strand == -1){
	    #reverse strand
	    if($start > $end){
	      print STDERR "\n>>$transcript_name strange 3' UTR (-) exon: ".$start." - ".$end;
	      next THREEUTR;
	    }
	    if($end == (($exons[$trans_end_exon{$transcript_name}]->start)-1)){
	      #combine exons, keep translation start
	      #print "\ncombine exons, keep current translation end....";
	      $translation_end = (($exons[$trans_end_exon{$transcript_name}]->end - $exons[$trans_end_exon{$transcript_name}]->start) +1);
	      #print $translation_end."(-) ";
	      $three_trans_end{$transcript_name} = $translation_end;
	      $exons[$trans_end_exon{$transcript_name}]->start($start);
	    }
	    elsif($end < (($exons[$trans_end_exon{$transcript_name}]->start)-1)){
	      #additional non-coding exon
	      #add to exon array
	      push(@exons, $three_prime_exon);
	      #print "\nadditional 3' non-coding exon (-)".$start." - ".$end;
	    }
	    else{
	      print STDERR "\n$transcript_name strange 3' UTR exon (-): $start - $end with 1.exons at ".$exons[$trans_end_exon{$transcript_name}]->start;
	      next THREEUTR;
	    }
	  }
	}	
	#print "\ntrans-exons: ".$trans_start_exon{$transcript_name}." - ".$trans_end_exon{$transcript_name}." (".scalar @exons.")";
      }
      #add exon data to transcript
      $transcripts{$transcript_name} = \@exons;

    }
    #print STDERR "\nCOUNT: ".keys %transcripts;
  }
#  my $c=0;
#  foreach my $tt (keys %transcripts){
#    print "\ntranscript $tt: ";
#    foreach my $ex (@{$transcripts{$tt}}){
#      print "..".$ex->start." -> ".$ex->end." (".$ex->strand."), ";
#    }
#  }

  return (\%transcripts, \%five_trans_start, \%three_trans_end, \%trans_start_exon, \%trans_end_exon);
}



=head2 create_transcripts

  Arg [1]   : hash ref from process transcripts
  Function  : creates actually transcript objects from the arrays of exons
  Returntype: hash ref keyed on gene id containg an array of transcripts
  Exceptions: 
  Caller    : 
  Example   : 

=cut


sub create_transcripts{
  my ($transcriptsRef, $five_start, $three_end, $trans_start_exon, $trans_end_exon) = @_;

  my @keys = keys(%$five_start);
  foreach my $key(@keys){
    print STDERR "have start of translation for ".$key." ".$five_start->{$key}."\n";
  }

  my %transcripts = %$transcriptsRef;
  my @non_translate;
  my %genes;
  my $gene_name;
  my $transcript_id;
  foreach my $transcript(keys(%transcripts)){
    my $time = time;
    my @exons = @{$transcripts{$transcript}};
    #print STDERR "\nWorking on $transcript.(".$exons[0]->strand.") "; 
    #get the gene-name
    if($transcript =~ /(\w+\.\d+)[a-z A-Z]*/){#if($transcript =~ /\w+\.\d+[a-z A-Z]*/){
     $gene_name = $1;#($gene_name) = $transcript =~ /(\w+\.\d+)[a-z A-Z]*/;
     $transcript_id = $transcript;
    }else{
      $gene_name = $transcript;
      $transcript_id = $transcript;
    }
    print STDERR "\nNote: Gene name= ".$gene_name." Transcript_id= ".$transcript_id." (for transcript ".$transcript.")";
    my $transcript = new Bio::EnsEMBL::Transcript;
    my $translation = new Bio::EnsEMBL::Translation;
    my @sorted_exons;
    if($exons[0]->strand == 1){
      @sorted_exons = sort{$a->start <=> $b->start} @exons;
    }else{
      @sorted_exons = sort{$b->start <=> $a->start} @exons;
    }
    my $exon_count = 1;
    my $phase = 0;
    foreach my $exon (@sorted_exons){
      #$exon->created($time);
      #$exon->modified($time);
      $exon->version(1);
      $exon->stable_id($transcript_id.".".$exon_count);
      $exon_count++;
      eval{
	$transcript->add_Exon($exon);
	#print "adding exon ".$exon->stable_id." \n";
      };
      if($@){ print STDERR "\n>>$transcript_id EXON ERROR: ".$@."\n"; }
    }
    my $start_exon_ind;
    if(exists($trans_start_exon->{$transcript_id})){
      print "\nadjusting coding exons to ".$trans_start_exon->{$transcript_id}." - ".$trans_end_exon->{$transcript_id};
      $translation->start_Exon($sorted_exons[$trans_start_exon->{$transcript_id}]);
      $translation->end_Exon  ($sorted_exons[$trans_end_exon->{$transcript_id}]);
      $start_exon_ind = $trans_start_exon->{$transcript_id};
    }
    else{
      print "not defined\n";
      $translation->start_Exon($sorted_exons[0]);
      $translation->end_Exon  ($sorted_exons[$#sorted_exons]);
      $start_exon_ind = 0;
    }
    print " creating translation for ".$transcript_id."\n";
    if (!exists ($trans_start_exon->{$transcript_id})){
      print "dne: no trans_start_exon for ".$transcript_id."\n";
    }
    
    if (exists($five_start->{$transcript_id})){
      #print "five_start->{transcript_id} is defined\n";
      print "1 setting translation start on transcript ".$transcript_id." to ".$five_start->{$transcript_id}."\n";
      $translation->start($five_start->{$transcript_id});
    } 
    elsif($sorted_exons[$start_exon_ind]->phase == 0) {
      print "case 0\n";
      $translation->start(1);
    }
    elsif ($sorted_exons[$start_exon_ind]->phase == 1) {
    print "case 1\n";
      $translation->start(3);
    }
    elsif ($sorted_exons[$start_exon_ind]->phase == 2) {
    print "case 2\n";
      $translation->start(2);
    }
    else{
      print "dies here ";
      die "Strange phase in $transcript_id ".$sorted_exons[0]->phase;
    }
	#print "done five prime\n";
    if((!defined($translation->start)) or ($translation->start <= 0) ){
      print STDERR ">> no translation start info for ".$transcript_id;
      print STDERR "..".$five_start->{$transcript_id}."\n";
      die();
    }
    
    if(exists($three_end->{$transcript_id})){
      print "2 setting translation end on transcript ".$transcript_id." to ".$three_end->{$transcript_id}." (1)\n";
      $translation->end($three_end->{$transcript_id});
    }else{
      if(defined($trans_end_exon->{$transcript_id})){
	$translation->end($sorted_exons[$trans_end_exon->{$transcript_id}]->end - $sorted_exons[$trans_end_exon->{$transcript_id}]->start +1);
	print "3 setting translation end on transcript ".$transcript_id." to ".
	  ($exons[$trans_end_exon->{$transcript_id}]->end - $exons[$trans_end_exon->{$transcript_id}]->start +1)." (2)\n";
      }
      else{
	$translation->end($sorted_exons[$#sorted_exons]->end - $sorted_exons[$#sorted_exons]->start +1);
	print "4 setting translation end on transcript ".$transcript_id." to ".
	  ($sorted_exons[$#sorted_exons]->end - $sorted_exons[$#sorted_exons]->start +1)." (2)\n";
      }
    }

    $translation->stable_id($transcript_id);
    $translation->version(1);
    $transcript->translation($translation);
    $transcript->version(1);
    $transcript->stable_id($transcript_id);
    if(!$genes{$gene_name}){
      $genes{$gene_name} = [];
      push(@{$genes{$gene_name}}, $transcript);
    }else{
      push(@{$genes{$gene_name}}, $transcript);
    }
    #print "\nstored: $gene_name / $transcript_id";
  }
  return \%genes;

}



=head2 create_gene

  Arg [1]   : array ref of Bio::EnsEMBL::Transcript
  Arg [2]   : name to be used as stable_id
  Function  : take an array of transcripts and create a gene
  Returntype: Bio::EnsEMBL::Gene
  Exceptions: 
  Caller    : 
  Example   : 

=cut


sub create_gene{
  my ($transcripts, $name) = @_;
  my $time = time;
  my $gene = new Bio::EnsEMBL::Gene; 
  my $exons = $transcripts->[0]->get_all_Exons;
  my $analysis = $exons->[0]->analysis;

  $gene->analysis($analysis);
  $gene->biotype($analysis->logic_name);
  $gene->version(1);
  $gene->stable_id($name);
  foreach my $transcript(@$transcripts){
    $gene->add_Transcript($transcript);
  }
  
  return $gene;
}



=head2 prune_Exons

  Arg [1]   : Bio::EnsEMBL::Gene
  Function  : remove duplicate exons between two transcripts
  Returntype: Bio::EnsEMBL::Gene
  Exceptions: 
  Caller    : 
  Example   : 

=cut

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




=head2 write_genes

  Arg [1]   : array ref of Bio::EnsEMBL::Genes
  Arg [2]   : Bio::EnsEMBL::DBSQL::DBAdaptor
  Function  : transforms genes into raw conti coords then writes them to the db provided
  Returntype: hash ref of genes keyed on clone name which wouldn't transform
  Exceptions: dies if a gene could't be stored
  Caller    : 
  Example   : 

=cut

sub write_genes{
  my ($genes, $db, $stable_id_check) = @_;
  my $e=0;
  my %stable_ids;
  if($stable_id_check){
    my $sql = 'select stable_id from gene_stable_id';
    my $sth = $db->dbc->prepare($sql);
    $sth->execute;
    while(my($stable_id) = $sth->fetchrow){
      $stable_ids{$stable_id} = 1;
    }
  }
  my %stored;

 GENE: foreach my $gene(@$genes){
   
    
    #print STDERR "BEFORE STORAGE \n";
    #&display_exons(@{$gene->get_all_Exons});
    if($stable_id_check){
      if($stable_ids{$gene->stable_id}){
        print STDERR $gene->stable_id." already exists\n";
        my $id = $gene->stable_id;
        $id .= '.pseudo';
        $gene->stable_id($id);
        foreach my $transcript(@{$gene->get_all_Transcripts}){
          my $trans_id = $transcript->stable_id;
          $trans_id .= '.pseudo';
          $transcript->stable_id($trans_id);
          foreach my $e(@{$transcript->get_all_Exons}){
            my $id = $e->stable_id;
            $id .= '.pseudo';
            $e->stable_id($id);
          }
        }
      }
    }
    if($stored{$gene->stable_id}){
      print STDERR "we have stored ".$gene->stable_id." already\n";
      next GENE;
    }
    my $gene_adaptor = $db->get_GeneAdaptor;
    eval{
      $stored{$gene->stable_id} = 1;
      $gene_adaptor->store($gene);
      $e++;
    };
    if($@){
      die "couldn't store ".$gene->stable_id." problems ".$@;
    }
  }
  print "\nStored gene: ".$e."\n";
}


=head2 translation_check

  Arg [1]   : Bio::EnsEMBL::Gene
  Function  : checks if the gene translates
  Returntype: Bio::EnsEMBL::Gene if translates undef if doesn't'
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub translation_check{
  my ($gene, $db) = @_;

  my @transcripts = @{$gene->get_all_Transcripts};
  foreach my $t (@transcripts){
    my $pep = $t->translate->seq;
    if($pep =~ /\*/g){
      if($gene->stable_id eq 'C06G3.7' and $db){
	#add Selenocysteine to translation. There seems to be only on Selenoc. in our worm...
	my $pos = pos($pep);
	print STDERR "transcript ".$t->stable_id." doesn't translate. Adding Selenocystein at position $pos.\n".
	  "Please beware of problems during further analysis steps because of this.";
	selenocysteine($t, $pos, $db);
      }
      else{
	print STDERR "transcript ".$t->stable_id." doesn't translate\n";
	print STDERR "translation start ".$t->translation->start." end ".$t->translation->end."\n";
	print STDERR "start exon coords ".$t->translation->start_Exon->start." ".$t->translation->start_Exon->end."\n";
	print STDERR "end exon coords ".$t->translation->end_Exon->start." ".$t->translation->end_Exon->end."\n";
	print STDERR "peptide ".$pep."\n";

	&display_exons(@{$t->get_all_Exons});
	&non_translate($t);
	return undef;
      }
    }
  }
  return $gene;
}

=head2 selenocysteine

  Arg [1]   : transcript object to be modified
  Arg [2]   : position (integer) in transcripts sequence to be modified
  Arg [3]   : Bio::EnsEMBL::DBSQL::DBAdaptor
  Function  : modifiy transcript/translation by adding a seleocysteine residue
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub selenocysteine{
  my ($transcript, $pos, $db) = @_;
  print "\nmodifying ".$transcript->stable_id;

  my $seq_edit = Bio::EnsEMBL::SeqEdit->new(
					    -CODE    => '_selenocysteine',
					    -NAME    => 'Selenocysteine',
					    -DESC    => 'Selenocysteine',
					    -START   => $pos,
					    -END     => $pos,
					    -ALT_SEQ => 'U' #the one-letter symbol for selenocysteine
					   );
  my $attribute = $seq_edit->get_Attribute();
  my $translation = $transcript->translation();
  my $attribute_adaptor = $db->get_AttributeAdaptor();
  $attribute_adaptor->store_on_Translation($translation, [$attribute]);
}


=head2 insert_agp_line

  Arg [1]   : the first 12 args are info for the assembly table
  Arg [2]   : Bio::EnsEMBL::DBSQL::DBAdaptor pointing to db where you want the assembly loaded
  Function  : load the provided info into the assembly table
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub insert_agp_line{
  my ($chr_id, $chr_start, $chr_end, $contig, $contig_start, $contig_end, $contig_ori, $db) = @_;

  if(!$contig){
    #print STDERR "trying to insert into ".$chr_id." ".$chr_start." ".$chr_end."\n";
    die "contig id must be defined for this to work\n";
  }
  my $sql = "insert into assembly(asm_seq_region_id, asm_start, asm_end, cmp_seq_region_id, cmp_start, cmp_end, ori) values(?, ?, ?, ?, ?, ?, ?)";
  
  my $sth = $db->dbc->prepare($sql);
  $sth->execute($chr_id, $chr_start, $chr_end, $contig, $contig_start, $contig_end, $contig_ori); 
}



=head2 display_exons

  Arg [1]   : array of Bio::EnsEMBL::Exons
  Function  : displays the array of exons provided for debug purposes put here for safe keeping
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub display_exons{
  my (@exons) = @_;

  @exons = sort{$a->start <=> $b->start || $a->end <=> $b->end} @exons if($exons[0]->strand == 1);
  
  @exons = sort{$b->start <=> $a->start || $b->end <=> $a->end} @exons if($exons[0]->strand == -1);
  
  foreach my $e(@exons){
       print STDERR $e->stable_id."\t ".$e->start."\t ".$e->end."\t ".$e->strand."\t ".$e->phase."\t ".$e->end_phase."\n";
    }
  
}


=head2 non_translate

  Arg [1]   : array of Bio::EnsEMBL::Transcripts
  Function  : displays the three frame translation of each exon here for safe keeping and debug purposes
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub non_translate{
  my (@transcripts) = @_;

  foreach my $t(@transcripts){

    my @exons = @{$t->get_all_Exons};
#    print "transcript sequence :\n".$t->seq."\n";
    foreach my $e(@exons){
      print "exon ".$e->stable_id." ".$e->start." ".$e->end." ".$e->strand."\n";
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


=head2 parse_operons

  Arg [1]   : gff file path
  Arg [2]   : Bio:Seq object
  Arg [3]   : analysis object
  Function  : parse operon information
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub parse_operons{
  my ($file, $seq, $analysis) = @_;

  #print STDERR "opening ".$file."\n";
  open(FH, $file) or die"couldn't open ".$file." $!";

  die " seq ".$seq." is not a Bio::Seq " unless($seq->isa("Bio::SeqI") || 
						$seq->isa("Bio::Seq")  || 
						$seq->isa("Bio::PrimarySeqI"));
 
  #CHROMOSOME_I	operon	transcription	13758388	13764178	.	-	.	Operon "CEOP1716"

  my @operons;
  LINE: while(<FH>){
      chomp;
      my($type, $start, $end, $strand, $id) = (split)[1,3,4,6,9];
      if($type ne 'operon'){
	next LINE;
      }
      if($strand eq '-'){
	$strand = -1;
      }else{
	$strand = 1;
      }
      my $simple_feature = Bio::EnsEMBL::SimpleFeature->new();
      $simple_feature->start($start);
      $simple_feature->strand($strand);
      $simple_feature->end($end);
      $id =~ s/\"//g;
      $simple_feature->display_label($id);
      $simple_feature->slice($seq);
      $simple_feature->analysis($analysis);
      push(@operons, $simple_feature);
  }

  return \@operons ;
}


sub parse_rnai{
  my ($file, $seq, $analysis) = @_;

  #print STDERR "opening ".$file."\n";
  open(FH, $file) or die "couldn't open ".$file." $!";

  die " seq ".$seq." is not a Bio::Seq " unless($seq->isa("Bio::SeqI") || 
						$seq->isa("Bio::Seq")  || 
						$seq->isa("Bio::PrimarySeqI"));

  print "ONLY FETCHES RNAi_PRIMARY\n";
  my @operons;
  LINE: while(<FH>){
      #CHROMOSOME_IV	RNAi	experimental	3195031	3196094	.	+	.	RNAi "JA:Y67D8A_375.b"
      #CHROMOSOME_I    RNAi_primary    RNAi_reagent    15040052        15041000        .       .       .       Target "RNAi:WBRNAi00027818" 1371 423
      #CHROMOSOME_I    RNAi_secondary  RNAi_reagent    15040723        15040982        .       .       .       Target "RNAi:WBRNAi00027816" 475 216      #print;
      chomp;
      my @values = split;
      if($values[1] ne 'cDNA_for_RNAi'){
	next LINE;
      }
      my ($start, $end, $strand, $id, $count);
      $count = 0;
#      if($values[2] ne 'experimental'){
#	#print "have no experimental tag\n";
#	$start = $values[2];
#	$end = $values[3];
#	$strand = $values[5];
#	$id = $values[8];
#      }else{
	$start = $values[3];
	$end = $values[4];
	$strand = $values[6];
	$id = $values[9];
#      }
#      print $_."\n";
      #print "have ".$start." ".$end." ".$strand." ".$id."\n";
      if($strand eq '+'){
	$strand = 1;
      }else{
	$strand = -1;
      }
      $id =~ s/\"//g;
      my $simple_feature = &create_simple_feature($start, $end, $strand, $id, $seq, $analysis);
      push(@operons, $simple_feature);
  }

  return \@operons ;
}


=head2 parse_operons

  Arg [1]   : gff file path
  Arg [2]   : Bio:Seq object
  Arg [3]   : analysis object
  Function  : parse expression information
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub parse_expr{
  my ($file, $seq, $analysis) = @_;

  #print STDERR "opening ".$file."\n";
  open(FH, $file) or die"couldn't open ".$file." $!";

  die " seq ".$seq." is not a Bio::Seq " unless($seq->isa("Bio::SeqI") || 
						$seq->isa("Bio::Seq")  || 
						$seq->isa("Bio::PrimarySeqI"));

  my @operons;
  LINE: while(<FH>){
      #CHROMOSOME_IV	RNAi	experimental	3195031	3196094	.	+	.	RNAi "JA:Y67D8A_375.b"
      #print;
      chomp;
      my @values = split;
      if($values[1] ne 'Expr_profile'){
	next LINE;
      }
      my ($start, $end, $strand, $id, $count);
      $count = 0;

      $start = $values[3];
      $end = $values[4];
      $strand = $values[6];
      $id = $values[9];

      #print $_."\n";
      #print "have ".$start." ".$end." ".$strand." ".$id."\n";
      if($strand eq '+'){
	$strand = 1;
      }else{
	$strand = -1;
      }
      $id =~ s/\"//g;
      my $simple_feature = &create_simple_feature($start, $end, $strand, $id, $seq, $analysis);
      push(@operons, $simple_feature);
  }

  return \@operons ;
}


sub parse_SL1{
  my ($file, $seq, $analysis) = @_;

  #print STDERR "opening ".$file."\n";
  open(FH, $file) or die"couldn't open ".$file." $!";

  die " seq ".$seq." is not a Bio::Seq " unless($seq->isa("Bio::SeqI") || 
						$seq->isa("Bio::Seq")  || 
						$seq->isa("Bio::PrimarySeqI"));

  my @operons;
  LINE: while(<FH>){
      #CHROMOSOME_I	SL1	trans-splice_acceptor	6473786	6473787	.	-	.	Note "SL1 trans-splice site; see yk719h1.5"
      #print;
      chomp;
      my @values = split;
      if($values[1] ne 'SL1'){
	next LINE;
      }
      if($_ =~ /personal\s+communication/){
	print STDERR "Can't put information in table about ".$_."\n";
	next LINE;
      }
      my ($start, $end, $strand, $id, $count);
      $count = 0;
  
      $start = $values[3];
      $end = $values[4];
      $strand = $values[6];
      $id = $values[14];
      if(!$id){
	$id = $values[13];
      }
      if(($id =~ /cDNA/) || ($id =~ /EST/)){
	$id = $values[15];
      }
      if($id eq '"'){
	$id = $values[13];
      }
      if((!$id) || ($id eq '"')){
	$id = $values[12];
      }
      #print $_."\n";
      #print "have ".$start." ".$end." ".$strand." ".$id."\n";
      if($strand eq '+'){
	$strand = 1;
      }else{
	$strand = -1;
      }
      #print "have id ".$id."\n";
      if((!$id) || ($id eq '.') || ($id =~ /EST/) || ($id =~ /cDNA/) || ($id eq '"')){
	#foreach my $v(@values){
	#  print $count." ".$v."\n";
	#  $count++;
	#}
	print "line ".$_." produced a weird id ".$id."\n";
	next LINE
      }
      $id =~ s/\"//g;
      my $simple_feature = &create_simple_feature($start, $end, $strand, $id, $seq, $analysis);
      push(@operons, $simple_feature);
      }

  return \@operons;
}

sub parse_SL2{
  my ($file, $seq, $analysis) = @_;

  #print STDERR "opening ".$file."\n";
  open(FH, $file) or die"couldn't open ".$file." $!";

  die " seq ".$seq." is not a Bio::Seq " unless($seq->isa("Bio::SeqI") || 
						$seq->isa("Bio::Seq")  || 
						$seq->isa("Bio::PrimarySeqI"));

  my @operons;
  LINE: while(<FH>){
      #CHROMOSOME_I	SL2	trans-splice_acceptor	6473786	6473787	.	-	.	Note "SL2 trans-splice site; see yk729h2.5"
      #print;
      chomp;
      my @values = split;
      if($values[1] ne 'SL2'){
	next LINE;
      }
      if($_ =~ /personal\s+communication/){
	print STDERR "Can't put information in table about ".$_."\n";
	next LINE;
      }
      my ($start, $end, $strand, $id, $count);
      $count = 0;
  
      $start = $values[3];
      $end = $values[4];
      $strand = $values[6];
      $id = $values[14];
      if(!$id){
	$id = $values[13];
      }
      if(($id =~ /cDNA/) || ($id =~ /EST/)){
	$id = $values[15];
      }
      if($id eq '"'){
	$id = $values[13];
	if($id eq 'ESTyk1004g06.5'){
	  $id =~ s/EST//g;
	}
      }
      if((!$id) || ($id eq '"')){
	$id = $values[12];
      }
      #print $_."\n";
      #print "have ".$start." ".$end." ".$strand." ".$id."\n";
      if($strand eq '+'){
	$strand = 1;
      }else{
	$strand = -1;
      }
      #print "have id ".$id."\n";
      if((!$id) || ($id eq '.') || ($id =~ /EST/) || ($id =~ /cDNA/) || ($id eq '"')){
	#foreach my $v(@values){
	#  print $count." ".$v."\n";
	#  $count++;
	#}
	print "line ".$_." produced a weird id ".$id."\n";
	next LINE
      }
      $id =~ s/\"//g;
      my $simple_feature = &create_simple_feature($start, $end, $strand, $id, $seq, $analysis);
      push(@operons, $simple_feature);
      }

  return \@operons ;
}


sub create_simple_feature{
  my ($start, $end, $strand, $id, $seq, $analysis) = @_;
  #warn "first: $start, $end, $strand, $id...";

  my $simple_feature = Bio::EnsEMBL::SimpleFeature->new();
  $simple_feature->start($start);
  $simple_feature->strand($strand);
  $simple_feature->end($end);
  $simple_feature->display_label($id);
  $simple_feature->slice($seq);
  $simple_feature->analysis($analysis);

  return $simple_feature;
}


sub write_simple_features{
  my ($operons, $db) = @_;
  eval{ print "\n check 1: ".$$operons[0]->display_label };
  eval{ print "\n check 2: ".$$operons[0]->start." - ".$$operons[0]->end."\n" };

  my $operon_adaptor = $db->get_SimpleFeatureAdaptor;

  eval{
    $operon_adaptor->store(@$operons);
  };
  if($@){
    die "couldn't store simple features problems ".$@;
  }
}


#not used at this time...
sub parse_tRNA_genes{
  my $type = "tRNA";
  my ($file, $seq, $analysis) = @_;
  &parse_pseudo_files($file, $seq, $analysis, $type);
}

sub parse_rRNA_genes{
  my $type = "rRNA";
  my ($file, $seq, $analysis) = @_;
  &parse_pseudo_files($file, $seq, $analysis, $type);
}
=head2 parse_pseudo_gff

  Arg [1]   : gff file path
  Arg [2]   : Bio:Seq object
  Arg [3]   : analysis object
  Function  : parse pseudogene information
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub parse_pseudo_gff{
  my $type = "Pseudogene";
  my ($file, $seq, $analysis) = @_;
  &parse_pseudo_files($file, $seq, $analysis, $type);
}


=head2 parse_pseudo_files

  Arg [1]   : gff file path
  Arg [2]   : Bio:Seq object
  Arg [3]   : analysis object
  Arg [4]   : type of feature to be parsed
  Function  : parse specific feature information from gff file
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub parse_pseudo_files{
  my ($file, $seq, $analysis, $types) = @_;

  open(FH, $file) or die "couldn't open ".$file." $!";

  die " seq ".$seq." is not a Bio::Seq " unless($seq->isa("Bio::SeqI") || 
						$seq->isa("Bio::Seq")  || 
						$seq->isa("Bio::PrimarySeqI"));
  my @genes;
  my ($transcripts) = &process_pseudo_files(\*FH, $types);
  print "there are ".keys(%$transcripts)." distinct special transcripts\n";

  my ($processed_transcripts) = &process_pseudo_transcripts($transcripts, $seq, $analysis);
  print "there are ".keys(%$processed_transcripts)." processed special transcripts\n";

  my $genes = undef;
  $genes = &create_pseudo_transcripts($processed_transcripts);
  print "PARSE GFF there are ".keys(%$genes)." special genes\n";

  foreach my $gene_id(keys(%$genes)){
    my $transcripts = $genes->{$gene_id};
    my $gene = &create_gene($transcripts, $gene_id);
    push(@genes, $gene);
  }
  close(FH);
  return \@genes;
}


#generic version
sub process_pseudo_files{
  my ($fh, $types) = @_;
  my %transcripts;

 LOOP: while(<$fh>){
    chomp;
    if(/^##/){
      next LOOP;
    }
    my($chr, $status, $type, $start, $end, $score, $strand, $frame, $sequence, $gene) = split;
    my $element = $_;
    if($chr =~ /sequence-region/){
      #print STDERR $_;
      next LOOP;
    }
    if(!$status && !$type){
      print "status and type no defined skipping\n";
      next LOOP;
    }
    my $line = $status." ".$type;
    #print "\n".$types." exon (".$line.")\n";
    if($line ne $types.' exon'){
      #print "ignoring ".$line."\n";
      next LOOP;
    }
    $gene =~ s/\"//g;
   # print "gene ".$gene."\n";
    if(!$transcripts{$gene}){
      $transcripts{$gene} = [];
      push(@{$transcripts{$gene}}, $element);
    }else{
      push(@{$transcripts{$gene}}, $element);
    }
  }
  return \%transcripts;
}


sub process_pseudo_transcripts{
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
      #print STDERR $line."\n";
      my($chr, $status, $type, $start, $end, $score, $strand, $frame, $sequence, $gene) = split /\s+/, $line;
      $chr =~ s/CHROMOSOME_//;
      if($start == $end){
	next;
      }

      my $exon = new Bio::EnsEMBL::Exon;
      if($frame eq '.'){
	$frame = 0;
      }
      my $phase = (3 - $frame)%3; # wormbase gff cotains frame which is effectively the opposite of phase 
                                  # for a good explaination of phase see the Bio::EnsEMBL::Exon documentation
      #print STDERR "phase calculated to be ".$phase."\n";
      $exon->start($start);
      $exon->end($end);
      $exon->analysis($analysis);
      $exon->slice($slice);
      $exon->phase($phase);
      my $end_phase = ($phase + ($exon->end-$exon->start) + 1)%3;
      #print STDERR "end phase calculated to be ".$end_phase."\n";
      $exon->end_phase($end_phase);
      if($strand eq '+'){
	$exon->strand(1);
      }else{
	$exon->strand(-1);
      }
      #$exon->score(100);
      push(@exons, $exon);
    }
    if($exons[0]->strand == -1){
      @exons = sort{$b->start <=> $a->start} @exons;
    }else{
      @exons = sort{$a->start <=> $b->start} @exons;
    }
    my $phase = 0;
    foreach my $e(@exons){
      push(@{$transcripts{$name}}, $e);
    }
  }

  return (\%transcripts);
}


sub create_pseudo_transcripts{
  my ($transcripts) = @_;

  my %transcripts = %$transcripts;
  my %genes;
  my $gene_name;
  my $transcript_id;
  foreach my $transcript(keys(%transcripts)){
    my $time = time;
    my @exons = @{$transcripts{$transcript}};
    if($transcript =~ /\w+\.\d+[a-z A-Z]/){
     ($gene_name) = $transcript =~ /(\w+\.\d+)[a-z A-Z]/;
     $transcript_id = $transcript;
    }else{
      $gene_name = $transcript;
      $transcript_id = $transcript;
    }
    my $transcript = new Bio::EnsEMBL::Transcript;
    my @sorted_exons;
    if($exons[0]->strand == 1){
      @sorted_exons = sort{$a->start <=> $b->start} @exons
    }else{
      @sorted_exons = sort{$b->start <=> $a->start} @exons  
    }
    my $exon_count = 1;
    my $phase = 0;
    foreach my $exon(@sorted_exons){
      $exon->version(1);
      $exon->stable_id($transcript_id.".".$exon_count);
      $exon_count++;
      $transcript->add_Exon($exon);
    }
    $transcript->version(1);
    $transcript->stable_id($transcript_id);
    if(!$genes{$gene_name}){
      $genes{$gene_name} = [];
      push(@{$genes{$gene_name}}, $transcript);
    }else{
      push(@{$genes{$gene_name}}, $transcript);
    }
  }
  return \%genes;

}

sub store_coord_system{
  my ($db, $name, $version, $sequence_level, $default, $rank) = @_;
  
  my $csa = $db->get_CoordSystemAdaptor();
  
  my $cs = Bio::EnsEMBL::CoordSystem->new
    (
     -NAME            => $name,
     -VERSION         => $version,
     -DEFAULT         => $default,
     -SEQUENCE_LEVEL  => $sequence_level,
     #-TOP_LEVEL       => $top_level,
     -RANK            => $rank
    );
  
  $csa->store($cs);

  return $cs;
}



sub store_slice{
  my ($db, $name, $start, $end, $strand, $coord_system, $sequence) = @_;
  
  my $sa  = $db->get_SliceAdaptor();

  my $slice = Bio::EnsEMBL::Slice->new
  (-seq_region_name  => $name,
   -start            => $start,
   -end              => $end,
   -strand           => $strand,
   -coord_system     => $coord_system);
  
  my $seq_ref;
  if($sequence){
    $seq_ref = \$sequence;
  }
  $sa->store($slice, $seq_ref);
  $slice->adaptor($sa);
  return $slice;
}


1;
