#
#
# Cared for by Ensembl  <ensembl-dev@ebi.ac.uk>
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::FPC_BlastMiniGenewise

=head1 SYNOPSIS

my $obj = Bio::EnsEMBL::Pipeline::RunnableDB::MiniGenewise->new(
					     -dbobj     => $db,
					     -input_id  => $id,
					     -golden_path => $gp,
					     -type      => $type,
                                             -threshold => $threshold		    
								    
                                             );
    $obj->fetch_input
    $obj->run

    my @newfeatures = $obj->output;


=head1 DESCRIPTION

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::RunnableDB::FPC_BlastMiniGenewise;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Getseqs;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;
use Bio::EnsEMBL::Pipeline::GeneConf qw (
					 GB_GOLDEN_PATH
					 GB_PROTEIN_INDEX
					 GB_SIMILARITY_TYPE
					 GB_SIMILARITY_THRESHOLD
					);

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB );

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);    
      
 print STDERR "checking: " . $::db_conf{'dbname'} ." : " . $::db_conf{'dbuser'} . " : " . $::db_conf{'dbpass'} . " : " . $::db_conf{'dbhost'} . "\n";

 
    if(!defined $self->seqfetcher) {
      my $seqfetcher =  $self->make_seqfetcher();
      $self->seqfetcher($seqfetcher);
    }
       
    my ($path, $type, $threshold) = $self->_rearrange([qw(GOLDEN_PATH TYPE THRESHOLD)], @args);

    if(!defined $path || $path eq ''){
      $path = $GB_GOLDEN_PATH;
    }

    if(!defined $type || $type eq ''){
      $type = $GB_SIMILARITY_TYPE;
    }
    
    if(!defined $threshold){
      $threshold = $GB_SIMILARITY_THRESHOLD
    }

    $path = 'UCSC' unless (defined $path && $path ne '');
    $type = 'sptr' unless (defined $type && $type ne '');
    $threshold = 200 unless (defined($threshold));

    $self->dbobj->static_golden_path_type($path);
    $self->type($type);
    $self->threshold($threshold);

    return $self; 
  }


sub type {
  my ($self,$type) = @_;

  if (defined($type)) {
    $self->{_type} = $type;
  }
  return $self->{_type};
}


=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output
    Function:   Writes output data to db
    Returns :   array of exons (with start and end)
    Args    :   none

=cut

sub write_output {
    my($self,@features) = @_;

    my $gene_adaptor = $self->dbobj->get_GeneAdaptor;

  GENE: foreach my $gene ($self->output) {	
      # do a per gene eval...
      eval {
	$gene_adaptor->store($gene);
	print STDERR "wrote gene " . $gene->dbID . "\n";
      }; 
      if( $@ ) {
	  print STDERR "UNABLE TO WRITE GENE\n\n$@\n\nSkipping this gene\n";
      }
	    
  }
   
}

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data from the database
    Returns :   nothing
    Args    :   none

=cut

  sub fetch_input {
    my( $self) = @_;
    
    print STDERR "Fetching input id : " . $self->input_id. " \n\n";

    $self->throw("No input id") unless defined($self->input_id);
    
    my $chrid  = $self->input_id;
       $chrid =~ s/\.(.*)-(.*)//;

    my $chrstart = $1;
    my $chrend   = $2;

    my $stadaptor = $self->dbobj->get_StaticGoldenPathAdaptor();
    my $contig    = $stadaptor->fetch_VirtualContig_by_chr_start_end($chrid,$chrstart,$chrend);
    my $genseq    = $contig->get_repeatmasked_seq;

    $contig->_chr_name($chrid);

    print STDERR "Chromosome id : $chrid\n";
    print STDERR "Range         : $chrstart - $chrend\n";
    print STDERR "Contig        : " . $contig->id . " \n";
    print STDERR "Length is     : " . $genseq->length . "\n\n";

    print STDERR "Fetching features \n\n";

    # need to pass in bp value of zero to prevent globbing on StaticContig.
    my @features  = $contig->get_all_SimilarityFeatures_above_score($self->type, $self->threshold);

    # lose version numbers - probably temporary till pfetch indices catch up

    foreach my $f(@features) {
      my $name = $f->hseqname;
      if ($name =~ /(\S+)\.\d+/) { 
	$f->hseqname($1);
      }
    }

    print STDERR "Number of features = " . scalar(@features) . "\n\n";

    my @genes     = $contig->get_Genes_by_Type('TGE_gw');

    print STDERR "Found " . scalar(@genes) . " genewise genes\n\n";

    my %redids;
    my $trancount = 1;

    # check which TargettedGenewise exons overlap similarity features
    foreach my $gene (@genes) {

#      print STDERR "Found genewise gene " . $gene->dbID . "\n";

      foreach my $tran ($gene->each_Transcript) {

	foreach my $exon ($tran->get_all_Exons) {
	  
#	  print STDERR "Exon " . $exon->dbID . " " . $exon->strand . "\n";

	  if ($exon->seqname eq $contig->id) {
	    
	  FEAT: foreach my $f (@features) {
	      if ($exon->overlaps($f)) {
		$redids{$f->hseqname} = 1;
#		print STDERR "ID " . $f->hseqname . " covered by genewise\n";
	      }
	    }
	  }
	}
	$trancount++;
      }
    }

    my %idhash;
    
    # collect those features which haven't been used by Targetted GeneWise
    foreach my $f (@features) {
#      print "Feature : " . $f->gffstring . "\n";
      
      if ($f->isa("Bio::EnsEMBL::FeaturePair") && 
	  defined($f->hseqname) &&
	  $redids{$f->hseqname} != 1) {
	$idhash{$f->hseqname} = 1;
	
      }
    }
    
    my @ids = keys %idhash;

#    print STDERR "Feature ids are @ids\n";

    my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise('-genomic'  => $genseq,
									   '-ids'      => \@ids,
									   '-seqfetcher' => $self->seqfetcher,
									   '-trim'     => 1);
    
    
    $self->runnable($runnable);
    # at present, we'll only ever have one ...
    $self->vc($contig);
}     


=head2 run

    Title   :   run
    Usage   :   $self->run
    Function:   calls the run method on each runnable, and then calls convert_output
    Returns :   nothing, but $self->output contains results
    Args    :   none

=cut

sub run {
    my ($self) = @_;

    # is there ever going to be more than one?
    foreach my $runnable ($self->runnable) {
	$runnable->run;
    }
    
    $self->convert_output;

}

=head2 convert_output

  Title   :   convert_output
  Usage   :   $self->convert_output
  Function:   converts output from each runnable into gene predictions
  Returns :   nothing, but $self->output contains results
  Args    :   none

=cut

sub convert_output {
  my ($self) =@_;
  
  my $trancount = 1;
  my $genetype;
  foreach my $runnable ($self->runnable) {
    if ($runnable->isa("Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise")){
      $genetype = "similarity_genewise";
    }
    else{
      $self->throw("I don't know what to do with $runnable");
    }

    my $anaAdaptor = $self->dbobj->get_AnalysisAdaptor;
    my @analyses = $anaAdaptor->fetch_by_logic_name($genetype);
    my $analysis_obj;
    if(scalar(@analyses) > 1){
      $self->throw("panic! > 1 analysis for $genetype\n");
    }
    elsif(scalar(@analyses) == 1){
      $analysis_obj = $analyses[0];
    }
    else{
      # make a new analysis object
      $analysis_obj = new Bio::EnsEMBL::Analysis
	(-db              => 'NULL',
	 -db_version      => 1,
	 -program         => $genetype,
	 -program_version => 1,
	 -gff_source      => $genetype,
	 -gff_feature     => 'gene',
	 -logic_name      => $genetype,
	 -module          => 'FPC_BlastMiniGenewise',
      );
    }

    my @results = $runnable->output;
    my $genes = $self->make_genes($genetype, $analysis_obj, \@results);

    my $remapped = $self->remap_genes($genes);

    $self->output(@$remapped);

  }
}


=head2 make_genes

  Title   :   make_genes
  Usage   :   $self->make_genes
  Function:   makes Bio::EnsEMBL::Genes out of the output from runnables
  Returns :   array of Bio::EnsEMBL::Gene  
  Args    :   $genetype: string
              $analysis_obj: Bio::EnsEMBL::Analysis
              $results: reference to an array of FeaturePairs

=cut

sub make_genes {
  my ($self, $genetype, $analysis_obj, $results) = @_;
  my @genes;

  foreach my $tmpf (@$results) {
    my $unchecked_transcript = $self->_make_transcript($tmpf, $self->vc, $genetype, $analysis_obj);
    
    next unless defined ($unchecked_transcript);

    # validate transcript
    my $valid_transcripts = $self->validate_transcript($unchecked_transcript);

    
    # make genes from valid transcripts
    foreach my $checked_transcript(@$valid_transcripts){
      # check coverage of parent protein
      my $covered = $self->check_coverage($checked_transcript);
      if($covered){
	my $gene = new Bio::EnsEMBL::Gene;
	$gene->type($genetype);
	$gene->analysis($analysis_obj);
	$gene->add_Transcript($checked_transcript);
	
	push (@genes, $gene);
      }
    }
  }

  return \@genes;

}

=head2 _make_transcript

 Title   : make_transcript
 Usage   : $self->make_transcript($gene, $contig, $genetype)
 Function: makes a Bio::EnsEMBL::Transcript from a SeqFeature representing a gene, 
           with sub_SeqFeatures representing exons.
 Example :
 Returns : Bio::EnsEMBL::Transcript with Bio::EnsEMBL:Exons(with supporting feature 
           data), and a Bio::EnsEMBL::translation
 Args    : $gene: Bio::EnsEMBL::SeqFeatureI, $contig: Bio::EnsEMBL::DB::ContigI,
  $genetype: string, $analysis_obj: Bio::EnsEMBL::Analysis


=cut

sub _make_transcript{
  my ($self, $gene, $contig, $genetype, $analysis_obj) = @_;
  $genetype = 'unspecified' unless defined ($genetype);

  unless ($gene->isa ("Bio::EnsEMBL::SeqFeatureI"))
    {print "$gene must be Bio::EnsEMBL::SeqFeatureI\n";}
  unless ($contig->isa ("Bio::EnsEMBL::DB::ContigI"))
    {print "$contig must be Bio::EnsEMBL::DB::ContigI\n";}

  my $transcript   = new Bio::EnsEMBL::Transcript;
  my $translation  = new Bio::EnsEMBL::Translation;    
  $transcript->translation($translation);

  my $excount = 1;
  my @exons;
    
  foreach my $exon_pred ($gene->sub_SeqFeature) {
    # make an exon
    my $exon = new Bio::EnsEMBL::Exon;
    
    $exon->contig_id($contig->id);
    $exon->start($exon_pred->start);
    $exon->end  ($exon_pred->end);
    $exon->strand($exon_pred->strand);
    
    $exon->phase($exon_pred->phase);
    $exon->attach_seq($contig);
    
    # sort out supporting evidence for this exon prediction
    foreach my $subf($exon_pred->sub_SeqFeature){
      $subf->feature1->source_tag($genetype);
      $subf->feature1->primary_tag('similarity');
      $subf->feature1->score(100);
      $subf->feature1->analysis($analysis_obj);
	
      $subf->feature2->source_tag($genetype);
      $subf->feature2->primary_tag('similarity');
      $subf->feature2->score(100);
      $subf->feature2->analysis($analysis_obj);
      
      $exon->add_Supporting_Feature($subf);
    }
    
    push(@exons,$exon);
    
    $excount++;
  }
  
  if ($#exons < 0) {
    print STDERR "Odd.  No exons found\n";
    return;
  } 
  else {
    
#    print STDERR "num exons: " . scalar(@exons) . "\n";

    if ($exons[0]->strand == -1) {
      @exons = sort {$b->start <=> $a->start} @exons;
    } else {
      @exons = sort {$a->start <=> $b->start} @exons;
    }
    
    foreach my $exon (@exons) {
      $transcript->add_Exon($exon);
    }
    
    $translation->start_exon($exons[0]);
    $translation->end_exon  ($exons[$#exons]);
    
    if ($exons[0]->phase == 0) {
      $translation->start(1);
    } elsif ($exons[0]->phase == 1) {
      $translation->start(3);
    } elsif ($exons[0]->phase == 2) {
      $translation->start(2);
    }
    
    $translation->end  ($exons[$#exons]->end - $exons[$#exons]->start + 1);
  }
  
  return $transcript;
}

=head2 validate_transcript

 Title   : validate_transcript 
 Usage   : my @valid = $self->validate_transcript($transcript)
 Function: Validates a transcript - rejects if mixed strands, splits if long introns
 Returns : Ref to @Bio::EnsEMBL::Transcript
 Args    : Bio::EnsEMBL::Transcript

=cut
sub validate_transcript {
  my ( $self, $transcript ) = @_;
  my @valid_transcripts;
  
  my $valid = 1;
  my $split = 0;

  my $previous_exon;
  foreach my $exon($transcript->get_all_Exons){
    if (defined($previous_exon)) {
      my $intron;
      
      if ($exon->strand == 1) {
	$intron = abs($exon->start - $previous_exon->end - 1);
      } else {
	$intron = abs($previous_exon->start - $exon->end - 1);
      }
      
      if ($intron > 100000) {
	print STDERR "Intron too long $intron  for transcript " . $transcript->{'temporary_id'} . "\n";
	$split = 1;
	$valid = 0;
      }
      
      if ($exon->strand != $previous_exon->strand) {
	print STDERR "Mixed strands for gene " . $transcript->{'temporary_id'} . "\n";
	$valid = 0;
	return;
      }
    }
    $previous_exon = $exon;
  }
  
  if ($valid) {
    # make a new transcript that's a copy of all the important parts of the old one
    # but without all the db specific gubbins
    my $newtranscript  = new Bio::EnsEMBL::Transcript;
    my $newtranslation = new Bio::EnsEMBL::Translation;

    $newtranscript->translation($newtranslation);
    $newtranscript->translation->start_exon($transcript->translation->start_exon);
    $newtranscript->translation->end_exon($transcript->translation->end_exon);
    $newtranscript->translation->start($transcript->translation->start);
    $newtranscript->translation->end($transcript->translation->end);
    foreach my $exon($transcript->get_all_Exons){
      $newtranscript->add_Exon($exon);
      foreach my $sf($exon->each_Supporting_Feature){
	  $sf->feature1->seqname($exon->contig_id);
      }
    }

    push(@valid_transcripts,$newtranscript);
  }
  elsif ($split){
    # split the transcript up.
    my $split_transcripts = $self->split_transcript($transcript);
    push(@valid_transcripts, @$split_transcripts);
  }

  return \@valid_transcripts;
}

=head2 split_transcript

 Title   : split_transcript 
 Usage   : my @splits = $self->split_transcript($transcript)
 Function: splits a transcript into multiple transcripts at long introns. Rejects single exon 
           transcripts that result. 
 Returns : Ref to @Bio::EnsEMBL::Transcript
 Args    : Bio::EnsEMBL::Transcript

=cut


sub split_transcript{
  my ($self, $transcript) = @_;
  $transcript->sort;

  my @split_transcripts   = ();

  if(!($transcript->isa("Bio::EnsEMBL::Transcript"))){
    $self->warn("[$transcript] is not a Bio::EnsEMBL::Transcript - cannot split");
    return (); # empty array
  }
  
  my $prev_exon;
  my $exon_added = 0;
  my $curr_transcript = new Bio::EnsEMBL::Transcript;
  my $translation     = new Bio::EnsEMBL::Translation;
  $curr_transcript->translation($translation);

EXON:   foreach my $exon($transcript->get_all_Exons){


    $exon_added = 0;
      # is this the very first exon?
    if($exon == $transcript->start_exon){

      $prev_exon = $exon;
      
      # set $curr_transcript->translation start and start_exon
      $curr_transcript->add_Exon($exon);
      $exon_added = 1;
      $curr_transcript->translation->start_exon($exon);
      $curr_transcript->translation->start($transcript->translation->start);
      push(@split_transcripts, $curr_transcript);
      next EXON;
    }
    
    if ($exon->strand != $prev_exon->strand){
      return (); # empty array
    }

    # We need to start a new transcript if the intron size between $exon and $prev_exon is too large
    my $intron = 0;
    if ($exon->strand == 1) {
      $intron = abs($exon->start - $prev_exon->end - 1);
    } else {
      $intron = abs($prev_exon->start - $exon->end - 1);
    }
    
    if ($intron > 100000) {
      # set translation end and end_exon of $curr_transcript->translation based on phase of $prev_exon
      $curr_transcript->translation->end_exon($prev_exon);
      $curr_transcript->translation->end($prev_exon->end - $prev_exon->start + 1);
      
      # start a new transcript 
      my $t  = new Bio::EnsEMBL::Transcript;
      my $tr = new Bio::EnsEMBL::Translation;
      $t->translation($tr);

      # add exon unless already added, and set translation start and start_exon
      $t->add_Exon($exon) unless $exon_added;
      $exon_added = 1;

      $t->translation->start_exon($exon);

      if ($exon->phase == 0) {
	$t->translation->start(1);
      } elsif ($exon->phase == 1) {
	$t->translation->start(3);
      } elsif ($exon->phase == 2) {
	$t->translation->start(2);
      }

      # this new transcript becomes the current transcript
      $curr_transcript = $t;

      push(@split_transcripts, $curr_transcript);
    }

    if($exon == $transcript->end_exon){
      # add it unless already added
      $curr_transcript->add_Exon($exon) unless $exon_added;
      $exon_added = 1;

      # set $curr_transcript end_exon and end
      $curr_transcript->translation->end_exon($exon);
      $curr_transcript->translation->end($transcript->translation->end);
    }

    else{
      # just add the exon
      $curr_transcript->add_Exon($exon) unless $exon_added;
    }
    foreach my $sf($exon->each_Supporting_Feature){
	  $sf->feature1->seqname($exon->contig_id);

      }
    # this exon becomes $prev_exon for the next one
    $prev_exon = $exon;

  }

  # discard any single exon transcripts
  my @final_transcripts = ();
  my $count = 1;
  
  foreach my $st(@split_transcripts){
    $st->sort;
    my @ex = $st->get_all_Exons;
    if(scalar(@ex) > 1){
      $st->{'temporary_id'} = $transcript->dbID . "." . $count;
      $count++;
      push(@final_transcripts, $st);
    }
  }

  return \@final_transcripts;

}


=head2 check_coverage

 Title   : check_coverage
 Usage   :
 Function: checks how much of the parent protein is covered by the genewise prediction
 Example :
 Returns : 1 if > 70% coverage, otherwise 0
 Args    :


=cut

sub check_coverage{
  my ($self, $transcript) = @_;
  my $threshold = 70;
  my $pstart = 0;
  my $pend = 0;
  my $protname;
  my $plength;

  my $matches = 0;

  foreach my $exon($transcript->get_all_Exons) {
    $pstart = 0;
    $pend   = 0;
    
    foreach my $f($exon->each_Supporting_Feature){
      
      if (!defined($protname)){
	$protname = $f->hseqname;
      }
      if($protname ne $f->hseqname){
	warn("$protname ne " . $f->hseqname . "\n");
      }
      
      if((!$pstart) || $pstart > $f->hstart){
	$pstart = $f->hstart;
      }
      
      if((!$pend) || $pend < $f->hend){
	$pend= $f->hend;
      }
    }
    $matches += ($pend - $pstart + 1);
  }
  
  my $seq; 
  eval{
    $seq = $self->seqfetcher->get_Seq_by_acc($protname);
  };
  if ($@) {
    $self->warn("Error fetching sequence for [$protname]\n");
  }
  
  if(!defined $seq){
    $self->warn("No sequence fetched for [$protname] - can't check coverage, letting gene through\n");
    return 1;
  }
  
  $plength = $seq->length;

  if(!defined($plength) || $plength == 0){
    warn("no sensible length for $protname - can't get coverage\n");
    return 0;
  }

  my $coverage = $matches/$plength;
  $coverage *= 100;
  if ($coverage < $threshold){
    $self->warn ("Coverage of $protname is only $coverage - will be rejected\n");
    return 0;
  }
  
  print STDERR "Coverage of $protname is $coverage - will be accepted\n";
  return 1;
}


=head2 remap_genes

 Title   : remap_genes
 Usage   : $self->remap_genes($runnable, @genes)
 Function: converts the coordinates of each Bio@EnsEMBL::Gene in @genes into RawContig
           coordinates for storage.
 Example : 
 Returns : array of Bio::EnsEMBL::Gene in RawContig coordinates
 Args    : @genes: array of Bio::EnsEMBL::Gene in virtual contig coordinates


=cut

sub remap_genes {
  my ($self, $genes) = @_;
  my $contig = $self->vc;

  my @newf;
  my $trancount=1;
  foreach my $gene (@$genes) {
    eval {
      my $newgene = $contig->convert_Gene_to_raw_contig($gene);
      # need to explicitly add back genetype and analysis.
      $newgene->type($gene->type);
      $newgene->analysis($gene->analysis);

      foreach my $tran ($newgene->each_Transcript) {
	foreach my $exon($tran->get_all_Exons) {
	  foreach my $sf($exon->each_Supporting_Feature) {
	    # this should be sorted out by the remapping to rawcontig ... strand is fine
	    if ($sf->start > $sf->end) {
	      my $tmp = $sf->start;
	      $sf->start($sf->end);
	      $sf->end($tmp);
	    }
	  }
	}
      }
      push(@newf,$newgene);

    };
    if ($@) {
      print STDERR "Couldn't reverse map gene " . $gene->id . " [$@]\n";
    }
    

  }

  return \@newf;
}

=head2 output

 Title   : output
 Usage   :
 Function: get/set for output array
 Example :
 Returns : array of Bio::EnsEMBL::Gene
 Args    :


=cut

sub output{
   my ($self,@genes) = @_;

   if (!defined($self->{'_output'})) {
     $self->{'_output'} = [];
   }
    
   if(defined @genes){
     push(@{$self->{'_output'}},@genes);
   }

   return @{$self->{'_output'}};
}

=head2 make_seqfetcher

 Title   : make_seqfetcher
 Usage   :
 Function: makes a Bio::EnsEMBL::SeqFetcher to be used for fetching protein sequences. If 
           GB_PROTEIN_INDEX is specified in GeneConf.pm, then a Getseqs 
           fetcher is made, otherwise it will be Pfetch
 Example :
 Returns : Bio::EnsEMBL::SeqFetcher
 Args    :


=cut

sub make_seqfetcher {
  my ( $self ) = @_;
  my $index   = $GB_PROTEIN_INDEX;

  my $seqfetcher;

  if(defined $index && $index ne ''){
    my @db = ( $index );
    $seqfetcher = new Bio::EnsEMBL::Pipeline::SeqFetcher::Getseqs('-db' => \@db,);
  }
  else{
    # default to Pfetch
    $seqfetcher = new Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;
  }

  return $seqfetcher;

}

1;
