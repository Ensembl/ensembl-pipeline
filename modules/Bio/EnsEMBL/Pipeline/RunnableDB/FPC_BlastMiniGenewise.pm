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
					     -db        => $db,
					     -input_id  => $id,
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

use Bio::EnsEMBL::DnaPepAlignFeature;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Pipeline::Runnable::Protein::Seg;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Databases qw (
							     GB_GW_DBNAME
							     GB_GW_DBHOST
							     GB_GW_DBUSER
							     GB_GW_DBPASS
							    );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Similarity qw (
							     GB_SIMILARITY_DATABASES
							     GB_SIMILARITY_COVERAGE
							     GB_SIMILARITY_MAX_INTRON
							     GB_SIMILARITY_MIN_SPLIT_COVERAGE
							     GB_SIMILARITY_GENETYPE
							     GB_SIMILARITY_MAX_LOW_COMPLEXITY
							    );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Targetted  qw (
							     GB_TARGETTED_GW_GENETYPE
							    );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::General    qw (
							     GB_INPUTID_REGEX
							     GB_REPEAT_MASKING
							    );

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB );

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);    
      
    # make sure at least one protein source database has been defined
    
    $self->throw("no protein source databases defined in Config::GeneBuild::Similarity::GB_SIMILARITY_DATABASES\n") 
      unless scalar(@{$GB_SIMILARITY_DATABASES});
    
    # make all seqfetchers
    foreach my $db(@{$GB_SIMILARITY_DATABASES}){
      my $seqfetcher =  $self->make_seqfetcher($db->{'index'}, $db->{seqfetcher});  
      $self->add_seqfetcher_by_type($db->{'type'}, $seqfetcher);
    }

    # IMPORTANT
    # SUPER creates db, which is a reference to GB_DBHOST@GB_DBNAME containing
    # features and dna
    # Here it is used as refdb only and we need to make a connection to GB_GW_DBNAME@GB_GW_DBHOST
    my $genewise_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
							 '-host'   => $GB_GW_DBHOST,
							 '-user'   => $GB_GW_DBUSER,
							 '-pass'   => $GB_GW_DBPASS,
							 '-dbname' => $GB_GW_DBNAME,
							);
    
    
    $genewise_db->dnadb($self->db);
    $self->output_db($genewise_db);
    
    return $self; 
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

    my $gene_adaptor = $self->output_db->get_GeneAdaptor;
    my @genes = $self->output;
    print STDERR "have ".@genes." genes\n";
  GENE: foreach my $gene (@genes) {	
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
    my $input_id = $self->input_id;
   
    $input_id =~ /$GB_INPUTID_REGEX/;
    my $chrid = $1;
    my $chrstart = $2;
    my $chrend   = $3;

    
    my $sla       = $self->db->get_SliceAdaptor();
    my $slice     = $sla->fetch_by_chr_start_end($chrid,$chrstart,$chrend);
    my $genseq    = $slice->get_repeatmasked_seq($GB_REPEAT_MASKING);

    
    $slice->chr_name($chrid);

    #print STDERR "Chromosome id : $chrid\n";
    #print STDERR "Range         : $chrstart - $chrend\n";
    #print STDERR "Contig        : " . $slice->id . " \n";
    #print STDERR "Length is     : " . $genseq->length . "\n\n";
    
    my @genes     = @{$slice->get_all_Genes_by_type($GB_TARGETTED_GW_GENETYPE)};
    #print STDERR "Found " . scalar(@genes) . " genewise genes\n\n";
    
    foreach my $database(@{$GB_SIMILARITY_DATABASES}){
      
      print STDERR "Fetching features for " . $database->{'type'} . 
	" with score above " . $database->{'threshold'}. "\n\n";
      my $pafa = $self->db->get_ProteinAlignFeatureAdaptor();
      my @features  = @{$pafa->fetch_all_by_Slice_and_score($slice, $database->{'threshold'}, $database->{'type'})};
      print STDERR "have ".@features." \n";
      # lose version numbers - probably temporary till pfetch indices catch up
      foreach my $f(@features) {
	#print STDERR "have ".$f." from the feature table\n";
	my $name = $f->hseqname;
	if ($name =~ /(\S+)\.\d+/) { 
	  $f->hseqname($1);
	}
      }
      
      #print STDERR "Number of features = " . scalar(@features) . "\n\n";
      my %redids;
      my $trancount = 1;
      
      # check which TargettedGenewise exons overlap similarity features
      foreach my $gene (@genes) {
	foreach my $tran (@{$gene->get_all_Transcripts}) {
	  foreach my $exon (@{$tran->get_all_Exons}) {
	    if ($exon->seqname eq $slice->id) {
	      
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
	      #print "Feature : " . $f->gffstring . "\n";
	
	if ($f->isa("Bio::EnsEMBL::FeaturePair") && 
	    defined($f->hseqname) &&
	    $redids{$f->hseqname} != 1) {
	  $idhash{$f->hseqname} = 1;
	  
	}
      }
      

      my @ids = keys %idhash;

      my $seqfetcher =  $self->get_seqfetcher_by_type($database->{'type'});
      #print STDERR "Feature ids are @ids\n";
      
      my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise('-genomic'  => $genseq,
									     '-ids'      => \@ids,
									     '-seqfetcher' => $seqfetcher,
									     '-trim'     => 1);
      
      $self->runnable($runnable);

      # at present, we'll only ever have one ...

      $self->query($slice);
    }
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
    #print STDERR "HAVE CONVERTED OUTPUT\n";
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
  my $genetype = $GB_SIMILARITY_GENETYPE;
  #print STDERR "in CONVERT OUTPUT\n";
  foreach my $runnable ($self->runnable) {
    $self->throw("I don't know what to do with $runnable") unless $runnable->isa("Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise");
										 
    if(!defined($genetype) || $genetype eq ''){
      $genetype = 'similarity_genewise';
      $self->warn("Setting genetype to $genetype\n");
    }

    my $anaAdaptor = $self->db->get_AnalysisAdaptor;
    print STDERR $anaAdaptor . "\n";

    my $analysis_obj = $anaAdaptor->fetch_by_logic_name($genetype);

    if ( !defined $analysis_obj ){
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

    my @remapped = @{$self->remap_genes($genes)};
    #print STDERR "HAVE ".@remapped." remapped genes\n";
    $self->output(@remapped);

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

  #print STDERR "making FPC_BlastMiniGenewise genes\n";

  $self->throw("[$analysis_obj] is not a Bio::EnsEMBL::Analysis\n") 
    unless defined($analysis_obj) && $analysis_obj->isa("Bio::EnsEMBL::Analysis");

  foreach my $tmpf (@{$results}) {
    my $unchecked_transcript = $self->_make_transcript($tmpf, $self->query, $genetype, $analysis_obj);
    
    next unless defined ($unchecked_transcript);

    # validate transcript
    my $valid_transcripts = $self->validate_transcript($unchecked_transcript);
    
    # make genes from valid transcripts
    foreach my $checked_transcript(@$valid_transcripts){
      #print "have ".$checked_transcript." as a transcript\n";
      my $gene = new Bio::EnsEMBL::Gene;
      $gene->type($genetype);
      $gene->analysis($analysis_obj);
      $gene->add_Transcript($checked_transcript);
      
      push (@genes, $gene);
    }
  }
  
  return \@genes;

}

=head2 _make_transcript

 Title   : make_transcript
 Usage   : $self->make_transcript($gene, $contig, $genetype, $analysis_obj)
 Function: makes a Bio::EnsEMBL::Transcript from a SeqFeature representing a gene, 
           with sub_SeqFeatures representing exons.
 Example :
 Returns : Bio::EnsEMBL::Transcript with Bio::EnsEMBL:Exons(with supporting feature 
           data), and a Bio::EnsEMBL::translation
 Args    : $gene: Bio::EnsEMBL::SeqFeatureI, $contig: Bio::EnsEMBL::RawContig,
  $genetype: string, $analysis_obj: Bio::EnsEMBL::Analysis


=cut

sub _make_transcript{
  my ($self, $gene, $contig, $genetype, $analysis_obj) = @_;
  $genetype = 'similarity_genewise' unless defined ($genetype);
  #print "have ".$gene." as a gene\n";
  unless ($gene->isa ("Bio::EnsEMBL::SeqFeatureI"))
    {print "$gene must be Bio::EnsEMBL::SeqFeatureI\n";}
#  unless ($contig->isa ("Bio::PrimarySeq"))
#    {print "$contig must be a RawContig or Slice\n";}

  my $transcript   = new Bio::EnsEMBL::Transcript;
  my $translation  = new Bio::EnsEMBL::Translation;    
  $transcript->translation($translation);

  my $excount = 1;
  my @exons;
    
  foreach my $exon_pred ($gene->sub_SeqFeature) {
    # make an exon
    my $exon = new Bio::EnsEMBL::Exon;
    
    $exon->start($exon_pred->start);
    $exon->end  ($exon_pred->end);
    $exon->strand($exon_pred->strand);
    $exon->adaptor($self->db->get_ExonAdaptor);
    $exon->phase($exon_pred->phase);
    $exon->end_phase($exon_pred->end_phase);
    $exon->contig($contig);  # replaces the attach_seq method

    # sort out supporting evidence for this exon prediction
    my @sf = $exon_pred->sub_SeqFeature;
    
    if(@sf){
      my $align = new Bio::EnsEMBL::DnaPepAlignFeature(-features => \@sf); 
      
      $align->seqname($contig->dbID);
      $align->contig($contig);
      my $prot_adp = $self->db->get_ProteinAlignFeatureAdaptor;
      $align->adaptor($prot_adp);
      $align->score(100);
      $align->analysis($analysis_obj);
      
      $exon->add_supporting_features($align);
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
    
    $translation->start_Exon($exons[0]);
    $translation->end_Exon  ($exons[$#exons]);
    
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
 Function: Validates a transcript - rejects if mixed strands, 
                                    rejects if low coverage, 
                                    rejects if stops in translation
                                    splits if long introns and insufficient coverage of parental protein
 Returns : Ref to @Bio::EnsEMBL::Transcript
 Args    : Bio::EnsEMBL::Transcript

=cut
sub validate_transcript {
  my ( $self, $transcript ) = @_;
  my @valid_transcripts;
  
  my $valid = 1;
  my $split = 0;

  

  # check coverage of parent protein
  my $coverage  = $self->check_coverage($transcript);

  if ($coverage < $GB_SIMILARITY_COVERAGE){
    $self->warn (" rejecting transcript for low coverage: $coverage\n");
    $valid = 0;
    return undef;
  }
  #  print STDERR "Coverage of $protname is $coverage - will be accepted\n";
  
  # check for stops in translation
  unless(Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_check_Translation($transcript)){
    $self->warn("discarding transcript - translation has stop codons\n");
    return undef;
  }
  
  # check for low complexity
  my $low_complexity = $self->check_low_complexity($transcript);
  #print STDERR "complexity computed by seg: $low_complexity, max_low_complexity we allow: $GB_SIMILARITY_MAX_LOW_COMPLEXITY\n";
  if($low_complexity > $GB_SIMILARITY_MAX_LOW_COMPLEXITY){
    $self->warn("discarding transcript - translation has $low_complexity% low complexity sequence\n");
    return undef;
  }

  my $previous_exon;
  my @exons = @{$transcript->get_all_Exons};
  foreach my $exon( @exons ){
    if (defined($previous_exon)) {
      
      # check phase consistency
      if ( $previous_exon ){
	my $end_phase = $previous_exon->end_phase;
	my $phase     = $exon->phase;
	if (  $end_phase != $phase ){
	  $self->warn("rejecting transcript with inconsistent phases ($phase - $end_phase) ");
	  return undef;
	}
      }
      
      # check intron size
      my $intron;
      if ($exon->strand == 1) {
	$intron = abs($exon->start - $previous_exon->end - 1);
      } 
      else {
	$intron = abs($previous_exon->start - $exon->end - 1);
      }
      
      # this isn't really working - it's letting through too many genes with long introns  - replace 
      # it with a simple split for the moment
      
      #      if ($intron > $GB_SIMILARITY_MAX_INTRON && $coverage < $GB_SIMILARITY_MIN_SPLIT_COVERAGE) {
      #	print STDERR "Intron too long $intron  for transcript " . $transcript->{'temporary_id'} . " with coverage $coverage\n";
      #	$split = 1;
      #	$valid = 0;
      #      }
      
      if ( $intron > $GB_SIMILARITY_MAX_INTRON ) {
	print STDERR "Intron too long $intron  for transcript " . $transcript->dbID . "\n";
	$split = 1;
	$valid = 0;
      }
      
     
      if ($exon->strand != $previous_exon->strand) {
	print STDERR "Mixed strands for gene " . $transcript->{'temporary_id'} . "\n";
	$valid = 0;
	return undef;
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
    $newtranscript->translation->start_Exon($transcript->translation->start_Exon);
    $newtranscript->translation->end_Exon($transcript->translation->end_Exon);
    $newtranscript->translation->start($transcript->translation->start);
    $newtranscript->translation->end($transcript->translation->end);
    foreach my $exon(@{$transcript->get_all_Exons}){
      $newtranscript->add_Exon($exon);
      foreach my $sf(@{$exon->get_all_supporting_features}){
	  $sf->seqname($exon->contig->dbID);
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

EXON:   foreach my $exon(@{$transcript->get_all_Exons}){


    $exon_added = 0;
      # is this the very first exon?
    if($exon == $transcript->start_Exon){

      $prev_exon = $exon;
      
      # set $curr_transcript->translation start and start_Exon
      $curr_transcript->add_Exon($exon);
      $exon_added = 1;
      $curr_transcript->translation->start_Exon($exon);
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
    
    if ($intron > $GB_SIMILARITY_MAX_INTRON) {
      $curr_transcript->translation->end_Exon($prev_exon);
      
      # need to account for end_phase of $prev_exon when setting translation->end
      $curr_transcript->translation->end($prev_exon->end - $prev_exon->start + 1 - $prev_exon->end_phase);
      
      # start a new transcript 
      my $t  = new Bio::EnsEMBL::Transcript;
      my $tr = new Bio::EnsEMBL::Translation;
      $t->translation($tr);

      # add exon unless already added, and set translation start and start_Exon
      $t->add_Exon($exon) unless $exon_added;
      $exon_added = 1;

      $t->translation->start_Exon($exon);

      if ($exon->phase == 0) {
	$t->translation->start(1);
      } elsif ($exon->phase == 1) {
	$t->translation->start(3);
      } elsif ($exon->phase == 2) {
	$t->translation->start(2);
      }
      
      # start exon always has phase 0
      $exon->phase(0);

      # this new transcript becomes the current transcript
      $curr_transcript = $t;

      push(@split_transcripts, $curr_transcript);
    }

    if($exon == $transcript->end_Exon){
      # add it unless already added
      $curr_transcript->add_Exon($exon) unless $exon_added;
      $exon_added = 1;

      # set $curr_transcript end_Exon and end
      $curr_transcript->translation->end_Exon($exon);
      $curr_transcript->translation->end($transcript->translation->end);
    }

    else{
      # just add the exon
      $curr_transcript->add_Exon($exon) unless $exon_added;
    }
    foreach my $sf(@{$exon->get_all_supporting_features}){
	  $sf->seqname($exon->contig->dbID);

      }
    # this exon becomes $prev_exon for the next one
    $prev_exon = $exon;

  }

  # discard any single exon transcripts
  my @final_transcripts = ();
  my $count = 1;
  
  foreach my $st(@split_transcripts){
    $st->sort;
    my @ex = @{$st->get_all_Exons};
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
 Function: returns how much of the parent protein is covered by the genewise prediction
 Example :
 Returns : percentage
 Args    :


=cut

sub check_coverage{
  my ($self, $transcript) = @_;
  print STDERR "entering check coverage\n";
  my $pstart = 0;
  my $pend = 0;
  my $protname;
  my $plength;

  my $matches = 0;

  foreach my $exon(@{$transcript->get_all_Exons}) {
    $pstart = 0;
    $pend   = 0;
    
    foreach my $f(@{$exon->get_all_supporting_features}){
      
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

# what to do if the same protein is found in >1 seq index? 

 SEQFETCHER:
  foreach my $seqfetcher( $self->each_seqfetcher){
    #print STDERR "FPC_BlastMiniGenewise: getting sequence for $protname\n";
    eval{
      $seq = $seqfetcher->get_Seq_by_acc($protname);
    };
    if ($@) {
      $self->warn("FPC_BMG:Error fetching sequence for [$protname] - trying next seqfetcher:[$@]\n");
    }

    if (defined $seq) {
      last SEQFETCHER;
    }

  }
  
  if(!defined $seq){
    $self->warn("FPC_BMG:No sequence fetched for [$protname] - can't check coverage, letting gene through\n");
    return 100;
  }
  
  $plength = $seq->length;

  if(!defined($plength) || $plength == 0){
    warn("no sensible length for $protname - can't get coverage\n");
    return 0;
  }

 # print STDERR "looking at coverage of $protname\n";

  my $coverage = $matches/$plength;
  $coverage *= 100;
  return $coverage;

}

=head2 check_low_complexity

 Title   : check_complexity
 Usage   :
 Function: uses seg to find low complexity regions in transcript->translate. 
           Calculates overall %low complexity of the translation
 Example :
 Returns : percentage low complexity sequence
 Args    :


=cut

sub check_low_complexity{
  my ($self, $transcript) = @_;
  my $low_complexity;
  eval{
    
    my $protseq = $transcript->translate;
    $protseq->display_id($transcript->{'temporary_id'} . ".translation");

    # ought to be got from analysisprocess table
    my $analysis = Bio::EnsEMBL::Analysis->new(
					       -db           => 'low_complexity',
					       -program      => '/usr/local/ensembl/bin/seg',
					       -program_file => '/usr/local/ensembl/bin/seg',
					       -gff_source   => 'Seg',
					       -gff_feature  => 'annot',
					       -module       => 'Seg',
					       -logic_name   => 'Seg'
	
					      );
    
    print STDERR "about to run seg with query $protseq\n";
    my $seg = new  Bio::EnsEMBL::Pipeline::Runnable::Protein::Seg(    
								  -query    => $protseq,
								  -analysis => $analysis,
								 );
    
    print STDERR "query = [".$seg->query."]\n";
    $seg->run;
    my $lc_length = 0;
    foreach my $feat($seg->output){
      if($feat->end < $feat->start){
	print STDERR "inverting feature\n";
	my $tmp = $feat->start;
	$feat->start($feat->end);
	$feat->end($tmp);
      }

      $lc_length += $feat->end - $feat->start + 1;
    }
    
    $low_complexity = (100*$lc_length)/($protseq->length)
    
  };
  
  if($@){
    print STDERR "problem running seg: \n[$@]\n";
    return 0; # let transcript through
  }

  return $low_complexity;

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
  my $contig = $self->query;

  my @newf;
  my $trancount=1;
  foreach my $gene (@$genes) {
    eval {
      $gene->transform;
      # need to explicitly add back genetype and analysis.
      #$newgene->type($gene->type);
      #$gene->analysis($gene->analysis);

      foreach my $tran (@{$gene->get_all_Transcripts}) {
	foreach my $exon(@{$tran->get_all_Exons}) {
	  foreach my $sf(@{$exon->get_all_supporting_features}) {
	    # this should be sorted out by the remapping to rawcontig ... strand is fine
	    if ($sf->start > $sf->end) {
	      my $tmp = $sf->start;
	      $sf->start($sf->end);
	      $sf->end($tmp);
	    }
	  }
	}
      }
      push(@newf,$gene);

    };
    if ($@) {
      print STDERR "Couldn't reverse map gene " . $gene . " [$@]\n";
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
    
   if(@genes){
     push(@{$self->{'_output'}},@genes);
   }
   #print STDERR "have ".@{$self->{'_output'}}." as output\n";
   return @{$self->{'_output'}};
}

=head2 make_seqfetcher

 Title   : make_seqfetcher
 Usage   :
 Function: makes a Bio::EnsEMBL::SeqFetcher to be used for fetching protein sequences. If 
           $index is specified, then a Getseqs fetcher is made, otherwise it throws
 Example :
 Returns : Bio::EnsEMBL::SeqFetcher
 Args    : string representing path to sequence index


=cut

sub make_seqfetcher {
  my ( $self, $index, $seqfetcher_class ) = @_;
  my $seqfetcher;

  (my $class = $seqfetcher_class) =~ s/::/\//g;
  require "$class.pm";

  if(defined $index && $index ne ''){
    my @db = ( $index );
    
    # make sure that your class is compatible with the index type
    $seqfetcher = "$seqfetcher_class"->new('-db' => \@db, );
    
  }
  else{
    $self->throw("Can't make seqfetcher\n");
  }

  return $seqfetcher;

}

=head2 each_seqfetcher

  Title   :   each_seqfetcher
  Usage   :   my @seqfetchers = $self->each_seqfetcher
  Function:   Returns an array of Bio::DB::RandomAccessI representing the various sequence indices 
              listed in Config::GeneBuild::Similarity::GB_SIMILARITY_DATABASES
  Returns :   Array of Bio::DB::RandomAccessI
  Args    :   none

=cut

sub each_seqfetcher {
  my ($self) = @_;
  
  if(!defined($self->{'_seqfetchers'})) {
     $self->{'_seqfetchers'} = {};
   }
    
  # NB array of seqfetchers
   return values(%{$self->{'_seqfetchers'}});
}

=head2 each_seqfetcher_by_type

  Title   :   each_seqfetcher_by_type
  Usage   :   my %seqfetchers_by_type = $self->each_seqfetcher_by_type
  Function:   Returns a hash of Bio::DB::RandomAccessI representing the various sequence indices 
              listed in Config::GeneBuild::Similarity::GB_SIMILARITY_DATABASES keyed by type listed therein.
  Returns :   Hash of all seqfetchers linking db_type to Bio::DB::RandomAccessI
  Args    :   none

=cut

sub each_seqfetcher_by_type {
  my ($self) = @_;
  
  if(!defined($self->{'_seqfetchers'})) {
     $self->{'_seqfetchers'} = {};
   }
    
  # NB hash of seqfetchers
   return %{$self->{'_seqfetchers'}};
}

=head2 add_seqfetcher_by_type

  Title   :   add_seqfetcher_by_type
  Usage   :   $self->add_seqfetcher_by_type('swall', $seqfetcher)
  Function:   Adds a Bio::DB::RandomAccessI into $self->{'_seqfetchers'} keyed by type
  Returns :   Nothing
  Args    :   $type - string representing db type
              $seqfetcher - Bio::DB::RandomAccesI

=cut

sub add_seqfetcher_by_type{
  my ($self, $type, $seqfetcher) = @_;
  $self->throw("no type specified\n") unless defined ($type); 
  $self->throw("no suitable seqfetcher specified: [$seqfetcher]\n") 
    unless defined ($seqfetcher) && $seqfetcher->isa("Bio::DB::RandomAccessI"); 
  $self->{'_seqfetchers'}{$type} = $seqfetcher;
}


=head2 get_seqfetcher_by_type

  Title   :   get_seqfetcher_by_type
  Usage   :   my $seqfetcher = $self->get_seqfetcher_by_type('swall')
  Function:   Fetches the seqfetcher associated with a particular db type as specified in Config::GeneBuild::Similarity::GB_SIMILARITY_DATABASES
  Returns :   Bio::DB::RandomAccessI
  Args    :   $type - string representing db type

=cut
sub get_seqfetcher_by_type{
  my ($self, $type) = @_;
  my %seqfetchers = $self->each_seqfetcher_by_type;
  foreach my $dbtype(keys %seqfetchers){
    if ($dbtype eq $type){
      return $seqfetchers{$dbtype};
    }
  }
}

=head2 output_db

 Title   : output_db
 Usage   : needs to be moved to a genebuild base class
 Function: 
           
 Returns : 
 Args    : 

=cut

sub output_db {
    my( $self, $output_db ) = @_;
    
    if ($output_db) 
    {
	$output_db->isa("Bio::EnsEMBL::DBSQL::DBAdaptor")
	    || $self->throw("Input [$output_db] isn't a Bio::EnsEMBL::DBSQL::DBAdaptor");
	$self->{_output_db} = $output_db;
    }
    return $self->{_output_db};
}

1;
