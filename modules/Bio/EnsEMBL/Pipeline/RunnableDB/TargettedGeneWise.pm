#
# Ensembl module for Bio::EnsEMBL::Pipeline::RunnableDB::TargettedGeneWise.pm
#
# Cared for by Ensembl <ensembl-dev@ebi.ac.uk>
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::TargettedGeneWise.pm - Targetted genewise Runnable DB

=head1 SYNOPSIS

my $tgw = new Bio::EnsEMBL::Pipeline::RunnableDB::TargettedGeneWise
    (  -dbobj => $dbobj,
       -input_id => $input_id);

  $tgw->fetch_input;
  $tgw->run();
  $tgw->output();
  $tgw->write_output(); # write to db

=head1 DESCRIPTION

This object manages the data fetching, running, output parsing, and data storing of Targetted Genewise in the Ensembl pipeline.

=head1 CONTACT

Ensembl - ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Pipeline::RunnableDB::TargettedGeneWise;

use vars qw(@ISA);
use strict;
# Object preamble - inherits from Bio::Root::RootI
use Bio::Root::RootI;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Getseqs;
use Bio::SeqIO;
use Bio::EnsEMBL::Pipeline::GeneConf qw (
					 GB_TARGETTED_PROTEIN_INDEX
					 GB_TARGETTED_SINGLE_EXON_COVERAGE
					 GB_TARGETTED_MULTI_EXON_COVERAGE
					 GB_TARGETTED_MAX_INTRON
					 GB_TARGETTED_MIN_SPLIT_COVERAGE
					 GB_TARGETTED_GW_GENETYPE
					 );

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  # protein sequence fetcher
  if(!defined $self->seqfetcher) {
    my $seqfetcher = $self->make_seqfetcher($GB_TARGETTED_PROTEIN_INDEX);
    $self->seqfetcher($seqfetcher);
  }

  return $self;
}

=head2 make_seqfetcher

 Title   : make_seqfetcher
 Usage   :
 Function: get/set
 Example :
 Returns : Bio::DB::RandomAccessI
 Args    :


=cut

sub make_seqfetcher{
  my ( $self, $index ) = @_;
  my $seqfetcher;

  if(defined $index && $index ne ''){
    my @db = ( $index );
    $seqfetcher = new Bio::EnsEMBL::Pipeline::SeqFetcher::Getseqs(
								  '-db' => \@db,
								 );
  }
  else{
    sssself->throw("Can't make seqfetcher\n");
  }

  return $seqfetcher;
}

=head2 protein_id

 Title   : protein_id
 Usage   :
 Function: get/set
 Example :
 Returns : 
 Args    :


=cut

sub protein_id {
    my( $self, $value ) = @_;    
    if ($value) {
        $self->{'_protein_id'} = $value;
    }
    return $self->{'_protein_id'};
}

=head2 fetch_input

 Title   : fetch_input
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub fetch_input{
  my ($self,@args) = @_;

  my $entry = $self->input_id;
  my $chrname;
  my $start;
  my $end;
  my $protein_id; 

  # chr12:10602496,10603128:Q9UGV6:
  if( !($entry =~ /(\S+):(\d+),(\d+):(\S+):/)) {
      $self->throw("Not a valid input id... $entry");
  }
  
  $chrname    = $1;
  $protein_id = $4;
  $start   = $2;
  $end     = $3;

  if ($2 > $3) { # let blast sort it out
      $start  = $3;
      $end    = $2;
  }

  
  my $sgpa = $self->dbobj->get_StaticGoldenPathAdaptor();
  my $vc = $sgpa->fetch_VirtualContig_by_chr_start_end($chrname,$start-10000,$end+10000);
  
  $self->vc($vc);
  $self->protein_id($protein_id);
  
  # genewise runnable
  # repmasking?
  my $r = Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise->new( '-genomic'    => $vc->primary_seq,
								    '-ids'        => [ $protein_id ] ,
								    '-seqfetcher' => $self->seqfetcher);
 
  $self->runnable($r);

}


=head2 run

 Title   : run
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub run {
   my ($self,@args) = @_;
   
   $self->runnable->run();
   $self->convert_gw_output;

   # clean up tmpfile
   my $tmpfile = $self->{'_tmpfile'};
   unlink $tmpfile;

   # remap genes to raw contig coords
   my @remapped = $self->remap_genes();
   $self->output(@remapped);
}

=head2 output

 Title   : output
 Usage   :
 Function:
 Example :
 Returns : 
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



=head2 runnable

 Title   : runnable
 Usage   : $obj->runnable($newval)
 Function: 
 Returns : value of runnable
 Args    : newvalue (optional)


=cut

sub runnable{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_runnable'} = $value;
    }
    return $obj->{'_runnable'};
}

=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output
    Function:   Writes output data to db
    Returns :   
    Args    :   none

=cut


sub write_output {
  my($self) = @_;
  
  my $gene_adaptor = $self->dbobj->get_GeneAdaptor;
  
 GENE: foreach my $gene ($self->output) {	
    # do a per gene eval...
    eval {
      $gene_adaptor->store($gene);
      print STDERR "wrote gene dbID " . $gene->dbID . "\n";
    }; 
    if( $@ ) {
      print STDERR "UNABLE TO WRITE GENE\n\n$@\n\nSkipping this gene\n";
    }
  }
}

=head2 gw_genes

 Title   : gw_genes
 Usage   :
 Function: get/set for genewise gene array
 Example :
 Returns : 
 Args    :


=cut

sub gw_genes {
  my ($self, @genes) = @_;
  if (!defined($self->{'_gw_genes'})) {
    $self->{'_gw_genes'} = [];
  }

  if (scalar(@genes)) {
    push(@{$self->{'_gw_genes'}},@genes);
  }
  
  return @{$self->{'_gw_genes'}};
}

=head2 convert_gw_output

 Title   : convert_gw_output
 Usage   :
 Function: converts output from Genewise into genes
 Example :
 Returns : 
 Args    :


=cut

sub convert_gw_output {
  my ($self) = @_;
  my $count = 1;
  my $genetype = $GB_TARGETTED_GW_GENETYPE;
  if(!defined $genetype || $genetype eq ''){
    $genetype = 'TGE_gw';
    $self->warn("Setting genetype to $genetype\n");
  }
  my @results  = $self->runnable->output;

  # get the appropriate analysis from the AnalysisAdaptor
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
       -module          => 'TargettedGeneWise',
      );
  }

  my @genes = $self->make_genes($count, $genetype, $analysis_obj, \@results);

  # check for stops?

  $self->gw_genes(@genes);
  
}

=head2 make_genes

 Title   : make_genes
 Usage   : $self->make_genes($count, $genetype, $analysis_obj, $runnable)
 Function: converts the output from $runnable into Bio::EnsEMBL::Genes in
           $contig(VirtualContig) coordinates. The genes have type $genetype, 
           and have $analysis_obj attached. Each Gene has a single Transcript, 
           which in turn has Exons(with supporting features) and a Translation
 Example : 
 Returns : array of Bio::EnsEMBL::Gene
 Args    : $count: integer, $genetype: string, $analysis_obj: Bio::EnsEMBL::Analysis, 
           $runnable: Bio::EnsEMBL::Pipeline::RunnableI


=cut

sub make_genes {
  my ($self, $count, $genetype, $analysis_obj, $results) = @_;
  my $contig = $self->vc;
  my @genes;
  
  $self->throw("[$analysis_obj] is not a Bio::EnsEMBL::Analysis\n") 
    unless defined($analysis_obj) && $analysis_obj->isa("Bio::EnsEMBL::Analysis");

 MAKE_GENE:  foreach my $tmpf (@$results) {
    my $transcript = $self->make_transcript($tmpf,$self->vc,$genetype,$count, $analysis_obj);
	
    # validate transcript - validate_transcript returns an array ref
    my $valid_transcripts = $self->validate_transcript($transcript);
    next MAKE_GENE unless defined $valid_transcripts;
      
    my $gene;
    # make one gene per valid transcript
    foreach my $valid (@$valid_transcripts){
      $gene   = new Bio::EnsEMBL::Gene;
      $gene->type($genetype);
      $gene->analysis($analysis_obj);
      $gene->add_Transcript($valid);

      push(@genes,$gene);
    }

  }
  return @genes;
}

=head2 validate_transcript

 Title   : validate_transcript 
 Usage   : my @valid = $self->validate_transcript($transcript)
 Function: Validates a transcript - rejects if mixed strands, 
                                    rejects if low coverage, 
                                    splits if long introns and insufficient coverage of parental protein
                                    rejects unless exon coordinates are sane
 Returns : Ref to @Bio::EnsEMBL::Transcript
 Args    : Bio::EnsEMBL::Transcript

=cut
sub validate_transcript {
  my ($self, $transcript) = @_;
  
  my @valid_transcripts;
  
  my $valid = 1;
  my $split = 0;
  
  # check coverage of parent protein
  my $threshold = $GB_TARGETTED_SINGLE_EXON_COVERAGE; 
  my @exons = $transcript->get_all_Exons;
    if(scalar(@exons) > 1){
      $threshold = $GB_TARGETTED_MULTI_EXON_COVERAGE;
    }

  if(!defined $threshold){
    print STDERR "You must define GB_TARGETTED_SINGLE_EXON_COVERAGE and GB_TARGETTED_MULTI_EXON_COVERAGE in GeneConf.pm\n";
    return undef;
  }

  my $coverage  = $self->check_coverage($transcript);
  if ($coverage < $threshold){
    $self->warn ("Coverage of ". $self->protein_id . " is only $coverage - will be rejected\n");
    return undef;
  }
  
  print STDERR "Coverage of ". $self->protein_id . " is $coverage%\n";

  my $previous_exon;
  foreach my $exon($transcript->get_all_Exons){
    if(!$self->validate_exon($exon)){
      print STDERR "Rejecting gene because of invalid exon\n";
      return undef;
    }
       
    # check intron size
    if (defined($previous_exon)) {
      my $intron;
      
      if ($exon->strand == 1) {
	$intron = abs($exon->start - $previous_exon->end - 1);
      } else {
	$intron = abs($previous_exon->start - $exon->end - 1);
      }
      
#      if ($intron > 250000 && $coverage < 95) {

      if ($intron > $GB_TARGETTED_MAX_INTRON && $coverage < $GB_TARGETTED_MIN_SPLIT_COVERAGE ) {
	print STDERR "Intron too long $intron  for transcript " . $transcript->{'temporary_id'} . " with coverage $coverage\n";
	$split = 1;
	$valid = 0;
      }
      
      # check sensible strands
      if ($exon->strand != $previous_exon->strand) {
	print STDERR "Mixed strands for gene " . $transcript->{'temporary_id'} . "\n";
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

  if(scalar(@valid_transcripts)){
    return \@valid_transcripts;
  }
  else { 
    return undef;
  }
}

=head2 remap_genes

 Title   : remap_genes
 Usage   :
 Function: 
 Example :
 Returns : 
 Args    :


=cut

sub remap_genes {
  my ($self) = @_;
  my @newf;  
  my $contig = $self->vc;

  my @genes = $self->gw_genes;

GENE:  foreach my $gene (@genes) {

    my @t = $gene->each_Transcript;
    my $tran = $t[0];

    # check that it translates
    if($gene->type eq 'TGE_gw'){
      
      my $translates = $self->check_translation($tran);
      if(!$translates){
	my $msg = "discarding gene - tranlation has stop codons\n";
	$self->warn($msg);
	next GENE;
      }
  }
    eval {
      my $genetype = $gene->type;
      my $newgene = $contig->convert_Gene_to_raw_contig($gene);
      # need to explicitly add back genetype and analysis.
      $newgene->type($genetype);
      $newgene->analysis($gene->analysis);
      push(@newf,$newgene);

      # sort out supporting feature coordinates
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
    };

    # did we throw exceptions?
    if ($@) {
      print STDERR "Couldn't reverse map gene:  [$@]\n";
    }
  }

  return @newf;
}

=head2 check_translation

 Title   : check_translation
 Usage   :
 Function: 
 Example :
 Returns : 1 if transcript translates with no stops, otherwise 0
 Args    :


=cut

sub check_translation {
  my ($self, $transcript) = @_;
  my $tseq;

  eval{
    $tseq = $transcript->translate;
  };

  if((!defined $tseq) || ($@)){
    my $msg = "problem translating :\n$@\n";
    $self->warn($msg);
    return 0;
  }

  if ($tseq->seq =~ /\*/ ) {
    return 0;
  }
  else{
    return 1;
  }
}

=head2 check_coverage

 Title   : check_coverage
 Usage   :
 Function: checks how much of the parent protein is covered by the genewise prediction
 Example :
 Returns : %coverage of parent protein
 Args    :


=cut

sub check_coverage{
  my ($self, $transcript) = @_;
  my $pstart = 0;
  my $pend = 0;
  my $protname = $self->protein_id;
  my $plength;
  my $fetcher = new Bio::EnsEMBL::Pipeline::SeqFetcher;

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
    $self->throw("Error fetching sequence for [$protname]: [$@]\n");
  }
  
  $self->throw("No sequence fetched for [$protname]\n") unless defined $seq;
  
  $plength = $seq->length;

  if(!defined($plength) || $plength == 0){
    warn("no sensible length for $protname - can't get coverage\n");
    return 0;
  }

  my $coverage = $matches/$plength;
  $coverage *= 100;
  return $coverage;
}

=head2 make_transcript

 Title   : make_transcript
 Usage   : $self->make_transcript($gene, $contig, $genetype, $count, $analysis_obj)
 Function: makes a Bio::EnsEMBL::Transcript from a SeqFeature representing a gene, 
           with sub_SeqFeatures representing exons.
 Example :
 Returns : Bio::EnsEMBL::Transcript with Bio::EnsEMBL:Exons(with supporting feature 
           data), and a Bio::EnsEMBL::translation
 Args    : $gene: Bio::EnsEMBL::SeqFeatureI, $contig: Bio::EnsEMBL::DB::ContigI,
           $genetype: string, $count: integer
           $analysis_obj: Bio::EnsEMBL::Analysis


=cut

sub make_transcript{
  my ($self, $gene, $contig, $genetype, $count, $analysis_obj)=@_;
  $genetype = 'TGE_gw' unless defined ($genetype);
  $count = 1 unless defined ($count);

  unless ($gene->isa ("Bio::EnsEMBL::SeqFeatureI"))
    { print "$gene must be Bio::EnsEMBL::SeqFeatureI\n"; }
  unless ($contig->isa ("Bio::EnsEMBL::DB::ContigI"))
    { print "$contig must be Bio::EnsEMBL::DB::ContigI\n"; }

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
  } 
  else {
    
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

=head2 validate_exon

 Title   : validate_exon
 Usage   : $self->validate_exon($exon)
 Function: checks start and end coordinates of exon are sane
 Example : 
 Returns : 1 if exon is valid, otherwise zero
 Args    : $exon: Bio::EnsEMBL::Exon


=cut

sub validate_exon{
  my ($self, $exon) = @_;

  if($exon->start < 0){
    my $msg = "rejecting exon, start < 0 : " . $exon->start . "\n";
    $self->warn($msg);
    return 0;
  }

  elsif($exon->start > $exon->end){
    my $msg = "rejecting exon, start > end : " . $exon->start . " > " . $exon->end . "\n";
    $self->warn($msg);
    return 0;
  }

  elsif($exon->start == $exon->end){
    my $msg = "naughty exon, start == end : " . $exon->start . " == " . $exon->end . " - letting it through\n";
    $self->warn($msg);
    return 1;
  }
  
  return 1;
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
      $intron = abs($exon->start - $prev_exon->end + 1);
    } else {
      $intron = abs($exon->end   - $prev_exon->start + 1);
    }
    
    if ($intron > $GB_TARGETTED_MAX_INTRON) {
      $curr_transcript->translation->end_exon($prev_exon);

      # need to account for end_phase of $prev_exon when setting translation->end
      $curr_transcript->translation->end($prev_exon->end - $prev_exon->start + 1 - $prev_exon->end_phase);
      
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

      # start exon always has phase 0
      $exon->phase(0);

      # this new transcript becomes the current transcript
      $curr_transcript = $t;

      push(@split_transcripts, $curr_transcript);
    }

    if($exon == $transcript->end_exon){
      $curr_transcript->add_Exon($exon) unless $exon_added;
      $exon_added = 1;

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


1;
