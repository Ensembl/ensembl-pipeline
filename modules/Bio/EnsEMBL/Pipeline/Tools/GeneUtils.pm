
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::GeneUtils - 

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 CONTACT

ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::Pipeline::Tools::GeneUtils;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Pipeline::Runnable::Protein::Seg;
use Bio::EnsEMBL::DnaPepAlignFeature;
use Bio::EnsEMBL::Pipeline::Tools::ExonUtils;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;

@ISA = qw(Bio::EnsEMBL::Root);

=head2 validate_Transcript

  Arg[1]    : Bio::EnsEMBL::Transcript $transcript
  Arg[2]    : Bio::EnsEMBML::Slice $slice
  Arg[3]    : int $coverage_threshold 
  Arg[4]    : ref to array of Bio::DB::RandomAccessI $seqfetchers
  Arg[5]    : int $max_intron
  Arg[6]    : int $min_split_coverage
  Arg[7]    : int $low_complexity
  Function  : evaluates $transcript - rejects if mixed strands, rejects if low coverage, 
              rejects if stops in translation, splits if long introns and insufficient 
              coverage of parental protein
  ReturnType: undef/ref to array of Bio::EnsEMBL::Transcript
  Exceptions: warns if transcript is to be rejected
  Caller    :
  Example   :

=cut

sub validate_Transcript {
  my ($self, $transcript, $slice, $coverage_threshold, $maxintron, $min_split_coverage, $seqfetchers, $low_complexity ) = @_;
  
  my @valid_transcripts;

  return undef unless Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_check_Transcript ($transcript, $slice);

  return undef unless Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_check_Translation($transcript);

  if(defined $low_complexity) {
    return undef unless Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_check_low_complexity($transcript, $low_complexity);
  }
  my $coverage     =  Bio::EnsEMBL::Pipeline::Tools::GeneUtils->_check_coverage($transcript, $coverage_threshold, $seqfetchers);

  return undef unless defined $coverage;

  # Do we really need to do this? Can we just take out the dbID adaptor stuff
  my $newtranscript  = Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_clone_Transcript($transcript);

  # split transcript if necessary - at present in practice all Similarity transcripts 
  # will be run through split_transcripts but very good Targetted predictions having 
  # long introns will not be
  if($coverage < $min_split_coverage){
    my $split_transcripts = Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->split_Transcript($newtranscript, $maxintron);
    push (@valid_transcripts, @$split_transcripts);
  }
  else{
    push (@valid_transcripts, $newtranscript);
  }

  return \@valid_transcripts;
}

=head2 _check_coverage

  Arg[1]    : Bio::EnsEMBL::Transcript $transcript
  Arg[2]    : int $coverage_threshold 
  Arg[3]    : ref to array of Bio::DB::RandomAccessI $seqfetchers
  Function  : calculates how much of the parent protein is covered by the predicted translation 
              and compares to a passed in coverage threshold.
  ReturnType: undef calculated coverage less than $coverage_threshold, otherwise int (coverage)
  Exceptions: warns if problems fetching sequence
  Caller    :
  Example   :

=cut

sub _check_coverage {
  my ( $self, $transcript, $coverage_threshold, $seqfetchers) = @_;
  
  my $matches = 0;
  my $pstart  = 0;
  my $pend    = 0;
  my $protname;
  my $plength;
    
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
 SEQFETCHER:
  foreach my $seqfetcher (@$seqfetchers){
    eval{
      $seq = $seqfetcher->get_Seq_by_acc($protname);
    };
    if ($@) {
      warn("GeneUtils:Error fetching sequence for [$protname] - trying next seqfetcher:[$@]\n");
    }
    
    if (defined $seq) {
      last SEQFETCHER;
    }
    
  }
  
  if(!defined $seq){
    warn("GeneUtils: No sequence fetched for [$protname] - can't check coverage, letting gene through\n");
    return 100;
  }
  
  $plength = $seq->length;
  
  if(!defined($plength) || $plength == 0){
    warn("GeneUtils: no sensible length for $protname - can't get coverage\n");
    return 0;
  }
  
  my $realcoverage = int(100 * $matches/$plength);
  
  if ($realcoverage < $coverage_threshold){
    warn("GeneUtils: Rejecting transcript for low coverage: $realcoverage\n");
    return undef;
  }
  
  return $realcoverage;
  
}


=head2 SeqFeature_to_Transcript

  Arg[1]    : Bio::EnsEMBL::SeqfeatureI $gene 
  Arg[2]    : Bio::EnsEMBL::RawContig $contig
  Arg[3]    : string $genetype 
  Arg[4]    : Bio::EnsEMBL::Analysis $analysis_obj
  Function  : makes a Bio::EnsEMBL::Transcript from a SeqFeature representing a gene, 
              with sub_SeqFeatures representing exons.
  ReturnType: Bio::EnsEMBL::Transcript
  Exceptions: none
  Caller    :
  Example   :

=cut
	    
sub SeqFeature_to_Transcript {
  my ($self, $gene, $contig, $analysis_obj, $db, $phase) = @_;
  unless ($gene->isa ("Bio::EnsEMBL::SeqFeatureI")){
    print "$gene must be Bio::EnsEMBL::SeqFeatureI\n";
  }

  my $curr_phase = 0;  

  my $excount = 1;
  my @exons;

  my $transcript   = new Bio::EnsEMBL::Transcript;
  my $translation  = new Bio::EnsEMBL::Translation;    
  $transcript->translation($translation);
  
  
  my $ea;
  if ($db ne "" && defined($db)) {
    $ea = $db->get_ExonAdaptor;
  }
  
  my @pred = $gene->sub_SeqFeature;
  if ($pred[0]->strand ==1 ) {
    @pred = sort {$a->start <=> $b->start} @pred;
  } else {
    @pred = sort {$b->start <=> $a->start} @pred;
  }

  foreach my $exon_pred (@pred) {
    my $exon = new Bio::EnsEMBL::Exon;
    
    $exon->start    ($exon_pred->start);
    $exon->end      ($exon_pred->end);
    $exon->strand   ($exon_pred->strand);
    
    if ($phase) {
      $exon->phase($curr_phase);# or $phase?
      my $exon_length =  $exon->end - $exon->start + 1;
      my $end_phase   = ( $exon_length + $exon->phase ) %3;
      
      $exon->end_phase($end_phase);
      
      $curr_phase = $end_phase;
    } else {
      $exon->phase    ($exon_pred->phase);
      $exon->end_phase($exon_pred->end_phase);
    }
    
    $exon->contig   ($contig); 
    $exon->adaptor  ($ea);
    
    # sort out supporting evidence for this exon prediction
    # Now these should be made in the correct type to start with
    
    my @sf = $exon_pred->sub_SeqFeature;
    my $prot_adp;
    
    if ($db ne "" && defined($db)) {
      $prot_adp = $db->get_ProteinAlignFeatureAdaptor;    
    }
    
    if(@sf){
      my $align = new Bio::EnsEMBL::DnaPepAlignFeature(-features => \@sf); 
      
      $align->seqname($contig->dbID);
      $align->contig($contig);
      $align->adaptor($prot_adp);
      $align->score(100); # Hmm!!!!
      $align->analysis($analysis_obj);
      
      $exon->add_supporting_features($align);
      
    }
    
    push(@exons,$exon);
    
    $excount++;
  }
  
  if ($#exons < 0) {
    print STDERR "Odd.  No exons found\n";
    return;
  } else {
    
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


=head2 prune_Exons

  Arg [1]   : Bio::EnsEMBL::Gene
  Function  : removes redundancy from a genes exon set 
  Returntype: Bio::EnsEMBL::Gene
  Exceptions: none 
  Caller    : 
  Example   : 

=cut


sub prune_Exons {
  my ($self,$gene) = @_;
  
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


=head2 _validate_Gene

  Arg [1]   : Bio::EnsEMBL::Gene
  Function  : checks sanity of Gene coords
  Returntype: 1/0
  Exceptions: warns if things aren't valid and return 0'
  Caller    : 
  Example   : 

=cut

sub _validate_gene{
  my ($self, $gene, $slice) = @_;

  foreach my $transcript(@{$gene->get_all_Transcripts}){
    # validate transcript
    if(!Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_check_Transcript($transcript, $slice)){
      $self->warn("Rejecting gene because of invalid transcript\n");
      return 0;
    }    
    foreach my $exon(@{$transcript->get_all_Exons}){
      if(!Bio::EnsEMBL::Pipeline::Tools::ExonUtils->validate_exon($exon)){
	$self->warn("Rejecting gene because of invalid exon\n");
	return 0;
      }
    }
  }
  
  return 1;
}

sub _print_Gene{
  my ($self,$gene) = @_;
   my $id;
  if ($gene->stable_id){
    $id = $gene->stable_id;
  }
  elsif ( $gene->dbID ){
    $id = $gene->dbID;
  }
  else{
    $id = "no id";
  }
  if ( defined( $gene->type ) ){
    $id .= " ".$gene->type;
  }
  foreach my $transcript (@{$gene->get_all_Transcripts} ){
    Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($transcript);
  }
}

############################################################

sub _clone_Gene{
  my ($self,$gene) = @_;
  
  my $newgene = new Bio::EnsEMBL::Gene;
  if ($gene->type){
      $newgene->type( $gene->type);
  }
  if ( defined $gene->dbID ){
      $newgene->dbID($gene->dbID);
  }
  if ( defined $gene->analysis ){
      $newgene->analysis($gene->analysis);
  }
  foreach my $transcript (@{$gene->get_all_Transcripts}){
    my $newtranscript = Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_clone_Transcript($transcript);
    $newgene->add_Transcript($newtranscript);
  }
  return $newgene;
}

############################################################


1;

