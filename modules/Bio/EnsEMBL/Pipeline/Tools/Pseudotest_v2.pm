#
# You may distribute this module under the same terms as perl itself
#

=pod

=head1 NAME

Bio::EnsEMBL::Pipeline::Tools::Pseudotest_v2

=head1 SYNOPSIS

    my $runname = "Bio::EnsEMBL::Pipeline::Tools::Pseudotest_v2";

    my $runnable = $runname->new
      ( '-query'  => $parameters{'-query'},
	'-max_intron_length' => $parameters{'-max_intron_length'}, #optional 
	'-max_intron_coverage' => $parameters{'-max_intron_coverage'},#optional 
	'-max_exon_coverage' => $parameters{'-max_exon_coverage'},#optional 
      );


=head1 DESCRIPTION

Pseudotest_v2 takes a Bio::EnsEMBL::Slice object and assesses genes and transcripts
for evidence of retrotransposition.
In the case of genes being identified as pseudogenes, the gene objects have their 
type set to pseudogene and all but the longest transcripts and translations are deleted.
If the gene has 1 or more pseudo transcripts but has other transcritps that are valid
the dubious transcripts are removed
The resulting gene objects are returned in an array.

=head1 CONTACT

Mail to B<ensembl-dev@ebi.ac.uk>

=cut

package Bio::EnsEMBL::Pipeline::Tools::Pseudotest_v2;

use strict;
use Carp;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::Root;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

# Object preamble - inherits from Bio::EnsEMBL::Root;



=head2 new

  Args       : various
  Description: Runnable constructor
  Returntype : Bio::EnsEMBL::Pipeline::Tools::Pseudotest_v2
  Caller     : general

=cut


sub new {
  my ($class,@args) = @_;

  my $self = $class->SUPER::new(@args);

  #SET UP ANY INSTANCE VARIABLES

  $self->{'_slice'} = undef; #Slice to run over
  $self->{'_max_intron_length'} = 50000; #default value for intron length cutoff
  $self->{'_max_intron_coverage'} = 80   ; #default value for intron length cutoff
  $self->{'_max_exon_coverage'} = 20   ; #default value for exon length cutoff
  $self->{'_modified_genes'} =[] ;# array ref to modified genes to write to new db
  $self->{'_discarded_transcripts'} = [];# array ref to discarded transcripts

  my( $slice,$max_intron_length, $max_intron_coverage, $max_exon_coverage) = $self->_rearrange([qw(
												   QUERY
												   MAX_INTRON_LENGTH
												   MAX_INTRON_COVERAGE
												   MAX_EXON_COVERAGE
												  )], @args);
  
  $self->_check_slice($slice) if ($slice);
  if ($max_intron_length){
    $self->max_intron_length($max_intron_length);
  }
  if ($max_intron_coverage){
    $self->max_intron_coverage($max_intron_coverage);
  }
  if ($max_exon_coverage){
    $self->max_exon_coverage($max_exon_coverage);
  }


  return $self;
}

=head2 max_intron_length

  Arg [1]    : scalar
  Description: get/set maximium allowed total intron length
  Returntype : scalar
  Exceptions : none
  Caller     : general

=cut

sub max_intron_length {
  my ($self, $max_intron_length) = @_;

  if ($max_intron_length > 0 ){
    $self->{'_max_intron_length'} = $max_intron_length;
  }
  return $self->{'_max_inton_length'};
}

=head2 max_intron_coverage

  Arg [1]    : scalar
  Description: get/set maximium allowed %coverage of introns with repeats (default 80%)
  Returntype : scalar
  Exceptions : none
  Caller     : general

=cut

sub max_intron_coverage {
  my ($self, $max_intron_coverage) = @_;

  if ($max_intron_coverage > 0 && $max_intron_coverage < 100  ){
    $self->{'_max_intron_coverage'} = $max_intron_coverage;
   }
  return $self->{'_max_intron_coverage'};
}

=head2 max_exon_coverage

  Arg [1]    : scalar
  Description: get/set maximium allowed %coverage of exons with repeats (default 20%)
  Returntype : scalar
  Exceptions : none
  Caller     : general

=cut

sub max_exon_coverage {
  my ($self, $max_exon_coverage) = @_;

  if ($max_exon_coverage > 0 && $max_exon_coverage < 100 ){
    $self->{'_max_exon_coverage'} = $max_exon_coverage;
     }
  return $self->{'_max_exon_coverage'};
}

=head2 check_slice

  Arg [1]    : Bio::EnsEMBL::Slice $query
  Description: accessor for query sequence
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : query not a Bio::EnsEMBL::Slice
  Caller     : general

=cut

sub _check_slice {
  my ($self, $slice) = @_;

  if ($slice) {
    unless ($slice->isa("Bio::EnsEMBL::Slice")) {
      $self->throw("Input isn't a Bio::EnsEMBL::Slice");
    }  
    $self->{_slice} = $slice ;
  }
  return $self->{_slice};
}

=head2 run

  Arg [none] :
  Description: runs the  runnable
  Returntype : none
  Exceptions : none
  Caller     : general

=cut

sub run {
  my ($self) = @_;
  #check sequence
  my $seq = $self->_check_slice || $self->throw("Clone required for dust\n");

  $self->test_genes;
  $self->summary;
  return 1;
}

=head2 summary

  Arg [none] :
  Description: prints out some data about the results
  Returntype : none
  Exceptions : none
  Caller     : general

=cut

sub summary {
  my ($self) = @_;
  print STDERR  $self->{'_real'}." real genes identified \n";
  print STDERR  $self->{'_pseudogenes'}." pseudogenes identified \n";  
  print STDERR  scalar(@{$self->{'_discarded_transcripts'}})." pseudotranscripts to be chucked \n";
  foreach my $transcript(@{$self->{'_discarded_transcripts'}}){
    print STDERR  $transcript->stable_id."\n";
  }
  return 1;
}


=head2 test_genes

  Arg [none] :
  Description: Check genes one transcript at at a time pushes test result for each 
  transcript onto an array, tests array to make final decision
  genes are classed as pseudogenes if the following criteria are met:

1. At least 80% of the introns are covered with repeats and the total intron length 
is smaller than 5kb (deFAULT VALUES).

2. All of the introns are short frameshifted introns
  Returntype : none
  Exceptions : none
  Caller     : general

=cut

sub test_genes{
  my $self = shift;
  my @evidence;
  my $num=0;
  my $pseudo= undef;
  my @genes = @{$self->{'_slice'}->get_all_Genes};


  foreach my $gene(@genes){
    my @pseudo_trans ;
    my @real_trans ;
    foreach my $transcript (@{$gene->get_all_Transcripts}){
      $num++;
      my $evidence = $self->transcript_evidence($transcript);
      $pseudo = undef;
      #transcript tests

      #CALL PSEUDOGENE IF AT LEAST 80% COVERAGE OF INTRONS BY REPEATS
      #TOTAL INTRON LENGTH < 5K

      if($evidence->{'total_intron_len'} < $self->{'_max_intron_length'} &&
	 #	 $evidence->{'frameshift_introns'} >= 1 &&
	 #	 $evidence->{'real_introns'} >= 1 &&
	 $evidence->{'covered_introns'} >=  $self->{'_max_intron_coverage'} ){
	$pseudo = 1;
      }
      #ALL FRAMESHIFTED
      
      if($evidence->{'num_introns'} && $evidence->{'frameshift_introns'} == $evidence->{'num_introns'}){$pseudo = 1;}
      
      #EXONS CONTAMINATED
      
      #    if($evidence->{'covered_exons'} >= $self->{'_max_exon_coverage'}){$pseudo = 1;}
      
      
      if ($pseudo){
	push (@pseudo_trans,$transcript);	
      }
      
      else{
	push (@real_trans,$transcript);	
      }
    }
    
    #########################################
    # gene tests

    #############################################
    # gene is pseudogene, set type to pseudogene
    # chuck away all but the longest transcript
    
    
    if (scalar(@pseudo_trans) > 0 && scalar(@real_trans) == 0){
      $gene->type('pseudogene');
      @pseudo_trans = sort {$a->length <=> $b->length} @pseudo_trans;
      my $only_transcript_to_keep = pop  @pseudo_trans;
      foreach my $pseudo_transcript (@pseudo_trans){
	push @{$self->{'_discarded_transcripts'}},$pseudo_transcript;
	$pseudo_transcript->translation(undef);
	$self->_remove_transcript_from_gene($gene,$pseudo_transcript);
      }
      push@{$self->{'_modified_genes'}}, $gene;   
      $self->{'_pseudogenes'}++;
    }


    ###############################################
    # gene is real but has some dodgy transcripts
    # delete the dodgy transcripts from the gene


    if (scalar(@pseudo_trans) > 0 && scalar(@real_trans) > 0){
      foreach my $trans (@pseudo_trans){
	$trans->translation(undef);
	$self->_remove_transcript_from_gene($gene,$trans); 
      }
      push @{$self->{'_modified_genes'}}, $gene;
      push @{$self->{'_discarded_transcripts'}},@pseudo_trans;
      $self->{'_real'}++;
    }

    ####################################
    # gene and transcripts are real real


    if (scalar(@pseudo_trans) == 0 && scalar(@real_trans) > 0){
      push (@{$self->{'_modified_genes'}},$gene); 
      $self->{'_real'}++;
    }
  }

return 1;
}

=head2 transcript_evidence

Arg [none] : Bio::EnsEMBL::Transcript
  Description: Test individual transcripts return a hash containing results evidence
  Returntype : hash
  Exceptions : none
  Caller     : general

=cut

##############################################################
#test individual transcripts return a hash containing evidence
##############################################################

sub transcript_evidence{

  my ($self,$transcript) =@_;
  my $repeat_blocks = $self->get_all_repeat_blocks($transcript);
  my $results;
  my  @exons =  @{$transcript->get_all_Exons};
  @exons = sort {$a->start <=> $b->start} @exons;
  my $prev_exon = undef;
  my $total_intron_len = 0;
  my $covered_intron_len = 0;
  my $total_exon_len = 0;
  my $covered_exon_len = 0;
  my $n_real_intron = 0;
  my $n_frameshift_intron = 0;
  my $covered_introns = 0 ;
  my $covered_exons = 0 ;
  
  foreach my $exon (@exons) {
    ###########################################
    #Need to convert exon object to seq feature
    #in order to use ranage methods

    my $seq_feature_exon = Bio::EnsEMBL::SeqFeature->new(
							 -START => $exon->start-$transcript->start ,
							 -END => $exon->end-$transcript->start,
							 -STRAND => $exon->strand
							);
    # Do intron
    if (defined($prev_exon)) {
      my $intron = Bio::EnsEMBL::SeqFeature->new(
						 -START => $prev_exon->end+1-$transcript->start,
						 -END => $exon->start-1-$transcript->start,
						 -STRAND => $exon->strand
						);
      if ($intron->length > 9) {
	$n_real_intron++;
      } else {
	$n_frameshift_intron++;
      }
      $total_intron_len+=$intron->length;
      $covered_intron_len+=$self->_len_covered($intron,$repeat_blocks);
    }
 
    # Do exon
    $total_exon_len+=$exon->length;
    $covered_exon_len+=$self->_len_covered($seq_feature_exon,$repeat_blocks);
    $prev_exon = $exon;
  }

  #calculate percentage coverage
  
  if ($total_intron_len > 0){
    $covered_introns = (($covered_intron_len/$total_intron_len)*100);
  }
  if ($total_exon_len >  0){
    $covered_exons =   (($covered_exon_len/$total_exon_len)*100)
  }
  
  $results = {
	      'num_introns' => $#exons,
	      'total_exon_len' =>  $total_exon_len,
	      'covered_exons' => $covered_exons,
	      'total_intron_len' =>  $total_intron_len,
	      'covered_introns' => $covered_introns,
	      'real_introns' => $n_real_intron,
	      'frameshift_introns' => $n_frameshift_intron
	     };
  return $results;
}


=head2 get_all_repeat_blocks

  Args       : none
  Description: gets slice from transcript and fetches all repeats and  
merges them into blocks
  Returntype : array of Seq_Feature blocks;

=cut 

sub get_all_repeat_blocks {
  my ($self,$transcript) = @_;
  my @repeat_blocks;
  my $rep_gene_slice = $transcript->feature_Slice;
  my @repeats = @{$rep_gene_slice->get_all_RepeatFeatures};
  @repeats = sort {$a->start <=> $b->start} @repeats;
  my $curblock = undef;

  REPLOOP: foreach my $repeat (@repeats) {
      my $rc = $repeat->repeat_consensus;
      if ($rc->repeat_class !~ /LINE/ && $rc->repeat_class !~ /LTR/ && $rc->repeat_class !~ /SINE/) { 
	next REPLOOP;
      }
      if ($repeat->start <= 0) { 
	$repeat->start(1); 
      }
      if (defined($curblock) && $curblock->end >= $repeat->start) {
	if ($repeat->end > $curblock->end) { 
	  $curblock->end($repeat->end); 
	}
      }
      else {
	$curblock = Bio::EnsEMBL::SeqFeature->new(
						  -START => $repeat->start,
						  -END => $repeat->end, 
						  -STRAND => $repeat->strand
						 );
	push (@repeat_blocks,$curblock);
      }
    }
  @repeat_blocks = sort {$a->start <=> $b->start} @repeat_blocks;
  return\@repeat_blocks;
}

=head2 len_covered

  Args       : Bio::Seq::Feature object, reference to an array of repeat blocks
  Description: measures how much of the seq feature (intron or exon) is covered by repeat blocks
  Returntype : scalar

=cut 

sub _len_covered {
  my ($self,$feat,$repeat_blocks_ref) = @_;

  my $covered_len = 0;
 RBLOOP: foreach my $repeat_block (@$repeat_blocks_ref) {
    if ($repeat_block->overlaps($feat, 'ignore')) {
      my $inter = $feat->intersection($repeat_block);
      $covered_len += $inter->length;
    } elsif ($repeat_block->start > $feat->end) {
      last RBLOOP;
    }
  }
  return  $covered_len;
}

=head2 remove_transcript_from_gene

  Args : Bio::EnsEMBL::Gene object , Bio::EnsEMBL::Transcript object
  Description: steves method for removing unwanted transcripts from genes
  Returntype : scalar

=cut 


sub _remove_transcript_from_gene {
  my ($self, $gene, $trans_to_del)  = @_;

  my @newtrans;
  foreach my $trans (@{$gene->get_all_Transcripts}) {
    if ($trans != $trans_to_del) {
      push @newtrans,$trans;
    }
  }

# The naughty bit!
  $gene->{_transcript_array} = [];

  foreach my $trans (@newtrans) {
    $gene->add_Transcript($trans);
  }

  return scalar(@newtrans);
}

=head2 output

  Arg [none] :
  Description: returns output of running Pseudotest_v2
  Returntype : @{Bio::EnsEMBL::Gene}
  Exceptions : none
  Caller     : general

=cut

sub output {
    my ($self) = @_;

    return @{$self->{'_modified_genes'}};
}

return 1;
