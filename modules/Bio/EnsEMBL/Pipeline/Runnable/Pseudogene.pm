#
# You may distribute this module under the same terms as perl itself
#

=pod

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::Pseudogene

=head1 SYNOPSIS

    my $runname = "Bio::EnsEMBL::Pipeline::Runnable::Pseudogene";

    my $runnable = $runname->new
      ( '-query'  => $parameters{'-query'},
	'-max_intron_length' => $parameters{'-max_intron_length'}, #optional 
	'-max_intron_coverage' => $parameters{'-max_intron_coverage'},#optional 
	'-max_exon_coverage' => $parameters{'-max_exon_coverage'},#optional 
      );


=head1 DESCRIPTION

Pseudogene takes a Bio::EnsEMBL::Slice object and assesses genes and transcripts
for evidence of retrotransposition.
In the case of genes being identified as pseudogenes, the gene objects have their 
type set to pseudogene and all but the longest transcripts and translations are deleted.
If the gene has 1 or more pseudo transcripts but has other transcritps that are valid
the dubious transcripts are removed
The resulting gene objects are returned in an array.

=head1 CONTACT

Mail to B<ensembl-dev@ebi.ac.uk>

=cut

package Bio::EnsEMBL::Pipeline::Runnable::Pseudogene;

use strict;
use Carp qw(cluck);
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
  Returntype : Bio::EnsEMBL::Pipeline::Tools::Pseudogene
  Caller     : general

=cut


sub new {
  my ($class,@args) = @_;

  my $self = $class->SUPER::new(@args);

  #SET UP ANY INSTANCE VARIABLES

  $self->{'_max_intron_length'} = 50000; #default value for intron length cutoff
  $self->{'_max_intron_coverage'} = 80   ; #default value for intron length cutoff
  $self->{'_max_exon_coverage'} = 20   ; #default value for exon length cutoff
  $self->{'_modified_genes'} =[] ;# array ref to modified genes to write to new db
  $self->{'_discarded_transcripts'} = [];# array ref to discarded transcripts
  $self->{'_genes'} = []; #array of genes to test;
  $self->{'_repeats'} = {}; # hash of repeat blocks corresponding to each gene;

  my( $genes,$repeats,$max_intron_length, $max_intron_coverage, $max_exon_coverage) = $self->_rearrange([qw(
												   GENES
												   REPEAT_FEATURES
												   MAX_INTRON_LENGTH
												   MAX_INTRON_COVERAGE
												   MAX_EXON_COVERAGE
												  )], @args);
  
#  $self->_check_slice($slice) if ($slice);
  if ($max_intron_length){
    $self->max_intron_length($max_intron_length);
  }
  if ($max_intron_coverage){
    $self->max_intron_coverage($max_intron_coverage);
  }
  if ($max_exon_coverage){
    $self->max_exon_coverage($max_exon_coverage);
  }
  if ($genes){
    $self->genes($genes);
  }
  if ($repeats){
    $self->repeats($repeats);
  }

#test for same number of repeats and genes?

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


sub genes {
  my ($self, $genes) = @_;
  foreach my $gene (@{$genes}){
    unless  ($gene->isa("Bio::EnsEMBL::Gene")){
      $self->throw("Input isn't a Bio::EnsEMBL::Gene, it is a $gene");
    }
  }
  $self->{'_genes'} = $genes;
  return $self->{'_genes'};
}

sub repeats {
  my ($self, $repeats) = @_;
  foreach my $repeat_array (values %{$repeats}){
    foreach my $repeat (@{$repeat_array}){
      unless ($repeat->isa("Bio::EnsEMBL::SeqFeature")){
        $self->throw("Input is not a Bio::EnsEMBL::SeqFeature, it is a $repeat");
      }
    }
  }
  $self->{'_repeats'} = $repeats;
  return $self->{'_repeats'};
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
  print STDERR   $self->{'_real'}." real genes identified \n";
  print STDERR   $self->{'_pseudogenes'}." pseudogenes identified \n";  
  print STDERR   scalar(@{$self->{'_discarded_transcripts'}})." pseudotranscripts to be chucked \n";
  foreach my $transcript(@{$self->{'_discarded_transcripts'}}){
    print STDERR   $transcript->stable_id."\n";
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
  my @genes = @{$self->{'_genes'}};


  foreach my $gene(@genes){
    my @pseudo_trans ;
    my @real_trans ;
    foreach my $transcript (@{$gene->get_all_Transcripts}){
      $num++;
      my $evidence = $self->transcript_evidence($transcript,$gene);
      $pseudo = undef;
      #transcript tests

      #CALL PSEUDOGENE IF AT LEAST 80% COVERAGE OF INTRONS BY REPEATS
      #AT LEAST 1 F/S EXON AND 1 REAL EXON (?)
      #TOTAL INTRON LENGTH < 5K

      if($evidence->{'total_intron_len'} < $self->{'_max_intron_length'} &&
	 $evidence->{'frameshift_introns'} >= 1 &&
	 $evidence->{'real_introns'} >= 1 &&
	 $evidence->{'covered_introns'} >=  $self->{'_max_intron_coverage'} ){
	$pseudo = 1;
	print STDERR $transcript->stable_id." - repeats in introns\n";
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
      if ($gene->type eq 'pseudogene'){$gene->type('changed');}
    }

    ####################################
    # gene and transcripts are real real


    if (scalar(@pseudo_trans) == 0 && scalar(@real_trans) > 0){
      push (@{$self->{'_modified_genes'}},$gene); 
      $self->{'_real'}++;
      if ($gene->type eq 'pseudogene'){$gene->type('changed');}
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

  my ($self,$transcript,$gene) =@_;
  my $repeat_blocks = $self->{'_repeats'}->{$gene};;
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
    #in order to use range methods

    my $seq_feature_exon = Bio::EnsEMBL::SeqFeature->new(
							 -START => $exon->start-$gene->start,
							 -END => $exon->end-$gene->start,
							 -STRAND => $exon->strand
							);
    # Do intron
    if (defined($prev_exon)) {
      my $intron = Bio::EnsEMBL::SeqFeature->new(
						 -START => $prev_exon->end+1-$gene->start,
						 -END => $exon->start-1-$gene->start,
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


=head2 len_covered

  Args       : Bio::Seq::Feature object, reference to an array of repeat blocks
  Description: measures how much of the seq feature (intron or exon) is covered by repeat blocks
  Returntype : scalar

=cut 

sub _len_covered {
  my ($self,$feat,$repeat_blocks_ref) = @_;

  my $covered_len = 0;
 RBLOOP: foreach my $repeat_block (@$repeat_blocks_ref) {
# print STDERR  "RB " . $repeat_block->start . " " . $repeat_block->end . "\n";
# print STDERR  "FT " . $feat->start . " " . $feat->end . "\n";
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
  Description: returns output of running Pseudogene
  Returntype : @{Bio::EnsEMBL::Gene}
  Exceptions : none
  Caller     : general

=cut

sub output {
    my ($self) = @_;
    return $self->{'_modified_genes'};
}



return 1;
