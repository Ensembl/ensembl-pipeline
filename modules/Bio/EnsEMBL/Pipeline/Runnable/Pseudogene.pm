#
# You may distribute this module under the same terms as perl itself
#

=pod

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::Pseudogene

=head1 SYNOPSIS

 my $runnable = Bio::EnsEMBL::Pipeline::Runnable::Pseudogene->new
      ( 
       '-genes' => \@_genes,
       '-repeat_features' => \%repeat_blocks, 
       '-homologs'        => \%homolog_hash,
      );
    $runnable->run;
    $runnable->output;

Where output returns an array of modified genes.
Repeat blocks is a hash of repeats covering each gene merged into blocks.
Homolog hash is a hash of transcript objects that are homologous the gene of interest

=head1 DESCRIPTION

Runnable for PseudogeneDB

Runs tests to identiy pseudogenes:
Calls it a pseudo gene if:

1. All of the introns are frameshifted.
2. Real introns are covered with repeats
3. Single exon and there is a copyelsewhere in the genome that is spliced

Pseudogene takes a Bio::EnsEMBL::Slice object and assesses genes and transcripts for evidence of retrotransposition.
In the case of genes being identified as pseudogenes, the gene objects have their type set to pseudogene and all but the longest transcripts and translations are deleted.
If the gene has 1 or more pseudo transcripts but has other transcritps that are valid the dubious transcripts are removed The resulting gene objects are returned in an array.

=head1 CONTACT

Mail to B<ensembl-dev@ebi.ac.uk>

=cut

package Bio::EnsEMBL::Pipeline::Runnable::Pseudogene;

use strict;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::Pipeline::Runnable::Bl2seq;
use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Pseudogene_config;

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


  $self->{'_modified_genes'} =[] ; # array ref to modified genes to write to new db
  $self->{'_discarded_transcripts'} = []; # array ref to discarded transcripts
  $self->{'_genes'} = [];	#array of genes to test;  
  $self->{'_homologs'} = {};	#2D hash ref of transcript homologs corresponding to single exon genes;
  $self->{'_repeats'} = {};	# hash of repeat blocks corresponding to each gene;
  $self->{'_real'} = 0;		# scalar number of real genes identified
  $self->{'_pseudogenes'} = 0;	#scalar number of pseudogenes identified
  
  my( $genes,$repeats,$homologs) = $self->_rearrange([qw(
							 GENES
							 REPEAT_FEATURES
							 HOMOLOGS
							)], @args);

  if ($genes) {
    $self->genes($genes);
  }
  if ($repeats) {
    $self->repeats($repeats);
  }
  if ($homologs) {
    $self->homologs($homologs);
  }

  #test for same number of repeats and genes?
  print "config things $PS_SPAN_RATIO, $PS_MIN_EXONS, $PS_PERCENT_ID_CUTOFF\n";
  print "more things $PS_MAX_INTRON_LENGTH, $PS_NUM_FRAMESHIFT_INTRONS, $PS_NUM_REAL_INTRONS, $PS_MAX_INTRON_COVERAGE\n";
  return $self;
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
  print STDERR   $self->real." real genes identified \n";
  print STDERR   $self->pseudogenes." pseudogenes identified \n";  
  print STDERR   scalar(@{$self->discarded_transcripts})." pseudotranscripts to be chucked \n";
  foreach my $transcript (@{$self->discarded_transcripts}) {
    print STDERR   $transcript->stable_id."\n";
  }
  return 1;
}


=head2 test_genes

Arg [none] :
  Description: Check genes one transcript at at a time pushes test result for each transcript onto an array, tests array to make final decision  genes are classed as pseudogenes if the following criteria are met:
1. At least 80% of the introns are covered with repeats and the total intron length is smaller than 5kb and the gene has a least 1 real and 1 frameshifted intron (default values).
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
  my @genes = @{$self->genes};
 
  foreach my $gene (@genes) {
    my @pseudo_trans ;
    my @real_trans ;

    ###################################################
    # if gene is single exon look for spliced elsewhere

    if (scalar(@{$gene->get_all_Exons()})==1){
	my $judgement = $self->spliced_elsewhere($gene);
	if ($judgement eq 'pseudogene') {
	  push (@pseudo_trans,@{$gene->get_all_Transcripts}[0]);
	} else {
	  push (@real_trans,@{$gene->get_all_Transcripts}[0]);
	}
      }
    else {

    ####################################
    # gene is multiexon, run other tests

	foreach my $transcript (@{$gene->get_all_Transcripts}) {
	  $num++;
	  my $evidence = $self->transcript_evidence($transcript,$gene);
	  
	$pseudo = undef;
	#transcript tests

	#CALL PSEUDOGENE IF AT LEAST 80% COVERAGE OF INTRONS BY REPEATS
	#AT LEAST 1 F/S EXON AND 1 REAL EXON (?)
	#TOTAL INTRON LENGTH < 5K
	


	if ($evidence->{'total_intron_len'} < $PS_MAX_INTRON_LENGTH &&
	    $evidence->{'frameshift_introns'} >= $PS_NUM_FRAMESHIFT_INTRONS &&
	    $evidence->{'real_introns'} >= $PS_NUM_REAL_INTRONS &&
	    $evidence->{'covered_introns'} >= $PS_MAX_INTRON_COVERAGE  ) {
	  $pseudo = 1;
	  print STDERR $gene->stable_id." - repeats in introns in transcript ".$transcript->stable_id."\n";
	  print STDERR join (', ',%{$evidence}),"\n"
	}
	#ALL FRAMESHIFTED
	
	if ($evidence->{'num_introns'} && $evidence->{'frameshift_introns'} == $evidence->{'num_introns'}) {
	  $pseudo = 1;
	}
	
	#EXONS CONTAMINATED
	
	#    if($evidence->{'covered_exons'} >= $PS_MAX_EXON_COVERAGE){$pseudo = 1;}
	
	
	if ($pseudo) {
	  push (@pseudo_trans,$transcript);	
	} else {
	  push (@real_trans,$transcript);	
	}
      }
    }
    
    #########################################
    # gene tests
      
    #############################################
    # gene is pseudogene, set type to pseudogene
    # chuck away all but the longest transcript
      
    if (scalar(@pseudo_trans) > 0 && scalar(@real_trans) == 0) {
      $gene->type('pseudogene');
      @pseudo_trans = sort {$a->length <=> $b->length} @pseudo_trans;
      my $only_transcript_to_keep = pop  @pseudo_trans;
      $only_transcript_to_keep->translation(undef);
      foreach my $pseudo_transcript (@pseudo_trans) {
	$self->discarded_transcripts($pseudo_transcript);
	$pseudo_transcript->translation(undef);
	$self->_remove_transcript_from_gene($gene,$pseudo_transcript);
      }
      $self->modified_genes($gene);
      $self->pseudogenes(1);
    }


    ###############################################
    # gene is real but has some dodgy transcripts
    # delete the dodgy transcripts from the gene


    if (scalar(@pseudo_trans) > 0 && scalar(@real_trans) > 0) {
      foreach my $trans (@pseudo_trans) {
	$trans->translation(undef);
	$self->_remove_transcript_from_gene($gene,$trans); 
      }
      $self->modified_genes($gene);
      $self->discarded_transcripts(@pseudo_trans);
      $self->real(1);
    }

    ####################################
    # gene and transcripts are real real


    if (scalar(@pseudo_trans) == 0 && scalar(@real_trans) > 0) {
      $self->modified_genes($gene); 
      $self->real(1);
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
  my $repeat_blocks = $self->get_repeats($gene);
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
  
  if ($total_intron_len > 0) {
    $covered_introns = (($covered_intron_len/$total_intron_len)*100);
  }
  if ($total_exon_len >  0) {
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

##############################################################
# Spliced elsewhere - tests for retrotransposition
##############################################################

=head2 spliced_elsewhere

Args : Bio::EnsEMBL::Gene object 
  Description: fetches all homologous transcripts for a given single exon  gene and runs a tblastx of the two transcripts, makes a judgement about retrotransposition   based on the stats returned. Looks at top scoring hit of each transcript homolog of the given gene and compares %ID and spans of the transcripts  looks for sequences > 80% id with a ratio of real transcript span over processed  transcript span of > 1.5
  Returntype : scalar

=cut 

sub spliced_elsewhere {
  my ($self,$gene) = @_;
  my $judgement = 'undecided';

  # If there is evidence to suggest that it may be spliced elsewhere

  if ($self->get_homologs_by_gene($gene)){

    my %homolog_hash = %{$self->get_homologs_by_gene($gene)};
    my $workdir = $self->workdir;
    my $seq1 =   @{$gene->get_all_Transcripts}[0]->seq;
    my $length = @{$gene->get_all_Transcripts}[0]->length;
    my $ratio;
    my @results;
    
    foreach my $homogene (keys %homolog_hash){
      my $seq2 = $homolog_hash{$homogene}->seq;
      my $span = $homolog_hash{$homogene}->end-$homolog_hash{$homogene}->start;
      my $n_exons = scalar(@{$homolog_hash{$homogene}->get_all_Exons});
      my $obj = Bio::EnsEMBL::Pipeline::Runnable::Bl2seq->new(
							      -seq1    => $seq1,
							      -seq2    => $seq2,
							      -alntype => 'tblastx',
							      -workdir => $workdir,
							     );
      $obj->run;
      my @output = $obj->output();
      @output = sort {$a->score <=> $b->score} @output;
      my $result = $output[0];
      $ratio = $span / $length;
      if ($result && $result->percent_id > $PS_PERCENT_ID_CUTOFF){
	push @results,$result;
	print  STDERR "\n".$gene->stable_id." "."matches transcript ".$homolog_hash{$homogene}->stable_id." with ".$result->percent_id."% ID,score - ".$result->score." s/l - ".$ratio." length - ".$length." span - ".$span."  exon no - ".$n_exons."\n";

	if ($ratio > $PS_SPAN_RATIO && $n_exons >= $PS_MIN_EXONS){
	  $judgement = 'pseudogene';
	  print STDERR "Calling it a pseudogene\n";
	}
      }
      print STDERR "\n";
    }
  }
 return $judgement;
}

=head2 remove_transcript_from_gene

  Args       : Bio::EnsEMBL::Gene object , Bio::EnsEMBL::Transcript object
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
  return $self->modified_genes;
}



#################################################################################
# Container methods



=head2 modified_genes

Arg [1]    : array ref
  Description: get/set modified gene set to return 
  Returntype : array ref to Bio::EnsEMBL::Gene objects
  Exceptions : throws if not a Bio::EnsEMBL::Gene
  Caller     : general

=cut

sub modified_genes {
  my ($self, $modified_genes) = @_;
  if ($modified_genes) {
    unless  ($modified_genes->isa("Bio::EnsEMBL::Gene")){
      $self->throw("Input isn't a Bio::EnsEMBL::Gene, it is a $modified_genes\n$@");
    }
    push @{$self->{'_modified_genes'}}, $modified_genes;
  }
  return $self->{'_modified_genes'};
}

=head2 discarded transcripts

Arg [1]    : array ref
  Description: get/set modified gene set to throw away
  Returntype : array ref to Bio::EnsEMBL::Gene objects
  Exceptions : throws if not a Bio::EnsEMBL::Gene
  Caller     : general

=cut

sub discarded_transcripts {
  my ($self, $discarded_transcripts) = @_;
  if ( $discarded_transcripts) {
    unless  ($discarded_transcripts->isa("Bio::EnsEMBL::Transcript")){
      $self->throw("Input isn't a Bio::EnsEMBL::Gene, it is a $discarded_transcripts\n$@");
    }
    push @{$self->{'_discarded_transcripts'}}, $discarded_transcripts;
  }
  return $self->{'_discarded_transcripts'};
}

=head2 genes

Arg [1]    : array ref
  Description: get/set gene set to run over
  Returntype : array ref to Bio::EnsEMBL::Gene objects
  Exceptions : throws if not a Bio::EnsEMBL::Gene
  Caller     : general

=cut

sub genes {
  my ($self, $genes) = @_;
  if ($genes) {
    foreach my $gene (@{$genes}) {
      unless  ($gene->isa("Bio::EnsEMBL::Gene")){
	$self->throw("Input isn't a Bio::EnsEMBL::Gene, it is a $gene\n$@");
      }
    }
    $self->{'_genes'} = $genes;
  }
  return $self->{'_genes'};
}

=head2 repeats

Arg [1]    : array ref
  Description: set repeat set to test genes against
  Returntype : none
  Exceptions : throws if not a Bio::EnsEMBL::SeqFeature
  Caller     : general

=cut

sub repeats {
  my ($self, $repeats) = @_;
  foreach my $repeat_array (values %{$repeats}) {
    foreach my $repeat (@{$repeat_array}) {
      unless ($repeat->isa("Bio::EnsEMBL::SeqFeature")){
        $self->throw("Input is not a Bio::EnsEMBL::SeqFeature, it is a $repeat");
      }
    }
  }
  $self->{'_repeats'} = $repeats;
  return  1;
}

=head2 get_repeats

Arg [1]    : array ref
 Description: get repeat array using a gene object as the key
  Returntype : array ref to Bio::EnsEMBL::SeqFeature objects
  Exceptions : warns if no values corresponding to key
  Caller     : general

=cut

sub get_repeats {
  my ($self, $gene) = @_;
  unless ( $self->{'_repeats'}->{$gene}) {
    warn ("repeat array not found for gene object $gene\n$@");
  }
  return  $self->{'_repeats'}->{$gene};
}

=head2 homologs

Arg [1]    : hash ref
  Description: ref to 2d hash containing single exon genes and homologous transcripts
  Returntype : 2d hash ref
  Exceptions : none
  Caller     : general

=cut

sub homologs {
  my ($self, $homologs) = @_;
  if ($homologs) {
    $self->{'_homologs'}= $homologs;
  } 
  return  $self->{'_homologs'};
}

=head2 get_homologs_by_gene

  Arg [1]    : hash ref
  Description: returns ref to hash of homologues given a gene
  Returntype : hash ref
  Exceptions : throws if _homologs not initialised
  Caller     : general

=cut

sub get_homologs_by_gene {
  my ($self, $gene) = @_;
  my %hash = %{$self->{'_homologs'}};
  unless ($hash{$gene}){
    warn ("homolog hash not found for gene $gene");
    }
  return $hash{$gene};
}

=head2 real

Arg [1]    : none
 Description: scalar number of real genes found
  Returntype : scalar integer
  Exceptions : none
  Caller     : general

=cut

sub real {
  my ($self,$num) = @_;
  if ($num){
    $self->{'_real'} += $num;
  }
  return $self->{'_real'};
}


=head2 pseudogenes

Arg [1]    : none
 Description: scalar number of pseudogenes found
  Returntype : scalar integer
  Exceptions : none
  Caller     : general

=cut

sub pseudogenes{
  my ($self,$num) = @_;
  if ($num){
    $self->{'_pseudogenes'}+= $num;
  }
  return $self->{'_pseudogenes'};
}


return 1;
