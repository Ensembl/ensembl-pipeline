#
# You may distribute this module under the same terms as perl itself
#

=pod

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::Pseudogene_2

=head1 SYNOPSIS

 my $runnable = Bio::EnsEMBL::Pipeline::Runnable::Pseudogene->new
      ( 
       '-genes' => \@_genes,
       '-repeat_features' => \%repeat_blocks, 
        );
    $runnable->run;
    $runnable->output;

Where output returns an array of modified genes.
Repeat blocks is a hash of repeats covering each gene merged into blocks.


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
  $self->{'_repeats'} = {};	# hash of repeat blocks corresponding to each gene;
  $self->{'_real'} = 0;		# scalar number of real genes identified
  $self->{'_pseudogenes'} = 0;	#scalar number of pseudogenes identified
  
  my( $genes,$repeats) = $self->_rearrange([qw(
							 GENES
							 REPEAT_FEATURES
							)], @args);

  if ($genes) {
    $self->genes($genes);
  }
  if ($repeats) {
    $self->repeats($repeats);
  }

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
    print STDERR   $transcript->dbID."\n";
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
  my $pseudo= 0;
  my $possible = 0;
  my @genes = @{$self->genes};
 
 GENE:  foreach my $gene (@genes) {
    my @pseudo_trans ;
    my @real_trans ;

    ###################################################
    # if gene is single exon assume its real for now...

    if (scalar(@{$gene->get_all_Exons()})==1) {
      $self->modified_genes($gene); 
      $self->real(1);
      next GENE;
    }


    ####################################
    # gene is multiexon, run other tests

  TRANS: foreach my $transcript (@{$gene->get_all_Transcripts}) {
      $num++;
      my $evidence = $self->transcript_evidence($transcript,$gene);

      #transcript tests

      #CALL PSEUDOGENE IF AT LEAST 80% COVERAGE OF INTRONS BY REPEATS
      #AT LEAST 1 F/S EXON AND 1 REAL EXON (?)
      #TOTAL INTRON LENGTH < 5K
	
      if ($evidence->{'total_intron_len'} < $PS_MAX_INTRON_LENGTH &&
	  $evidence->{'frameshift_introns'} >= $PS_NUM_FRAMESHIFT_INTRONS &&
	  $evidence->{'real_introns'} >= $PS_NUM_REAL_INTRONS &&
	  $evidence->{'covered_introns'} >= $PS_MAX_INTRON_COVERAGE  ) {
	push @pseudo_trans, $transcript;
	print STDERR $gene->dbID." - repeats in introns in transcript ".$transcript->dbID."\n";
	print STDERR join (', ',%{$evidence}),"\n";
	next TRANS;
      }


      #ALL FRAMESHIFTED - it is a pseudogene
	
      if ($evidence->{'num_introns'} && 
	  $evidence->{'frameshift_introns'} == $evidence->{'num_introns'}) {
	push @pseudo_trans, $transcript;
	next TRANS;
      }

      # Tests for situation where 2 exon gene has a protein feasture
      # covering intron, ie it has spliced around somthing bad...

      if ($evidence->{'real_introns'} == 1) {
	my $judgement = $self->protein_covered_intron($transcript,$gene);
	if ($judgement eq "dodgy"){
	  push @pseudo_trans, $transcript;
	  next TRANS;
	}	
      }

      # transcript passes all tests, it is real

      push @real_trans, $transcript;
    }

    #########################################
    # gene tests

    #############################################
    # gene is pseudogene, set type to pseudogene
    # chuck away all but the longest transcript

    if (scalar(@pseudo_trans) > 0 && 
	scalar(@real_trans) == 0) {
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
      next GENE;
    }

    ###############################################
    # gene is real but has some dodgy transcripts
    # delete the dodgy transcripts from the gene

    if (scalar(@pseudo_trans) > 0 && 
	scalar(@real_trans) > 0) {
      foreach my $trans (@pseudo_trans) {
	$trans->translation(undef);
	$self->_remove_transcript_from_gene($gene,$trans); 
      }
      $self->modified_genes($gene);
      $self->discarded_transcripts(@pseudo_trans);
      $self->real(1);
      next GENE;
    }

    ####################################
    # gene and transcripts are real


    if (scalar(@pseudo_trans) == 0 && 
	scalar(@real_trans) > 0) {
      $self->modified_genes($gene); 
      $self->real(1);
      next GENE;
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
    ##############################################################
    # Need to make intron and exon Features to compare with repeats
    # features are relative to gene as are the repeats

    my $seq_feature_exon = Bio::EnsEMBL::Feature->new(
							 -START => $exon->start-$gene->start,
							 -END => $exon->end-$gene->start,
							 -STRAND => $exon->strand
							);
    # Do intron
    if (defined($prev_exon)) {
      my $intron = Bio::EnsEMBL::Feature->new(
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
#   print STDERR  "FT " . $feat->start . " " . $feat->end . "\n";
  my $covered_len = 0;
 RBLOOP: foreach my $repeat_block (@$repeat_blocks_ref) {
#    print STDERR  "RB " . $repeat_block->start . " " . $repeat_block->end . "\n";
    if ($repeat_block->overlaps($feat, 'ignore')) {
      my $inter = $self->intersection($feat,$repeat_block);
      $covered_len += $inter->length;
    } elsif ($repeat_block->start > $feat->end) {
      last RBLOOP;
    }
  }
  return  $covered_len;
}

=head2 protein_covered_intron

  Args       : Bio::EnsEMBL::Transcript object, Bio::EnsEMBL::Gene object 
  Description: decides if 'real' intron in transcript is covered with a protein feature
  Returntype : scalar

=cut 


sub protein_covered_intron{
  my ($self,$transcript,$gene) =@_;
  my %seq_features;
  my $identified;
  my @all_exons  = @{$transcript->get_all_Exons};
  @all_exons = sort {$a->start <=> $b->start} @all_exons;  

  my @exons;
  # find real intron
 EXON: for (my $i = 1 ; $i <= $#all_exons ; $i++){
    my $intron_length  = $all_exons[$i]->start - $all_exons[$i-1]->end;
    if ($intron_length > 9) {
      # real intron
      push @exons , $all_exons[$i-1];
      push @exons , $all_exons[$i];
      last EXON;
    }
  }
  $self->throw("real intron not found for gene " . $gene->dbID . " exons : @all_exons\n")  unless (scalar(@exons) == 2); 
  my @exon_features = @{$exons[0]->get_all_supporting_features};
  push @exon_features,@{$exons[1]->get_all_supporting_features};

  
  ########################################
  # make a seq feature represening the intron
  
  my $intron = Bio::EnsEMBL::Feature->new(
					  -START => $exons[0]->end+1,
					  -END => $exons[1]->start-1,
					  -STRAND => $exons[0]->strand
					 );
  if (@exon_features) {
    ###########################################################
    # get featues off both exons, split them into ungapped
    # sections and test them against the intron
    # feature see if there is an overlap (80% by default)
    # need to group features by sequence name

  FEATURES:   foreach my $feat (@exon_features) {
      my @sub_features = $feat->ungapped_features;
      foreach my $subfeature (@sub_features) {
	my $seq_feature = Bio::EnsEMBL::Feature->new(
						     -START => $subfeature->start,
						     -END => $subfeature->end,
						     -STRAND => $subfeature->strand
						       );
	push @{$seq_features{$feat->hseqname}}, $seq_feature;
      }
    }
  }
  foreach my $key (keys %seq_features){
    if ($intron && scalar(@{$seq_features{$key}}) > 0) {
      my @features = @{$seq_features{$key}};
      my @feature_blocks;

      #########################################################
      # merge overlapping features together for protein features
      # grouped by name

      my $curblock = undef;
      @features = sort {$a->start <=> $b->start} @features;
      foreach my $feature (@features){
	if (defined($curblock) && $curblock->end >= $feature->start) {
	  if ($feature->end > $curblock->end) { 
	    $curblock->end($feature->end); 
	  }
	} else {
	  $curblock = Bio::EnsEMBL::Feature->new(
						 -START => $feature->start,
						 -END => $feature->end, 
						 -STRAND => $feature->strand
						);
	  push (@feature_blocks,$curblock);
	}
      }
      my $coverage =  $self->_len_covered($intron,\@feature_blocks)."\n";
      if ($coverage/$intron->length*100 > $PS_MAX_INTRON_COVERAGE) {
	$identified++;
	print STDERR $transcript->dbID." two exon with $key covering intron ".$coverage/$intron->length*100 . "%.\t";

	# need more than one peice of protein evidence to make the call

	if ($identified >1){
	  print STDERR "\ncalling ".$transcript->dbID." as having a protein covered intron\n";
	  return "dodgy";
	}
	else{
	  print "need another piece of evidence though....\n";
	  }
      }
    }
  }
  return 1;
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

=head2 intersection

Arg [none] :
  Description: returns length of intersecion of two features
  Returntype : scalar
  Exceptions : throws if not given two Bio::EsEMBL::Feature objects
  Caller     : len_covered

=cut

sub intersection {
  my ($self,$feat1,$feat2) = @_;
    unless ($feat1->isa("Bio::EnsEMBL::Feature") && $feat2->isa("Bio::EnsEMBL::Feature")){
    $self->throw("object is not a Bio::EnsEMBL::Feature they are feat1 $feat1, feat2 $feat2\n$@\n");
  }
  my @start = sort {$a<=>$b}
    ($feat1->start(), $feat2->start());
    my @end   = sort {$a<=>$b}
    ($feat1->end(),   $feat2->end());

    my $start = pop @start;
    my $end = shift @end;

    if($start > $end) {
	return undef;
    } else {
	return Bio::EnsEMBL::Feature->new('-start'  => $start,
			                  '-end'    => $end,
			                  '-strand' => $feat1->strand
                       			  );
    }
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
      unless ($repeat->isa("Bio::EnsEMBL::Feature")){
        $self->throw("Input is not a Bio::EnsEMBL::Feature, it is a $repeat");
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
