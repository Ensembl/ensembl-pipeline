#
# Written by Eduardo Eyras
#
# Copyright GRL/EBI 2002
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::PseudoGeneFinder;

=head1 SYNOPSIS

my $exonerate2genes = Bio::EnsEMBL::Pipeline::RunnableDB::PseudoGeneFinder->new(
                              -db         => $refdb,
			      -input_id   => \@sequences,
			      -rna_seqs   => \@sequences,
			      -analysis   => $analysis_obj,
			      -database   => $GENOMIC,
			      -options    => $EXONERATE_OPTIONS,
			     );
    

$exonerate2genes->fetch_input();
$exonerate2genes->run();
$exonerate2genes->output();
$exonerate2genes->write_output(); #writes to DB

=head1 DESCRIPTION

This object contains the logic for (processed) pseudogene finding.
It is based on the alignment of cDNAs on the genome.
For every cDNA sequence, it will try to detect the locations
for processed pseudogenes with similarity to the cDNA sequence.
For the alignments we use Exonerate (G. Slater).
Potential processed pseudogene loci are test further for
lack of homology with other genomes. The tests carried out are the following:

* whether the supporting evidence is found spliced elsewhere in the genome - processed pseudogenes 
  represent an unspliced copy of a functional transcript.
  This is the first test performed on the alignments.
 
* lack of introns, i.e. single exon transcripts
  Using the ensembl representation, this means we are looking for
  transcripts with >=0 introns such that all the introns are frameshifts.

* presence of a poly A tail downstream of the disrupted open reading frame
 
* absence of methionine at the start of the predicted translation
  This is tricky to test with cDNAs - we will find an alternative for this

* whether there is sequence similarity in homologous regions in other species - 
  most of the detectable processed pseudogenes have appeared 
  after speciation, hence they are independently integrated in 
  the genome and therefore unlikely to have any sequence similarity in homologous regions.

  For this we will need a compara database with dna-dna alignments (only)

We found that the strongest signals for processed pseudogenes are for
single-exon predictions with frameshifts, which are based on protein
evidence that is spliced elsewhere in the genome and that have no
sequence similarity in the homologous region in other genomes.

=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::PseudoGeneFinder;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::NewExonerate;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Bio::EnsEMBL::Pipeline::Config::PseudoGenes::PseudoGenes;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;

use vars qw(@ISA);

@ISA = qw (Bio::EnsEMBL::Pipeline::RunnableDB);

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  
  ############################################################
  # SUPER::new(@args) has put the refdb in $self->db()
  #

  my ($database, $rna_seqs, $query_type, $target_type, $exonerate, $options) =  
      $self->_rearrange([qw(
			    DATABASE 
			    RNA_SEQS
			    QUERY_TYPE
			    TARGET_TYPE
			    EXONERATE
			    OPTIONS
			    )], @args);
  
  # must have a query sequence
  unless( @{$rna_seqs} ){
      $self->throw("ExonerateToGenes needs a query: @{$rna_seqs}");
  }
  $self->rna_seqs(@{$rna_seqs});
  
  # you can pass a sequence object for the target or a database (multiple fasta file);
  if( $database ){
      $self->database( $database );
  }
  else{
      $self->throw("ExonerateToGenes needs a target - database: $database");
  }
  
  # Target type: dna  - DNA sequence
  #              protein - protein sequence
  if ($target_type){
      $self->target_type($target_type);
  }
  else{
    print STDERR "Defaulting target type to dna\n";
    $self->target_type('dna');
  }
  
  # Query type: dna  - DNA sequence
  #             protein - protein sequence
  if ($query_type){
    $self->query_type($query_type);
  }
  else{
    print STDERR "Defaulting query type to dna\n";
    $self->query_type('dna');
  }


  # can choose which exonerate to use
  $self->exonerate($EXONERATE);
  
  # can add extra options as a string
  if ($options){
    $self->options($options);
  }
  return $self;
}

############################################################

sub fetch_input {
  my( $self) = @_;
  
  my @sequences = $self->rna_seqs;
  my $target;
  if ($self->database ){
    $target = $self->database;
  }
  else{
    $self->throw("sorry, it can only run with a database");
  }
  
  my @chr_names = $self->get_chr_names;

  foreach my $chr_name ( @chr_names ){
    
    next unless $chr_name eq '6';
    print STDERR "Running only with chr6\n";
    
    my $database = $target."/".$chr_name.".fa";
    
    # check that the file exists:
    if ( -s $database){
      
      #print STDERR "creating runnable for target: $database\n";
      my $runnable = Bio::EnsEMBL::Pipeline::Runnable::NewExonerate
	->new(
	      -database    => $database,
	      -query_seqs  => \@sequences,
	      -exonerate   => $self->exonerate,
	      -options     => $self->options,
	      -target_type => $self->target_type,
	      -query_type  => $self->query_type,
	     );
      $self->runnables($runnable);
    }
    else{
      $self->warn("file $database not found. Skipping");
    }
  }
}

############################################################

sub run{
  my ($self) = @_;
  my @results;
  
  ############################################################
  # align the cDNAs with Exonerate
  $self->throw("Can't run - no funnable objects") unless ($self->runnables);
  foreach my $runnable ($self->runnables){
    
    # run the funnable
    $runnable->run;  
    
    # store the results
    my @results_here = $runnable->output;
    push ( @results, @results_here );
    print STDERR scalar(@results_here)." matches found in ".$runnable->database."\n";
    
  }
  
  ############################################################
  # filter the output to get the candidate pseudogenes
  my @filtered_results = $self->filter_output(@results);
  
  ############################################################
  # make genes/transcripts
  my @genes = $self->make_genes(@filtered_results);
  
  ############################################################
  # print out the results:
  foreach my $gene (@genes){
    foreach my $trans (@{$gene->get_all_Transcripts}){
      Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Evidence($trans);
    }
  }
  
  ############################################################
  # test for homology
  my @processed_pseudogenes = $self->_test_homology( @genes );
 
  ############################################################
  # map candidates into one-block structures
  my @candidate_genes = $self->_convert_to_block( @processed_pseudogenes );

  
  ############################################################
  # check that the 'block' has in-frame stops
  my @candidate_genes2 = $self->_check_inframe_stops( @candidate_genes );
  
  ############################################################
  # need to convert coordinates
  my @mapped_genes = $self->convert_coordinates( @candidate_genes2 );
  
  ############################################################
  # store the results
  if( !@mapped_genes ){
    print STDERR "No genes stored - exiting\n";
  }
  else{
    $self->output(@mapped_genes);
  }
}

############################################################
# this method will collect all the alignments for each
# cDNA sequence and will try to get the potential processed
# pseudogenes out of them. The logic is as follows: 
# If there is an alignment for which all the introns are in fact frameshifts
# (or there is no introns) and there is a equal or better spliced alignment
# elsewhere in the genome, we consider that a potential processed pseudogene.
# We also consider a potential processed pseudogene when a real intron
# contains a repeat over at least 80% of its length.

sub filter_output{
  my ($self,@results) = @_;
  
  # results are Bio::EnsEMBL::Transcripts with exons and supp_features  
  my @potential_processed_pseudogenes;
  
  ############################################################
  # collect the alignments by cDNA identifier
  my %matches;
  foreach my $transcript (@results ){
    my $id = $self->_evidence_id($transcript);
    push ( @{$matches{$id}}, $transcript );
  }
  
  ############################################################
  # sort the alignments by coverage - number of exons - percentage identity
  my %matches_sorted;
  my %selected_matches;
 RNA:
  foreach my $rna_id ( keys( %matches ) ){
    
    my @selected;
  TRAN:
    foreach my $tran ( @{$matches{$rna_id}} ){ 
      my $score    = $self->_coverage($tran);
      my $perc_id  = $self->_percent_id($tran);
      
      ############################################################
      # we allow the rna to align with perc_id below threshold if
      # if the coverage is much better than just above threshold
      next TRAN unless ( ( $score >= $MIN_COVERAGE && $perc_id >= $MIN_PERCENT_ID )
			    ||
			    ( $score >= (1 + 5/100)*$MIN_COVERAGE && $perc_id >= ( 1 - 3/100)*$MIN_PERCENT_ID )
			  );
      push(@selected,$tran);
    }
    
    @{$matches_sorted{$rna_id}} = 
      sort { my $result = ( $self->_radial_score($b) <=> $self->_radial_score($a) );
	     if ( $result){
	       return $result;
	     }
	     else{
	       return ( scalar(@{$b->get_all_Exons}) <=> scalar(@{$a->get_all_Exons}) );
	     }
	   } @selected;
    
    my $count = 0;
    my $is_spliced = 0;
    my $max_score;
    my $perc_id_of_best;
    my $best_has_been_seen = 0;
    
    print STDERR "####################\n";
    print STDERR "Matches for $rna_id:\n";
    
  TRANSCRIPT:
    foreach my $transcript ( @{$matches_sorted{$rna_id}} ){
      $count++;
      unless ( $max_score ){
	$max_score = $self->_radial_score($transcript);
      }
      unless ( $perc_id_of_best ){
	$perc_id_of_best = $self->_percent_id($transcript);
      }
      
      Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($transcript);
      
      my $score   = $self->_coverage($transcript);
      my $perc_id = $self->_percent_id($transcript);
      
      my @exons  = sort { $a->start <=> $b->start } @{$transcript->get_all_Exons};
      my $start  = $exons[0]->start;
      my $end    = $exons[$#exons]->end;
      my $strand = $exons[0]->strand;
      my $seqname= $exons[0]->seqname;
      $seqname   =~ s/\.\d+-\d+$//;
      my $extent = $seqname.".".$start."-".$end;
      
      my $label;
      if ( $count == 1 ){
	$label = 'best_match';
      }
      elsif ( $count > 1 
	      && $is_spliced 
	      && !Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->is_spliced( $transcript )
	    ){
	$label = 'potential_processed_pseudogene';
      }
      else{
	$label = $count;
      }
      
      ############################################################
      # put flag if the first one is spliced
      if ( $count == 1 && Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->is_spliced( $transcript ) ){
	$is_spliced = 1;
      }
      
      my $accept;
      
      ############################################################
      # we keep anything which is 
      # within the 2% of the best score
      # with score >= $MIN_COVERAGE and percent_id >= $MIN_PERCENT_ID
      if ( $score >= (0.98*$max_score) ){
	
	# we want to keep the unspliced cases (only frameshifts) where the
	# best match is spliced 
	# and the spliced cases with repeats in the introns
	if ( $count > 1 
	     && $is_spliced 
	     && !Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->is_spliced( $transcript )
	   ){
	  $accept = 'YES';
	  push( @potential_processed_pseudogenes, $transcript);
	}
	############################################################
	# if it has real introns, but the intron falls in a repeat
	# it is also a candidate
	elsif( $self->_has_repeat_in_intron( $transcript ) ){
	  $accept = 'YES';
	  push( @potential_processed_pseudogenes, $transcript);
	}
	else{
	  $accept = 'NO';
	}
      }
      else{
	$accept = 'NO';
      }
      print STDERR "match:$rna_id coverage:$score perc_id:$perc_id extent:$extent strand:$strand comment:$label accept:$accept\n";
      print STDERR "--------------------\n";
    }
  }
  
  return @potential_processed_pseudogenes;
}

############################################################
# this method check whether the found gene
# has any homology in a region in another species
# homologous to the regions it is contained within

sub _test_homology{
  my ($self,@genes) = @_;
  my @selected;
  
  my $focus_species = $FOCUS_SPECIES;
  
  # compara database
  my $compara_db = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(
								-user      => 'ensro',
								-dbname    => $COMPARA_DBNAME,
								-host      => $COMPARA_DBHOST,
							       );

  # database where the focus species dna is
  my $focus_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
						   '-host'   => $REF_DBHOST,
						   '-user'   => 'ensro',
						   '-dbname' => $REF_DBNAME,
						  );
  
  
  my ($target_db, $target_db2);
  my ($target_species, $target_species2);

  if ( @$COMPARATIVE_DBS ){
    if ( $COMPARATIVE_DBS->[0] ){
      $target_species = $COMPARATIVE_DBS->[0]->{SPECIES};
      
      $target_db = $compara_db->get_db_adaptor($target_species,$COMPARATIVE_DBS->[0]->{PATH});
      
      unless ( $target_db ){
	$target_db =  
	  Bio::EnsEMBL::DBSQL::DBAdaptor
	    ->new(
		  -user      => 'ensro',
		  -dbname    => $COMPARATIVE_DBS->[0]->{DBNAME},
		  -host      => $COMPARATIVE_DBS->[0]->{DBHOST},
		 );
	$target_db->assembly_type( $COMPARATIVE_DBS->[0]->{PATH} );
	$compara_db->add_db_adaptor( $target_db );
      }
    }
    if ( $COMPARATIVE_DBS->[1] ){
      $target_species2 = $COMPARATIVE_DBS->[1]->{SPECIES};
      $target_db2 = $compara_db->get_db_adaptor($target_species2,$COMPARATIVE_DBS->[1]->{PATH});
      
      unless( $target_db2){
	$target_db2 =  
	  Bio::EnsEMBL::DBSQL::DBAdaptor
	    ->new(
		  -user      => 'ensro',
		  -dbname    => $COMPARATIVE_DBS->[1]->{DBNAME},
		  -host      => $COMPARATIVE_DBS->[1]->{DBHOST},
		 );
	$target_db2->assembly_type( $COMPARATIVE_DBS->[1]->{PATH} );
	$compara_db->add_db_adaptor( $target_db2 );
      }
    }
  }
 
  my %gene_ref;
  my $threshold = 40;
  foreach my $gene (@genes){
    
    my ($homology1,$homology2) = (0,0);
    my @transcripts = @{$gene->get_all_Transcripts};
    my $transcript = $transcripts[0];

    $gene_ref{$transcript} = $self->_evidence_id($transcript);
    
    ############################################################
    # has it got homology in the first species?
    if ( $target_species ){
      print STDERR "\n--- testing for homology in $target_species ---\n";
      $homology1 = Bio::EnsEMBL::Pipeline::GeneComparison::ComparativeTools
	->test_for_orthology_with_tblastx($transcript, $focus_db, $focus_db, $focus_species, $compara_db, $target_db, $target_species, $threshold, \%gene_ref );
    }
    
    ############################################################
    # has it got homology in the second species?
    if ( $target_species2 ){
      print STDERR "\n--- testing for homology in $target_species2 ---\n";
      $homology2 = Bio::EnsEMBL::Pipeline::GeneComparison::ComparativeTools
	->test_for_orthology_with_tblastx($transcript, $focus_db, $focus_db, $focus_species, $compara_db, $target_db2, $target_species2, $threshold, \%gene_ref );
    }
    unless ( $homology1 || $homology2 ){
      push (@selected, $gene);
    }
  }
  return @selected;
}

############################################################
# method to check whether real introns in a transcript
# ( if intron <=9 it is considered a frameshift )
# overlap with a repeat at least over 80% of its length.
# That's another tell-tale sign of a processed pseudogene

sub _has_repeat_in_intron{
  my ($self, $tran) = @_;
  
  # instantiate a db with repeats
  my $repeat_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
						     -host   => $REPEAT_DBHOST,
						     -dbname => $REPEAT_DBNAME,
						     -user   => 'ensro',  
						    );
  
  ############################################################
  # pseudogene is a transcript
  my @exons  = sort { $a->start <=> $b->start} @{$tran->get_all_Exons};
  my $start  = $exons[0]->start;
  my $end    = $exons[$#exons]->end;
  my $strand = $exons[0]->strand;
  my $length = $end - $start + 1;

  ############################################################
  # need to get a slice where pseudogene sits:
  my $chr_name  = $exons[0]->contig->chr_name;
  my $chr_start = $exons[0]->contig->chr_start;
  my $chr_end   = $exons[0]->contig->chr_end;
  
  my $slice_start = $chr_start + $start - 1 - 1000;
  my $slice_end   = $chr_start + $end   - 1 + 1000;

  ############################################################
  # get slice from the same db where our pseudogene is
  my $slice = $repeat_db->get_SliceAdaptor->fetch_by_chr_start_end( $chr_name, $slice_start, $slice_end );
  my @features = @{$slice->get_all_RepeatFeatures( 'RepeatMask' )};
  print STDERR "found ".scalar(@features)." repeats\n";
  unless ( @features ){
    return 0;
  }

  my @clusters = @{$self->_cluster_Features(@features)};

  my $introns_with_repeats = 0;
 INTRON:
  for (my $i=0; $i<$#exons; $i++ ){
    my $intron_start  = $exons[$i]->end + 1;
    my $intron_end    = $exons[$i+1]->start - 1;
    my $intron_strand = $exons[$i]->strand;
    my $intron_length = $intron_end - $intron_start + 1;

    #print STDERR "intron: $intron_start-$intron_end $intron_strand\n";
    next INTRON unless ( $intron_length > 9 );

    ############################################################
    # should we skip really long introns?

    my $overlap_length = 0;
  FEAT:
    foreach my $cluster ( @clusters ){
      my $repeat_start  = $cluster->start + $slice_start - 1;
      my $repeat_end    = $cluster->end   + $slice_start - 1;
      my $repeat_strand = $cluster->strand;
      
      #print STDERR "repeat: $repeat_start-$repeat_end $repeat_strand\n";
      next FEAT unless ( $repeat_strand == $intron_strand );
      next FEAT if ( $repeat_start > $intron_end || $repeat_end < $intron_start );
      
      my $overlap_end   = $self->_min( $repeat_end  , $intron_end );
      my $overlap_start = $self->_max( $repeat_start, $intron_start );
      
      $overlap_length += $overlap_end - $overlap_start + 1;
    }
    
    my $overlap_proportion = 100.00*$overlap_length/$intron_length;
    print STDERR "overlap proportion $overlap_proportion\n";

    if ( $overlap_proportion >= 80 ){
      $introns_with_repeats++;
    }
  }
  return $introns_with_repeats;
} 

############################################################

sub _cluster_Features{
  my ($self, @features) = @_;

  # no point if there are no exons!
  return unless ( scalar( @features) > 0 );   

  # main cluster feature - holds all clusters
  my $cluster_list = [];
  
  # sort exons by start coordinate
  @features = sort { $a->start <=> $b->start } @features;

  # Create the first exon_cluster
  my $cluster = new Bio::EnsEMBL::SeqFeature;
  
  # Start off the cluster with the first exon
  $cluster->add_sub_SeqFeature($features[0],'EXPAND');

  $cluster->strand($features[0]->strand);    
  push( @$cluster_list, $cluster );
  
  # Loop over the rest of the exons
  my $count = 0;
  
 EXON:
  foreach my $feature (@features) {
    if ($count > 0) {
      
      # Add to cluster if overlap AND if strand matches
      if ( $cluster->overlaps($feature) && ( $feature->strand == $cluster->strand) ) { 
	$cluster->add_sub_SeqFeature($feature,'EXPAND');
      }  
      else {
	# Start a new cluster
	$cluster = new Bio::EnsEMBL::SeqFeature;
	$cluster->add_sub_SeqFeature($feature,'EXPAND');
	$cluster->strand($feature->strand);
		
	# and add it to the main_cluster feature
	push( @$cluster_list, $cluster );
      }
    }
    $count++;
  }
  return $cluster_list;
}

############################################################

sub _max{
  my ($self,$max,$other) = @_;
  if ($other > $max){
    $max = $other;
  }
  return $max;
}

############################################################

sub _min{
  my ($self,$min,$other) = @_;
  if ($other < $min){
    $min = $other;
  }
  return $min;
}

############################################################
# method to convert a transcript with exons into
# a 1-exon transcript. Processed pseudogenes are described
# as just 1-block

sub _convert_to_block{
  my ( $self, @genes ) = @_;
  
  my @new_genes;
  foreach my $gene ( @genes ){
    my $new_gene = Bio::EnsEMBL::Gene->new();
    $new_gene->analysis($gene->analysis);
    $new_gene->type($gene->analysis->logic_name);
    foreach my $tran ( @{$gene->get_all_Transcripts} ){
      my $new_transcript  = new Bio::EnsEMBL::Transcript;
      
      my $newexon = Bio::EnsEMBL::Exon->new();
      my @exons = sort{ $a->start <=> $b->start } @{$tran->get_all_Exons};
      my $strand = $exons[0]->strand;
      
      ############################################################
      # make all exons into one
      $newexon->start    ($exons[0]->start);
      $newexon->end      ($exons[-1]->end);
      $newexon->phase    (0);
      $newexon->end_phase(0);
      $newexon->strand   ($strand);
      $newexon->contig   ($exons[0]->contig);
      $newexon->seqname  ($exons[0]->seqname);
  
      ############################################################
      # transfer supporting evidence
      my %evidence_hash;
      foreach my $exon ( @exons ){
	foreach my $sf ( @{$exon->get_all_supporting_features} ){
	  if ( $evidence_hash{$sf->hseqname}{$sf->hstart}{$sf->hend}{$sf->start}{$sf->end} ){
	    next;
	  }
	  $evidence_hash{$sf->hseqname}{$sf->hstart}{$sf->hend}{$sf->start}{$sf->end} = 1;
	  $newexon->add_supporting_features( $sf );
	}
      }
      
      $new_transcript->add_Exon($newexon);
      $new_gene->add_Transcript($new_transcript);
    }
    push(@new_genes, $new_gene);
  }
  return @new_genes;
}

############################################################
# method to check the translation of the transcript-block
# for inframe stops which fall midway in the transcript-block

sub _check_inframe_stops{
  my ($self,@genes) = @_;

  my @selected;
  foreach my $gene ( @genes ){
    my @trans      = @{$gene->get_all_Transcripts};
    my $tran       = $trans[0];
    my $mrna       = $tran->seq->seq;
    my $display_id = $self->_evidence_id($tran);
    
    my $peptide    = Bio::Seq->new( -seq      => $mrna,
				    -moltype  => "dna",
				    -alphabet => 'dna',
				    -id       => $display_id );
    
    my $pep0 = $peptide->translate(undef,undef,0);
    my $pep1 = $peptide->translate(undef,undef,1);
    my $pep2 = $peptide->translate(undef,undef,2);
  }
  
  
  return @genes;

}


############################################################

sub _evidence_id{
    my ($self,$tran) = @_;
    my @exons = @{$tran->get_all_Exons};
    my @evi = @{$exons[0]->get_all_supporting_features};
    return $evi[0]->hseqname;
}

############################################################

sub _coverage{
  my ($self,$tran) = @_;
  if ( $self->{_coverage}{$tran} ){
    return  $self->{_coverage}{$tran};
  }
  my @exons = @{$tran->get_all_Exons};
  my @evi = @{$exons[0]->get_all_supporting_features};
  $self->{_coverage}{$tran} = $evi[0]->score;
  return  $self->{_coverage}{$tran};
}

############################################################

sub _percent_id{
  my ($self,$tran) = @_;
  if ( $self->{_percent_id}{$tran} ){
    return  $self->{_percent_id}{$tran};
  }
  my @exons = @{$tran->get_all_Exons};
  my @evi = @{$exons[0]->get_all_supporting_features};
  $self->{_percent_id}{$tran} = $evi[0]->percent_id;
  return  $self->{_percent_id}{$tran};
}

############################################################

sub _radial_score{
  my ($self,$tran) = @_;
  if ( $self->{_radial_score}{$tran} ){
    return  $self->{_radial_score}{$tran};
  }
  my $p  = $self->_percent_id($tran);
  my $c = $self->_coverage($tran);
  $self->{_radial_score}{$tran} = sqrt( $p**2 + $c**2 );         
  return  $self->{_radial_score}{$tran};
}

############################################################

sub write_output{
  my ($self,@output) = @_;

  ############################################################
  # here is the only place where we need to create a db adaptor
  # for the database where we want to write the genes
  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					      -host             => $PSEUDO_DBHOST,
					      -user             => $PSEUDO_DBUSER,
					      -dbname           => $PSEUDO_DBNAME,
					      -pass             => $PSEUDO_DBPASS,
					      );

  # Get our gene adaptor, depending on the type of gene tables
  # that we are working with.
  my $gene_adaptor;


  $gene_adaptor = $db->get_GeneAdaptor;  
  
  unless (@output){
      @output = $self->output;
  }
  foreach my $gene (@output){
    print STDERR "gene is a $gene\n";
    print STDERR "about to store gene ".$gene->type." $gene\n";
    foreach my $tran (@{$gene->get_all_Transcripts}){
      Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($tran);
    }
    
    eval{
      $gene_adaptor->store($gene);
    };
    
    if ($@){
      $self->warn("Unable to store gene!!\n$@");
    }
  }
}

############################################################

sub make_genes{
  my ($self,@transcripts) = @_;
  
  my @genes;
  my $slice_adaptor = $self->db->get_SliceAdaptor;
  my $gene;
  my $checked_transcript;
  foreach my $tran ( @transcripts ){
    $gene = Bio::EnsEMBL::Gene->new();
    $gene->analysis($self->analysis);
    $gene->type($self->analysis->logic_name);
    
    ############################################################
    # put a slice on the transcript
    my $slice_id      = $tran->start_Exon->seqname;
    my $chr_name;
    my $chr_start;
    my $chr_end;
    #print STDERR " slice_id = $slice_id\n";
    if ($slice_id =~/$INPUTID_REGEX/){
      $chr_name  = $1;
      $chr_start = $2;
      $chr_end   = $3;
      #print STDERR "chr: $chr_name start: $chr_start end: $chr_end\n";
    }
    else{
      $self->warn("It cannot read a slice id from exon->seqname. Please check.");
    }
    my $slice = $slice_adaptor->fetch_by_chr_start_end($chr_name,$chr_start,$chr_end);
    foreach my $exon (@{$tran->get_all_Exons}){
      $exon->contig($slice);
      foreach my $evi (@{$exon->get_all_supporting_features}){
	$evi->contig($slice);
	$evi->analysis($self->analysis);
      }
    }
    
    $checked_transcript = $self->check_strand( $tran );
    $gene->add_Transcript($checked_transcript);
    push( @genes, $gene);
  }
  return @genes;
}

############################################################

sub convert_coordinates{
  my ($self,@genes) = @_;
  
  my $rawcontig_adaptor = $self->db->get_RawContigAdaptor;
  my $slice_adaptor     = $self->db->get_SliceAdaptor;
  
  my @transformed_genes;
 GENE:
  foreach my $gene (@genes){
  TRANSCRIPT:
    foreach my $transcript ( @{$gene->get_all_Transcripts} ){
      
      # is it a slice or a rawcontig?
      my $rawcontig = $rawcontig_adaptor->fetch_by_name($transcript->start_Exon->seqname);
      if ( $rawcontig ){
	foreach my $exon (@{$transcript->get_all_Exons}){
	  $exon->contig( $rawcontig);
	}
      }
      
      my $contig = $transcript->start_Exon->contig;
      
      if ( $contig && $contig->isa("Bio::EnsEMBL::RawContig") ){
	print STDERR "transcript already in raw contig, no need to transform:\n";
	Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($transcript);
	next TRANSCRIPT;
      }
      my $slice_id      = $transcript->start_Exon->seqname;
      my $chr_name;
      my $chr_start;
      my $chr_end;
      if ($slice_id =~/$INPUTID_REGEX/){
	$chr_name  = $1;
	$chr_start = $2;
	$chr_end   = $3;
      }
      else{
	$self->warn("It looks that you haven't run on slices. Please check.");
	next TRANSCRIPT;
      }
      
      my $slice = $slice_adaptor->fetch_by_chr_start_end($chr_name,$chr_start,$chr_end);
      foreach my $exon (@{$transcript->get_all_Exons}){
	$exon->contig($slice);
	foreach my $evi (@{$exon->get_all_supporting_features}){
	  $evi->contig($slice);
	}
      }
    }
    
    my $transformed_gene;
    eval{
      $transformed_gene = $gene->transform;
    };
    if ( !$transformed_gene || $@ ){
      my @t        = @{$gene->get_all_Transcripts};
      my $id       = $self->_evidence_id( $t[0] );
      my $coverage = $self->_coverage( $t[0] );
      print STDERR "gene $id with coverage $coverage falls on a gap\n";
    }
    else{
      push( @transformed_genes, $transformed_gene);
    }
  }
  return @transformed_genes;
}

############################################################

sub get_chr_names{
  my ($self) = @_;
  my @chr_names;
  
  print STDERR "fetching chromosomes info\n";
  my $chr_adaptor = $self->db->get_ChromosomeAdaptor;
  my @chromosomes = @{$chr_adaptor->fetch_all};
  
  foreach my $chromosome ( @chromosomes ){
    push( @chr_names, $chromosome->chr_name );
  }
  print STDERR "retrieved ".scalar(@chr_names)." chromosome names\n";
  return @chr_names;
}
############################################################

=head2 check_splice_sites

processed pseudogenes are single-exon objects. Using dna-2-dna alignments
remains a questions of what is the correct strand for the alignment.
We try to estimate that.

=cut

sub check_strand{
  my ($self,$transcript) = @_;

  my $verbose = 0;
    
  #print STDERR "checking splice sites in transcript:\n";
  #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($transcript);
  #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_TranscriptEvidence($transcript);
  
  my $strand = $transcript->start_Exon->strand;
  my @exons;
  if ( $strand == 1 ){
    @exons = sort { $a->start <=> $b->start } @{$transcript->get_all_Exons};
  }
  else{
    @exons = sort { $b->start <=> $a->start } @{$transcript->get_all_Exons};
  }
  my $introns  = scalar(@exons) - 1 ; 
  if ( $introns <= 0 ){
    return $transcript;
  }
  
  my $correct  = 0;
  my $wrong    = 0;
  my $other    = 0;
  
  # all exons in the transcripts are in the same seqname coordinate system:
  my $slice = $transcript->start_Exon->contig;
  
  if ($strand == 1 ){
    
  INTRON:
    for (my $i=0; $i<$#exons; $i++ ){
      my $upstream_exon   = $exons[$i];
      my $downstream_exon = $exons[$i+1];
      
      my $upstream_start = ($upstream_exon->end     + 1);
      my $upstream_end   = ($upstream_exon->end     + 2);      
      my $downstream_start = $downstream_exon->start - 2;
      my $downstream_end   = $downstream_exon->start - 1;

      #eval{
#	$upstream_site = 
#	  $slice->subseq( ($upstream_exon->end     + 1), ($upstream_exon->end     + 2 ) );
#	$downstream_site = 
#	  $slice->subseq( ($downstream_exon->start - 2), ($downstream_exon->start - 1 ) );
#      };
            
      #print STDERR "upstream $upstream_site, downstream: $downstream_site\n";
      ## good pairs of upstream-downstream intron sites:
      ## ..###GT...AG###...   ...###AT...AC###...   ...###GC...AG###.
      
      ## bad  pairs of upstream-downstream intron sites (they imply wrong strand)
      ##...###CT...AC###...   ...###GT...AT###...   ...###CT...GC###...
      
      # print STDERR "From database:\n";
#      print STDERR "strand: + upstream (".
#	($upstream_exon->end + 1)."-".($upstream_exon->end + 2 ).") = $upstream_site, downstream ".
#	  ($downstream_exon->start - 2)."-".($downstream_exon->start - 1).") = $downstream_site\n";
      
      ############################################################
      # use chr_subseq from Tim Cutts
      
      my $upstream_site   = 
	$self->get_chr_subseq($slice->chr_name, $upstream_start, $upstream_end, $strand );
      my $downstream_site = 
	$self->get_chr_subseq($slice->chr_name, $downstream_start, $downstream_end, $strand );
      
      unless ( $upstream_site && $downstream_site ){
	print STDERR "problems retrieving sequence for splice sites\n";
	next INTRON;
      }

      print STDERR "Tim's script:\n" if $verbose;
      print STDERR "strand: + upstream (".
      	($upstream_start)."-".($upstream_end).") = $upstream_site, downstream ".
        ($downstream_start)."-".($downstream_end).") = $downstream_site\n" if $verbose;

      if (  ($upstream_site eq 'GT' && $downstream_site eq 'AG') ||
	    ($upstream_site eq 'AT' && $downstream_site eq 'AC') ||
	    ($upstream_site eq 'GC' && $downstream_site eq 'AG') ){
	$correct++;
      }
      elsif (  ($upstream_site eq 'CT' && $downstream_site eq 'AC') ||
	       ($upstream_site eq 'GT' && $downstream_site eq 'AT') ||
	       ($upstream_site eq 'CT' && $downstream_site eq 'GC') ){
	$wrong++;
      }
      else{
	$other++;
      }
	} # end of INTRON
  }
  elsif ( $strand == -1 ){
    
    #  example:
    #                                  ------CT...AC---... 
    #  transcript in reverse strand -> ######GA...TG###... 
    # we calculate AC in the slice and the revcomp to get GT == good site
    
  INTRON:
    for (my $i=0; $i<$#exons; $i++ ){
      my $upstream_exon   = $exons[$i];
      my $downstream_exon = $exons[$i+1];
      my $up_site;
      my $down_site;
      
      my $up_start   = $upstream_exon->start - 2;
      my $up_end     = $upstream_exon->start - 1;
      my $down_start = $downstream_exon->end + 1;
      my $down_end   = $downstream_exon->end + 2;

      #eval{
#	$up_site = 
#	  $slice->subseq( ($upstream_exon->start - 2), ($upstream_exon->start - 1) );
#	$down_site = 
#	  $slice->subseq( ($downstream_exon->end + 1), ($downstream_exon->end + 2 ) );
#      };
       
#      ( $upstream_site   = reverse(  $up_site  ) ) =~ tr/ACGTacgt/TGCAtgca/;
#      ( $downstream_site = reverse( $down_site ) ) =~ tr/ACGTacgt/TGCAtgca/;
      	  
#      print STDERR "From database:\n";
#      print STDERR "strand: + ".
#	"upstream ($up_start-$up_end) = $upstream_site,".
#	  "downstream ($down_start-$down_end) = $downstream_site\n";
      
      ############################################################
      # use chr_subseq from Tim Cutts
      my $upstream_site   = 
	$self->get_chr_subseq($slice->chr_name, $up_start, $up_end, $strand );
      my $downstream_site = 
	$self->get_chr_subseq($slice->chr_name, $down_start, $down_end, $strand );
 
      unless ( $upstream_site && $downstream_site ){
	print STDERR "problems retrieving sequence for splice sites\n";
	next INTRON;
      }
      
      print STDERR "Tim's script:\n" if $verbose;
      print STDERR "strand: + upstream (".
	($up_start)."-".($up_end).") = $upstream_site, downstream ".
	  ($down_start)."-".($down_end).") = $downstream_site\n" if $verbose;

      
      #print STDERR "strand: - upstream $upstream_site, downstream: $downstream_site\n";
      if (  ($upstream_site eq 'GT' && $downstream_site eq 'AG') ||
	    ($upstream_site eq 'AT' && $downstream_site eq 'AC') ||
	    ($upstream_site eq 'GC' && $downstream_site eq 'AG') ){
	$correct++;
      }
      elsif (  ($upstream_site eq 'CT' && $downstream_site eq 'AC') ||
	       ($upstream_site eq 'GT' && $downstream_site eq 'AT') ||
	       ($upstream_site eq 'CT' && $downstream_site eq 'GC') ){
	$wrong++;
      }
      else{
	$other++;
      }
      
    } # end of INTRON
  }
  unless ( $introns == $other + $correct + $wrong ){
    print STDERR "STRANGE: introns:  $introns, correct: $correct, wrong: $wrong, other: $other\n";
  }
  if ( $wrong > $correct ){
    print STDERR "changing strand\n" if $verbose;
    return  $self->change_strand($transcript);
  }
  else{
    return $transcript;
  }
}

############################################################

sub get_chr_subseq{
  my ( $self, $chr_name, $start, $end, $strand ) = @_;

  my $chr_file = $GENOMIC."/".$chr_name.".fa";
  my $command = "chr_subseq $chr_file $start $end |";
 
  #print STDERR "command: $command\n";
  open( SEQ, $command ) || $self->throw("Error running chr_subseq within ExonerateToGenes");
  my $seq = uc <SEQ>;
  chomp $seq;
  close( SEQ );
  
  if ( length($seq) != 2 ){
    print STDERR "WRONG: asking for chr_subseq $chr_file $start $end and got = $seq\n";
  }
  if ( $strand == 1 ){
    return $seq;
  }
  else{
    ( my $revcomp_seq = reverse( $seq ) ) =~ tr/ACGTacgt/TGCAtgca/;
    return $revcomp_seq;
  }
}


############################################################
=head2 change_strand

    this method changes the strand of the exons

=cut

sub change_strand{
    my ($self,$transcript) = @_;
    my $original_strand = $transcript->start_Exon->strand;
    my $new_strand      = (-1)*$original_strand;
    foreach my $exon (@{$transcript->get_all_Exons}){
      $exon->strand($new_strand);
      foreach my $evi ( @{$exon->get_all_supporting_features} ){
	$evi->strand($new_strand);
	$evi->hstrand( $evi->hstrand*(-1) );
      }
    }
    $transcript->sort;
    return $transcript;
}

############################################################


=head2 _get_SubseqFetcher

Prototype fast SubseqFetcher to get splice sites.

=cut

sub _get_SubseqFetcher {
    
  my ($self, $slice) = @_;
  
  if (defined $self->{'_current_chr_name'}  && 
      $slice->chr_name eq $self->{'_current_chr_name'} ){
    return $self->{'_chr_subseqfetcher'};
  } else {
    
    $self->{'_current_chr_name'} = $slice->chr_name;
    my $chr_filename = $GENOMIC . "/" . $slice->chr_name . "\.fa";

    $self->{'_chr_subseqfetcher'} = Bio::EnsEMBL::Pipeline::SubseqFetcher->new($chr_filename);
  }
}



############################################################
#
# get/set methods
#
############################################################

############################################################

# must override RunnableDB::output() which is an eveil thing reading out from the Runnable objects

sub output {
  my ($self, @output) = @_;
  if (@output){
    push( @{$self->{_output} }, @output);
  }
  return @{$self->{_output}};
}

############################################################

sub runnables {
  my ($self, $runnable) = @_;
  if (defined($runnable) ){
    unless( $runnable->isa("Bio::EnsEMBL::Pipeline::RunnableI") ){
      $self->throw("$runnable is not a Bio::EnsEMBL::Pipeline::RunnableI");
    }
    push( @{$self->{_runnable}}, $runnable);
  }
  return @{$self->{_runnable}};
}

############################################################

sub rna_seqs {
  my ($self, @seqs) = @_;
  if( @seqs ) {
    unless ($seqs[0]->isa("Bio::PrimarySeqI") || $seqs[0]->isa("Bio::SeqI")){
      $self->throw("query seq must be a Bio::SeqI or Bio::PrimarySeqI");
    }
    push( @{$self->{_rna_seqs}}, @seqs);
  }
  return @{$self->{_rna_seqs}};
}

############################################################

sub genomic {
  my ($self, $seq) = @_;
  if ($seq){
    unless ($seq->isa("Bio::PrimarySeqI") || $seq->isa("Bio::SeqI")){
      $self->throw("query seq must be a Bio::SeqI or Bio::PrimarySeqI");
    }
    $self->{_genomic} = $seq ;
  }
  return $self->{_genomic};
}

############################################################

sub exonerate {
  my ($self, $location) = @_;
  if ($location) {
    $self->throw("Exonerate not found at $location: $!\n") unless (-e $location);
    $self->{_exonerate} = $location ;
  }
  return $self->{_exonerate};
}

############################################################

sub options {
  my ($self, $options) = @_;
  if ($options) {
    $self->{_options} = $options ;
  }
  return $self->{_options};
}

############################################################

sub database {
  my ($self, $database) = @_;
  if ($database) {
    $self->{_database} = $database;
  }
  return $self->{_database};
}
############################################################

sub query_type {
  my ($self, $mytype) = @_;
  if (defined($mytype) ){
    my $type = lc($mytype);
    unless( $type eq 'dna' || $type eq 'protein' ){
      $self->throw("not the right query type: $type");
    }
    $self->{_query_type} = $type;
  }
  return $self->{_query_type};
}

############################################################

sub target_type {
  my ($self, $mytype) = @_;
  if (defined($mytype) ){
    my $type = lc($mytype);
    unless( $type eq 'dna' || $type eq 'protein' ){
      $self->throw("not the right target type: $type");
    }
    $self->{_target_type} = $type ;
  }
  return $self->{_target_type};
}

############################################################





1;
