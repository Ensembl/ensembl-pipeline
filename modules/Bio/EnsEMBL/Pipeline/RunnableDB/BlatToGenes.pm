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

Bio::EnsEMBL::Pipeline::RunnableDB::BlatToGenes;


=head1 SYNOPSIS

my $blat2genes = Bio::EnsEMBL::Pipeline::RunnableDB::BlatToGenes->new(
                              -db         => $db,
			      -input_id   => \@sequences,
			      -rna_seqs   => \@sequences,
			      -analysis   => $analysis_obj,
			      -database   => $EST_BLAT_GENOMIC,
			      -query_type => 'rna',
			      -target_type=> 'dna',
			      -options    => $EST_BLAT_OPTIONS,
			     );
    

$blat2genes->fetch_input();
$blat2genes->run();
$blat2genes->output();
$blat2genes->write_output(); #writes to DB

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Pipeline::Runnable::Blat
It is meant to provide the interface for mapping ESTs to the genome
sequence and writing the results as genes. By the way Blat is run
(similar to the way Exonerate is run) we do not cluster transcripts into
genes and only write one transcript per gene.
A Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor is required for database storage.

=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::BlatToGenes;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::Blat;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Bio::EnsEMBL::Pipeline::ESTConf;

use vars qw(@ISA);

@ISA = qw (Bio::EnsEMBL::Pipeline::RunnableDB);

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  
  my ($database,$rna_seqs, $query_type, $target_type, $blat, $options) =  $self->_rearrange([qw(
												DATABASE 
												RNA_SEQS
												QUERY_TYPE
												TARGET_TYPE
												BLAT
												OPTIONS
											       )], @args);
  
  # must have a query sequence
  unless( @{$rna_seqs} ){
    $self->throw("BlatToGenes needs a query: @{$rna_seqs}");
  }
  $self->rna_seqs(@{$rna_seqs});
  
  # you can pass a sequence object for the target or a database (multiple fasta file);
  if( $database ){
    $self->database( $database );
  }
  else{
    $self->throw("BlatToGenes needs a target - database: $database");
  }
  
  # Target type: dna  - DNA sequence
  #              prot - protein sequence
  #              dnax - DNA sequence translated in six frames to protein
  #              The default is dna
  if ($target_type){
    $self->target_type($target_type);
  }
  else{
    print STDERR "Defaulting target type to dna\n";
    $self->target_type('dna');
  }
  
  # Query type: dna  - DNA sequence
  #             rna  - RNA sequence
  #             prot - protein sequence
  #             dnax - DNA sequence translated in six frames to protein
  #             rnax - DNA sequence translated in three frames to protein
  #             The default is rna
  if ($query_type){
    $self->query_type($query_type);
  }
  else{
    print STDERR "Defaulting query type to rna\n";
    $self->query_type('rna');
  }

  # can choose which blat to use
  $self->blat('/usr/local/ensembl/bin/blat');
  #$self->blat($self->find_executable($blat));
  
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
    my $database = $target."/".$chr_name.".fa";
    
    # check that the file exists:
    if ( -s $database){
      
      print STDERR "creating runnable for target: $database\n";
      my $runnable = Bio::EnsEMBL::Pipeline::Runnable::Blat->new(
								 -database    => $database,
								 -query_seqs  => \@sequences,
								 -query_type  => $self->query_type,
								 -target_type => $self->target_type,
								 -blat        => $self->blat,
								 -options     => $self->options,
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

  # get the funnable
  $self->throw("Can't run - no runnable objects") unless ($self->runnables);
  
  foreach my $runnable ($self->runnables){

    # run the funnable
    $runnable->run;  
    
    # store the results
    print STDERR scalar(@results)." matches found\n";
    push ( @results, $runnable->output );
    
  }
  
  #filter the output
  my @filtered_results = $self->filter_output(@results);
  
  # make genes out of the features
  my @genes = $self->make_genes(@filtered_results);
  
  # print out the results:
  foreach my $gene (@genes){
    foreach my $trans (@{$gene->get_all_Transcripts}){
      Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($trans);
      foreach my $exon (@{$trans->get_all_Exons}){
	foreach my $evi (@{$exon->get_all_supporting_features}){
	  print STDERR $evi->hseqname." ".
	    $evi->start." ".$evi->end." ".
	      $evi->hstart." ".$evi->hend." ".
		$evi->primary_tag." ".$evi->source_tag."\n";
	}
      }
    }
  }
  
  # need to convert coordinates?
  my @mapped_genes = $self->convert_coordinates( @genes );
  print STDERR "mapped_gene is a $mapped_genes[0]\n";

  $self->output(@mapped_genes);
}

############################################################

sub filter_output{
  my ($self,@results) = @_;

  # recall that the results are Bio::EnsEMBL::SeqFeatures
  # where each one contains a set of sub_SeqFeatures representing the exons

  my @good_matches;

  my %matches;
  foreach my $result (@results ){
    push ( @{$matches{$result->seqname}}, $result );
  }
  
  my %matches_sorted_by_coverage;
  my %selected_matches;
  foreach my $rna_id ( keys( %matches ) ){
    @{$matches_sorted_by_coverage{$rna_id}} = sort { $b->score <=> $a->score  } @{$matches{$rna_id}};
    
    my $max_score;
    print STDERR "matches for $rna_id:\n";
    foreach my $match ( @{$matches_sorted_by_coverage{$rna_id}} ){
      unless ($max_score){
	$max_score = $match->score;
      }
      foreach my $sub_feat ( $match->sub_SeqFeature ){
	print STDERR $sub_feat->gffstring." ".$sub_feat->percent_id."\n";
      }
      my $score = $match->score;
      
      # we keep anything which is 
      # within the 2% fo the best score
      # with score >= $EST_MIN_COVERAGE and percent_id >= $EST_MIN_PERCENT_ID
      if ( $score >= (0.98*$max_score) && 
	   $score >= $EST_MIN_COVERAGE && 
	   $match->percent_id >= $EST_MIN_PERCENT_ID ){

	print STDERR "Accept!\n";
	push( @good_matches, $match);
      }
      else{
	print STDERR "Reject!\n";
      }
    }
  }
  
  return @good_matches;
}

############################################################

sub write_output{
  my ($self,@output) = @_;
  
  my $gene_adaptor = $self->db->get_GeneAdaptor;

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
  my ($self,@features) = @_;
  
  my @genes;
  foreach my $feature ( @features ){
    my $transcript = Bio::EnsEMBL::Transcript->new();
    my $gene       = Bio::EnsEMBL::Gene->new();
    $gene->analysis($self->analysis);
    $gene->type($self->analysis->gff_feature);
    
    $gene->add_Transcript($transcript);
    
    foreach my $sub_feature ($feature->sub_SeqFeature){
      # each sub_feature is a feature pair
      
      # make the exon out of the feature1 (the genomic feature)
      my $exon = Bio::EnsEMBL::Exon->new();
      $exon->seqname($sub_feature->feature1->seqname);
      $exon->contig($sub_feature->feature1->contig);
      $exon->start  ($sub_feature->feature1->start);
      $exon->end    ($sub_feature->feature1->end);
      $exon->strand ($sub_feature->feature1->strand);
      
      # we haven't set any translations here!!
      $exon->phase    (0);
      $exon->end_phase(0);
      # score is actually the coverage for the entire rna/est transcript
      $exon->score      ($sub_feature->feature1->score);
      $exon->adaptor    ($self->db->get_ExonAdaptor);
      
      # what about the supporting evidence?
      my @supp_features = ($sub_feature);
      my $supp_feature = Bio::EnsEMBL::DnaDnaAlignFeature->new( -features => \@supp_features);
      $supp_feature->contig     ($exon->contig);
      $supp_feature->seqname    ($sub_feature->feature1->seqname);
      $supp_feature->score      ($sub_feature->feature2->score);
      $supp_feature->percent_id ($sub_feature->feature2->percent_id);
      $supp_feature->analysis   ($self->analysis );
      
      $exon->add_supporting_features($supp_feature);
      $transcript->add_Exon($exon);
    }
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
      if ($slice_id =~/$EST_INPUTID_REGEX/){
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
    

    my $transformed_gene = $gene->transform;
    push( @transformed_genes, $transformed_gene);
  }
  return @transformed_genes;
}

############################################################

sub get_chr_names{
  my ($self) = @_;
  my @chr_names;
  
  print STDERR "featching chromosomes info\n";
  my $chr_adaptor = $self->db->get_ChromosomeAdaptor;
  my @chromosomes = @{$chr_adaptor->fetch_all};
  
  foreach my $chromosome ( @chromosomes ){
    push( @chr_names, $chromosome->chr_name );
  }
  print STDERR "retrieved ".scalar(@chr_names)." chromosome names\n";
  return @chr_names;
}

############################################################



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

sub blat {
  my ($self, $location) = @_;
  if ($location) {
    $self->throw("Blat not found at $location: $!\n") unless (-e $location);
    $self->{_blat} = $location ;
  }
  return $self->{_blat};
}

############################################################

sub query_type {
  my ($self, $mytype) = @_;
  if (defined($mytype) ){
    my $type = lc($mytype);
    unless( $type eq 'dna' || $type eq 'rna' || $type eq 'prot' || $type eq 'dnax' || $type eq 'rnax' ){
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
    unless( $type eq 'dna' || $type eq 'prot' || $type eq 'dnax' ){
      $self->throw("not the right target type: $type");
    }
    $self->{_target_type} = $type ;
  }
  return $self->{_target_type};
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





1;
