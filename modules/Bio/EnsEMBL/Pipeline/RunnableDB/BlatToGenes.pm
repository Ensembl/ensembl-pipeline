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

my $blat2genes = Bio::EnsEMBL::Pipeline::RunnableDB::BlatSSAHA->new(
    -dbobj      => $db,
    -input_id   => $input_id
    -analysis   => $analysis
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
  
  my ($genomic,$database,$rna_seqs, $query_type, $target_type, $blat, $options) =  $self->_rearrange([qw(GENOMIC
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
  elsif ($genomic){
    $self->genomic($genomic);
  }
  else{
    $self->throw("BlatToGenes needs a target - genomic: $genomic or  database: $database");
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
    $self->($options);
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
  elsif( $self->genomic ){
    $target = $self->genomic;
  }

  my $runnable = Bio::EnsEMBL::Pipeline::Runnable::Blat->new(
							     -database    => $target,
							     -query_seqs  => \@sequences,
							     -query_type  => $self->query_type,
							     -target_type => $self->target_type,
							     -blat        => $self->blat,
							     -options     => $self->options,
							    );
  
  $self->runnable($runnable);
}

############################################################

sub run{
  my ($self) = @_;
  my @results;

  # get the funnable
  $self->throw("Can't run - no runnable objects") unless defined($self->runnable);
  my $runnable = $self->runnable;

  # run the funnable
  $runnable->run;  
  push ( @results, $runnable->output );
  
  # make genes out of the features
  my @genes = $self->make_genes(@results);
  
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
  #my @mapped_genes = $self->convert_to_raw_contig( @genes );
  print STDERR "transforming genes here\n";

  $self->output(@genes);
}

############################################################

sub write_output{
  my ($self,@output) = @_;
  unless (@output){
    @output = $self->output;
  }

}

############################################################

sub make_genes{
  my ($self,@features) = @_;
  
  my @genes;
  foreach my $feature ( @features ){
    my $transcript = Bio::EnsEMBL::Transcript->new();
    my $gene       = Bio::EnsEMBL::Gene->new();
    $gene->add_Transcript($transcript);
    
    foreach my $sub_feature ($feature->sub_SeqFeature){
      # each sub_feature is a feature pair
      
      # make the exon out of the feature1 (the genomic feature)
      my $exon = Bio::EnsEMBL::Exon->new();
      $exon->contig ($sub_feature->feature1->contig);
      $exon->seqname($sub_feature->feature1->seqname);
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
#
# get/set methods
#
############################################################

sub runnable {
  my ($self, $runnable) = @_;
  if (defined($runnable) ){
    unless( $runnable->isa("Bio::EnsEMBL::Pipeline::RunnableI") ){
      $self->throw("$runnable is not a Bio::EnsEMBL::Pipeline::RunnableI");
    }
    $self->{_runnable} = $runnable;
  }
  return $self->{_runnable};
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
