#
# Cared for by Eduardo Eyras  <eae@sanger.ac.uk>
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documenta

=pod 

=head1 NAME

  Bio::EnsEMBL::Pipeline::Runnable::GeneGraphGenerator

=head1 SYNOPSIS

    my $graph = Bio::EnsEMBL::Pipeline::Runnable::GeneGraph->new();
    my @transcripts_in_graph = $graph->_find_connected_graphs(@transcripts);
    
=head1 DESCRIPTION

It creates a graph where the nodes are the transcripts and the vertices is a relation between the
transcripts. In general this relation is 'having one exon in common'. This will most likely be modified 
in the future. For instance, to include the realtion 'having one intron in common'.

The method  _find_connected_graphs will retrieve from the graph the connected components
with more than one element, and returns all the transcripts that are in these connected components.
This will be soon modified to be able to retrieve either the components separately or the total list of
transcripts.

=head1 CONTACT

eae@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...
package Bio::EnsEMBL::Pipeline::Runnable::GeneGraphGenerator;

use diagnostics;
use vars qw(@ISA);
use strict;

use Bio::Range;
use Bio::EnsEMBL::Pipeline::GeneComparison::GeneCluster;
use Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster;
use Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptComparator;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Bio::EnsEMBL::Pipeline::Tools::GeneUtils;
use Bio::EnsEMBL::Pipeline::Tools::ExonUtils;
use Bio::EnsEMBL::Pipeline::GeneCombinerConf;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

############################################################

sub new{
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  
  return $self;
  
}


############################################################
# this function takes transcripts 
# and tries to figure out whether they
# make an acceptable set of isoforms

sub _check_est_Cluster{
  my ($self,@est_transcripts) = @_;
  my %color;

  #if ( scalar(@est_transcripts) == 1 ){
  print STDERR "cluster with ".scalar(@est_transcripts)." transcripts\n";
  #}

  # adjacency lists:
  my %adj;
  my %seen;
  my @linked;

  for(my $i=0;$i<scalar(@est_transcripts);$i++){
    for(my $j=0;$j<scalar(@est_transcripts);$j++){
      
      next if $j==$i;
      print STDERR "Comparing transcripts:\n";
      Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($est_transcripts[$i]);
      Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($est_transcripts[$j]);
      
      # we only check on coincident exon:
      if ( $self->_check_exact_intron_Match( $est_transcripts[$i], $est_transcripts[$j]) 
	   && $self->_check_exact_exon_Match( $est_transcripts[$i], $est_transcripts[$j]) 
	   ){
	  print STDERR "they are linked\n";
	  push ( @{ $adj{$est_transcripts[$i]} } , $est_transcripts[$j] );
	  unless ( defined $seen{ $est_transcripts[$i] } &&  $seen{ $est_transcripts[$i] } ){
	      push ( @linked, $est_transcripts[$i] );
	      $seen{ $est_transcripts[$i] } =1 ;
	  }
	  unless ( defined $seen{ $est_transcripts[$j] } &&  $seen{ $est_transcripts[$j] } ){
	      push ( @linked, $est_transcripts[$j] );
	      $seen{ $est_transcripts[$j] } = 1;
	   }
       }
    }
  }
  
  my @accepted_transcripts;
  ############################################################
  # find the connected components doing a depth-first search
  # only with the transcripts that have been linked
  ############################################################
  if ( @linked ){
      print STDERR scalar(@linked). " linked transcripts\n";
      
      print STDERR "adjacency lists:\n";
      foreach my $tran (@linked){
	  print STDERR $tran->dbID." -> ";
	  foreach my $link ( @{ $adj{ $tran } } ){
	      print STDERR $link->dbID.",";
	  }
	  print STDERR "\n";
      }
      
      foreach my $tran ( @linked ){
	  $color{$tran} = "white";
      }
      
      my @potential_genes;
      
      
      foreach my $tran ( @linked ){
	  if ( $color{$tran} eq 'white' ){
	      my @potential_gene;
	      $self->_visit( $tran, \%color, \%adj, \@potential_gene);
	      push ( @potential_genes, \@potential_gene );
	  }
      }
      print STDERR scalar(@potential_genes)." potential genes created\n";
      my @accepted_transcripts;
      
      # we accept 
      foreach my $gene (@potential_genes){
	  push ( @accepted_transcripts, @$gene );
      }
  }
  elsif ( @est_transcripts ){
      ############################################################
      # if there are no linked transcripts, but we have some transcripts to work with,
      # it is worthwhile to try to keep at least one, every gene counts
      ############################################################
      
      print STDERR "no linked est-transcripts, taking the one with the best score/perc_id\n";
      my @sorted_transcripts = sort { my $result = ( $self->_get_score($b) <=> 
						     $self->_get_score($a) );
				      unless($result){
					  return ( $self->_get_precent_id($b) <=>
						   $self->_get_percent_id($a) );
				      }
				      return $result;
				  } @est_transcripts;
      
      push( @accepted_transcripts, shift @sorted_transcripts );
  }
  
  print STDERR "returning ".scalar( @accepted_transcripts)." transcripts\n";
  return @accepted_transcripts;
      
}

#########################################################################

sub _visit{
  my ($self, $node, $color, $adj, $potential_gene) = @_;
  
  # node is a transcript object;
  $color->{ $node } = 'gray';

  foreach my $trans ( @{ $adj->{$node} } ){
    if ( $color->{ $trans } eq 'white' ){
      $self->_visit( $trans, $color, $adj, $potential_gene );
    }
  }
  unless ( $color->{$node} eq 'black'){
    push( @{ $potential_gene }, $node);
  }
  $color->{ $node } = 'black';    
  return;
}

############################################################
#
# METHODS FOR DEFINING THE VERTEX RELATION
#
############################################################

#########################################################################
# having a set of est_genes only, we have no reference transcript (ensembl one),
# so to determine whether two transcripts are two alternative forms of the same gene
# we check whether they share at least an exon

sub _check_exact_exon_Match{
 my ($self, $tran1, $tran2 ) = @_;
 my @exons1 = @{$tran1->get_all_Exons};
 my @exons2 = @{$tran2->get_all_Exons};
 my $exact_match = 0;
 
 # how many exact matches we need (maybe 1 is enough)
 foreach my $exon1 (@exons1){
   foreach my $exon2 (@exons2){
     return 1 if  ( $exon1->start == $exon2->start && $exon1->end == $exon2->end );
   }
 }
 return 0;
}

############################################################

sub _check_exact_intron_Match{
  my ($self,$tran1,$tran2) = @_;
  my @introns1 = $self->_get_introns_from_Transcript($tran1);
  my @introns2 = $self->_get_introns_from_Transcript($tran2);
  my %range_metric;
  foreach my $intron1 ( @introns1 ){
    $range_metric{ $intron1->start }{ $intron1->end } = 1;
  }
  foreach my $intron2 ( @introns2 ){
    if ( defined  $range_metric{ $intron2->start }{ $intron2->end } 
	 &&  $range_metric{ $intron2->start }{ $intron2->end } == 1 ){
      return 1;
    }
  }
  return 0;
}


############################################################

sub _get_introns_from_Transcript{
  my ($self,$trans) = @_;
  my @exons = sort { $a->start <=> $b->start } @{$trans->get_all_Exons};
  
  my @introns;
  
  for(my $i=0; $i<scalar(@exons);$i++){
    if ( $i > 0 ){
      my $intron_range = Bio::Range->new();
      $intron_range->start($exons[$i-1]->end + 1);
      $intron_range->end(  $exons[$i]->start - 1);
      push( @introns, $intron_range);
    }
  }
  return @introns;
}

#########################################################################
#
# having a set of est_genes only, we have no reference transcript (ensembl one),
# so to determine whether two transcripts are two alternative forms of the same gene
# we check whether they have a similar protein product

sub _check_protein_Match{
 my ($self, $tran1, $tran2 ) = @_;

 my $seq1;
 my $seq2;

 my $compatible_proteins = 0;
 eval{
   $seq1 = $tran1->translate;
   $seq2 = $tran2->translate;
 };
 if ( $seq1 && $seq2 ){
   
   if ( $seq1 =~/\*/ || $seq2 =~/\*/ ){ 
     print STDERR "On of the peptides has a stop codon\n";
     return 0;
   }
   elsif ( $seq1 eq $seq2 ){
     print STDERR "Identical translation\n";
     return 1;
   }
   elsif( $seq1 =~/$seq2/ || $seq2 =~/$seq1/ ){
     return 1;
   }
 }
 else{
   print STDERR "unable to compare translations\n";
   return 0;
 }
}

############################################################

# we average over the number of exons with supporting evidence
# we take the largest value from all the evidences in each exon

sub _get_score{
    my ( $self, $tran ) = @_;
    my $score = 0;
    my $count = 0;
    foreach my $exon ( @{$tran->get_all_Exons} ){
	$count++;
	my @evidence = sort { $b->score <=> $a->score } @{$exon->get_all_supporting_features};
	$score += $evidence[0];
    }
    my $final_score = $score/$count;
    return $final_score;
}

############################################################

# we average over the number of exons with supporting evidence
# we take the largest value from all the evidences in each exon

sub _get_percent_id{
    my ( $self, $tran ) = @_;
    my $percent_id = 0;
    my $count = 0;
    foreach my $exon ( @{$tran->get_all_Exons} ){
	$count++;
	my @evidence = sort { $b->percent_id <=> $a->percent_id } @{$exon->get_all_supporting_features};
	$percent_id += $evidence[0];
    }
    my $final_percent_id = $percent_id/$count;
    return $final_percent_id;
}


1;
