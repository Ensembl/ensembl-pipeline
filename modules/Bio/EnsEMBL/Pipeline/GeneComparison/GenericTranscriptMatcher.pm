#
# Written by Eduardo Eyras
#
# Cared for by EnsEMBL  <ensembl-dev@ebi.ac.uk>
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::GenericTranscriptMatcher

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::Runnable::GenericTranscriptMatcher->new(
									      -reference_set => \@transcript_set1,
									      -match_set => \@transcript_set2,
									      );

    $obj->run;
   
    # the output is a Bio::EnsEMBL::Pipeline::GeneComparison::ObjectMap object
    my $matching_map = $obj->output;

    # see the documentation in Bio::EnsEMBL::Pipeline::GeneComparison::ObjectMap for details

 
=head1 DESCRIPTION

Class to map one set of transcripts (match_set) to a second set (reference_set) 
according to certain coordinate-based comparison rules

=head1 CONTACT

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::Runnable::GenericTranscriptMatcher;

use diagnostics;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::Pipeline::GeneComparison::ObjectMap;
use Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptComparator;


@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

######################################################################

sub new{
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  my( $reference_set, $match_set ) = 
      $self->_rearrange([qw(REFERENCE_SET MATCH_SET)], @args);
  
  unless( $reference_set && $match_set ){
      $self->warn("lacking of the inputs. Nothing to match");
      exit(0);
  }
  unless ( $reference_set->[0]->isa('Bio::EnsEMBL::Transcript') &&
	   $match_set->[0]->isa('Bio::EnsEMBL::Transcript') ){
      $self->throw("wrong input type, needs two sets of transcripts");
  }
  
  my $matching_map = Bio::EnsEMBL::Pipeline::GeneComparison::ObjectMap->new();
  
  $self->matching_map( $matching_map );	
  $self->_reference_set($reference_set);
  $self->_match_set($match_set);
  return $self;
  
}

#########################################################################
#
# GET/SET METHODS 
#
#########################################################################

sub reference_set{
  my ( $self, $listref ) = @_;
  if ( $listref ){
      $self->{_reference_set} = $listref;
  }
  return $self->{_reference_set};
}

############################################################

sub match_set{
  my ( $self, $listref ) = @_;
  if ( $listref ){
      $self->{_matach_set} = $listref;
  }
  return $self->{_match_set};
}

############################################################

sub gene_Clusters {
  my ($self, @clusters) = @_;
  if (@clusters){
    push ( @{$self->{'_gene_clusters'} }, @clusters);
  }
  return @{ $self->{'_gene_clusters'} };
}

############################################################

sub matching_map{
    my ($self,$object_map) = @_;
    if ( $object_map ){
	$self->{_matching_map} = $object_map;
    }
    return 	$self->{_matching_map};
} 

#########################################################################

sub output{
  my ($self)= @_;
  
  return $self->matching_Map
}

#########################################################################

############################################################
#
# RUN METHOD
#
############################################################

sub run{
  my ($self,@args) = @_;

  my @reference_trans = @{$self->reference_set};
  my @match_trans     = @{$self->match_set};

  print STDERR "Trying to match ".scalar( @match_trans )." transcripts to ".
      scalar(@reference_trans)." reference transcripts\n";
  
  foreach my $transcript ( @reference_trans ){
      $self->_map_Transcripts( $transcript, \@match_trans );
  }
  
  # before returning, check that we have written anything
  unless( $self->matching_Map ){
      print STDERR "not matches found\n";
      exit(0);
  }
  return;
}

############################################################
#
# METHODS CALLED FROM RUN METHOD... DOING ALL THE MAGIC
#
############################################################
  
############################################################
#
#  Method for mapping ESTs or any transcripts to genes:
#
#  For a given transcript,
#
# link to it all the ests that have exons which
# are consecutively included in the transcript
# (this will avoid linking to alternative variants which may in fact be expressed differently)
# allow the 5' and 3' exons of the ests to extend beyond the transcript
# (they may contain UTRs that we failed to annotate)

# ests           (1)###------###
#                          (2)#####--------#######
#           (3)########------###
#                (4)#####----#####
#             (5)####--------#####

#transcript      ######------######--------######

# (1),(2) and (3) would be linked. 
# cases like (4) and (5) would not be linked, unless one uses a different
# method in the transcript comparator 

# case (2) could be a hint for alternative polyA site if we have
# already annotated an UTR for that transcript. We could check for this case
# and only add (2) if there is no 3'UTR, to be sure.

# Case (3) could be also related to an alternative start of transcription,
# we could add it only for cases that a 5'UTR is not annotated.

# Part of the alternative polyA sites and start of transcription seems to be 
# correlated with alternative splicing so maybe this 'ambiguity' cases will
# not cause too many problems.

# there are some methods to include a check of UTRs for (2) and (3). ESTs that are not linked
# would be rejected in principle, I cannot predict yet how much data will
# remain unused by doing this.

sub _map_Transcripts{
  my ($self,$transcript,$ests) = @_;

  # a map from transcript to est/transcript
  my %expression_map;
  
  # a comparison tool
  my $transcript_comparator = Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptComparator->new(
												-comparison_level         => 3,
												-exon_match               => 0,
												-splice_mismatch          => 1,
												-intron_mismatch          => 0,
											       );
  my $sublist;new();
  
 EST:
  foreach my $est (@$ests){
      
      # compare this est
      my ($merge, $overlaps) = $transcript_comparator->compare($transcript,$est);
      
      # (this method checks exact exon boundary matches but
      # allows mismatches at outer end of the 5' and 3' exons)
      # check 5' and 3' ends in case ESTs give an alternative transcription 
      
      # if match, put est in $expression_map{ $transcript }
      if ($merge){
	  $self->matching_map->match($transcript,$est);
      }
  }
}



############################################################
#
# METHODS TO CHECK THE UTRs (Not in use yet)
#
############################################################

sub _check_5prime{
  my ($self,$transcript,$est) = @_;
  my $alt_start = 0;
  
  # first find out whether the transcript has 5' UTR
  my $utr5;
  eval{
    $utr5 = $transcript->five_prime_utr;
  };
  unless( $utr5 ){
    return 0;
  }
  
  $transcript->sort;
  $est->sort;
  
  my $start_exon = $transcript->start_exon;
  #my $start_exon = $transcript->translation->start_exon;
  my $strand     = $start_exon->strand;
  foreach my $exon ( @{$transcript->get_all_Exons} ){
    my $est_exon_count = 0;
    
    foreach my $est_exon ( @{$est->get_all_Exons} ){
      $est_exon_count++;
      if ( $exon == $start_exon ){
	if ( $exon->overlaps( $est_exon ) ){
	  if ($strand == 1){
	    if ( $est_exon->start < $exon->start ){
	      print STDERR "potential alternative transcription start in forward strand\n";
	      $alt_start = 1;
	    }
	  }
	  if ($strand == -1){
	    if ( $est_exon->end > $exon->end ){
	      print STDERR "potential alternative transcription start in reverse strand\n";
	      $alt_start = 1;
	    }
	  }
	  if ($est_exon_count > 1){
	    print STDERR "There are more est exons upstream\n";
	    if ( $alt_start == 1){
	      return 1;
	    }
	  }
	}
      }
    }
  }
  return 0;
}


#########################################################################

sub _check_3prime{
  my ($self,$transcript,$est) = @_;
  my $alt_polyA = 0;
  
  # first find out whether the transcript has 5' UTR
  my $utr3;
  eval{
    $utr3 = $transcript->three_prime_utr;
  };
  unless( $utr3 ){
    return 0;
  }
  
  $transcript->sort;
  $est->sort;

  my $end_exon = $transcript->end_exon;
  #my $end_exon = $transcript->translation->end_exon;
  my $strand   = $end_exon->strand;
  
  foreach my $exon ( @{$transcript->get_all_Exons} ){
    my $est_exon_count = 0;
    my @est_exons = @{$est->get_all_Exons};
    
    foreach my $est_exon ( @est_exons ){
      $est_exon_count++;
      if ( $exon == $end_exon ){
	if ( $exon->overlaps( $est_exon ) ){
	  if ($strand == 1){
	    if ( $est_exon->end > $exon->end ){
	      print STDERR "potential alternative polyA site in forward strand\n";
	      $alt_polyA = 1;
	    }
	  }
	  if ($strand == -1){
	    if ( $est_exon->start < $exon->start ){
	      print STDERR "potential alternative polyA site in reverse strand\n";
	      print STDERR "looking at : exon:".$exon->start."-".$exon->end." and est_exon:".$est_exon->start."-".$est_exon->end."\n";
	      $alt_polyA = 1;
	    }
	  }
	  if ($est_exon_count != scalar(@est_exons) ){
	    print STDERR "There are more est exons downstream\n";
	    print STDERR "est exon count = $est_exon_count, exons = ".scalar(@est_exons)."\n";
	    
	    if ( $alt_polyA ==1 ){
	      return 1;
	    }
	  }
	}
      }
    }
  }
  return 0;
}

#########################################################################
 
1;
