# Cared for by Ensembl <ensembl-dev@ebi.ac.uk>
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::ESTTranscriptFilter 

=head1 SYNOPSIS

    my $est_filter = Bio::EnsEMBL::Pipeline::Runnable::ESTTranscriptFilter
                       ->new( -coverage => 95,
			      -perc_id  => 99,
			      -depth    => 10,
			     );
    
    my @accepted_ests = $est_filter->filter(@ests);


=head1 DESCRIPTION

    Filters ESTs (as transcripts) according to a coverage, perc_id and depth thresholds.
    Depth 'd' means that only 'd' number of completely containing higher scores will be
    permitted for this EST. 
    
=head1 CONTACT

Ensembl - ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::Runnable::ESTTranscriptFilter;

use Bio::EnsEMBL::SeqFeature;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Root;
@ISA = qw(Bio::EnsEMBL::Root);

sub new {
    my($class,@args) = @_;
    my $self = $class->SUPER::new(@args);  
    my($coverage,$perc_id,$depth) = $self->_rearrange([qw(COVERAGE
							  PERC_ID
							  DEPTH
							 )],
						      @args);
    
    $coverage  = 90  unless $coverage;
    $perc_id   = 97  unless $perc_id;
    $depth     = 100 unless $depth;
    
    $self->coverage_threshold($coverage);
    $self->perc_id_threshold($perc_id);
    $self->depth_threshold($depth);
    
    return $self;
}

############################################################
# coverage means here means how much of the EST is
# represented in the (transcript) alignment

sub _coverage{
  my ($self,$tran) = @_;
  my @exons = @{$tran->get_all_Exons};
  my $score;
  foreach my $exon ( @exons ){
    my @evi = @{$exons[0]->get_all_supporting_features};
    $score = $evi[0]->score;
    last if defined($score);
  }
  return $score;
}

############################################################

sub _exon_score{
  my ($self,$exon) = @_;
  my @evi = @{$exon->get_all_supporting_features};
  my $score;
  foreach my $evi ( @evi){
    $score = $evi[0]->score;
    last if $score;
  }
  return $score;
}

############################################################

sub _exon_perc_id{
  my ($self,$exon) = @_;
  my @evi = @{$exon->get_all_supporting_features};
  my $score;
  foreach my $evi ( @evi){
    $score = $evi[0]->percent_id;
    last if $score;
  }
  return $score;
}

############################################################

sub _perc_id{
  my ($self,$tran) = @_;
  my @exons = @{$tran->get_all_Exons};
  my @evi = @{$exons[0]->get_all_supporting_features};
  return $evi[0]->percent_id;
}

############################################################

sub _id{
  my ($self,$tran) = @_;
  my @exons = @{$tran->get_all_Exons};
  my @evi = @{$exons[0]->get_all_supporting_features};
  return $evi[0]->hseqname;
} 

############################################################

sub _start{
  my ($self,$tran) = @_;
  my @exons = sort {$a->start <=> $b->start} @{$tran->get_all_Exons};
  return $exons[0]->start;
}  

############################################################

sub _end{
  my ($self,$tran) = @_;
  my @exons = sort {$a->start <=> $b->start} @{$tran->get_all_Exons};
  return $exons[-1]->end;
}  

############################################################

sub filter {
  my ($self,$transcripts) = @_;  
  my $min_coverage = $self->coverage_threshold;
  my $min_perc_id  = $self->perc_id_threshold;
  my $depth        = $self->depth_threshold;
  
  my %exon2est;
 
  # valid hits are stored in a hash of arrays
  # we sort by score to know that the first score for a hseqname is its best  
  @$transcripts = sort { $self->_coverage($b) <=> $self->_coverage($a) } @$transcripts;

  #print STDERR "before filtering ".scalar(@$transcripts)."\n";
  #foreach my $est ( @$transcripts ){
  #  print STDERR $self->_id($est)." coverage:".$self->_coverage($est)." perc_id:".$self->_perc_id($est)."\n";
  #  Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript($est);
  #}
  
  my @accepted_exons;
  my @tmp_ests;
  foreach my $t ( @$transcripts ) {
    
    my $c  = $self->_coverage($t);
    my $p  = $self->_perc_id($t);
    
    next unless ( $c > $min_coverage );
    next unless ( $p > $min_perc_id );
    
    if ( $self->depth_threshold ){
      foreach my $exon ( @{$t->get_all_Exons} ){
	
	$exon2est{$exon} = $t;
	push( @accepted_exons, $exon );
	
      }
    }
    else{
      push ( @tmp_ests, $t );
    }
  }
  
  # empty input array - saves on memory!
  $transcripts = [];
  
  if ( $self->depth_threshold ){
    my @accepted_ests = $self->depth_filter( \@accepted_exons, \%exon2est );
    return @accepted_ests;
  }
  else{
    return @tmp_ests;
  }
}

############################################################

sub depth_filter {
  my ($self,$exons,$exon2est) = @_;

  print STDERR "filtering ".scalar(@$exons)." exons\n";
  # no point if there are no exons!
  return unless ( scalar( @$exons) > 0 );   

  # keep track about in which cluster is each exon
  my %exon2cluster;
  my %exon2est = %$exon2est;
  
  # main cluster feature - holds all clusters
  my $cluster_list = [];
  
  # sort exons by start coordinate
  my @exons = sort { $a->start <=> $b->start } @$exons;

  # Create the first exon_cluster
  my $exon_cluster = new Bio::EnsEMBL::SeqFeature;
  
  # Start off the cluster with the first exon
  $exon_cluster->add_sub_SeqFeature($exons[0],'EXPAND');

  $exon_cluster->strand($exons[0]->strand);    
  push( @$cluster_list, $exon_cluster);
  
  # Loop over the rest of the exons
  my $count = 0;
  
 EXON:
  foreach my $exon (@exons) {
    if ($count > 0) {
      
      # Add to cluster if overlap AND if strand matches
      if ( $exon_cluster->overlaps($exon) 
	   && 
	   ( $exon->strand == $exon_cluster->strand) 
	 ) { 
	$exon_cluster->add_sub_SeqFeature($exon,'EXPAND');
      }  
      else {
	# Start a new cluster
	$exon_cluster = new Bio::EnsEMBL::SeqFeature;
	$exon_cluster->add_sub_SeqFeature($exon,'EXPAND');
	$exon_cluster->strand($exon->strand);
	
	# and add it to the main_cluster feature
	push( @$cluster_list,$exon_cluster);
      }
    }
    $count++;
  }
  
  # sort the clusters by the number of exons and other things...
  my @sorted_clusters = 
    sort { $self->_order($b) <=> $self->_order($a) 
	     or
	       $self->_average_score($b) <=> $self->_average_score($a)
		 or
		   $self->_average_perc_id($b) <=> $self->_average_perc_id($a)
		 } 
  @$cluster_list;
  
  my $depth = $self->depth_threshold;
  
  my %banned;
  my %taken;
  my @accepted_ests;
  
  print STDERR "looking at depth of ".scalar(@sorted_clusters)." clusters\n";
  foreach my $cluster ( @sorted_clusters ){
      my @exons = sort { $self->_exon_score($b) <=> $self->_exon_score($a)
			   or
			     $self->_exon_perc_id($b) <=> $self->_exon_perc_id($a) 
			   } 
      $cluster->sub_SeqFeature;
      
      my $count = 0;
      while ( @exons ){
	my $exon = shift @exons;
	my $est  = $exon2est{$exon};
	if ( $count >= $depth ){
	  $banned{$est} = 1 ;
	}
	
	unless ( $banned{$est} ){
	  unless ( $taken{$est} ){
	    push( @accepted_ests, $est );
	    $taken{$est} = 1;
	  }
	  $count++;
	}
      }
  }

  print STDERR "returning ".scalar(@accepted_ests)." ests after scores and depth filtering\n";
  foreach my $est ( @accepted_ests ){
    print STDERR $self->_id($est)." coverage:".$self->_coverage($est)." perc_id:".$self->_perc_id($est)."\n";
    Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript($est);
  }
  return @accepted_ests;
}

############################################################

sub _order{
    my ($self, $cluster) = @_;
    return scalar( $cluster->sub_SeqFeature );
}

############################################################

sub _average_score{
  my ($self, $cluster) = @_;
  my $coverage = 0;
  my $count = 0;
  foreach my $exon ( $cluster->sub_SeqFeature ){
    my @evi = @{$exon->get_all_supporting_features};
    $coverage += $evi[0]->score;
    $count++;
  }
  my $av = sprintf "%.2f", ( $coverage/$count );
  return $av;
}

############################################################

sub _average_perc_id{
  my ($self, $cluster) = @_;
  my $perc_id = 0;
  my $count   = 0;
  foreach my $exon ( $cluster->sub_SeqFeature ){
    my @evi = @{$exon->get_all_supporting_features};
    $perc_id += $evi[0]->percent_id;
    $count++;
  }
  my $av = sprintf "%.2f", ( $perc_id/$count );
  return $av;
}

############################################################

sub coverage_threshold{
  my $obj = shift;
  if( @_ ) {
    my $value = shift;
    $obj->{'_min_coverage'} = $value;
  }
  return $obj->{'_min_coverage'};    
}

############################################################

sub perc_id_threshold{
    my $obj = shift;
    if( @_ ) {
	my $value = shift;
	$obj->{'_perc_id'} = $value;
    }
    return $obj->{'_perc_id'};
}

############################################################

sub depth_threshold{
   my $obj = shift;
   if( @_ ) {
       my $value = shift;
       $obj->{'_depth'} = $value;
   }
   return $obj->{'_depth'};
}

############################################################

1;
