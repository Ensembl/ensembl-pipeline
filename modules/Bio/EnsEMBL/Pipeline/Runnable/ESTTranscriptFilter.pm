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

use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
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
    
    $self->verbose(0);

    return $self;
}

############################################################

sub verbose{
  my ($self,$boolean) = @_;
  if (defined $boolean ){
    $self->{_verbose} = $boolean;
  }
  return $self->{_verbose};
}

############################################################
# coverage means here means how much of the EST is
# represented in the (transcript) alignment

sub _coverage{
    my ($self,$tran) = @_;
    my @exons = @{$tran->get_all_Exons};
    my $score;
    my $hseqname;
  EXON:
    foreach my $exon ( @exons ){
	my @evi = @{$exons[0]->get_all_supporting_features};
	foreach my $evi ( @evi ){
	    $score = $evi->score;
	    $hseqname = $evi->hseqname;
	    last EXON if defined($score);
	}
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
  my $perc_id;
  my @exons = @{$tran->get_all_Exons};
  EXON:
  foreach my $exon ( @exons ){
      my @evi = @{$exon->get_all_supporting_features};
      foreach my $evi ( @evi ){
	  $perc_id = $evi->percent_id;
	  last EXON if defined $perc_id;
      }
  }
  return $perc_id;
}

############################################################

sub _id{
  my ($self,$tran) = @_;
  my @exons = @{$tran->get_all_Exons};
  my $id;
  EXON:
  foreach my $exon ( @exons ){
      my @evi = @{$exon->get_all_supporting_features};
      foreach my $evi ( @evi ){
	  $id = $evi->hseqname;
	  last EXON if defined $id;
      }
  }
  return $id;
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
  
  my $verbose = $self->verbose();

  print STDERR "filter input: ".scalar( @$transcripts )."\n" if $verbose;
  my @filtered_by_length = $self->filter_by_length( $transcripts );
  print STDERR "after filter by length: ".scalar( @filtered_by_length )."\n" if $verbose;


  my $min_coverage = $self->coverage_threshold;
  my $min_perc_id  = $self->perc_id_threshold;
  my $depth        = $self->depth_threshold;
  
  my %exon2est;
 
  if ( $verbose ){
    foreach my $est ( @filtered_by_length ){
      #print STDERR $self->_id($est)." coverage:".$self->_coverage($est)." perc_id:".$self->_perc_id($est)."\n" ;
      #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript($est);
    }
  }
  
  my @accepted_exons;
  my @tmp_ests;
  #print STDERR "Depth threshold = ".$self->depth_threshold."\n";
  foreach my $t ( @filtered_by_length ) {
    
    my $c  = $self->_coverage($t);
    my $p  = $self->_perc_id($t);
        
    next unless ( $c > $min_coverage );
    next unless ( $p > $min_perc_id );
    
    if ( $self->depth_threshold ){
      
      foreach my $exon ( @{$t->get_all_Exons} ){
        if ($verbose){
          #print STDERR "taking: ";
          #$self->_print_EST($t);
        }
        $exon2est{$exon} = $t;
        push( @accepted_exons, $exon );
      }
    }else{
      push ( @tmp_ests, $t );
    }
  }
  
  # empty input array - saves on memory!
  $transcripts = [];
  
  if ( $self->depth_threshold ){
    my @accepted_ests = $self->depth_filter( \@accepted_exons, \%exon2est );
    
    if ($verbose){
      #print STDERR "final list:\n";
      foreach my $est ( @accepted_ests ){ 
	#print STDERR $self->_id($est)." coverage:".$self->_coverage($est)." perc_id:".$self->_perc_id($est)."\n"; 
	#Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript($est); 
      } 
    }
    
    return @accepted_ests;
  }
  else{
    return @tmp_ests;
  }
}

############################################################
# this method consider the distribution of lengths of ests P(L)
# which is close to normal, and selects those
# with lengths L>=Lm where Lm is the maximum for P(L)

sub filter_by_length{
    my ( $self, $ests ) = @_;
    
    my $verbose = $self->verbose();
    my %lengths;
    my @length_list;

    # the genomic distances of ESTs in a small genomic region ( ~1Mb )
    # do not seem to be normally distributed ( unlike the estlengths, which
    # is normally distributed )
    # there is a approximate normal distribution for shorter genomic lengths
    # and then there is a low long tail, i.e. ests are distributed
    # homogeneously for long genomic extents - which also means there are fewer
    # of them.
    # The mean value usually falls outside the normal distribution, but the
    # median seems to correlate quite well with the peak of the normal part.
    # We therefore take the median as threshold. This will still keep some
    # short ESTs, as it is less restrictive than the mean, but it will be more representative.
    my $mean = 0;
    my $median = 0;
    my $strand = 1;
    foreach my $est ( @$ests ){
      #print STDERR "est is a $est\n";
      #print STDERR "exons: ".scalar( @{$est->get_all_Exons} )."\n";
      my @exons = sort { $a->start <=> $b->start } @{$est->get_all_Exons};
      $strand = -1 if $exons[0]->strand == -1;
      #print STDERR $exons[0]->strand." ".$exons[0]->slice->strand."\n";
      # we consider the genomic extension of the transcripts
      # this will not favour big unspliced ESTs over spliced ones with small exons
      #print STDERR "end ".$exons[-1]->end." start ".$exons[0]->start."\n";
      my $length = $exons[-1]->end - $exons[0]->start + 1 ;  
      push( @length_list, $length );
      push( @{ $lengths{ 10*int($length/10 + 1) } }, $est );
      $mean += $length;
    }
    $mean /= scalar( @$ests );
    print STDERR "mean length = $mean\n" if $verbose;
    print STDERR "number of exon lengths ".@length_list."\n";
    my @sorted_lengths = sort{$a <=> $b} @length_list;
    if ( scalar(@sorted_lengths)>2 ){
      print @sorted_lengths."\n";
      if ( scalar(@sorted_lengths)%2 == 0 ){
        my $i1 = scalar(@sorted_lengths)/2;
        my $i2 = $i1 - 1;
        $median = ($sorted_lengths[$i1] + $sorted_lengths[$i2])/2;
      }else{
        my $i = ( scalar(@sorted_lengths) -1 )/2;
        $median = $sorted_lengths[$i];
      }
      print STDERR $median."\n";
    }else{
      $median = $mean;
    }
    #print STDERR "median length = $median\n" if $verbose;
    
    #test
    if ($verbose && $strand == 1){
      foreach my $key ( sort { $a <=> $b } keys %lengths ){
        #print STDERR "$key\t".scalar( @{ $lengths{$key} })."\n";
      }
    }
  
    my @selected;
    foreach my $key ( keys %lengths ){
      if($key >= int($median)){
        push( @selected, @{ $lengths{$key} });
      }
    }
    return @selected;
}

############################################################

sub depth_filter {
  my ($self,$exons,$exon2est) = @_;
  
  my $verbose = $self->verbose;
  
  print STDERR "filtering ".scalar(@$exons)." exons\n" if $verbose;
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
  my $exon_cluster = new Bio::EnsEMBL::Feature;
  
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
           && ( $exon->strand == $exon_cluster->strand)) { 
        $exon_cluster->add_sub_SeqFeature($exon,'EXPAND');
      } else {
        # Start a new cluster
        $exon_cluster = new Bio::EnsEMBL::Feature;
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
    sort { $self->_average_score($b) <=> $self->_average_score($a) or
             $self->_average_perc_id($b) <=> $self->_average_perc_id($a) or
               $self->_order($b) <=> $self->_order($a) } 
      @$cluster_list;
  
  my $depth = $self->depth_threshold;
  
  my %banned;
  my %taken;
  my @accepted_ests;
  
  print STDERR "looking at depth of ".scalar(@sorted_clusters)." clusters\n" if $verbose;
  foreach my $cluster ( @sorted_clusters ){
    my @exons = sort {($self->_exon_score($b) <=> $self->_exon_score($a)) 
                        or ($self->_exon_perc_id($b) <=> 
                            $self->_exon_perc_id($a)) or
                              ($self->_length( $exon2est{$a} ) 
                               <=> $self->_length( $exon2est{$b}))} 
      $cluster->sub_SeqFeature;
    
    my $count = 0;
    while ( @exons ){
      my $exon = shift @exons;
      my $est  = $exon2est{$exon};
      if ( $count >= $depth ){
        $banned{$est} = 1 ;
        if ($verbose){
          #print STDERR "rejecting: ";
          #$self->_print_EST($est);
        }
      }
      
      unless ( $banned{$est} ){
        unless ( $taken{$est} ){
          push( @accepted_ests, $est );
          $taken{$est} = 1;
          if ($verbose){
            #print STDERR "taking: ";
            #$self->_print_EST($est);
          }
        }
        $count++;
      }
    }
  }
  
  if ($verbose){
    print STDERR "returning ".scalar(@accepted_ests)." ests after scores and depth filtering\n";
    foreach my $est ( @accepted_ests ){
      # print STDERR $self->_id($est)." coverage:".$self->_coverage($est)." perc_id:".$self->_perc_id($est)."\n";
      # Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript($est);
    }
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

sub _print_EST{
  my ($self,$t) = @_;
  my $id      = $self->_id($t);
  my $score   = $self->_coverage($t);
  my $perc_id = $self->_perc_id($t);
  my $length  = $self->_length($t);
  print STDERR "$id coverage:$score perc_id:$perc_id length:$length\n";
}

############################################################

sub _length{
  my ($self,$t) = @_;
  my @exons = sort{ $a->start <=> $b->start } @{$t->get_all_Exons};
  my $length = $exons[-1]->end - $exons[0]->start + 1;
  return $length;
}

############################################################

1;
