=head1 NAME - Bio::EnsEMBL::Pipeline::GeneComparison::ExonComparison;

=head1 DESCRIPTION

Perl Class for comparison of two sets of genes at the exon level.

=head1 SYNOPSIS

  my $exon_comparison = Bio::EnsEMBL::Pipeline::GeneComparison::ExonComparison
                        ->new(
                               -annotation => \@annotation_genes,
                               -prediction => \@prediction_genes,
                               -annotation_types => ['type1',...],
                               -prediction_types => ['type2',...],
                             );

  $exon_comparison->compare_exons_base_pair_level();

=head1 CONTACT

eae@sanger.ac.uk

=cut

# Let the code begin ...

package Bio::EnsEMBL::Pipeline::GeneComparison::ExonComparison;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::GeneComparison::GeneCluster;
use Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster;
use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::GeneComparison::GeneCompConf;

@ISA = qw(Bio::EnsEMBL::Root);

####################################################################################

=head2 new()

the new() method accepts two array references

=cut

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);
  
    my( $annotation, $prediction, $annotation_types, $prediction_types, $parameter ) 
	= $self->_rearrange([qw(
				ANNOTATION
				PREDICTION
				ANNOTATION_TYPES
				PREDICTION_TYPES
				PARAMETER
				)], 
			    @args);
    
    unless( $annotation && $prediction && $annotation_types && $prediction_types ){
	$self->throw("need anotation and prediction genes and types");
    }
    
    $self->annotation_Genes(@$annotation);
    $self->prediction_Genes(@$prediction);
    $self->annotation_Types($annotation_types);
    $self->prediction_Types($prediction_types);
    

    if ( $parameter ){
	$self->parameter( $parameter );
    }
    
    return $self;
}

############################################################

sub parameter{
    my ($self,$value) = @_; 
	if ( $value ){
	    $self->{'_parameter'} = $value;
	}
    return $self->{'_parameter'};
}

######################################################################################

sub annotation_Genes {
  my ($self,@genes) = @_;
  if ( @genes ){
    push ( @{ $self->{'_annotation_genes'} }, @genes );
  }
  return $self->{'_annotation_genes'};
}

############################################################

sub annotation_Types{
    my ($self,$types) = @_;
    if ($types){
	$self->{'_annotation_Types'} = $types;
    }
    return $self->{'_annotation_Types'};
}

sub prediction_Types{
    my ($self,$types) = @_;
    if ($types){
	$self->{'_annotation_Types'} = $types;
    }
    return $self->{'_prediction_Types'};
}

######################################################################################

sub prediction_Genes {
  my ($self,@genes) = @_;
  if ( @genes ){
    push ( @{ $self->{'_prediction_genes'} }, @genes );
  }
  return $self->{'_prediction_genes'};
}

######################################################################################

sub compare_exons_base_pair_level{
    my ($self,$coding) = @_;
    
    my @annotated_genes = @{$self->annotation_Genes};
    my @predicted_genes = @{$self->prediction_Genes};
    my @annotated_exons;
    my @predicted_exons;

    #my $gene_comparison = 
    #Bio::EnsEMBL::Pipeline::GeneComparison::GeneComparison
    #    ->new(
    #	    '-annotation_genes' => $self->annotation_Genes,
    #	    '-prediction_genes' => $self->prediction_Genes,
    #	    );
    
    ############################################################
    ## cluster the genes we have passed to $gene_comparison
    #   my @gene_clusters    = $gene_comparison->cluster_Genes;
    
    ############################################################
    my $total_ann_length;
    my $total_pred_length;
    my $covered;
    my $missed;
    my $overpredicted;
    ############################################################
    
    # CLUSTER:
    #   foreach my $cluster ( @gene_clusters ){
    #	my @annotated_exons;
    #	my @predicted_exons;
    #
    #	my @annotated_genes;
    #	my @predicted_genes;
    #
    #	foreach my $type ( @{$self->annotated_Types} ){
    #	    push( @annotated_genes , $cluster->get_Genes_by_Type( $type ) );
    #	}
    #	foreach my $type ( @{$self->predicted_Types} ){
    #	    push( @predicted_genes , $cluster->get_Genes_by_Type( $type ) );
    #	}
    
    unless ( $coding ){
      
      ############################################################
      # get all annotated exons
      foreach my $gene ( @annotated_genes ){
	push( @annotated_exons, @{$gene->get_all_Exons} );
      }
      
      ############################################################
      # get all predicted exons and cluster them
      foreach my $gene ( @predicted_genes ){
	push( @predicted_exons, @{$gene->get_all_Exons} );
      }
      
    }
    if ($coding ){
      
      ############################################################
      # get all annotated exons
      my %ann_seen;
      foreach my $gene ( @annotated_genes ){
	foreach my $tran ( @{$gene->get_all_Transcripts} ){
	  if ( defined $tran->translation ){
	    foreach my $exon ( @{$tran->get_all_translateable_Exons} ){
	      next if $ann_seen{$exon->start}{$exon->end}{$exon->strand};
	      push( @annotated_exons, $exon );
	      $ann_seen{$exon->start}{$exon->end}{$exon->strand} = 1;
	    }
	  }
	}
      }
      
      ############################################################
      # get all predicted exons and cluster them
      my %pred_seen;
      foreach my $gene ( @predicted_genes ){
	foreach my $tran ( @{$gene->get_all_Transcripts} ){
	  if ( defined $tran->translation ){
	    foreach my $exon ( @{$tran->get_all_translateable_Exons} ){
	      next if $pred_seen{$exon->start}{$exon->end}{$exon->strand};
	      push( @predicted_exons, $exon );
	      $pred_seen{$exon->start}{$exon->end}{$exon->strand} = 1;
	    }
	  }
	}
      }
    }
    
    
    my @ann_clusters  = @{$self->_cluster_Exons(@annotated_exons)};
    my @pred_clusters = @{$self->_cluster_Exons(@predicted_exons)};
    
    ############################################################
    # for each annotated exon, find the overlapping exon_cluster
    # and calculate how much we cover it and how much we exceed it
    
    my %seen_pred;
    foreach my $ann_cluster ( @ann_clusters ){
      
      $total_ann_length  += ( $ann_cluster->end - $ann_cluster->start + 1 );
      foreach my $pred_cluster ( @pred_clusters ){
	
	unless( $seen_pred{$pred_cluster} ){
	  $total_pred_length += ( $pred_cluster->end - $pred_cluster->start + 1 );
	}
	$seen_pred{$pred_cluster} = 1;
	
	next if ( $pred_cluster->start > $ann_cluster->end 
		  || 
		  $pred_cluster->end < $ann_cluster->start );
	
	next if ( $pred_cluster->strand != $ann_cluster->strand );
	
	$covered += $self->min($ann_cluster->end,$pred_cluster->end) - $self->max($ann_cluster->start,$pred_cluster->start) + 1;
	
	$missed += $self->max(0,$pred_cluster->start - $ann_cluster->start);
	$missed += $self->max(0,$ann_cluster->end - $ann_cluster->end);
	
	$overpredicted += $self->max(0,$ann_cluster->start - $pred_cluster->start);
	$overpredicted += $self->max(0,$pred_cluster->end - $ann_cluster->end);
      }
    }
    #} # end of CLUSTER
	
	
    print STDERR "total_annotated: $total_ann_length bp\n";
    print STDERR "total_predicted: $total_pred_length bp\n";
    print STDERR "covered        : $covered\n";
    print STDERR "missed         : $missed\n";
    print STDERR "overpredicted  : $overpredicted\n";
    print STDERR "\n";
    
    my $Sn = sprintf "%.2f", ( $covered/$total_ann_length );
    my $Sp = sprintf "%.2f", ( $covered/$total_pred_length );

    print STDERR "At the base pair level:\n";
    print STDERR "Sn = $Sn\t Sp = $Sp\n";

    ############################################################
    # calculate the unclustered exons:
    
    #my @unclustered = $gene_comparison->unclustered_Genes;
    #my @ann_unclustered;
    #my @pred_unclustered;
    
  #UNCLUSTER:
  #  foreach my $uncluster ( @unclustered ){
#	my @gene = $uncluster->get_Genes;
#	#if ( scalar(@gene)>1 ){
#	#    print STDERR "genes @gene are wrongly unclustered\n";
#	#}
#	my $this_type = $gene[0]->type;
#	foreach my $type ( @{ $self->annotation_Types } ){
#	    if ( $this_type eq $type ){
#		push( @ann_unclustered, $uncluster );
#		next UNCLUSTER;
#	    }
#	}
#	foreach my $type ( @{ $self->prediction_Types } ){
#	    if ( $this_type eq $type ){
#		push( @prediction_unclustered, $uncluster );
#		next UNCLUSTER;
#	    }
#	}
#    }
    

}

############################################################

sub max{
    my ($self, $max, @rest ) = @_;
    for ( my $i=0; $i<=$#rest; $i++ ){
	$max = $rest[$i] if $rest[$i]>$max;
    }
    return $max;
}

############################################################

sub min{
    my ($self, $min, @rest ) = @_;
    for ( my $i=0; $i<=$#rest; $i++ ){
	$min = $rest[$i] if $rest[$i]<$min;
    }
    return $min;
}

############################################################

sub _cluster_Exons{
    my ($self, @exons) = @_;
    
    # no point if there are no exons!
    return unless ( scalar( @exons) > 0 );   
    
    # keep track about in which cluster is each exon
    my %exon2cluster;
    
    # main cluster feature - holds all clusters
    my $cluster_list;
    
    # sort exons by start coordinate
    @exons = sort { $a->start <=> $b->start } @exons;
    
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
	    if ( $exon_cluster->overlaps($exon) && ( $exon->strand == $exon_cluster->strand) ) { 
		$exon_cluster->add_sub_SeqFeature($exon,'EXPAND');
	    }  
	    else {
		# Start a new cluster
		$exon_cluster = new Bio::EnsEMBL::SeqFeature;
		$exon_cluster->add_sub_SeqFeature($exon,'EXPAND');
		$exon_cluster->strand($exon->strand);
		
		# and add it to the main_cluster feature
		push( @$cluster_list, $exon_cluster);
	    }
	}
	$count++;
    }
    return $cluster_list;
}

1;
