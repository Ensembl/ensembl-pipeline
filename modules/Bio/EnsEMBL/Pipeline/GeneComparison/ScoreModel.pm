#
# Written by Eduardo Eyras
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::ScoreModel

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 AUTHOR

eae@sanger.ac.uk

=head1 CONTACT

ensembl-dev@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::GeneComparison::ScoreModel;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster;
use Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptComparator;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;

@ISA = qw(Bio::EnsEMBL::Root);

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);
    
    my( $hold_list, $transcripts, $label ) = $self->_rearrange([qw(
								   HOLD_LIST
								   TRANSCRIPTS
								   LABEL
								  )], 
							       @args);
    
    ############################################################
    # hold_list is a hashref where we hold for each transcript
    # the list of ESTs/cDNAs used to build it
    unless( $hold_list ){
	$self->throw("Need the hashref hold_list, cannot work without the lists");
    }
    $self->{_hold_list} = $hold_list;
    
    unless( $transcripts ){
	$self->throw("Need the list of transcript predictions");
    }
    $self->transcripts($transcripts);
    
    if ( defined $label ){
      $self->_label($label);
    }

    $self->verbose(0);

    return $self;
  }

############################################################

sub _label{
  my ($self,$label) = @_;
  if ( defined $label ){
    $self->{_label} = $label;
  }
  return $self->{_label};
}

############################################################

sub verbose{
  my ($self, $boolean) = @_;
  if ( defined $boolean ){
    $self->{_verbose} = $boolean;
  }
  return $self->{_verbose};
}

############################################################
#
# $list is an arrayref with the transcripts that make up the transcript $tran
#
sub hold_list{
  my ($self,$tran,$list) = @_;
  
  if ($tran){
    unless ( $self->{_hold_list}{$tran} ){
      $self->{_hold_list}{$tran} = [];
    }
    if ( $list ){
      $self->{_hold_list}{$tran} = $list;
    }
    return $self->{_hold_list}{$tran};
  }
  else{
    return $self->{_hold_list};
  }
}

############################################################
# it holds an arrayref with the predictions
sub transcripts{
  my ($self, $trans ) = @_;
  if ( $trans ){
    $self->{_predictions} = $trans;
  }
  return $self->{_predictions};
}

############################################################

sub score_Transcripts{
  my ($self) = @_;

  print STDERR "**************** label : ".$self->_label."\n";
  my $verbose = $self->verbose;
  
  ############################################################
  # first need  to estimate how many sites of 
  # alternative splicing there are per transcript
  
  ############################################################
  # cluster transcripts into genes according to exon overlap
  # but without make exons unique:
  my @clusters = @{$self->cluster_Transcripts( $self->transcripts )};
  
  ############################################################
  # get the sites of alt-splicing on each transcript cluster
  my $cluster_count = 0;
  CLUSTER:
  foreach my $cluster ( @clusters ){

    $cluster_count++;
    my $label;
    
    my @trans = @{ $cluster->get_Transcripts };
    print STDERR "Looking at cluster:\n" if $verbose;
    if ( $verbose ){
      foreach my $tran ( @trans ){
	Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript($tran);
      }
    }
    
    if (scalar(@trans)==1 ){
      print STDERR "1 transcript-> all exons get 100 as score\n" if $verbose;
      $label = '';
      if ( $self->_label ){
	$label = $self->_label;
      }
      my $tran_id = $label."_".$cluster_count."_1";
      print STDERR "transcript: $tran_id (single transcript)\n";
      $trans[0]->stable_id($tran_id);
      
      foreach my $exon ( @{$trans[0]->get_all_Exons} ){
	$exon->score( 100 );
      }
      next CLUSTER;
    }
    
    # my $exon_clusters_count; # counts howmany exon clusters are in this transcript clusters

    my ($sites, $exon_clusters_count ) = $self->get_alternative_sites( $cluster );
    my $average_score = 0;
    my $average_missed_sites;
    my %trans_with_site; # counts how many transcripts have this site
    
    ############################################################
    # now get the sites of alternative splicing
    # contained in each transcript
    # and then calculate the distances between 
    # those sites and the lengths
    # of the list of ESTs making up the transcript
    # and compare them!
    my $tran_count = 0;
  TRAN:
    foreach my $tran ( @trans ){

      # invent a tran stable id with the slice id, cluster_id and a number
      $tran_count++;
      $label = '';
      if ( $self->_label ){
	$label = $self->_label;
      }
      my $tran_id = $label."_".$cluster_count."_".$tran_count;
      $tran->stable_id($tran_id);
      print STDERR "transcript: $tran_id\n";
      # list of ESTs:
      my @list = @{ $self->hold_list($tran) };
      
      print STDERR "finding sites in transcript: " if $verbose;
	if ($verbose ){
	  Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript( $tran );
	}
	############################################################
	# which sites does this transcript have?
	# $site is a SeqFeature with exons as sub_SeqFeatures 
	my @these_sites;
	foreach my $site ( @$sites ){
	    my ($start,$end,$strand) = $self->get_start_end_strand_of_transcript($tran);
	    if ( !( $site->start > $end) && !( $site->end < $start ) ){
		push( @these_sites, $site );
		$trans_with_site{$site}++;
	    }
	}
	
	
	############################################################
	# @these_sites contains now all the sites
	# present in the transcript.
	# Now we follow a greedy approach to find out how many
	# potential 'splits' ( pairs of sites unlinked ) in the evidence there are.
	
	my $n = scalar( @these_sites );
	my $inv_score = 1;
	
	# we consider arrays of consecutive sites
	# from the largest possible (n) to the smallest (2)
	my $covered_sites = $n;
      LEVEL:
	for (my $i=0; $i<$n-1; $i++ ){
	    
	  SUBLEVEL:
	    for (my $diff = 0; $diff<=$i; $diff++){
		
		print STDERR "sites: " if $verbose;
		my @site_combination;             
		for ( my $j=0; $j<$n-$i; $j++ ){
		  print STDERR ($j+$diff)." " if $verbose;
		    push( @site_combination, $these_sites[$j+$diff] );
		}
		print STDERR "\n" if $verbose;
		
		############################################################
		# check whether this transcript contains ESTs
		# covering all the sites in @site_combination
		my $covered = 0;
		
		my @covered_sites;
		
	      EST:
		foreach my $est ( @list ){
		    my ($est_start, $est_end, $est_strand) = $self->get_start_end_strand_of_transcript( $est );
		    
		    ############################################################
		    # by construction, it is enough to check that the EST 
		    # covers both extreme sites to know that it covers all in the middle,
		    if ( $est_start <= $site_combination[0]->end 
			 &&
			 $est_end   >  $site_combination[0]->end
			 &&
			 $est_end   >= $site_combination[-1]->start
			 &&
			 $est_start <  $site_combination[-1]->start
			 ){
			$covered = 1;
			if ($verbose ){
			    print STDERR "covered by est: $est_start-$est_end\n";
			  Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils
			      ->_print_SimpleTranscript( $est );
			}
			last LEVEL;
		    }
		    
		} # end of EST
	    }     # end of SUBLEVEL
	    
	    ############################################################
	    # everytime we leave a level unsuccessful, we 
	    # have a certainty loss factor of 1/2 and the
	    # (maximum) number of covered sites is reduced by one
	    $inv_score *= 2;
	    $covered_sites--;

	} # end of LEVEL
	
	my $score = sprintf "%.2f", ( 100/$inv_score );
	$average_score += $score;
	$average_missed_sites += ( $n - $covered_sites );
	my $exons = scalar( @{$tran->get_all_Exons} );

	############################################################
	# TRAN number_ests number_sites max_num_sites_covered:
	
	print STDERR "TRAN\t".
	  $tran_id."\t".
	    "exons:".$exons."\t".
	      "ests:".scalar(@list)."\t".
		"sites:".$n."\t".
		  "covered:".$covered_sites."\t".
		    "score:".$score."\t".
		      "s-c:".($n - $covered_sites)."\n";
      
      ############################################################
      # put the transcript score in the exons:
      foreach my $exon ( @{$tran->get_all_Exons} ){
	$exon->score( $score );
      }
    } # end of TRAN
    
    my $trans_number  = scalar( @trans );
    
    # we calculate on each site whether it is covered
    # by the same est in different transcripts
    my %site_coverage;
  SITE:
    foreach my $site ( @$sites ){
      my %used_est;
      foreach my $tran ( @trans ){
	my %seen_est;
	foreach my $est ( @{ $self->hold_list($tran) } ){
	  next if $seen_est{$est};
	  my ($est_start, $est_end, $est_strand) = $self->get_start_end_strand_of_transcript( $est ); 
	  unless (  $est_start > $site->end || $est_end   <  $site->start ){
	    $used_est{$est}++;
	    $seen_est{$est} = 1;
	    $site_coverage{$site}++ unless $used_est{$est}>1;
	  }
	}
      }
      my %bin_used_est;
      foreach my $key ( keys %used_est ){
	$bin_used_est{ $used_est{$key} }++;
      }
      my $site_string;
      my @keys = sort { $a <=> $b } keys %bin_used_est;
      for (my $i=1; $i<=$keys[-1]; $i++ ){
	my $string = "ests_used$i:0\t";
	if ( $bin_used_est{$i} ){
	  $string = "ests_used$i:$bin_used_est{$i}\t";
	}
	$site_string .= $string;
      }
      
      print STDERR "SITE\tcoverage:$site_coverage{$site}\t".
	"trans:$trans_number\t".
	  "trans_with_site:$trans_with_site{$site}\t".
	    "$site_string\n";
   
    } # end of SITE

    ############################################################
    # cluster info:
    my $cluster_sites = scalar( @$sites );
    my $max_sites     = 2 ** $cluster_sites;
    $average_score   /= $trans_number;
    $average_missed_sites /= $trans_number;

    my $gene_id = $label."_".$cluster_count;

    print STDERR "GENE\t".
      $gene_id."\t".
	"sites:".$cluster_sites."\t".
	  "trans:".$trans_number."\t".
	    "exon_clust:".$exon_clusters_count."\t".
	      "2^N:".$max_sites."\t".
		"av_score:".$average_score."\t".
		  "av_missed_sites:".$average_missed_sites."\n";
    
    
  }   # end of CLUSTER
  return @{$self->transcripts};
}

############################################################
#
# this method takes a cluster of transcripts which
# can make up a gene. It cluster the exons and
# walks along the exon clusters to find out where are
# the sites of alternative splicing. Each of these
# sites is modelled as an exon-cluster which contains
# the exons in that position (contitutive and alternative exons )
# and has a position in the genomic coordinates, so that we
# do not lose information of the relative position between
# each of them in the genomic

# Note: this method will not be able to
# see a site of alternative splicing which it
# has been actually predicted not to have alternative splicing.
# Thus this is not the ed of the story.

sub get_alternative_sites{
  my ($self, $cluster ) = @_;

  my $verbose = $self->verbose;

  print STDERR "finding alternative sites\n" if $verbose;
  my @trans = @{ $cluster->get_Transcripts };
  my %exon2transcript;
  my @all_exons;
  foreach my $tran ( @trans ){
    foreach my $exon ( @{$tran->get_all_Exons} ){
      $exon2transcript{ $exon } = $tran;
      push( @all_exons, $exon );
    }
  }
  
  ############################################################
  # cluster the exons according to overlap
  my $exon_cluster_list = $self->_cluster_Exons( @all_exons );
  my @clusters = sort { $a->start <=> $b->start } $exon_cluster_list->sub_SeqFeature;
  my $exon_clusters_count = scalar(@clusters);

  ############################################################
  # get the sites of alternative splicing:
  my @sites;
  my %added;
  my $exon_cluster_pos = 0;
 EXON_CLUSTER:
  foreach my $exon_cluster ( @clusters ){
    $exon_cluster_pos++;
    my @exons = $exon_cluster->sub_SeqFeature;
    my $previous_exon;
    my %seen_transcript;
    
    print STDERR "cluster with ".scalar(@exons)." exons\n" if $verbose;
  EXON:
    while ( @exons ){
      my $exon = shift @exons;
      my $found = 0;
      
    TRAN:
      foreach my $t ( @trans ){
	#print STDERR "exon2trans: ".$exon2transcript{ $exon }." looking at trans: $t\n";
	if ( $t == $exon2transcript{ $exon } ){
	  $seen_transcript{$t} = 1;
	  #print STDERR "transcript found\n";
	  last TRAN;
	}
      }
      
      ############################################################
      # unless all exons have the same splice coordinates
      # we have a case of alternative splicing: alterantive 3'/5' site, intron retention
      if ( $previous_exon ){
	unless ( $previous_exon->start == $exon->start && $previous_exon->end == $exon->end ){
	  unless ( $added{$exon_cluster} ){
	    if ( $exon_cluster_pos == 1 || $exon_cluster_pos == $#clusters ){
	      print STDERR "Alternative terminal exon: " if $verbose;
	    }
	    else{
	      print STDERR "Alternative site: previous_exon: " if $verbose;
	    }
	    print STDERR $previous_exon->start."-".$previous_exon->end." ".
	      "exon: ".
		$exon->start."-".$exon->end."\n" if $verbose;
	    $added{ $exon_cluster }++;
	    push( @sites, $exon_cluster );
	  }
	}
      }
      $previous_exon = $exon;
      
    } # end of EXON
    ############################################################
    # if there are transcripts that had no exons in this cluster
    # we have a case of exon skipping:
    if ( scalar( @trans ) != scalar ( keys %seen_transcript ) ){
      unless ( $added{$exon_cluster} ){
	print STDERR "Exon skipping: ".$exon_cluster->start."-".$exon_cluster->end."\n" if $verbose;
	print STDERR "cluster with ".scalar(@trans)." transcripts\n" if $verbose;
	foreach my $tran ( @trans ){
	  Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript($tran);
	}	
	push( @sites, $exon_cluster );
      }
    }
  }   # end of EXON_CLUSTER
  
  ############################################################
  # sites of alternative splicing are described by a cluster of exons 
  # which has genomic coordinates
  print STDERR scalar(@sites)." sites found\n";
  return (\@sites, $exon_clusters_count);
}
  
############################################################

=head2 _cluster_Exons
 
 Function: it cluster exons according to overlap,
           it returns a Bio::EnsEMBL::SeqFeature, where the sub_SeqFeatures
           are exon_clusters, which are at the same time Bio::EnsEMBL::SeqFeatures,
           whose sub_SeqFeatures are exons
=cut

sub _cluster_Exons{
  my ($self, @exons) = @_;
  
  # no point if there are no exons!
  return unless ( scalar( @exons) > 0 );   

  # keep track about in which cluster is each exon
  my %exon2cluster;
  
  # main cluster feature - holds all clusters
  my $cluster_list = new Bio::EnsEMBL::SeqFeature; 
  
  # sort exons by start coordinate
  @exons = sort { $a->start <=> $b->start } @exons;

  # Create the first exon_cluster
  my $exon_cluster = new Bio::EnsEMBL::SeqFeature;
  
  # Start off the cluster with the first exon
  $exon_cluster->add_sub_SeqFeature($exons[0],'EXPAND');
  $exon_cluster->strand($exons[0]->strand);    
  $cluster_list->add_sub_SeqFeature($exon_cluster,'EXPAND');
  
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
	$cluster_list->add_sub_SeqFeature($exon_cluster,'EXPAND');	
      }
    }
    $count++;
  }
  return $cluster_list;
}


############################################################
#
# this method gets the start and end of transcript meaning:
# start: lowest coordinate
# end  : highest coordinate

sub get_start_end_strand_of_transcript {
  my ($self,$t) = @_;
  my @exons = sort { $a->start <=> $b->start } @{$t->get_all_Exons}; 
  my $start = $exons[0]->start;
  my $end   = $exons[-1]->end;
  return ($start, $end, $exons[0]->strand);
}

############################################################

sub by_transcript_low {
    my @aexons = sort { $a->start <=> $b->start } @{$a->get_all_Exons};
    my @bexons = sort { $a->start <=> $b->start } @{$b->get_all_Exons};
    my $alow = $aexons[0]->start;
    my $blow = $bexons[0]->start;
    return $alow <=> $blow;
}

############################################################

sub cluster_Transcripts {
    my ($self,$trans) = @_;
    my $forward_trans;
    my $reverse_trans;
    foreach my $t ( @$trans ){
	if ( $t->start_Exon->strand == 1 ){
	    push ( @{$forward_trans}, $t );
	}
	else{
	    push ( @{$reverse_trans}, $t );
	}
    }
    my @clusters;
    if ( $forward_trans ){
	my $f_clusters = $self->cluster_Transcripts_by_strand($forward_trans);
	push( @clusters, @{$f_clusters} );
	print STDERR scalar( @{$f_clusters} ). " clusters on forward strand\n";
    }
    if ( $reverse_trans ){
	my $r_clusters = $self->cluster_Transcripts_by_strand($reverse_trans);
	push( @clusters, @{$r_clusters} );
	print STDERR scalar( @{$r_clusters} ). " clusters on reverse strand\n";
    }
    return \@clusters;
}

############################################################

sub cluster_Transcripts_by_strand {
    my ($self, $transcripts_unsorted) = @_;
    my @transcripts = sort by_transcript_low @$transcripts_unsorted;
    my @clusters = ();
    my $verbose = 0;
    
    # clusters transcripts by whether or not any exon overlaps with an exon in
    foreach my $tran (@transcripts) {
	my @matching_clusters;
	my ($trans_start, $trans_end, $trans_strand) = 
	    $self->get_start_end_strand_of_transcript($tran);
	
	print STDERR "transcript limits: $trans_start $trans_end \n" if $verbose;
	
      CLUSTER: 
	foreach my $cluster (@clusters) {
	    
	    print STDERR "Testing against cluster with limits " . $cluster->start ." to " . $cluster->end . "\n" if $verbose;
	    
	    if (!($trans_start > $cluster->end || $trans_end < $cluster->start) &&
		$trans_strand == $cluster->strand) {
		
		print STDERR "In range\n" if $verbose;
		foreach my $cluster_transcript (@{$cluster->get_Transcripts()}) {
		    foreach my $exon1 (@{$tran->get_all_Exons}) {
			
			foreach my $cluster_exon (@{$cluster_transcript->get_all_Exons}) {
			    if ($exon1->overlaps($cluster_exon) &&
				$exon1->strand == $cluster_exon->strand) {
				push (@matching_clusters, $cluster);
				next CLUSTER;
			    }
			}
		    }
		}
	    }
	}
	
	if (scalar(@matching_clusters) == 0) {
	    print STDERR "Found new cluster for " . $tran->dbID . "\n" if $verbose;
	    my $newcluster = new Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster;
	    $newcluster->put_Transcripts($tran);
	    push(@clusters,$newcluster);
	    
	} 
	elsif (scalar(@matching_clusters) == 1) {
	    print STDERR "Adding to cluster for " . $tran->dbID . "\n" if $verbose;
	    $matching_clusters[0]->put_Transcripts($tran);
	    
	} 
	else {
	    # Merge the matching clusters into a single cluster
	    print STDERR "Merging clusters for " . $tran->dbID . "\n" if $verbose;
	    my @new_clusters;
	    my $merged_cluster = new Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster;
	    foreach my $clust (@matching_clusters) {
		$merged_cluster->put_Transcripts(@{$clust->get_Transcripts});
	    }
	    $merged_cluster->put_Transcripts($tran);
	    push @new_clusters,$merged_cluster;
	    # Add back non matching clusters
	    foreach my $clust (@clusters) {
		my $found = 0;
	      MATCHING: 
		foreach my $m_clust (@matching_clusters) {
		    if ($clust == $m_clust) {
			$found = 1;
			last MATCHING;
		    }
		}
		if (!$found) {
		    push @new_clusters,$clust;
		}
	    }
	    @clusters = @new_clusters;
	}
    }
    
    #Safety checks
    my $ntrans = 0;
    my %trans_check_hash;
    foreach my $cluster (@clusters) {
	$ntrans += scalar(@{$cluster->get_Transcripts});
	foreach my $trans (@{$cluster->get_Transcripts}) {
	    if (defined($trans_check_hash{"$trans"})) {
		print STDERR ("Transcript " . $trans->dbID . " added twice to clusters\n") ;
	    }
	    $trans_check_hash{"$trans"} = 1;
	}
	if (!scalar(@{$cluster->get_Transcripts})) {
	    print STDERR ("Empty cluster");
	}
    }
    
    if ($ntrans != scalar(@transcripts)) {
	print STDERR ("Not all transcripts have been added into clusters $ntrans and " . scalar(@transcripts). " \n");
    }
    #end safety checks
    
    #print "Started with " . scalar(@transcripts) . " transcripts and ended with " . scalar(@clusters) . " clusters\n";
    
    return \@clusters;
}

############################################################


#############################################################
## Method to score each transcript
##
## If the transcript only has 1 site of alternative splicing
## we give a 100 score
##
## else the score depends on the length of the ESTs relative
## to the length between the two or more sites of alternative splicing
##

#sub _score_Transcript{
#    my ($tran, @other ) = @_;
    
#    my @list = @{ $self->hold_list($tran) };
    
#    ############################################################
#    # calculate the lengths of the ESTs used:
#    my $average = 0;
#    foreach my $est ( @list ){
	
#	# est is a transcript object
#	my $length = $est->length;
	
#	$average += $length;
#    }
#    if ( @list ){
#	$average = int( $average/scalar(@list) );
#    }
    
#}
	
#############################################################

1;
