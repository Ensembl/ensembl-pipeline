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
    
    my( $hold_list, $transcripts ) = $self->_rearrange([qw(
							   HOLD_LIST
							   TRANSCRIPTS
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
    return $self;
}

############################################################
# $list is an arrayref with the transcripts that make up the transcript $tran
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
    if ( $tran ){
	$self->{_predictions} = $trans;
    }
    return $self->{_predictions};
}

############################################################

sub score_Transcripts{
    my ($self) = @_;
    
    ############################################################
    # first need  to estimate how many sites of 
    # alternative splicing there are per transcript

    ############################################################
    # cluster transcripts into genes according to exon overlap
    # but without make exons unique:
    my @clusters = $self->cluster_Transcripts( $self->transcripts );
    

}

############################################################

sub get_transcript_start_end_strand {
  my ($t) = @_;
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
	    $self->get_transcript_start_end_strand($tran);
	
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


############################################################
# Method to score each transcript
#
# If the transcript only has 1 site of alternative splicing
# we give a 100 score
#
# else the score depends on the length of the ESTs relative
# to the length between the two or more sites of alternative splicing
#

sub _score_Transcript{
    my ($tran, @other ) = @_;
    
    my @list = @{ $self->hold_list($tran) };
    
    ############################################################
    # calculate the lengths of the ESTs used:
    my $average = 0;
    foreach my $est ( @list ){
	
	# est is a transcript object
	my $length = $est->length;
	
	$average += $length;
    }
    if ( @list ){
	$average = int( $average/scalar(@list) );
    }
    
}
	
############################################################

1;
