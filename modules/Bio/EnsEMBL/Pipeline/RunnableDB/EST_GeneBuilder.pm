#
#
# Cared for by EnsEMBL  <ensembl-dev@ebi.ac.uk>
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::EST_GeneBuilder

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::RunnableDB::EST_GeneBuilder->new(
								       -dbobj     => $db,
								       -input_id  => $id
								      );
    $obj->fetch_input
    $obj->run

    my @newfeatures = $obj->output;


=head1 DESCRIPTION

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::RunnableDB::EST_GeneBuilder;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::RootI;
use Bio::EnsEMBL::Pipeline::RunnableDBI;
use Bio::EnsEMBL::Pipeline::Runnable::Genomewise;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Getz;

use Bio::EnsEMBL::Pipeline::GeneConf qw (EXON_ID_SUBSCRIPT
					 TRANSCRIPT_ID_SUBSCRIPT
					 GENE_ID_SUBSCRIPT
					 PROTEIN_ID_SUBSCRIPT
					 );

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);

    # dbobj input_id mandatory and read in by BlastableDB
    if (!defined $self->seqfetcher) {
      my $seqfetcher = new Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;
      $self->seqfetcher($seqfetcher);
    }

    return $self; 
}

=head2 RunnableDB methods

=head2 analysis

    Title   :   analysis
    Usage   :   $self->analysis($analysis);
    Function:   Gets or sets the stored Analusis object
    Returns :   Bio::EnsEMBLAnalysis object
    Args    :   Bio::EnsEMBL::Analysis object

=head2 dbobj

    Title   :   dbobj
    Usage   :   $self->dbobj($obj);
    Function:   Gets or sets the value of dbobj
    Returns :   A Bio::EnsEMBL::Pipeline::DB::ObjI compliant object
                (which extends Bio::EnsEMBL::DB::ObjI)
    Args    :   A Bio::EnsEMBL::Pipeline::DB::ObjI compliant object

=head2 input_id

    Title   :   input_id
    Usage   :   $self->input_id($input_id);
    Function:   Gets or sets the value of input_id
    Returns :   valid input id for this analysis (if set) 
    Args    :   input id for this analysis 

=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output
    Function:   Writes output data to db
    Returns :   array of exons (with start and end)
    Args    :   none

=cut

sub write_output {
    my($self,@features) = @_;

    my $db = $self->dbobj;
   
    if( !defined $db ) {
      $self->throw("unable to make write db");
    }
    
    my %contighash;
    my $gene_obj = $db->gene_Obj;

    my @newgenes = $self->output;
    return unless ($#newgenes >= 0);

    # get new ids
    eval {

	my $genecount  = 0;
	my $transcount = 0;
	my $translcount = 0;
	my $exoncount  = 0;

	# get counts of each type of ID we need.

	foreach my $gene ( @newgenes ) {
	    $genecount++;
	    foreach my $trans ( $gene->each_Transcript ) {
		$transcount++;
		$translcount++;
	    }
	    foreach my $exon ( $gene->each_unique_Exon() ) {
		$exoncount++;
	    }
	}

	# get that number of ids. This locks the database

	my @geneids  =  $gene_obj->get_New_external_id('gene',$GENE_ID_SUBSCRIPT,$genecount);
	my @transids =  $gene_obj->get_New_external_id('transcript',$TRANSCRIPT_ID_SUBSCRIPT,$transcount);
	my @translids =  $gene_obj->get_New_external_id('translation',$PROTEIN_ID_SUBSCRIPT,$translcount);
	my @exonsid  =  $gene_obj->get_New_external_id('exon',$EXON_ID_SUBSCRIPT,$exoncount);

	# database locks are over.

	# now assign ids. gene and transcripts are easy. Exons are harder.
	# the code currently assummes that there is one Exon object per unique
	# exon id. This might not always be the case.

	foreach my $gene ( @newgenes ) {
	    $gene->id(shift(@geneids));
	    my %exonhash;
	    foreach my $exon ( $gene->each_unique_Exon() ) {
		my $tempid = $exon->id;
		$exon->id(shift(@exonsid));
		$exonhash{$tempid} = $exon->id;
	    }
	    foreach my $trans ( $gene->each_Transcript ) {
		$trans->id(shift(@transids));
		$trans->translation->id(shift(@translids));
		$trans->translation->start_exon_id($exonhash{$trans->translation->start_exon_id});
		$trans->translation->end_exon_id($exonhash{$trans->translation->end_exon_id});
	    }
	    
	}

	# paranoia!
	if( scalar(@geneids) != 0 || scalar(@exonsid) != 0 || scalar(@transids) != 0 || scalar (@translids) != 0 ) {
	    $self->throw("In id assignment, left with unassigned ids ".scalar(@geneids)." ".scalar(@transids)." ".scalar(@translids)." ".scalar(@exonsid));
	}

    };
    if( $@ ) {
	$self->throw("Exception in getting new ids. Exiting befor write\n\n$@" );
    }


    # this now assummes that we are building on a single VC.



  GENE: foreach my $gene (@newgenes) {	
      # do a per gene eval...
      eval {
	
	  $gene_obj->write($gene);
      }; 
      if( $@ ) {
	  print STDERR "UNABLE TO WRITE GENE\n\n$@\n\nSkipping this gene\n";
      }
	    
  }
   
}

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data from the database
    Returns :   nothing
    Args    :    string: chr1.1-10000

=cut

sub fetch_input {
    my( $self) = @_;

    my $genetype = 'bmeg';

    print STDERR "Fetching input: " . $self->input_id. " \n";
    $self->throw("No input id") unless defined($self->input_id);

    # get genomic 
    my $chrid  = $self->input_id;
       $chrid =~ s/\.(.*)-(.*)//;

    my $chrstart = $1;
    my $chrend   = $2;

    print STDERR "Chromosome id = $chrid , range $chrstart $chrend\n";

    $self->dbobj->static_golden_path_type('UCSC');

    my $stadaptor = $self->dbobj->get_StaticGoldenPathAdaptor();
    my $contig    = $stadaptor->fetch_VirtualContig_by_chr_start_end($chrid,$chrstart,$chrend);

    $contig->_chr_name($chrid);
    $self->vc($contig);

    print STDERR "got vc\n";

    # forward strand
    print STDERR "*** forward strand\n";

    # get genes
    my @genes  = $contig->get_Genes_by_Type($genetype, 'evidence');
    
    print STDERR "Number of genes = " . scalar(@genes) . "\n";
    if(!scalar(@genes)){
      $self->warn("No forward strand genes found");
    }

    my @plus_transcripts;

    my $single = 0;

    # split by strand
GENE:    
    foreach my $gene(@genes) {
      my @transcripts = $gene->each_Transcript;
      if(scalar(@transcripts) > 1 ){
	$self->warn($gene->id . " has more than one transcript - skipping it\n");
	next GENE;
      }

      my @exons = $transcripts[0]->each_Exon;

      if(scalar(@exons) == 1){
	$single++;
	next GENE;
      }

      if($exons[0]->strand == 1){
	push (@plus_transcripts, $transcripts[0]);
	next GENE;
      }

    }
    # cluster each strand
    if(scalar(@plus_transcripts)){
      my @plus_clusters  = $self->cluster_transcripts(@plus_transcripts);

      
      print STDERR scalar(@plus_transcripts) . " forward strand genes, and $single single exon genes\n";
      
      # make a genomewise runnable for each cluster of transcripts
      foreach my $cluster(@plus_clusters){
	print STDERR "new genomewise " . scalar(@$cluster) . "\n";
	
	my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::Genomewise();
	$runnable->seq($contig->primary_seq);
	#    my $genseq    = $contig->get_repeatmasked_seq; use repmasked seq instead of primary seq?
	$self->add_runnable($runnable);
	foreach my $t(@$cluster){
	  $runnable->add_Transcript($t);
	}
      }
    }

    # minus strand - flip the vc and hope it copes ...
    # but SLOW - get the same genes twice ...
    print STDERR "*** reverse strand\n";
    my $reverse = 1;
    my $revcontig = $contig->invert;
    my @revgenes  = $revcontig->get_Genes_by_Type($genetype, 'evidence');
    my @minus_transcripts;
    
    print STDERR "Number of genes = " . scalar(@revgenes) . "\n";
    if(!scalar(@genes)){
      $self->warn("No reverse strand genes found");
    }

    REVGENE:    
    foreach my $gene(@revgenes) {
      my @transcripts = $gene->each_Transcript;
      if(scalar(@transcripts) > 1 ){
	$self->warn($gene->id . " has more than one transcript - skipping it\n");
	next REVGENE;
      }

      my @exons = $transcripts[0]->each_Exon;

      if(scalar(@exons) == -1){
	next REVGENE;
      }

      # these are really - strand, but the VC is reversed, so they are realtively + strand
      # owwwww my brain hurts
      elsif($exons[0]->strand == 1){
	push (@minus_transcripts, $transcripts[0]);
	next REVGENE;
      }
    }
    
    if(scalar(@minus_transcripts)){

      my @minus_clusters = $self->cluster_transcripts(@minus_transcripts);  
      # minus clusters need some fancy schmancy feature and contig inversion. 
      foreach my $cluster(@minus_clusters){
	print STDERR "new genomewise " . scalar(@$cluster) . "\n";
	
	my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::Genomewise();
	$runnable->seq($revcontig->primary_seq);
	#    my $genseq    = $contig->get_repeatmasked_seq; use repmasked seq instead of primary seq?
	$self->add_runnable($runnable, $reverse);
	foreach my $t(@$cluster){
	  $runnable->add_Transcript($t);
	}
      }
    }

}
 
=head2 cluster_transcripts

    Title   :   cluster_transcripts
    Usage   :   $self->cluster_transcripts
    Function:   clusters input array of transcripts
    Returns :   
    Args    :   @Bio::EnsEMBL::Transcript

=cut

# borrows very heavily from Spangle
sub cluster_transcripts {
  my ($self, @alltranscripts) = @_;

  # need all the exons - we're going to cluster them; while we're at it, check consistency of hid & strand
  my %est2transcript; # relate transcript to est id
  my %exon2est;        # relate exon to est id
  my @exons = $self->check_transcripts(\%est2transcript, \%exon2est, @alltranscripts);
 
  # store references to the exon and transcript hashes
  $self->{'_est2transcript'} = \%est2transcript;
  $self->{'_exon2est'}       = \%exon2est;
 
  print STDERR scalar(@exons) . " exons found\n";
  
  my $main_cluster    = $self->make_clusters(@exons);
 
  $main_cluster = $self->find_common_ends($main_cluster);

  my @linked_clusters = $self->link_clusters($main_cluster);

  my @est_sets        = $self->process_clusters(@linked_clusters);

  my @transcript_clusters;

  foreach my $set(@est_sets){
    print STDERR "new transcript set with " .scalar(@$set) .  " ests\n";
    my %transcripts;

    # putting Transcripts as keys into hashes stringifies them  - key by transcript->id instead
    foreach my $est(@$set){
      my $trans = $est2transcript{$est};
#      print STDERR "transcript: " . $trans->id . "\n";
      $transcripts{$trans->id} = $trans;
    }
    push (@transcript_clusters, [ values %transcripts ])
  }

 return @transcript_clusters;

}

=head2 check_transcripts

    Title   :   check_transcripts
    Usage   :   $self->check_transcripts
    Function:   checks transcripts for strand consistency, hid consistency & exon content
    Returns :   @Bio::EnsEMBL::Exon
    Args    :   @Bio::EnsEMBL::Transcript, ref to hash for linking hid to transcript, ref to hash for linking exon to hid, 

=cut

sub check_transcripts {
  my ($self, $hid_trans, $exon_hid, @transcripts) = @_;
  my @allexons;

  TRANS: 
  foreach my $transcript(@transcripts){
    my @exons = $transcript->each_Exon;
    my $hid;
    my $strand;
    
    foreach my $exon(@exons){
      my $hstart;
      my $hend;

      # check strand consistency
      if(!defined $strand) { 
	$strand = $exon->strand; 
      }
      
      if($strand ne $exon->strand){
	$self->warn("strand not consistent among exons for " . $transcript->id . " - skipping it\n");
	next TRANS;
      }
      
      # check supporting_feature consistency
      my @sf = $exon->each_Supporting_Feature;      
    SF:
      foreach my $feature(@sf) {
	# because we get all the supporting features indiscriminately
	next SF unless $feature->source_tag eq 'bmeg'; 
	
	if(!defined $hid) { $hid = $feature->hseqname; }
	
	if($hid ne $feature->hseqname){
	  $self->warn("hid not consistent among exons for " . $transcript->id . " - skipping it\n");
	  next TRANS;
	}

      }

      $hstart = $sf[0]->hstart;
      $hend   = $sf[$#sf]->hend;
      my %exhash = {
		    'hid'    => $hid,
		    'hstart' => $hstart,
		    'hend'   => $hend
		   };


      $$exon_hid{$exon}{"hid"} = $hid;
      $$exon_hid{$exon}{"hstart"} = $hid;
      $$exon_hid{$exon}{"hend"} = $hid;
      # make sure we can get out the hid corresponding to this exon
    }

    push(@allexons, @exons);
    if(defined $hid && defined $$hid_trans{$hid}) { 
      $self->warn("$hid is being used by more than one transcript!\n"); 
    }
    $$hid_trans{$hid} = $transcript;
  }

  return @allexons;

}



sub make_clusters {
  my ($self, @exons) = @_;

  my $main_cluster = new Bio::EnsEMBL::SeqFeature; # main cluster feature - holds all subclusters
  @exons = sort { $a->start <=> $b->start } @exons;

  return unless @exons >= 0; # no point if there are no exons!
  
  # Create the first cluster object
  my $subcluster = new Bio::EnsEMBL::SeqFeature;
  
  # Start off the cluster with the first exon
  $subcluster->add_sub_SeqFeature($exons[0],'EXPAND');
  $subcluster->strand($exons[0]->strand);  

  $main_cluster->add_sub_SeqFeature($subcluster,'EXPAND');
  
  # Loop over the rest of the exons
  my $count = 0;

 EXON:
  foreach my $e(@exons) {
#    print STDERR "Clustering " . $e->id . "\n";

    if ($count > 0) {
      my @overlap = $self->match($e, $subcluster);    
 #     print STDERR "Overlap :@overlap:$#overlap\n";

      # Add to cluster if overlap AND if strand matches
      if ($overlap[0]) {
	$subcluster->add_sub_SeqFeature($e,'EXPAND');
      }  
      else {
	# Start a new cluster
	  $subcluster = new Bio::EnsEMBL::SeqFeature;
	  $subcluster->add_sub_SeqFeature($e,'EXPAND');
	  $subcluster->strand($e->strand);
	  
	  # And add to the top cluster feature
	  $main_cluster->add_sub_SeqFeature($subcluster,'EXPAND');	
      }
    }
    $count++;
  }
  
  # remove any clusters that have only 1 exon - we don't want to build genes from single est hits.
  # At least, I don't think we do.
  # This may cause fragmentation problems ... if only 1 est spans a gap between 2 clusters that 
  # really are related .
  my @exon_clusters = $main_cluster->sub_SeqFeature;
  $main_cluster->flush_sub_SeqFeature;

  my $count = 0;
  foreach my $ec(@exon_clusters) {
    $count++;
    my @exons    = $ec->sub_SeqFeature;
    $main_cluster->add_sub_SeqFeature($ec,'EXPAND') unless scalar(@exons) <= 1;
  }

  @exon_clusters = $main_cluster->sub_SeqFeature;

  return $main_cluster;
}

=head2 find_common_ends

    Title   :   find_common_ends
    Usage   :   $cluster = $self->find_common_ends($precluster)
    Function:   Resets the ends of the clusters to the most frequent coordinate
    Returns :   
    Args    :   

=cut

# should we introduce some sort of score weighting?
sub find_common_ends {
  my ($self, $main_cluster) = @_;
  my $exon_hid = $self->{'exon2est'};
  my @clusters = $main_cluster->sub_SeqFeature;
  @clusters    = sort { $a->start <=> $b->start  } @clusters;
  my $count    =  0 ;

  foreach my $cluster(@clusters) {
    my %start;
    my %end;
    my $newstart;
    my $newend;
    my $maxend = 0;
    my $maxstart = 0;
    my $newstart = $cluster->start;
    my $newend   = $cluster->end;

    # count frequency of starts & ends
    foreach my $exon ($cluster->sub_SeqFeature) {
      my $est = $$exon_hid{$exon}{'hid'};
      $start{$exon->start}++;
      $end{$exon->end}++;
      
    }

    print STDERR "cluster $count\nstarts\n";
    while ( my ($key, $value) = each %start) {
      if ($value > $maxstart){
	$newstart = $key;
	$maxstart = $value;
      } elsif ($value == $maxstart){
	print STDERR "$key is just as common as $newstart - ignoring it\n";
      }
      print STDERR "$key => $value\n";
    }

    print STDERR "ends\n";

    while ( my ($key, $value) = each %end) {
      if ($value > $maxend){
	$newend = $key;
	$maxend = $value;
      } elsif ($value == $maxend){
	print STDERR "$key is just as common as $newend - ignoring it\n";
      }
      print STDERR "$key => $value\n";
    }

    # if we haven't got a clear winner, we might as well stick with what we had
    if( $maxstart <= 1 ) {
      $newstart = $cluster->start;
    }
    if( $maxend <= 1 ) {
      $newend = $cluster->end;
    }

    # first and last clusters are special cases - potential UTRs, take the longest one.
    if( $count == 0) {
      print STDERR "first cluster\n";
      $newstart = $cluster->start;
    } elsif ( $count == $#clusters ) {
      print STDERR "last cluster\n";
      $newend = $cluster->end;
    }


    print STDERR "we used to have: " . $cluster->start . " - " . $cluster->end . "\n";
    print STDERR "and the winners are ... $newstart - $newend\n\n";

    $count++;

    # reset the ends of all the exons in this cluster to $newstart - $newend
    foreach my $exon($cluster->sub_SeqFeature) {
      $exon->start($newstart);
      $exon->end($newend);
    }
  }

  return $main_cluster;

}

=head2 is_in_cluster

  Title   : is_in_cluster
  Usage   : my $in = $self->is_in_cluster($id, $cluster)
  Function: Checks whether a feature with id $id is in a cluster
  Returns : 0,1
  Args    : String,Bio::SeqFeature::Generic

=cut

sub is_in_cluster {
  my ($self, $id, $cluster) = @_;

  my $exon2est = $self->{'_exon2est'};

  foreach my $exon ( $cluster->sub_SeqFeature ) {
    return (1) if $$exon2est{$exon}{'hid'} eq $id;
  }

  return 0;

}

=head2 link_clusters

    Title   :   link_clusters
    Usage   :   $self->link_clusters
    Function:   Links clusters together to find groups of transcripts predicting the same gene, for input to genomewise
    Returns :   
    Args    :   

=cut

sub link_clusters{
  my ($self, $main_cluster, $tol) = @_;
  my $exon_hid = $self->{'_exon2est'};
  $tol = 10 unless $tol; # max allowed discontinuity in sequence

  my @clusters = $main_cluster->sub_SeqFeature;
  print STDERR "Created " . scalar (@clusters) . " clusters\n";

  # now need to link up sub_clusters

  # Sort the clusters by start position
  @clusters = sort { $a->start <=> $b->start } @clusters;

  # Loop over all clusters
 CLUSTER1:  
  for (my $c = 0; $c < $#clusters; $c++) {
    my $p1    = $clusters[$c];
    my @sub   = $p1->sub_SeqFeature;
    
    print("Finding links in cluster number $c \n");
    
    # Sort the supporting features by score
    @sub = sort  { $b->score <=> $a->score } @sub;
    
    my $found = 0;
    my $c2    = $c+1;

    # Now look in all the other clusters to find a link
  CLUSTER2:    
    while ($c2 <= $#clusters) {
      
      my $p2   = $clusters[$c2];
      my @sub2 = $p2->sub_SeqFeature;
      
      # Only look at a cluster if it starts after the end of
      # the current cluster
      if ($p2->start > $p1->end) {
	# Sort by score
	@sub2 = sort { $a->score <=> $b->score } @sub2;
	
	my $foundlink = 0;

	# Now loop over the supporting features in the first linking cluster
	foreach my $s1 (@sub) {
	  # Does the second cluster contain this database id?
	  my $hid1 = $$exon_hid{$s1}{'hid'}; 

	  my $s2;
	  foreach my $exon(@sub2){
	    if ($$exon_hid{$exon}{'hid'} eq $hid1){
	      print STDERR "cluster $c " . $s1->id . " links to cluster $c2 " . $exon->id . "\n";
	      $s2 = $exon;
	    }
	  }

	  # If we have the same id in two clusters then...
	  if ( defined $s2 ) {
	    print("Found ids in two clusters: " . $s1->id . " " . $s2->id  . "\n");
#	    print($s2->toString . "\n");
	    # Check the strand is the same
	    if ($s1->strand eq $s2->strand) {
	      # Only allow a $tol discontinuity in the sequence

	      if ($s1->strand == $s2->strand && 
		  (abs($$exon_hid{$s1}{'hend'} - $$exon_hid{$s2}{'start'}) < $tol || 
		   abs($$exon_hid{$s1}{'hstart'} - $$exon_hid{$s2}{'hend'}) < $tol)) {

		$found = 1;
		$foundlink = 1;
	      }
	      else {
		print STDERR "not going to be able to make a link\n";
		print STDERR "s1 : " . $s1->start . " - " .$s1->end 
		  . " s2 : " . $s2->start . " - " . $s2->end . " tol " . $tol . "\n";
	      }
	    }

	    else {
	      print STDERR " strands differ, or no overlap\n";
	    }

	    if ($foundlink == 1) {
	      # We have found a link between p1 and p2
#	      print("Creating link from " . $p1->start . "\t" . $p1->end . "\n");
#	      print("to                 " . $p2->start . "\t" . $p2->end . "\n");
	      push(@{$p1->{_forward}} ,$p2);
	      push(@{$p2->{_backward}},$p1);
	      next CLUSTER1;
	    }
	  }

	} # end of loop over p2 sub seqfeatures

      } # end of p2->start > p1->end
      $c2++;
    }  # END of loop over p2

  } # End of loop over p1

  return @clusters;

}

sub process_clusters{
  my ($self, @clusters) = @_;

  my $exon_hid = $self->{'_exon2est'};

  # pullout sets of linked transcripts
  my @transcript_sets;
  
  for ( my $c = 0; $c < $#clusters; $c++ ){
    my $exon_cluster = $clusters[$c];
    my %t_hids;
    
    # if this the first of a new transcript cluster
    if (!defined($exon_cluster->{_backward}) && defined($exon_cluster->{_forward})) {
      # start collecting hids - we can maop back to transcripts using hids.

      my @exons = $exon_cluster->sub_SeqFeature;
      
      # collect the hids for all the exons in this cluster
      foreach my $ex(@exons){
	my $e_hid = $$exon_hid{$ex}{'hid'};
	$t_hids{$e_hid} = 1;
      }
      
      # look at all the exon_clusters that link on from this one
      while (ref($exon_cluster->{_forward}) ne "") {
	$exon_cluster = $exon_cluster->{_forward}[0];
	@exons = $exon_cluster->sub_SeqFeature;
	# collect the hids for all the exons in this cluster
	foreach my $ex(@exons){
	  my $e_hid = $$exon_hid{$ex}{'hid'};
		  $t_hids{$e_hid} = 1;
	}
      } 
    }

    # we have reached the end of the first group of exon_clusters.
    if(scalar(keys %t_hids)){
      push (@transcript_sets, [ keys %t_hids  ]);
    }
  }
  
  return @transcript_sets;
}

sub add_runnable{
  my ($self, $value, $reverse) = @_;

  if (!defined($self->{'_forward_runnables'})) {
    $self->{'_forward_runnables'} = [];
  }
  if (!defined($self->{'_reverse_runnables'})) {
    $self->{'_reverse_runnables'} = [];
  } 
  if (defined($value)) {
    if ($value->isa("Bio::EnsEMBL::Pipeline::RunnableI")) {
      if(defined $reverse){
	push(@{$self->{'_reverse_runnables'}},$value);
      }
      else {
	push(@{$self->{'_forward_runnables'}},$value);
      }
    } else {
      $self->throw("[$value] is not a Bio::EnsEMBL::Pipeline::RunnableI");
    }
  }
  
} 

sub each_runnable{
  my ($self, $reverse) = @_;
  
  if (!defined($self->{'_forward_runnables'})) {
    $self->{'_forward_runnables'} = [];
  }
  
  if (!defined($self->{'_reverse_runnables'})) {
    $self->{'_reverse_runnables'} = [];
  } 

  if(defined $reverse){
    return @{$self->{'_reverse_runnables'}};
  }
  else {
    return @{$self->{'_forward_runnables'}};
  }
  
}

sub run {
  my ($self) = @_;

  # run genomewise 

  # sort out analysis & genetype here or we will get into trouble with duplicate analyses
#  my $genetype = 'genomewise'; 
  my $genetype = 'new_genomewise'; 
  my $anaAdaptor = $self->dbobj->get_AnalysisAdaptor;
  my @analyses = $anaAdaptor->fetch_by_logic_name($genetype);
  my $analysis_obj;
  if(scalar(@analyses) > 1){
    $self->throw("panic! > 1 analysis for $genetype\n");
  }
  elsif(scalar(@analyses) == 1){
      $analysis_obj = $analyses[0];
  }
  else{
      # make a new analysis object
      $analysis_obj = new Bio::EnsEMBL::Analysis
	(-db              => 'NULL',
	 -db_version      => 1,
	 -program         => $genetype,
	 -program_version => 1,
	 -gff_source      => $genetype,
	 -gff_feature     => 'gene',
	 -logic_name      => $genetype,
	 -module          => 'EST_GeneBuilder',
      );
  }

  $self->genetype($genetype);
  $self->analysis($analysis_obj);

  # plus strand
  foreach my $gw_runnable( $self->each_runnable) {
    print STDERR "about to run genomewise\n";
    $gw_runnable->run;

    # convert_output
    $self->convert_output($gw_runnable);
  }

  # minus strand
  my $reverse = 1;
  foreach my $gw_runnable( $self->each_runnable($reverse)) {
    print STDERR "about to run genomewise\n";
    $gw_runnable->run;

    # convert_output
    $self->convert_output($gw_runnable, $reverse);
  }
}

sub genetype {
  my ($self, $genetype) = @_;

  if(defined $genetype){
    $self->{'_genetype'} = $genetype;
  }

  return $self->{'_genetype'};
}

# override method from RunnableDB.pm
sub analysis {
  my ($self, $analysis) = @_;

  if(defined $analysis){
    $self->throw("$analysis is not a Bio::EnsEMBL::Analysis") unless $analysis->isa("Bio::EnsEMBL::Analysis");
    $self->{'_analysis'} = $analysis;
  }

  return $self->{'_analysis'};
}

# convert genomewise output into genes
sub convert_output {
  my ($self, $gwr, $reverse) = @_;

  my @genes = $self->make_genes($gwr, $reverse);

  my @remapped = $self->remap_genes(\@genes, $reverse);

  # store genes
 
  $self->output(@remapped);
 
}



sub make_genes {

  my ($self, $runnable, $reverse) = @_;
  my $genetype = $self->genetype;
  my $analysis_obj = $self->analysis;
  my $count = 0;
  my @genes;

  my $time  = time; chomp($time);
  my $contig = $self->vc;

  # are we working on the reverse strand?
  if(defined $reverse){
    $contig = $contig->invert;
  }
  my @trans = $runnable->output;
  print "transcripts: " . scalar(@trans) . "\n";

  foreach my $transcript ($runnable->output) {
    $count++;
    my $gene   = new Bio::EnsEMBL::Gene;
    $gene->type($genetype);
    $gene->id($contig->id . ".$genetype.$count");
    $gene->version(1);
    $gene->created($time);
    $gene->modified($time);
    $gene->analysis($analysis_obj);

    # add transcript to gene
    $transcript->id($contig->id . ".$genetype.$count");
    $transcript->version(1);
    $gene->add_Transcript($transcript);


    # and store it
    push(@genes,$gene);

    # sort the exons 
    $transcript->sort;
    my $excount = 1;
    my @exons = $transcript->each_Exon;

    foreach my $exon(@exons){
      $exon->id($contig->id . ".$genetype.$count.$excount");
      $exon->contig_id($contig->id);
      $exon->created($time);
      $exon->modified($time);
      $exon->version(1);
      $exon->attach_seq($contig->primary_seq);
      $excount++;
      }
    
    print STDERR "exons: " . scalar(@exons) . "\n";

    # sort out translation
    my $translation  = new Bio::EnsEMBL::Translation;    
    $translation->id($contig->id . ".$genetype.$count");
    $translation->version(1);    
    
    $translation->start_exon_id($exons[0]->id);

    $translation->end_exon_id  ($exons[$#exons]->id);
    if ($exons[0]->phase == 0) {
      $translation->start(1);
    } elsif ($exons[0]->phase == 1) {
      $translation->start(3);
    } elsif ($exons[0]->phase == 2) {
      $translation->start(2);
    }
    $translation->end  ($exons[$#exons]->end - $exons[$#exons]->start + 1);

    $transcript->translation($translation);
  }
  return @genes;
}

=head2 _remap_genes

    Title   :   _remap_genes
    Usage   :   $self->_remap_genes(@genes)
    Function:   Remaps predicted genes into genomic coordinates
    Returns :   array of Bio::EnsEMBL::Gene
    Args    :   Bio::EnsEMBL::Virtual::Contig, array of Bio::EnsEMBL::Gene

=cut

sub remap_genes {
  my ($self, $genes, $reverse) = @_;
  my $contig = $self->vc;
  if (defined $reverse){
    $contig = $contig->invert;
  }

  print STDERR "genes before remap: " . scalar(@$genes) . "\n";
  
  my @newf;
  my $trancount=1;
  foreach my $gene (@$genes) {
    eval {
      my $newgene = $contig->convert_Gene_to_raw_contig($gene);
      # need to explicitly add back genetype and analysis.
      $newgene->type($gene->type);
      $newgene->analysis($gene->analysis);

      push(@newf,$newgene);
      
    };
    if ($@) {
      print STDERR "Couldn't reverse map gene " . $gene->id . " [$@]\n";
    }
    
  }
  
  return @newf;
  
}


=head2 vc

 Title   : vc
 Usage   : $obj->vc($newval)
 Function: 
 Returns : value of vc
 Args    : newvalue (optional)


=cut

sub vc{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'vc'} = $value;
    }
    return $obj->{'vc'};

}

# Adapted from Spangle

=head2 match

 Title   : match
 Usage   : @overlap = $feat->match($f)
 Function: Returns an array of 3 numbers detailing
           the overlap between the two features.
           $overlap[0]  = 1 if any overlap, 0 otherwise
           The remaining elements of the array are only set if $overlap[0] = 1
           $overlap[1]  = left hand overlap (-ve if $f starts within $self, +ve if outside)
           $overlap[2]  = right hand overlap (-ve if $f ends within $self, $+ve if outside)
 Returns : @int
 Args    : Bio::SeqFeature::Generic

=cut

sub match {
  my ($self, $f1,$f2) = @_;
  
  my ($start1,
      $start2,
      $end1,
      $end2,
      $rev1,
      $rev2,
     );

  # Swap the coords round if necessary
  if ($f1->start > $f1->end) {
    $start1 = $f1->end;
    $end1   = $f1->start;
    $rev1   = 1;
  } else {
    $start1 = $f1->start;
    $end1   = $f1->end;
  }

  if ($f2->start > $f2->end) {
    $start2 = $f2->end;
    $end2   = $f2->start;
    $rev2   = 1;
  } else {
    $start2 = $f2->start;
    $end2   = $f2->end;
  }

  # Now check for an overlap
  if (($end2 > $start1 && $start2 < $end1) ||
	($start2 < $end1 && $end2 > $start1)) {
	
	#  we have an overlap so we now need to return 
	#  two numbers reflecting how accurate the span 
	#  is. 
	#  0,0 means an exact match with the exon
	# a positive number means an over match to the exon
	# a negative number means not all the exon bases were matched

	my $left = ($start2 - $start1);
	my $right = ($end1 - $end2);
	
	if ($rev1) {
	    my $tmp = $left;
	    $left = $right;
	    $right = $tmp;
	}
	
	my @overlap;

	push (@overlap,1);

	push (@overlap,$left);
	push (@overlap,$right);

	return @overlap;
      }
    
}

=head2 output

 Title   : output
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub output{
   my ($self,@genes) = @_;

   if (!defined($self->{'_output'})) {
     $self->{'_output'} = [];
   }
    
   if(defined @genes){
     push(@{$self->{'_output'}},@genes);
   }

   return @{$self->{'_output'}};
}


1;


