#
# BioPerl module for GeneBuilder
#
# Cared for by EnsEMBL <ensembl-dev@ebi.ac.uk>
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::GeneBuilder

=head1 SYNOPSIS

# This is the main analysis database

my $db = new Bio::EnsEMBL::DBSQL::Obj(-host   => 'obi-wan',
				      -user   => 'ensro',
				      -dbname => 'ens500',
				      );

# Fetch a clone and its contigs from the database
my $clone       = $db   ->get_Clone($clone);
my @contigs     = $clone->get_all_Contigs;

# The genebuilder object will fetch all the features from the contigs
# and use them to first construct exons, then join those exons into
# exon pairs.  These exon apris are then made into transcripts and
# finally all overlapping transcripts are put together into one gene.


my $genebuilder = new Bio::EnsEMBL::Pipeline::GeneBuilder
    (-contigs => \@contigs);

my @genes       = $genebuilder->build_Genes;

# After the genes are built they can be used to order the contigs they
# are on.

my @contigs     = $genebuilder->order_Contigs;


=head1 DESCRIPTION

Takes in contigs and returns genes.  The procedure is currently
reimplementing the TimDB method of building genes where genscan exons
are confirmed by similarity features which are then joined together
into exon pairs.  An exon pair is constructed as follows :

  ---------          --------    genscan exons
    -------          ----->      blast hit which spans an intron
    1     10        11    22        

For an exon pair to make it into a gene there must be at least 2 blast
hits (features) that span across an intron.  This is called the
coverage of the exon pair.

After all exon pairs have been generated for all the genscan exons
there is a recursive routine (_recurseTranscripts) that looks for all
exons that are the start of an exon pair with no preceding exons.  The
exon pairs are followed recursively (including alternative splices) to
build up full set of transcripts.

To generate the genes the transcripts are grouped together into sets
with overlapping exons.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::GeneBuilder;

use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::Runnable::PredictionGeneBuilder;
use Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster;
use Bio::EnsEMBL::Pipeline::GeneConf qw (
					 TRANSCRIPT_ID_SUBSCRIPT
					 GB_MIN_GENSCAN_EXONS
					 GB_GENSCAN_MAX_INTRON
					 GB_TARGETTED_GW_GENETYPE
					 GB_SIMILARITY_GENETYPE
					 GB_COMBINED_GENETYPE
					 GB_MIN_FEATURE_SCORE
					 GB_MIN_FEATURE_LENGTH
					 GB_INPUTID_REGEX
					 GB_ABINITIO_TYPE
					 GB_ABINITIO_SUPPORTED_TYPE
					 GB_MAXSHORTINTRONLEN
					 GB_MINSHORTINTRONLEN
					);

use vars qw(@ISA);
use strict;

@ISA = qw(Bio::EnsEMBL::Root);


############################################################

sub new {
    my ($class,@args) = @_;

    my $self = $class->SUPER::new(@args);

    my ($slice,$input_id) = $self->_rearrange([qw(SLICE INPUT_ID)],
					      @args);

    $self->throw("Must input a slice to GeneBuilder") unless defined($slice);
    $self->slice($slice);
    $self->{_final_genes}          = [];
    $self->{_gene_types}           = [];
    $self->gene_types($GB_COMBINED_GENETYPE);
    $self->gene_types($GB_TARGETTED_GW_GENETYPE);
    $self->gene_types($GB_SIMILARITY_GENETYPE);

    unless ( $input_id =~ /$GB_INPUTID_REGEX/ ){
      $self->throw("format of the input is not defined in GeneConf::GB_INPUTID_REGEX = $GB_INPUTID_REGEX");
    }
    $self->input_id($input_id);

    return $self;
}

############################################################

=head2 input_id

 Function: get/set for input id
 Returns : string
 Args    : string (it expects a string of the format chr_name.start_coord-end_coord

=cut
  
sub input_id {
  my ($self,$id) = @_;
  
  if (defined($id)) {
    $self->{_input_id} = $id;
  }
  return $self->{_input_id};
}

############################################################

=head2 build_Genes

 Example    : my @genes = $self->build_Genes
 Description: builds genes. It is like the run method in Runnables. It calls everything that needs to be done.
 Returns    : none
 Args       : none
 Caller     : Bio::EnsEMBL::Pipeline::RunnableDB::Gene_Builder

=cut

sub build_Genes{
  my ($self) = @_;
  
  print STDERR "Building genes...\n";
  
  # get all genes of type defined in gene_types() on this slice
  $self->get_Genes;
  print STDERR "Number of genewise and combined transcripts " . scalar($self->genewise_combined_Transcripts) . "\n";
  
  # get all Genscan predictions on this slice
  $self->get_Predictions;
  print STDERR "Number of ab initio predictions ". scalar($self->predictions)  . "\n";
  
  # get all the dna/protein align features from the pre-computes pipeline on this slice
  $self->get_Similarities;
  print STDERR "Number of similarity features ". scalar($self->features) . "\n";
  
  # process PredictionTranscripts using the features and the annotations:
  my @predictions = $self->predictions;
  my @features    = $self->features;
  my @annotations = $self->genewise_combined_Transcripts;
  
  my $genecooker = Bio::EnsEMBL::Pipeline::Runnable::PredictionGeneBuilder->new(
									       -predictions => \@predictions,
									       -features    => \@features,
									       -annotations => \@annotations,
									      );
  
  my @supported_predictions = $genecooker->run;
  
  # cluster all the transcripts according to mere genomic overlap
  my @all_transcripts     = ( @annotations, @supported_predictions );
  unless( @all_transcripts ){
      print STDERR "no transcripts left to cook. Exiting...\n";
      exit(0);
  }

  print STDERR "clustering transcripts...\n";
  my @transcript_clusters = $self->cluster_Transcripts(@all_transcripts);
  print STDERR scalar(@transcript_clusters)." clusters formed\n";

  # prune the redundant transcripts for each cluster
  print STDERR "pruning transcripts...\n";
  my @pruned_transcripts = $self->prune_Transcripts(@transcript_clusters);
  print STDERR scalar(@pruned_transcripts)." transcripts obtained\n";

  # cluster transcripts into genes
  print STDERR "clustering into genes...\n";
  my @genes = $self->cluster_into_Genes(@pruned_transcripts);
  print STDERR scalar(@genes)." genes built\n";
  
  print STDERR "Final_results:\n";
  my $count = 0;
  foreach my $gene ( @genes ){
      $count++;
      print STDERR "Gene $count:\n";
      foreach my $tran ( @{$gene->get_all_Transcripts} ){
	  $self->_print_Transcript($tran);
      }
  }

  # final_genes is not working, check it!!
  print STDERR scalar($self->final_genes)." final genes\n";
    
}



############################################################


=head2 get_Genes

 Description: retrieves genewise and combined gene annotations with supporting evidence. 
              Splits transcripts with very long introns, discards transcripts with strand problems etc.
 ReturnType : none, but $self->genewise_combined_Transcripts is filled
 Args       : none

=cut

sub get_Genes {
  my ($self) = @_;
  my @transcripts;
  my $db = $self->genes_db;
  my $sa = $db->get_SliceAdaptor;
  
  my $input_id = $self->input_id;
  $input_id =~/$GB_INPUTID_REGEX/;
  my $chr   = $1;
  my $start = $2;
  my $end   = $3;
  my $slice = $sa->fetch_by_chr_start_end($chr,$start,$end);
  
  my @unchecked_genes;
  
  foreach my $type ($self->gene_types) {
    foreach my $gene ( @{$slice->get_all_Genes_by_type($type, 'evidence')} ){
    
    TRANSCRIPT:
      foreach my $tran (@{$gene->get_all_Transcripts}) {
	
	# set temporary_id to be dbID
	$tran->{'temporary_id'} = ($tran->dbID) unless (defined $tran->{'temporary_id'} && $tran->{'temporary_id'} ne '');
	
	# my @valid_transcripts = $self->validate_transcript($t);
	# next TRANSCRIPT unless scalar(@valid_transcripts);
	
	unless ( $self->_check_Transcript( $tran ) && $self->_check_Translation($tran)){
	  next TRANSCRIPT;
	}
	$tran->type($type);
	push(@transcripts, $tran);
      }
    }
  }
  $self->genewise_combined_Transcripts(@transcripts);
}


###########################################################c

=head2 cluster_Transcripts

 Description : It separates transcripts according to strand and then clusters 
               each set of transcripts by calling _cluster_Transcripts_by_genomic_range()
  Args       : Array of Bio::EnsEMBL::Transcript
  Return     : Array of Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster

=cut

sub cluster_Transcripts {
  my ($self,@transcripts) = @_;
 
  my @forward_transcripts;
  my @reverse_transcripts;
 
  foreach my $transcript (@transcripts){
    my @exons = @{ $transcript->get_all_Exons };
    if ( $exons[0]->strand == 1 ){
      push( @forward_transcripts, $transcript );
    }
    else{
      push( @reverse_transcripts, $transcript );
    }
  }
  
  my @forward_clusters = $self->_cluster_Transcripts_by_genomic_range( @forward_transcripts );
  my @reverse_clusters = $self->_cluster_Transcripts_by_genomic_range( @reverse_transcripts );
  
  return ( @forward_clusters, @reverse_clusters );
}

############################################################

=head2 _cluster_Transcripts_by_genomic_range

 Description : It clusters transcripts according to genomic overlap
  Args       : Array of Bio::EnsEMBL::Transcript
  Return     : Array of Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster

=cut

sub _cluster_Transcripts_by_genomic_range{
  my ($self,@mytranscripts) = @_;

  # first sort the transcripts by their start position coordinate
  my %start_table;
  my $i=0;
  foreach my $transcript (@mytranscripts){
    my $start;
    my $seqname;
    my @exons = @{$transcript->get_all_Exons};
    @exons = sort { $a->start <=> $b->start } @exons;
    if ( $exons[0]->start > $exons[0]->end){
      $start = $exons[0]->end;
    }
    else{
      $start = $exons[0]->start;
    }
    $start = $transcript->start_Exon->start;
    $start_table{$i} = $start;
    $i++;
  }
  
  my @transcripts;
  foreach my $pos ( sort { $start_table{$a} <=> $start_table{$b} } keys %start_table ){
    push (@transcripts, $mytranscripts[$pos]);
  }

  # create a new cluster 
  my $cluster=Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster->new();
  my $count = 0;
  my @cluster_starts;
  my @cluster_ends;
  my @clusters;
  
  # put the first transcript into these cluster
  $cluster->put_Transcripts( $transcripts[0] );
  $cluster_starts[$count] = $transcripts[0]->start;
  $cluster_ends[$count]   = $transcripts[0]->end;
  
  # store the list of clusters
  push( @clusters, $cluster );
  
  # loop over the rest of the transcripts
 LOOP1:
  for (my $c=1; $c<=$#transcripts; $c++){
    my $found=0;
    
    if ( !( $transcripts[$c]->end < $cluster_starts[$count] ||
	    $transcripts[$c]->start > $cluster_ends[$count] ) ){
      $cluster->put_Transcripts( $transcripts[$c] );
      
      # re-adjust size of cluster
      if ($transcripts[$c]->start < $cluster_starts[$count]) {
	$cluster_starts[$count] = $transcripts[$c]->start;
      }
      if ( $transcripts[$c]->end  > $cluster_ends[$count]) {
	$cluster_ends[$count] =  $transcripts[$c]->end;
      }
    }
    else{
      # else, create a new cluster with this feature
      $count++;
      my $cluster=Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster->new();
      $cluster->put_Transcripts( $transcripts[$c] );
      $cluster_starts[$count] = $transcripts[$c]->start;
      $cluster_ends[$count]   = $transcripts[$c]->end;
      
      # store it in the list of clusters
      push(@clusters,$cluster);
    }
  }
  return @clusters;
}

############################################################


=head2 prune_Transcripts

 Example    : my @pruned_transcripts = $self->prune_Transcripts(@transcript_clusters)
 Description: rejects duplicate transcripts, transfers supporting feature data from the rejected transcripts
               to the accepted ones
 Returns    : array of Bio::EnsEMBL::Transcript
  Args      : array of Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster

=cut

sub prune_Transcripts {
  my ($self, @transcript_clusters) = @_;
  my @newtran;
  

 CLUSTER:
  foreach my $transcript_cluster ( @transcript_clusters ){
    my @transcripts = $transcript_cluster->get_Transcripts;
    
    my @clusters;
    my %lengths;
    
    my @newgenes;
    my $max_num_exons = 0;
    
    # sizehash holds transcript length - based on sum of exon lengths
    my %sizehash;
    
    # orfhash holds orf length - based on sum of translateable exon lengths
    my %orfhash;
    
    my %tran2orf;
    my %tran2length;
    
    # keeps track of to which transcript(s) this exon belong
    my %exon2transcript;
    
    foreach my $tran (@transcripts) {
      
      # keep track of number of exons in multiexon transcripts
      my @exons = @{ $tran->get_all_Exons };
      if(scalar(@exons) > $max_num_exons){ $max_num_exons = scalar(@exons); }
      
      # total exon length
      my $length = 0;
      foreach my $e ( @{ $tran->get_all_Exons} ){
	$length += $e->end - $e->start + 1;
	
	push ( @{ $exon2transcript{ $e } }, $tran );
	
      }
      $sizehash{$tran->{'temporary_id'}} = $length;
      $tran2length{ $tran } = $length;
      
      # now for ORF length
      $length = 0;
      foreach my $e(@{$tran->get_all_translateable_Exons}){
	$length += $e->end - $e->start + 1;
      }
      $tran2orf{ $tran } = $length;
      push(@{$orfhash{$length}}, $tran);
   }  
    

    # VAC 15/02/2002 sort transcripts based on total exon length - this
    # introduces a problem - we can (and have) masked good transcripts
    # with long translations in favour of transcripts with shorter
    # translations and long UTRs that are overall slightly longer. This
    # is not good.
    # better way? hold both total exon length and length of translateable exons. Then sort:
    # long translation + UTR > long translation no UTR > short translation + UTR > short translation no UTR
    
    
    # eae: Notice that this sorting:
    my @sordid_transcripts =  sort { my $result = ( 
						   $tran2orf{ $b }
						   <=> 
						   $tran2orf{ $a }
						  );
				     if ($result){
				       return $result;
				     }
				     else{
				       return ( $tran2length{ $b } <=> $tran2length{ $a } )
				     }
				   } @transcripts;
    
    # is only equivalent to the one used below when all transcripts have UTRs. 
    # The one below is actually the desired behaviour.
    # since we want long UTRs and long ORFs but the sorting must be fuzzy in the sense that we want to give priority 
    # to a long ORF with UTR over a long ORF without UTR which could only slightly longer.
    
    #test
    #print STDERR "1.- sordid transcripts:\n";
    #foreach my $tran (@sordid_transcripts){
    #  print STDERR $tran." orf_length: $tran2orf{$tran}, total_length: $tran2length{$tran}\n";
    #}
    
    @transcripts = ();
    
    # sort first by orfhash{'length'}
    my @orflengths = sort {$b <=> $a} (keys %orfhash);
    
    # strict sort by translation length is just as wrong as strict sort by UTR length
    # bin translation lengths - 4 bins (based on 25% length diff)? 10 bins (based on 10%)?
    my %orflength_bin;
    my $numbins = 4;
    my $currbin = 1;
    
    foreach my $orflength(@orflengths){
      last if $currbin > $numbins;
      my $percid = ($orflength*100)/$orflengths[0];
      if ($percid > 100) { $percid = 100; }
      my $currthreshold = $currbin * (100/$numbins);
      $currthreshold = 100 - $currthreshold;
      
      if($percid <$currthreshold) { $currbin++; }
      my @tmp = @{$orfhash{$orflength}};
      push(@{$orflength_bin{$currbin}}, @{$orfhash{$orflength}});
    }
    
    # now, foreach bin in %orflengthbin, sort by exonlength
    $currbin = 1;
  EXONLENGTH_SORT:
    while( $currbin <= $numbins){
      if(!defined $orflength_bin{$currbin} ){
	$currbin++;
	next EXONLENGTH_SORT;
      }
      
      my @sorted_transcripts = sort {$sizehash{$b->{'temporary_id'}} <=> $sizehash{$a->{'temporary_id'}}} @{$orflength_bin{$currbin}};
      push(@transcripts, @sorted_transcripts);
      $currbin++;
  }
    
    #test
    my $debug = 0;
    if ($debug == 1 ){
      
      print STDERR "2.- sorted transcripts:\n";
      
      
      foreach my $tran (@transcripts){
	if ( $tran->dbID ){
	  print STDERR $tran->dbID." ";
	}
	print STDERR "orf_length: $tran2orf{$tran}, total_length: $tran2length{$tran}\n";
	my @exons = @{$tran->get_all_Exons};
	if ( $exons[0]->strand == 1 ){
	  @exons = sort { $a->start <=> $b->start } @exons;
	}
	else{
	  @exons = sort { $b->start <=> $a->start } @exons;
	}
	foreach my $exon ( @{$tran->get_all_Exons} ){
	  print "  ".$exon->start."-".$exon->end." ".( $exon->end - $exon->start + 1)." phase: ".$exon->phase." end_phase ".$exon->end_phase." strand: ".$exon->strand."\n";
	}
      }
   }
    # old way - sort strictly on exon length
    #    @transcripts = sort {$sizehash{$b->{'temporary_id'}} <=> $sizehash{$a->{'temporary_id'}}} @transcripts;
    
    ########################################
    # deal with single exon genes
    my @maxexon = @{$transcripts[0]->get_all_Exons};
    
    # do we really just want to take the first transcript only? What about supporting evidence from other transcripts?
    # also, if there's a very long single exon gene we will lose any underlying multi-exon transcripts
    # this may increase problems with the loss of valid single exon genes as mentioned below. 
    # it's a balance between keeping multi exon transcripts and losing single exon ones
    if ($#maxexon == 0 && $max_num_exons == 1) {
      push(@newtran, $transcripts[0] );
	
      # we are done with this cluster
      next CLUSTER;
    }
    
    # otherwise we need to deal with multi exon transcripts and reject duplicates.
    
    # links each exon in the transcripts of this cluster with a hash of other exons it is paired with
    my %pairhash;
    
    # allows retrieval of exon objects by exon->id - convenience
    my %exonhash;
    
    foreach my $tran (@transcripts) {
      my @exons = @{$tran->get_all_Exons};
      $tran->sort;
      
      #print STDERR "\ntranscript: ".$tran->{'temporary_id'}."\n";
      #foreach my $exon ( @exons ){
      #  print STDERR $exon->start."-".$exon->end." ";
      #}
      #print STDERR "\n";
      
      
      my $i     = 0;
      my $found = 1;
      
      # if this transcript has already been seen, this
      # will be used to transfer supporting evidence
      my @evidence_pairs;
      
      # 10.1.2002 VAC we know there's a potential problem here - single exon transcripts which are in a 
      # cluster where the longest transcriopt has > 1 exon are not going to be considered in 
      # this loop, so they'll always be marked "transcript already seen"
      # How to sort them out? If the single exon overlaps an exon in a multi exon transcript then 
      # by our rules it probably ought to be rejected the same way transcripts with shared exon-pairs are.
      # Tough one.
      # more of a problem is that if the transcript with the largest number of exons is really a
      # single exon with frameshifts, it will get rejected here based on intron size but in addition
      # any valid non-frameshifted single exon transcripts will get rejected - which is definitely not right.
      # We need code to represent frameshifted exons more sensibly so the frameshifted one doesn't 
      # get through the check for single exon genes above.
      
    EXONS:
      for ($i = 0; $i < $#exons; $i++) {
	my $foundpair = 0;
	my $exon1 = $exons[$i];
	my $exon2 = $exons[$i+1];
	
	# Only count introns > 50 bp as real introns
	my $intron;
	if ($exon1->strand == 1) {
	  $intron = abs($exon2->start - $exon1->end - 1);
	} else {
	  $intron = abs($exon1->start - $exon2->end - 1);
	}
	
	#	GB_MINSHORTINTRONLEN
	# GB_MAXSHORTINTRONLEN 
	# print STDERR "Intron size $intron\n";
	if ($intron < $GB_MAXSHORTINTRONLEN && $intron > $GB_MINSHORTINTRONLEN ) {
	  $foundpair = 1;	# this pair will not be compared with other transcripts
	} 
	else {
	  
	  # go through the exon pairs already stored in %pairhash. 
	  # If there is a pair whose exon1 overlaps this exon1, and 
	  # whose exon2 overlaps this exon2, then these two transcripts are paired
	  
	  foreach my $exon1id (keys %pairhash) {
	    my $exon1a = $exonhash{$exon1id};
	    
	    foreach my $exon2id (keys %{$pairhash{$exon1id}}) {
	      my $exon2a = $exonhash{$exon2id};
	      
	      if (($exon1->overlaps($exon1a) && $exon2->overlaps($exon2a))) {
		$foundpair = 1;
		
		# eae: this method allows a transcript to be covered by exon pairs
		# from different transcripts, rejecting possible
		# splicing variants
		
		# we put first the exon from the transcript being tested:
		push( @evidence_pairs, [ $exon1 , $exon1a ] );
		push( @evidence_pairs, [ $exon2 , $exon2a ] );
		
		# transfer evidence between exons, assuming the suppfeat coordinates are OK.
		# currently not working as the supporting evidence is not there - 
		# can get it for genewsies, but why not there for genscans?
		#	      $self->transfer_supporting_evidence($exon1, $exon1a);
		#	      $self->transfer_supporting_evidence($exon1a, $exon1);
		#	      $self->transfer_supporting_evidence($exon2, $exon2a);
		#	      $self->transfer_supporting_evidence($exon2a, $exon2);
	      }
	    }
	  }
	}
	
	if ($foundpair == 0) {	# ie this exon pair does not overlap with a pair yet found in another transcript
	  
	  $found = 0;		# ie currently this transcript is not paired with another
	  
	  # store the exons so they can be retrieved by id
	  $exonhash{$exon1->{'temporary_id'}} = $exon1;
	  $exonhash{$exon2->{'temporary_id'}} = $exon2;
	  
	  # store the pairing between these 2 exons
	  $pairhash{$exon1->{'temporary_id'}}{$exon2->{'temporary_id'}} = 1;
	}
      }				# end of EXONS
      
      # decide whether this is a new transcript or whether it has already been seen
      if ($found == 0) {
	  #print STDERR "found new transcript " . $tran->{'temporary_id'} . "\n";
	  push(@newtran,$tran);
	  @evidence_pairs = ();
      } 
      else {
	#print STDERR "\n\nTranscript already seen " . $tran->{'temporary_id'} . "\n";
	
	## transfer supporting feature data. We transfer it to exons
	foreach my $pair ( @evidence_pairs ){
	  my @pair = @$pair;
	  
	  # first in the pair is the 'already seen' exon
	  my $source_exon = $pair[0];
	  my $target_exon = $pair[1];
	  
	  #print STDERR "\n";
	  $self->transfer_supporting_evidence($source_exon, $target_exon)
	}
      }
    } # end of this transcript
  }
  
  return @newtran;
}
 

############################################################

=head2 cluster_into_Genes

    Example :   my @genes = $self->cluster_into_Genes(@transcripts);
Description :   it clusters transcripts into genes according to exon overlap.
                It will take care of difficult cases like transcripts within introns.
                It also unify exons that are shared among transcripts.
    Returns :   a beautiful list of geen objects
    Args    :   a list of transcript objects

=cut

sub cluster_into_Genes{
  my ($self, @transcripts_unsorted) = @_;
  
  my $num_trans = scalar(@transcripts_unsorted);
  print STDERR "clustering $num_trans transcripts into genes\n";

  
  # flush old genes
  #$self->flush_Genes;
  
  my @transcripts = sort by_transcript_high @transcripts_unsorted;
  my @clusters;
  
  # clusters transcripts by whether or not any exon overlaps with an exon in 
  # another transcript (came from original prune in GeneBuilder)
  foreach my $tran (@transcripts) {
    my @matching_clusters;
  CLUSTER: 
    foreach my $cluster (@clusters) {
      foreach my $cluster_transcript (@$cluster) {
        
        foreach my $exon1 (@{$tran->get_all_Exons}) {
	  foreach my $cluster_exon (@{$cluster_transcript->get_all_Exons}) {
            if ($exon1->overlaps($cluster_exon) && $exon1->strand == $cluster_exon->strand) {
              push (@matching_clusters, $cluster);
              next CLUSTER;
            }
          }
        }
      }
    }
    
    if (scalar(@matching_clusters) == 0) {
      my @newcluster;
      push(@newcluster,$tran);
      push(@clusters,\@newcluster);
    } 
    elsif (scalar(@matching_clusters) == 1) {
      push @{$matching_clusters[0]}, $tran;
      
    } 
    else {
      # Merge the matching clusters into a single cluster
      my @new_clusters;
      my @merged_cluster;
      foreach my $clust (@matching_clusters) {
        push @merged_cluster, @$clust;
      }
      push @merged_cluster, $tran;
      push @new_clusters,\@merged_cluster;
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
      @clusters =  @new_clusters;
    }
  }
  
  # safety and sanity checks
  $self->check_Clusters(scalar(@transcripts), \@clusters);

  # make and store genes
  
  print STDERR scalar(@clusters)." created, turning them into genes...\n";
  my @genes;
  foreach my $cluster(@clusters){
    my $count = 0;
    my $gene = new Bio::EnsEMBL::Gene;
    
    # sort them, get the longest CDS + UTR first (like in prune_Transcripts() )
    my @sorted_transcripts = $self->_bin_sort_transcripts( @{$cluster} );
    foreach my $transcript( @sorted_transcripts ){
	if ($count < 10) {
	    $gene->add_Transcript($transcript);
	    print STDERR "accepting:\n";
	    print STDERR "check this, prune_Transcripts may not be working properly!\n";
	    $self->_print_Transcript($transcript);
	}
	$count++;
    }
    
    # prune out duplicate exons
    my $new_gene = $self->prune_Exons($gene);
    push( @genes, $new_gene );
   }
  
  $self->final_genes(@genes);
  return @genes;
}

############################################################

sub check_Clusters{
  my ($self, $num_transcripts, $clusters) = @_;
  #Safety checks
  my $ntrans = 0;
  my %trans_check_hash;
  foreach my $cluster (@$clusters) {
    $ntrans += scalar(@$cluster);
    foreach my $trans (@$cluster) {
      if (defined($trans_check_hash{$trans})) {
        $self->throw("Transcript " . $trans->dbID . " added twice to clusters\n");
      }
      $trans_check_hash{$trans} = 1;
    }
    if (!scalar(@$cluster)) {
      $self->throw("Empty cluster");
    }
  }
  if ($ntrans != $num_transcripts) {
    $self->throw("Not all transcripts have been added into clusters $ntrans and " . $num_transcripts. " \n");
  } 
  #end safety checks
  return;
}


############################################################

sub by_transcript_high {
  my $alow;
  my $blow;
  my $ahigh;
  my $bhigh;

  if ($a->start_Exon->strand == 1) {
    $alow = $a->start_Exon->start;
    $ahigh = $a->end_Exon->end;
  } else {
    $alow = $a->end_Exon->start;
    $ahigh = $a->start_Exon->end;
  }

  if ($b->start_Exon->strand == 1) {
    $blow = $b->start_Exon->start;
    $bhigh = $b->end_Exon->end;
  } else {
    $blow = $b->end_Exon->start;
    $bhigh = $b->start_Exon->end;
  }

  if ($ahigh != $bhigh) {
    return $ahigh <=> $bhigh;
  } else {
    return $alow <=> $blow;
  }
}



############################################################

sub prune_Exons {
  my ($self,$gene) = @_;
  
  my @unique_Exons; 
  
  # keep track of all unique exons found so far to avoid making duplicates
  # need to be very careful about translation->start_Exon and translation->end_Exon
  
  foreach my $tran (@{$gene->get_all_Transcripts}) {
    my @newexons;
    foreach my $exon (@{$tran->get_all_Exons}) {
      my $found;
      #always empty
    UNI:foreach my $uni (@unique_Exons) {
	if ($uni->start  == $exon->start  &&
	    $uni->end    == $exon->end    &&
	    $uni->strand == $exon->strand &&
	    $uni->phase  == $exon->phase  &&
	    $uni->end_phase == $exon->end_phase
	   ) {
	  $found = $uni;
	  last UNI;
	}
      }
      if (defined($found)) {
	push(@newexons,$found);
	if ($exon == $tran->translation->start_Exon){
	  $tran->translation->start_Exon($found);
	}
	if ($exon == $tran->translation->end_Exon){
	  $tran->translation->end_Exon($found);
	}
      } else {
	push(@newexons,$exon);
	push(@unique_Exons, $exon);
      }
    }          
    $tran->flush_Exon;
    foreach my $exon (@newexons) {
      $tran->add_Exon($exon);
    }
  }
  return $gene;
}

############################################################

=head2 get_Predictions

Description:  gets your favourite ab initio predictions (genscan,genefinder,fgenesh,genecooker...8-)
Returns : none, but $self->predictions is filled
Args    : none

=cut
  
sub get_Predictions {
  my ($self) = @_;
  my @checked_predictions;
  foreach my $prediction ( @{ $self->slice->get_all_PredictionTranscripts } ){
    #unless ( $self->_check_Transcript( $prediction ) ){
    #  next;
    #}
    push ( @checked_predictions, $prediction );
  }
  $self->predictions(@checked_predictions);
}

############################################################

=head2 get_Similarities

 Title   : get_Similarities
 Usage   : $self->get_Similarities
 Function: gets similarity features for this region
 Returns : none, but $self->feature is filled
 Args    : none

=cut

sub get_Similarities {
  my ($self) = @_;

  my @features = @{ $self->slice->get_all_SimilarityFeatures('',$GB_MIN_FEATURE_SCORE) };
  
  my %idhash;
  my @other_features;
  
  foreach my $feature (@features) {
    unless ( $idhash{ $feature->hseqname } ){
      $idhash{ $feature->hseqname } = [];
    }
    if ($feature->length > $GB_MIN_FEATURE_LENGTH) {
      if ($feature->isa("Bio::EnsEMBL::BaseAlignFeature")) {
	push( @{ $idhash{ $feature->hseqname } }, $feature );
      }
    }
    else {
      push(@other_features,$feature);
    }
  }
  
  my @merged_features = $self->merge(\%idhash);

  my @newfeatures;
  push(@newfeatures,@merged_features);
  push(@newfeatures,@other_features);
  
  $self->features(@newfeatures);
}
   

############################################################


=head2 merge

 Description: wicked meethod that merges two or more homol features into one if they are close enough together
  Returns   : nothing
  Args      : none

=cut

sub merge {
  my ($self,$feature_hash,$overlap,$query_gap,$homol_gap) = @_;
  
  $overlap   = 20  unless $overlap;
  $query_gap = 15  unless $query_gap;
  $homol_gap = 15  unless $homol_gap;
  
  my @mergedfeatures;
  
  foreach my $id (keys %{ $feature_hash }) {
    
    my $count = 0;
    my @newfeatures;
    my @features = @{$feature_hash->{$id}};
    
    @features = sort { $a->start <=> $b->start} @features;
    
    # put the first feature in the new array;
    push(@newfeatures,$features[0]);
    
    for (my $i=0; $i < $#features; $i++) {
      my $id  = $features[$i]  ->id;
      my $id2 = $features[$i+1]->id;
      
      # First case is if start of next hit is < end of previous
      if ( $features[$i]->end > $features[$i+1]->start && 
	  ($features[$i]->end - $features[$i+1]->start + 1) < $overlap) {
	
	if ($features[$i]->strand == 1) {
	  $newfeatures[$count]-> end($features[$i+1]->end);
	  $newfeatures[$count]->hend($features[$i+1]->hend);
	} 
	else {
	  $newfeatures[$count]-> end($features[$i+1]->end);
	  $newfeatures[$count]->hend($features[$i+1]->hstart);
	}
	
	# Take the max score
	if ($features[$i+1]->score > $newfeatures[$count]->score) {
	  $newfeatures[$count]->score($features[$i+1]->score);
	}
	
	if ($features[$i+1]->hstart == $features[$i+1]->hend) {
	  $features[$i+1]->strand($features[$i]->strand);
	}
	
	# Allow a small gap if < $query_gap, $homol_gap
      } elsif (($features[$i]->end < $features[$i+1]->start) &&
	       abs($features[$i+1]->start - $features[$i]->end) <= $query_gap) {
	
	if ($features[$i]->strand eq "1") {
	  $newfeatures[$count]->end($features[$i+1]->end);
	  $newfeatures[$count]->hend($features[$i+1]->hend);
	} else {
	  $newfeatures[$count]->end($features[$i+1]->end);
	  $newfeatures[$count]->hstart($features[$i+1]->hstart);
	}
	
	if ($features[$i+1]->score > $newfeatures[$count]->score) {
	  $newfeatures[$count]->score($features[$i+1]->score);
	}
	
	if ($features[$i+1]->hstart == $features[$i+1]->hend) {
	  $features[$i+1]->strand($features[$i]->strand);
	}
	
      } else {
	# we can't extend the merged homologies so start a
	# new homology feature
	
	# first do the coords on the old feature
	if ($newfeatures[$count]->hstart > $newfeatures[$count]->hend) {
	  my $tmp = $newfeatures[$count]->hstart;
	  $newfeatures[$count]->hstart($newfeatures[$count]->hend);
	  $newfeatures[$count]->hend($tmp);
	}
	
	$count++;
	$i++;
	
	push(@newfeatures,$features[$i]);
	$i--;
      }
    }
    
    # Adjust the last new feature coords
    if ($newfeatures[$#newfeatures]->hstart > $newfeatures[$#newfeatures]->hend) {
      my $tmp = $newfeatures[$#newfeatures]->hstart;
      $newfeatures[$#newfeatures]->hstart($newfeatures[$#newfeatures]->hend);
      $newfeatures[$#newfeatures]->hend($tmp);
    }
    
    my @pruned = $self->prune_features(@newfeatures);
    
    push(@mergedfeatures,@pruned);
  }
  return @mergedfeatures;
}

############################################################

=head2 prune_features

 Description: prunes out duplicated features
 Returntype : array of Bio::EnsEMBL::SeqFeature
 Args       : array of Bio::EnsEMBL::SeqFeature

=cut

sub prune_features {
  my ($self,@features)  = @_;
    
  my @pruned;

  @features = sort {$a->start <=> $b->start} @features;

  my $prev = -1;

  F: 
  foreach  my $f (@features) {
    if ($prev != -1 && $f->hseqname eq $prev->hseqname &&
	$f->start   == $prev->start &&
	$f->end     == $prev->end   &&
	$f->hstart  == $prev->hstart &&
	$f->hend    == $prev->hend   &&
	$f->strand  == $prev->strand &&
	$f->hstrand == $prev->hstrand) {
    } 
    else {
      push(@pruned,$f);
      $prev = $f;
    }
  }
  return @pruned;
}

############################################################

sub _bin_sort_transcripts{
    my ($self,@transcripts) = @_;

    my %lengths;    
    my @newgenes;
    my $max_num_exons = 0;
    
    # sizehash holds transcript length - based on sum of exon lengths
    my %sizehash;
    
    # orfhash holds orf length - based on sum of translateable exon lengths
    my %orfhash;
    
    my %tran2orf;
    my %tran2length;
    
    # keeps track of to which transcript(s) this exon belong
    my %exon2transcript;
    
    foreach my $tran (@transcripts) {
      
      # keep track of number of exons in multiexon transcripts
      my @exons = @{ $tran->get_all_Exons };
      if(scalar(@exons) > $max_num_exons){ 
	  $max_num_exons = scalar(@exons); 
      }
      
      # total exon length
      my $length = 0;
      foreach my $e ( @exons ){
	  $length += $e->end - $e->start + 1;
	  push ( @{ $exon2transcript{ $e } }, $tran );
      }
      $sizehash{$tran} = $length;
      $tran2length{ $tran } = $length;
      
      # now for ORF length
      $length = 0;
      foreach my $e(@{$tran->get_all_translateable_Exons}){
	$length += $e->end - $e->start + 1;
      }
      $tran2orf{ $tran } = $length;
      push(@{$orfhash{$length}}, $tran);
   }  

    # VAC 15/02/2002 sort transcripts based on total exon length - this
    # introduces a problem - we can (and have) masked good transcripts
    # with long translations in favour of transcripts with shorter
    # translations and long UTRs that are overall slightly longer. This
    # is not good.
    # better way? hold both total exon length and length of translateable exons. Then sort:
    # long translation + UTR > long translation no UTR > short translation + UTR > short translation no UTR
    ##########
    @transcripts = ();
    ##########

    # sort first by orfhash{'length'}
    my @orflengths = sort {$b <=> $a} (keys %orfhash);
    
    # strict sort by translation length is just as wrong as strict sort by UTR length
    # bin translation lengths - 4 bins (based on 25% length diff)? 10 bins (based on 10%)?
    my %orflength_bin;
    my $numbins = 4;
    my $currbin = 1;
    
    foreach my $orflength(@orflengths){
      last if $currbin > $numbins;
      my $percid = ($orflength*100)/$orflengths[0];
      if ($percid > 100) { 
	  $percid = 100; 
      }
      my $currthreshold = $currbin * (100/$numbins);
      $currthreshold    = 100 - $currthreshold;
      
      if($percid <$currthreshold) { 
	  $currbin++; 
      }
      my @tmp = @{$orfhash{$orflength}};
      push(@{$orflength_bin{$currbin}}, @{$orfhash{$orflength}});
    }
    
    # now, foreach bin in %orflengthbin, sort by exonlength
    $currbin = 1;
  EXONLENGTH_SORT:
    while( $currbin <= $numbins){
	if(!defined $orflength_bin{$currbin} ){
	    $currbin++;
	    next EXONLENGTH_SORT;
	}	
	my @sorted_transcripts = sort { $sizehash{$b} <=> $sizehash{$a} } @{$orflength_bin{$currbin}};
	push(@transcripts, @sorted_transcripts);
	$currbin++;
    }
    
    #test
    my $debug = 0;
    if ($debug == 1 ){
      
      print STDERR "2.- sorted transcripts:\n";
      
      
      foreach my $tran (@transcripts){
	if ( $tran->dbID ){
	  print STDERR $tran->dbID." ";
	}
	print STDERR "orf_length: $tran2orf{$tran}, total_length: $tran2length{$tran}\n";
	my @exons = @{$tran->get_all_Exons};
	if ( $exons[0]->strand == 1 ){
	  @exons = sort { $a->start <=> $b->start } @exons;
	}
	else{
	  @exons = sort { $b->start <=> $a->start } @exons;
	}
	foreach my $exon ( @{$tran->get_all_Exons} ){
	  print "  ".$exon->start."-".$exon->end." ".( $exon->end - $exon->start + 1)." phase: ".$exon->phase." end_phase ".$exon->end_phase." strand: ".$exon->strand."\n";
      }
    }
  }
    return @transcripts;
}


############################################################
#
# GETSET METHODS
#
############################################################

# get/set method holding a reference to the db with genewise and combined genes
# this reference is set in Bio::EnsEMBL::Pipeline::RunnableDB::Gene_Builder

sub genes_db{
 my ($self,$genes_db) = @_;
 if ( $genes_db ){
   $self->{_genes_db} = $genes_db;
 }
 return $self->{_genes_db};
}

############################################################

sub genewise_combined_Transcripts {
    my ($self,@transcripts) = @_;

    if (!defined($self->{_genewise_andthelike_transcripts})) {
        $self->{_genewise_andthelike_transcripts} = [];
    }

    if (scalar @transcripts > 0) {
	push(@{$self->{_genewise_andthelike_transcripts}},@transcripts);
    }

    return @{$self->{_genewise_andthelike_transcripts}};
}

############################################################

=head2 my_genes

 Description: this holds and returns the genes that are produced after putting together genewise, combined and
              processed_supporte_ab_initio predictions and removing the redundant set, giving priority
              to long CDSs + UTR

=cut


sub my_genes {
  my ($self,@genes) = @_;
  
  unless($self->{_my_genes}){
    $self->{_my_genes} = [];
  }

  if (@genes){
    push(@{$self->{_my_genes}},@genes);
  }
  return @{$self->{_my_genes}};
}

############################################################

=head2 final_genes

 Descripton: this holds/returns the final genes produced after clustering transcripts and sharing common exons

=cut

sub final_genes{
  my ($self, @genes);
  unless ( $self->{_final_genes} ){
      $self->{_final_genes} = [];
  }
  if ( @genes ){
    push( @{$self->{_final_genes}}, @genes );
  }
  return @{$self->{_final_genes}};
}

############################################################

=head2 gene_types

 Description: get/set for the type(s) of genes (usually TGE_gw, similarity_genewise and combined_e2g genes) 
              to be used in the genebuilder they get set in new()
              Does not include the ab inition predictions
=cut

sub gene_types {
  my ($self,$type) = @_;

  if (defined($type)) {
     push(@{$self->{_gene_types}},$type);
  }

  return @{$self->{_gene_types}};
}
############################################################

=head2 predictions

 Description: get/set for the PredictionTranscripts. It is  set in new()

=cut

sub predictions {
  my ($self,@predictions) = @_;

  if ( @predictions ) {
     push(@{$self->{_predictions}},@predictions);
  }
  return @{$self->{_predictions}};
}

############################################################

sub features {
  my ($self,@features) = @_;
  
  if (!defined($self->{_feature})) {
    $self->{_feature} = [];
  }
  if ( scalar @features ) {
    push(@{$self->{_feature}},@features);
  }
  return @{$self->{_feature}};
}

############################################################

sub slice {
    my ($self,$slice) = @_;
    
    if (defined($slice)) {
      $self->{_slice} = $slice;
    }
    return $self->{_slice};
}

#############################################################################
# 
# Printing routines
#
#############################################################################

sub print_Exon {
  my ($self,$exon) = @_;
  
  print STDERR 
    $exon->seqname." ".$exon->start ."-".$exon->end." ".$exon->strand." [".$exon->phase.",".$exon->end_phase."]\n";
}

############################################################
  
sub _print_Transcript{
  my ($self,$transcript) = @_;
  my @exons = @{$transcript->get_all_Exons};
  my $id;
  if ( $transcript->dbID ){
    $id = $transcript->dbID;
  }
  else{
    $id = "no id";
  }
  print STDERR "transcript id: ".$id."\n";
  foreach my $exon ( @exons){
    $self->print_Exon($exon);
  }
  #print STDERR "Translation : ".$transcript->translation."\n";
  if ( $transcript->translation ){
    print STDERR "translation start exon: ".
      $transcript->translation->start_Exon->start."-".$transcript->translation->start_Exon->end.
	" start: ".$transcript->translation->start."\n";
    print STDERR "translation end exon: ".
      $transcript->translation->end_Exon->start."-".$transcript->translation->end_Exon->end.
	  " end: ".$transcript->translation->end."\n";
    print STDERR "\n";
  }
}

############################################################

sub print_ExonPairs {
  my ($self) = @_;
  
  foreach my $pair ($self->get_all_ExonPairs) {
    $self->print_ExonPair($pair);
  }
}

############################################################

sub print_ExonPair {
  my ($self,$pair) = @_;
  
  $self->print_Exon($pair->exon1);
  $self->print_Exon($pair->exon2);
  
  print(STDERR "\nExon Pair (splice - " . $pair->splice_seq->seq . ")\n");
  
  foreach my $ev ($pair->get_all_Evidence) {
    print(STDERR "   -  " . $ev->hseqname . "\t" . $ev->hstart . "\t" . $ev->hend . "\t" . $ev->strand . "\n");
  }
}

############################################################


=head2 split_transcript

 Title   : split_transcript 
 Usage   : my @splits = $self->split_transcript($transcript)
 Function: splits a transcript into multiple transcripts at long introns. Rejects single exon 
           transcripts that result. 
 Returns : @Bio::EnsEMBL::Transcript
 Args    : Bio::EnsEMBL::Transcript

=cut


sub split_transcript{
  my ($self, $transcript) = @_;
  $transcript->sort;
  my @split_transcripts   = ();

  if(!($transcript->isa("Bio::EnsEMBL::Transcript"))){
    $self->warn("[$transcript] is not a Bio::EnsEMBL::Transcript - cannot split");
    return (); # empty array
  }
  
  my $prev_exon;
  my $exon_added = 0;
  my $curr_transcript = new Bio::EnsEMBL::Transcript;
  my $translation     = new Bio::EnsEMBL::Translation;
  $curr_transcript->type($transcript->type);
  $curr_transcript->translation($translation);


EXON:   
  foreach my $exon( @{$transcript->get_all_Exons} ){
    
    $exon_added = 0;
      # is this the very first exon?
    if($exon == $transcript->start_Exon){
      $prev_exon = $exon;      
      # set $curr_transcript->translation start and start_exon
      $curr_transcript->add_Exon($exon);
      $exon_added = 1;
      $curr_transcript->translation->start_Exon($exon);
      $curr_transcript->translation->start($transcript->translation->start);
      push(@split_transcripts, $curr_transcript);
      next EXON;
    }
    
    if ($exon->strand != $prev_exon->strand){
      return (); # empty array
    }

    # We need to start a new transcript if the intron size between $exon and $prev_exon is too large
    my $intron = 0;
    if ($exon->strand == 1) {
      $intron = abs($exon->start - $prev_exon->end + 1);
    } else {
      $intron = abs($exon->end   - $prev_exon->start + 1);
    }
    
    if ($intron > $GB_GENSCAN_MAX_INTRON) {
      $curr_transcript->translation->end_Exon($prev_exon);
      # need to account for end_phase of $prev_exon when setting translation->end
      $curr_transcript->translation->end($prev_exon->end - $prev_exon->start + 1 - $prev_exon->end_phase);
      
      # start a new transcript 
      my $t  = new Bio::EnsEMBL::Transcript;
      my $tr = new Bio::EnsEMBL::Translation;
      $t->type($transcript->type);
      $t->translation($tr);

      # add exon unless already added, and set translation start and start_exon
      $t->add_Exon($exon) unless $exon_added;
      $exon_added = 1;

      $t->translation->start_Exon($exon);

      if ($exon->phase == 0) {
	$t->translation->start(1);
      } elsif ($exon->phase == 1) {
	$t->translation->start(3);
      } elsif ($exon->phase == 2) {
	$t->translation->start(2);
      }

      # start exon always has phase 0
      $exon->phase(0);      

      # this new transcript becomes the current transcript
      $curr_transcript = $t;

      push(@split_transcripts, $curr_transcript);
    }

    if($exon == $transcript->end_Exon){
      # add it unless already added
      $curr_transcript->add_Exon($exon) unless $exon_added;
      $exon_added = 1;

      # set $curr_transcript end_exon and end
      $curr_transcript->translation->end_Exon($exon);
      $curr_transcript->translation->end($transcript->translation->end);
    }

    else{
      # just add the exon
      $curr_transcript->add_Exon($exon) unless $exon_added;
    }
    
    # this exon becomes $prev_exon for the next one
    $prev_exon = $exon;

  }

  # discard any single exon transcripts
  my @t = ();
  my $count = 1;
  
  foreach my $st(@split_transcripts){
    $st->sort;
    my @ex = @{$st->get_all_Exons};
    if(scalar(@ex) > 1){
      $st->{'temporary_id'} = $transcript->dbID . "." . $count;
      $count++;
      push(@t, $st);
    }
  }

  return @t;

}

############################################################

=head2 validate_transcript

 Title   : validate_transcript 
 Usage   : my @valid = $self->validate_transcript($transcript)
 Function: Validates a transcript - rejects if mixed strands, splits if long introns
 Returns : @Bio::EnsEMBL::Transcript
 Args    : Bio::EnsEMBL::Transcript

=cut

sub validate_transcript{
  my ($self, $transcript) = @_;
  my @valid_transcripts;

  my $valid = 1;
  my $split = 0;

  # check exon phases:
  my @exons = @{$transcript->get_all_Exons};
  $transcript->sort;
  for (my $i=0; $i<(scalar(@exons-1)); $i++){
    my $end_phase = $exons[$i]->end_phase;
    my $phase    = $exons[$i+1]->phase;
    if ( $phase != $end_phase ){
      $self->warn("rejecting transcript with inconsistent phases( $phase <-> $end_phase) ");
      return undef;
    }
  }
  

  my $previous_exon;
  foreach my $exon (@exons){
    if (defined($previous_exon)) {
      my $intron;
      
      if ($exon->strand == 1) {
	$intron = abs($exon->start - $previous_exon->end + 1);
      } else {
	$intron = abs($exon->end   - $previous_exon->start + 1);
      }
      
      if ($intron > $GB_GENSCAN_MAX_INTRON) {
	print STDERR "Intron too long $intron  for transcript " . $transcript->{'temporary_id'} . "\n";
	$split = 1;
	$valid = 0;
      }
      
      if ($exon->strand != $previous_exon->strand) {
	print STDERR "Mixed strands for gene " . $transcript->{'temporary_id'} . "\n";
	$valid = 0;
	return;
      }
    }
    $previous_exon = $exon;
  }
  
       if ($valid) {
	 push(@valid_transcripts,$transcript);
       }
       elsif ($split){
	 # split the transcript up.
	 my @split_transcripts = $self->split_transcript($transcript);
    push(@valid_transcripts, @split_transcripts);
  }
       return @valid_transcripts;
}

############################################################

=head2 transfer_supporting_evidence

 Title   : transfer_supporting_evidence
 Usage   : $self->transfer_supporting_evidence($source_exon, $target_exon)
 Function: Transfers supporting evidence from source_exon to target_exon, 
           after checking the coordinates are sane and that the evidence is not already in place.
 Returns : nothing, but $target_exon has additional supporting evidence

=cut

sub transfer_supporting_evidence{
  my ($self, $source_exon, $target_exon) = @_;
  
  my @target_sf = @{$target_exon->get_all_supporting_features};
  #  print "target exon sf: \n";
  #  foreach my $tsf(@target_sf){ print STDERR $tsf; $self->print_FeaturePair($tsf); }
  
  #  print "source exon: \n";
 
  # keep track of features already transferred, so that we do not duplicate
  my %unique_evidence;
  my %hold_evidence;

 SOURCE_FEAT:
  foreach my $feat ( @{$source_exon->get_all_supporting_features}){
    next SOURCE_FEAT unless $feat->isa("Bio::EnsEMBL::FeaturePair");
    
    # skip duplicated evidence objects
    next SOURCE_FEAT if ( $unique_evidence{ $feat } );
    
    # skip duplicated evidence 
    if ( $hold_evidence{ $feat->hseqname }{ $feat->start }{ $feat->end }{ $feat->hstart }{ $feat->hend } ){
      #print STDERR "Skipping duplicated evidence\n";
      next SOURCE_FEAT;
    }

    #$self->print_FeaturePair($feat);
    
  TARGET_FEAT:
    foreach my $tsf (@target_sf){
      next TARGET_FEAT unless $tsf->isa("Bio::EnsEMBL::FeaturePair");
      
      if($feat->start    == $tsf->start &&
	 $feat->end      == $tsf->end &&
	 $feat->strand   == $tsf->strand &&
	 $feat->hseqname eq $tsf->hseqname &&
	 $feat->hstart   == $tsf->hstart &&
	 $feat->hend     == $tsf->hend){
	
	#print STDERR "feature already in target exon\n";
	next SOURCE_FEAT;
      }
    }
    #print STDERR "from ".$source_exon->{'temporary_id'}." to ".$target_exon->{'temporary_id'}."\n";
    #$self->print_FeaturePair($feat);
    $target_exon->add_supporting_features($feat);
    $unique_evidence{ $feat } = 1;
    $hold_evidence{ $feat->hseqname }{ $feat->start }{ $feat->end }{ $feat->hstart }{ $feat->hend } = 1;
  }
}

############################################################

sub print_FeaturePair{
  my ($self, $fp) = @_;
  return unless $fp->isa("Bio::EnsEMBL::FeaturePair");
  print STDERR $fp;
  print STDERR $fp->seqname . " " .
    $fp->start . " " .
      $fp->end . " " .
	$fp->strand . " " .
	  $fp->hseqname . " " .
	    $fp->hstart . " " .
	      $fp->hend . "\n";
}

############################################################

sub _check_Transcript{
  my ($self,$transcript) = @_;
  my $slice = $self->slice;
  
  my $id = $self->transcript_id( $transcript );
  
  my $valid = 1;

  # check that transcripts are not completely outside the slice
  if ( $transcript->start > $slice->length || $transcript->end < 1 ){
    print STDERR "transcript $id outside the slice\n";
    $valid = 0;
  }
  # allow transcripts that fall partially off the slice only at one end, the 'higher' end of the slice
  elsif ( $transcript->start < 1 && $transcript->end > 1 ){
      print STDERR "transcript $id falls off the slice by its lower end\n";
    $valid = 0;
  }
  
  # sort the exons 
  $transcript->sort;
  my @exons = @{$transcript->get_all_Exons};
  
  if ($#exons > 0) {
    for (my $i = 1; $i <= $#exons; $i++) {
      
      # check phase consistency:
      if ( $exons[$i-1]->end_phase != $exons[$i]->phase  ){
	print STDERR "transcript $id has phase inconsistency\n";
	$valid = 0;
	last;
      }
      
      # check for folded transcripts
      if ($exons[0]->strand == 1) {
	if ($exons[$i]->start < $exons[$i-1]->end) {
	  print STDERR "transcript $id folds back on itself\n";
	  $valid = 0;
	} 
      } 
      elsif ($exons[0]->strand == -1) {
	if ($exons[$i]->end > $exons[$i-1]->start) {
	  print STDERR "transcript $id folds back on itself\n";
	  $valid = 0;
	} 
      }
    }
  }
  if ($valid == 0 ){
    $self->_print_Transcript($transcript);
  }
  return $valid;
}


############################################################

sub _check_Translation{
  my ($self,$transcript) = @_;
  
  my $id = $self->transcript_id( $transcript );
  
  my $valid = 1;
  
  # check that they have a translation
  my $translation = $transcript->translation;
  my $sequence;
  eval{
    $sequence = $transcript->translate;
  };
  unless ( $sequence ){
    print STDERR "transcript $id has no translation\n";
    return 0;
  }
  if ( $sequence ){
    my $peptide = $sequence->seq;
    if ( $peptide =~ /\*/ ){
      print STDERR "translation of transcript $id has STOP codons\n";
      $valid = 0;
    }
  }
  if ($valid == 0 ){
    $self->_print_Transcript($transcript);
  }
  return $valid;
}


############################################################

sub transcript_id {
  my ( $self, $t ) = @_;
  my $id;
  if ( $t->stable_id ){
    $id = $t->stable_id;
  }
  elsif( $t->dbID ){
    $id = $t->dbID;
  }
  elsif( $t->temporary_id ){
    $id = $t->temporary_id;
  }
  else{
    $id = 'no-id';
  }
  return $id;
}

############################################################

1;
