
# Cared for by Eduardo Eyras  <eae@sanger.ac.uk>
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
								       -db        => $db,
								       -input_id  => $id
								      );
    $obj->fetch_input
    $obj->run

    my @newfeatures = $obj->output;


=head1 DESCRIPTION

EST_GeneBuilder processes est2genome gene predictions and feed them
to genomewise to create transcripts with translations and UTRs.

=head1 CONTACT

eae@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::RunnableDB::EST_GeneBuilder;

use diagnostics;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::MiniGenomewise;
use Bio::EnsEMBL::Pipeline::Runnable::Genomewise;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::DenormGeneAdaptor;
use Bio::EnsEMBL::Pipeline::Runnable::ClusterMerge;
use Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster;
use Bio::EnsEMBL::Pipeline::Tools::TranslationUtils;

use Bio::EnsEMBL::Pipeline::Config::cDNAs_ESTs::EST_GeneBuilder_Conf qw (
									 EST_INPUTID_REGEX
									 EST_REFDBHOST
									 EST_REFDBUSER
									 EST_REFDBNAME
									 EST_REFDBPASS
									 EST_DBNAME
									 EST_DBHOST
									 EST_DBUSER
									 EST_DBPASS     
									 EST_GENOMIC
									 EST_GENEBUILDER_INPUT_GENETYPE
									 EST_USE_DENORM_GENES
									 EST_MIN_INTRON_SIZE
									 BRIDGE_OVER_SMALL_INTRONS
									 EST_MAX_INTRON_SIZE
									 EST_MAX_EVIDENCE_DISCONTINUITY
									 EST_GENEBUILDER_INTRON_MISMATCH
									 ESTGENE_TYPE
									 USE_cDNA_DB
									 cDNA_DBNAME
									 cDNA_DBHOST
									 cDNA_DBUSER
									 cDNA_DBPASS
									 cDNA_GENETYPE
									 REJECT_SINGLE_EXON_TRANSCRIPTS
									 USE_GENOMEWISE
									 GENOMEWISE_SMELL
									 EST_MIN_EXON_SIZE
									 EST_GENEBUILDER_COMPARISON_LEVEL
									 EST_GENEBUILDER_SPLICE_MISMATCH
									 EST_GENEBUILDER_INTRON_MISMATCH
									 EST_GENEBUILDER_EXON_MATCH
									 CHECK_SPLICE_SITES
									 RAISE_SINGLETON_COVERAGE
									 FILTER_ON_SINGLETON_SIZE
									 MAX_TRANSCRIPTS_PER_GENE
									 CLUSTERMERGE_MIN_EVIDENCE_NUMBER
									 EST_USE_DENORM_GENES
);


@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);

    my $refdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
						   -host             => $EST_REFDBHOST,
						   -user             => $EST_REFDBUSER,
						   -dbname           => $EST_REFDBNAME,
						   -pass             => $EST_REFDBPASS,
						  );
    
    my $est_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
						   -host             => $EST_DBHOST,
						   -user             => $EST_DBUSER,
						   -dbname           => $EST_DBNAME,
						    -pass             => $EST_DBPASS,
						   );
    
    $est_db->dnadb($refdb);
    $self->est_db($est_db);
    
    $self->db->dnadb($refdb);
   
   
    if ( $USE_cDNA_DB ){
      my $cdna_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
						       -host             => $cDNA_DBHOST,
						       -user             => $cDNA_DBUSER,
						       -dbname           => $cDNA_DBNAME,
						       -pass             => $cDNA_DBPASS,
						       -dnadb            => $refdb,
						      );
      
      $self->cdna_db($cdna_db);
    }
    
    $self->genetype($ESTGENE_TYPE);



    return $self; 
}

############################################################

sub est_db{
    my ($self, $est_db) = @_;
    if ($est_db){
	$self->{_est_db} = $est_db;
    }
    return $self->{_est_db};
}
############################################################

sub cdna_db{
 my ($self, $cdna_db);
 if ($cdna_db){
  $self->{_cdna_db} = $cdna_db;
 }
 return $self->cdna_db;
}

sub revcomp_query{
    my ($self,$slice) = @_;
    if ($slice){
	$self->{_revcomp_query} = $slice;
    }
    return $self->{_revcomp_query};
}

############################################################

=head2 write_output

    Function:   Writes output data to db

=cut

sub write_output {
  my ($self) = @_;
  
  my $gene_adaptor = $self->db->get_GeneAdaptor;
  
 GENE: 
  foreach my $gene ($self->output) {	
    my $num = scalar(@{$gene->get_all_Transcripts});
    
    eval {
	  $gene_adaptor->store($gene);
	  print STDERR "wrote gene " . $gene->dbID . " with $num transcripts\n";
      }; 
      if( $@ ) {
	  print STDERR "UNABLE TO WRITE GENE\n\n$@\n\nSkipping this gene\n";
	  next;
      }
      my @transcripts = @{ $gene->get_all_Transcripts};
      #foreach my $tran (@transcripts){
      #	Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Evidence($tran);
      #} 
      
  }
}


############################################################

=head2 fetch_input

    Function:   Fetches input data from the database

=cut

sub fetch_input {
    my( $self) = @_;
    my $strand;

    # the type of the genes being read is specified in Bio/EnsEMBL/Pipeline/ESTConf.pm
    my $genetype =  $EST_GENEBUILDER_INPUT_GENETYPE;

    # make sure you have an analysis
    $self->throw("No analysis") unless defined( $self->analysis );

    #print STDERR "Fetching input: " . $self->input_id. " \n";
    $self->throw("No input id") unless defined($self->input_id);
    
    # get genomic region 
    my $input_id    = $self->input_id;
    unless ($input_id =~ /$EST_INPUTID_REGEX/ ){
      $self->throw("input $input_id not compatible with EST_INPUTID_REGEX $EST_INPUTID_REGEX");
    }
    my $chrname  = $1;
    my $chrstart = $2;
    my $chrend   = $3;

    print STDERR "Chromosome id = $chrname , range $chrstart $chrend\n";

    my $slice = $self->est_db->get_SliceAdaptor->fetch_by_chr_start_end($chrname,$chrstart,$chrend);    
    $slice->chr_name($chrname);
    $self->query($slice);
    
    ############################################################
    # forward strand
    ############################################################

    $strand = 1;
    print STDERR "\n****** forward strand ******\n\n";

    # get genes

    my $genes;
    if ($EST_USE_DENORM_GENES){
      my $dga = Bio::EnsEMBL::Pipeline::DBSQL::DenormGeneAdaptor->new($self->est_db);
      print STDERR "getting denormalised genes from slice: ".
	$slice->chr_name.".".$slice->chr_start."-".$slice->chr_end." with genetype $genetype\n";
      $genes = $dga->get_genes_by_Slice_and_type($slice, $genetype);
    } 
    else {
      $genes  = $slice->get_all_Genes_by_type($genetype);
    } 
    
    print STDERR "Number of genes of type $genetype   = " . scalar(@$genes) . "\n";

    my $cdna_slice;
    if ( $USE_cDNA_DB ){
	my $cdna_db = $self->cdna_db;
	
	$cdna_slice = $cdna_db->get_SliceAdaptor->fetch_by_chr_start_end($chrname,$chrstart,$chrend);
	my $cdna_genes  = $cdna_slice->get_all_Genes_by_type($cDNA_GENETYPE);
	print STDERR "Number of genes from cdnas = " . scalar(@$cdna_genes) . "\n";
	push (@$genes, @$cdna_genes);
	
    }
    
    my @forward_transcripts;
    my $single = 0;
    
    # split by strand
  GENE:    
    foreach my $gene (@$genes) {
	foreach my $transcript ( @{$gene->get_all_Transcripts} ){
	    
	    # Don't skip genes with one exon, potential info for UTRs
	    my $exons = $transcript->get_all_Exons;
	    if(scalar(@$exons) == 1){
		$single++;
	    }
	    
	    # keep only genes in the forward strand
	    if ($exons->[0]->strand == 1){
		push (@forward_transcripts, $transcript );
	    }
	}
    }
    
    print STDERR "In EST_GeneBuilder.fetch_input(): ".scalar(@forward_transcripts) . " forward strand genes\n";
    
    # process transcripts in the forward strand
    
    
    if( scalar(@forward_transcripts) ){
      $self->_forward_transcripts( @forward_transcripts );
    }

    ############################################################
    # reverse strand
    ############################################################    
    print STDERR "\n****** reverse strand ******\n\n";
    $strand = -1;
    
    # this will return a slice which corresponds to the reversed complement of $slice:
    my $rev_slice = $slice->invert;
    $self->revcomp_query($rev_slice);

    my $revgenes;
    if ($EST_USE_DENORM_GENES){
      my $dga = Bio::EnsEMBL::Pipeline::DBSQL::DenormGeneAdaptor->new($self->est_db);
      $revgenes = $dga->get_genes_by_Slice_and_type($rev_slice, $genetype);
    } else {
      $revgenes  = $rev_slice->get_all_Genes_by_type($genetype);
    } 

    my @reverse_transcripts;
    
    if ( $USE_cDNA_DB ){
	my $cdna_revslice = $cdna_slice->invert;
	my $cdna_revgenes  = $cdna_revslice->get_all_Genes_by_type($cDNA_GENETYPE);
	print STDERR "Number of genes from cdnas = " . scalar(@$cdna_revgenes) . "\n";
	push ( @$revgenes, @$cdna_revgenes ); 
    }
    
    $single=0;
  REVGENE:    
    foreach my $gene (@$revgenes) {
      foreach my $transcript ( @{$gene->get_all_Transcripts} ){
	
	my @exons = @{$transcript->get_all_Exons};
	
	# DON'T throw away single-exon genes
	if(scalar(@exons) == 1){
	  $single++;
	}
	
	# these are really - strand, but the Slice is reversed, so they are relatively + strand
	if( $exons[0]->strand == 1){
	  push (@reverse_transcripts, $transcript);
	}
      }
    }
    print STDERR "In EST_GeneBuilfer.fetch_input(): ".scalar(@reverse_transcripts) . " reverse strand genes\n";
    
    if(scalar(@reverse_transcripts)){
      $self->_reverse_transcripts(@reverse_transcripts);
    }
  }

############################################################

sub _forward_transcripts{
  my ($self, @transcripts) = @_;
  if ( @transcripts ){
    $self->{_forward_transcripts} = \@transcripts;
  }
  return @{$self->{_forward_transcripts}};
}

############################################################

sub _reverse_transcripts{
  my ($self, @transcripts) = @_;
  if ( @transcripts ){
    $self->{_reverse_transcripts} = \@transcripts;
  }
  return @{$self->{_reverse_transcripts}};
}

############################################################

=head2 _process_Transcripts

    Usage   :   @new_transcripts= $self->_process_Transcripts(@read_transcripts)
    Function:   main magic and witchcraft on the transcripts. 
                It checks, clusters and  merges an input array of transcripts
    Returns :   @Bio::EnsEMBL::Transcript
    Args    :   @Bio::EnsEMBL::Transcript

=cut

sub _process_Transcripts {
  my ($self, $alltranscripts, $strand) = @_;

  print STDERR "EST_GeneBuilder: processing input transcripts...\n";

  # first check transcripts and hold info about est_evidence, etc...
  my @transcripts = $self->_check_Transcripts($alltranscripts,$strand);
  print STDERR scalar(@transcripts)." transcripts returned from _check_Transcripts\n";
  
  # reject ests/cdnas if they have more than one non-standard intron splice site consensus sequence
  # or if the only intron they have is non standard.
  # the standard introns are taken to be:  (GT-AG, AT-AC, GC-AG)
  
  my @checked_transcripts = @transcripts;
  
  
  if ( scalar(@checked_transcripts) == 0 ){
    print STDERR "No transcripts left\n";
    return;
  }
  my $merge_object 
    = Bio::EnsEMBL::Pipeline::Runnable::ClusterMerge
      ->new(
	    -transcripts      => \@checked_transcripts,
	    -comparison_level => $EST_GENEBUILDER_COMPARISON_LEVEL,
	    -splice_mismatch  => $EST_GENEBUILDER_SPLICE_MISMATCH,
	    -intron_mismatch  => $EST_GENEBUILDER_INTRON_MISMATCH,
	    -exon_match       => $EST_GENEBUILDER_EXON_MATCH,
	    -minimum_order    => $CLUSTERMERGE_MIN_EVIDENCE_NUMBER,
	   );
  
  $merge_object->run;
  my @merged_transcripts = $merge_object->output;

  # reject the single exon transcripts
  my @filtered_transcripts;
  if ( $REJECT_SINGLE_EXON_TRANSCRIPTS ){
    @filtered_transcripts = @{$self->_reject_single_exon_Transcripts(@merged_transcripts)};      
    print STDERR scalar(@filtered_transcripts)
      ." transcripts left after rejecting single-exon transcripts\n";
  }
  else{
    @filtered_transcripts = @merged_transcripts;
    print STDERR scalar(@filtered_transcripts)
      ." transcripts obtained ( Not rejecting single-exon transcripts)\n";
  }
  
  return @filtered_transcripts;
}

############################################################

sub _reject_single_exon_Transcripts{
    my ($self,@transcripts) = @_;
    my @filtered_transcripts;
    foreach my $tran (@transcripts){
	unless ( scalar(@{$tran->get_all_Exons}) <= 1 ){
	    push( @filtered_transcripts, $tran );
	}
    }
    return \@filtered_transcripts;
}

############################################################

=head2 _check_Transcripts
    
    Function:   checks transcripts (representing ESTs or cDNAs) for consistency:
                the maximum allowed discontinuity in the evidence is $EST_MAX_EVIDENCE_DISCONTINUITY;
    Returns :   @Bio::EnsEMBL::Transcript (only those that get through the checks)
    Args    :   @Bio::EnsEMBL::Transcript

=cut

sub _check_Transcripts {
  my ($self, $ref_transcripts, $strand) = @_;

  # the maximum allowed discontinuity in EST hits = $EST_MAX_EVIDENCE_DISCONTINUITY;
  # reject ests with introns larger than $EST_MAX_INTRON_SIZE;
  # reject exons that are smaller or equal to $EST_MIN_EXON_SIZE;


  print STDERR "EST_GeneBuilder: checking consistency of transcripts...\n";

  my @allexons;       # here we'll put all exons that pass the check
  my @alltranscripts; # here we'll put all the transcripts that pass the check
  my %hid_trans;
  #my $exon_adaptor    = $self->db->get_ExonAdaptor;
  my $total_rejected        = 0;

  my $slice;
  if ( $strand == +1 ){
    $slice = $self->query;
  }
  else{
    $slice = $self->revcomp_query;
  }
  
  #print STDERR "transcripts:\n";
  #foreach my $t (@$ref_transcripts){
  #    print $t->dbID."\n";
  #}
  
 TRANSCRIPT: 
  foreach my $transcript (@$ref_transcripts){
      
      my $new_transcript = Bio::EnsEMBL::Transcript->new();
      if ( defined $transcript->dbID){
	  $new_transcript->dbID($transcript->dbID);
      }
      
      ############################################################
      # reject the transcripts that fall off the slice at the lower end
      ############################################################
      unless ( $self->_check_Transcript_Location($transcript,$slice,$strand) ){
	  next TRANSCRIPT;
      }
      my @exons = sort { $a->start <=> $b->end } @{$transcript->get_all_Exons};
      my $exon_count = 0;
      my $previous_exon;
    
      ############################################################
      # for single exon ests, take only those that are >= 200bp and have coverage >= 95%
      if ( scalar(@exons) == 1 ){
	
	if ( $FILTER_ON_SINGLETON_SIZE && ( $exons[0]->end - $exons[0]->start + 1) < $FILTER_ON_SINGLETON_SIZE ){
	  next TRANSCRIPT;
	}
	if ( $RAISE_SINGLETON_COVERAGE ){
	  my @evidence = @{$exons[0]->get_all_supporting_features};
	  my $coverage = $evidence[0]->score;
	  if ( $coverage < $RAISE_SINGLETON_COVERAGE ){
	    next TRANSCRIPT;
	  }
	}
      }
      
    EXON:
      foreach my $exon (@exons){
	  
	  my $hstart;
	  my $hend;
	  # get the supporting_evidence for each exon
	  my @sf = sort { $a->hstart <=> $b->hstart } @{$exon->get_all_supporting_features};
	  
	  # check that you get suporting_evidence at all
	  if ( scalar( @sf ) == 0 ){
	      print STDERR "exon ".$exon->dbID." in transcript ".
		  $transcript->dbID."has no supporting evidence - ";
	      if ( $exon->is_sticky ){
		  print STDERR "it is a sticky exon\n";
	      }
	      else{
		  print STDERR "it is NOT a sticky exon\n";
	      }
	  }
	  ############################################################
	  # reject transcripts with too-small exons:
	  ############################################################
	  my $size = $exon->end - $exon->start + 1;
	  if ( $size  < $EST_MIN_EXON_SIZE ){
	    print STDERR "rejecting transcript with very small exon size = $size\n";
	    Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($transcript);
	    #$self->warn("rejecting very small exon: size=$size");
	    next TRANSCRIPT;
	  }
	  ############################################################
	  # check the gap with the evidence of the next exon
	  # if the ESTs are of good quality, this should not reject any
	  ############################################################
	  if ( $exon_count > 1 ){
	      my $est_gap = 0;
	      
	      my @previous_sf = sort { $a->hstart <=> $b->hstart } @{$previous_exon->get_all_supporting_features};
	      
	      if ( @previous_sf ){
		  
		  # if the hstart increases per exon, the EST runs in the same direction of the gene 
		  if ( $previous_sf[0]->hstart < $sf[0]->hstart ){
		      $est_gap = abs( $sf[0]->hstart - $previous_sf[$#previous_sf]->hend ) - 1;
		  }
		  # if hstart decreases that means that the EST runs in the opposite direction
		  elsif (  $previous_sf[0]->hstart > $sf[0]->hstart ){
		      $est_gap = abs( $previous_sf[0]->hstart - $sf[$#sf]->hend) - 1;
		  }
		  # else, same EST piece is hitting two exons, not good!
		  else{
		      print STDERR "same bit of evidence is hitting two exons!\n";
		  }
		  
                  # check the evidence gap between both exons
		  if ( $est_gap > $EST_MAX_EVIDENCE_DISCONTINUITY ){
		      print STDERR "Rejecting transcript: EST evidence with gap too large: $est_gap\n";
		    Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Evidence($transcript);
		      next TRANSCRIPT;
		  }
	      }
	  }
	  ############################################################
	  # reject transcript with too large intron length
	  ############################################################
	  my $intron_length;
	  if ($exon_count > 1 ){
	    my ( $s, $e, $intron_length);
	    $s             = $previous_exon->end;
	    $e             = $exon->start;
	    $intron_length = $e - $s - 1;
	    if ( $intron_length > $EST_MAX_INTRON_SIZE ){
	      print STDERR "Rejecting transcript for having too long intron: $intron_length\n";
	      Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($transcript);
	      next TRANSCRIPT;
	    }
	  }
	  
	  ############################################################
	  # check tiny introns
	  ############################################################
	  if ( $EST_MIN_INTRON_SIZE && $exon_count>1 ){
	    if ( $exon->start - $previous_exon->end - 1 <=  $EST_MIN_INTRON_SIZE ){
	      if ( $BRIDGE_OVER_SMALL_INTRONS ){
		$previous_exon->end( $exon->end );
		next EXON;
	      }
	      ############################################################
	      # if not bridged - reject it
	      ############################################################
	      else{
		print STDERR "Rejecting transcript with small intron size =  ".
		  ($exon->start - $previous_exon->end - 1)."\n";
		Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($transcript);
		next TRANSCRIPT;
	      }
	    }
	  }
      	  
	  $previous_exon = $exon;
	  $new_transcript->add_Exon( $previous_exon );
	  $exon_count++;
	  
      } # end of EXON
      
      ############################################################
      # check the splice sites
      ############################################################
      if ( $CHECK_SPLICE_SITES ){
	next TRANSCRIPT unless $self->check_splice_sites($new_transcript,$strand);
      }
      # if the transcript made it to this point, keep it
      push (@alltranscripts, $new_transcript);
      
    }    # end of TRANSCRIPT
  return @alltranscripts;
  
}

############################################################

sub _check_Transcript_Location{
    my ($self,$transcript, $slice, $strand) = @_;
      
    my $valid = 1;
    
    my $id = Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->transcript_id( $transcript );
    
    # never ever trus the order in which the exons come out! 
    $transcript->sort;

    # check that transcripts are not completely outside the slice
    if ( $transcript->start > $slice->length || $transcript->end < 1 ){
	print STDERR "transcript $id outside the slice\n";
	$valid = 0;
    }
    ############################################################
    # if the transcript is in the forward strand,
    # allow transcripts that fall partially off the slice only at one end, the 'higher' end of the slice
    ############################################################
    elsif ( $strand == 1 && $transcript->start < 1 && $transcript->end > 1 ){
      print STDERR "transcript $id falls off the slice by its lower end\n";
      $valid = 0;
    }
    ############################################################
    # else, if the transcripts are in the reverse strand, since we revcomp
    # the slice, the transcript falling off the lower end would be now fallign off the higher end, and
    # vice versa, so  
    ############################################################
    elsif ( $strand == -1 && $transcript->start <=  $slice->length && $transcript->end >  $slice->length ){
      print STDERR "transcript $id in reverse strand falls off the revcomp-slice by its higher end\n";
      $valid = 0;
    }

    ############################################################
    # check for possible coordinate nonsense in the exons
    ############################################################    
    my @exons = @{$transcript->get_all_Exons};
    if ($#exons > 0) {
	for (my $i = 1; $i <= $#exons; $i++) {
	    
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
       Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($transcript);
    }
    return $valid;
}



############################################################

=head2 _cluster_Transcripts

    Usage   :   @clusters = $self->_cluster_Transcripts(\@transcripts)
    Function:   it clusters transcripts, if run on a long piece of sequence it can be very slow
                since it checks through all previous clusters until it finds the matching one.
                It can be speeded up by only checking a given number of clusters. Rather than doing that,
                I would suggest that short pieces are used: 1MB seems to be o.k.
    
=cut 
  
sub _cluster_Transcripts{
  my ($self,$ref_transcripts) = @_;
  my @transcripts = @{ $ref_transcripts };
  my @clusters;
  print STDERR "EST_GeneBuilder: clustering transcripts...\n";
			 
  # first sort the transcripts by their start position coordinate
  my %start_table;
  my $i=0;
  foreach my $transcript (@transcripts){
    $start_table{$i} = $transcript->start_Exon->start;
    $i++;
  }
  my @sorted_transcripts=();
  foreach my $pos ( sort { $start_table{$a} <=> $start_table{$b} } keys %start_table ){
    push (@sorted_transcripts, $transcripts[$pos]);
  }
  @transcripts = @sorted_transcripts;
  
  # create a new cluster 
  my $cluster = Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster->new();
  my $cluster_count = 1;

  # put the first transcript into this cluster
  $cluster->put_Transcripts( $sorted_transcripts[0] );
  push( @clusters, $cluster );
    
  # keep track of the edges of the cluster, useful for a negative check
  my %start;
  my %end;
  $sorted_transcripts[0]->sort;
  $start{ $cluster } = $sorted_transcripts[0]->start_Exon->start;
  $end{ $cluster }   = $sorted_transcripts[0]->end_Exon->end;

  ## just a test to see whether we can trust start_Exon() and end_Exon()
  #print STDERR "start from transcript: ".$start{ $cluster }."\n";
  #print STDERR "      from method    : ".$self->_get_start_of_Transcript( $sorted_transcripts[0] )."\n";
  #print STDERR "end from transcript  : ". $end{ $cluster }."\n";
  #print STDERR "    from method      : ".$self->_get_end_of_Transcript( $sorted_transcripts[0] )."\n";

  # loop over the rest of the genes
 LOOP1:
  for (my $c=1; $c<=$#sorted_transcripts; $c++){
    my $found=0;

    # first do a negative-check on this cluster
    $sorted_transcripts[$c]->sort;
    my $this_start = $sorted_transcripts[$c]->start_Exon->start;
    my $this_end   = $sorted_transcripts[$c]->end_Exon->end;
    
    # only look if they potentially overlap
    #print STDERR "1:comparing transcript ".$sorted_transcripts[$c]->dbID." [ $this_start , $this_end ] with cluster $cluster [ $start{$cluster} , $end{ $cluster } ]\n";
    if ( !( $this_start > $end{ $cluster } || $this_end < $start{ $cluster } ) ){
      #print STDERR "  inside!\n";

      # compare with the transcripts in this cluster
    LOOP2:
      foreach my $t_in_cluster ( $cluster->get_Transcripts ){       
	if ( $self->_compare_Transcripts( $sorted_transcripts[$c], $t_in_cluster ) ){	
	  $cluster->put_Transcripts( $sorted_transcripts[$c] );                       
	  $found=1;
	  
	  # reset start/end if necessary
	  if ( $this_start < $start{$cluster} ){
	    $start{ $cluster } = $this_start;
	  }
	  if ( $this_end   > $end{ $cluster }  ){
	    $end{ $cluster } = $this_end;
	  }
	  next LOOP1;
	}
      }
    }
    # if not in this cluster compare to the previous clusters:

    # to restrict this to the ($limit) previous clusters
    # set my $limit = 6; (for example) and include in the while the following condition
    # while ( !(...)  && !($lookup > $limit) )

    #print STDERR "  found = $found\n";
    if ( $found == 0 && $cluster_count > 1 ) {
      my $lookup = 1;
      
      # loop through the clusters backwards
      while ( !($cluster_count <= $lookup ) ){ 
	#print STDERR "cluster_count: $cluster_count, looking at ".($cluster_count - $lookup)."\n";
	my $previous_cluster = $clusters[ $cluster_count - 1 - $lookup ];
	
	# only look if it is potentially overlapping
	#print STDERR "2:comparing transcript ".$sorted_transcripts[$c]->dbID." [ $this_start , $this_end ] with cluster $previous_cluster [ $start{$previous_cluster} , $end{ $previous_cluster } ]\n";
	if ( !(  $this_start > $end{ $previous_cluster } || $this_end < $start{ $previous_cluster } ) ){
	  #print STDERR "  inside!\n";
	  # loop over the transcripts in this previous cluster
	  foreach my $t_in_cluster ( $previous_cluster->get_Transcripts ){
	    if ( $self->_compare_Transcripts( $sorted_transcripts[$c], $t_in_cluster ) ){	
	      $previous_cluster->put_Transcripts( $sorted_transcripts[$c] );                       
	      $found=1;
	      
	      # reset start/end if necessary
	      if ( $this_start < $start{ $previous_cluster} ){
		$start{ $previous_cluster } = $this_start;
	      }
	      if ( $this_end   > $end{ $previous_cluster }  ){
		$end{ $previous_cluster } = $this_end;
	      }
	      next LOOP1;
	    }
	  }
	}
	$lookup++;
      }
    }
    # if not-clustered create a new TranscriptCluster
    #print STDERR "  found = $found\n";
    if ( $found == 0 ){  
      $cluster = new Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster; 
      $cluster->put_Transcripts( $sorted_transcripts[$c] );
      $start{ $cluster } = $sorted_transcripts[$c]->start_Exon->start;
      $end{ $cluster }   = $sorted_transcripts[$c]->end_Exon->end;
      #print STDERR "  creating a new cluster $cluster [ $start{ $cluster }, $end{ $cluster } ]\n";
      push( @clusters, $cluster );
      $cluster_count++;
    }
  }

#  ## print out the clusters
#  my $number  = 1;
#  foreach my $cluster (@clusters){
#    my $count = 1;
#    print STDERR "cluster $number :\n";
#    foreach my $tran ($cluster->get_Transcripts){
#      print STDERR "$count:\n";
#      print STDERR $tran->dbID." >";
#      foreach my $exon ( $tran->get_all_Exons ){
#  	print STDERR $exon->start.":".$exon->end." ";
#      }
#      print STDERR "\n";
#      $count++;
#    }
#    $number++;
#  }		
  
  return @clusters;
}
############################################################

sub _get_start_of_Transcript{        
  my ($self,$transcript) = @_;
  my @exons = @{$transcript->get_all_Exons};
  my @sorted_exons = sort { $a->start <=> $b->start } @exons;
  my $start = $sorted_exons[0]->start;
  return $start;
}    

############################################################
sub _get_end_of_Transcript {        
  my ($self,$transcript) = @_;
  my @exons = @{$transcript->get_all_Exons};
  my $end = 0;
  my $this_end;
  foreach my $exon (@exons){
   $this_end = $exon->end;
   if ( $this_end > $end ){
     $end = $this_end;
   }
  }
  return $this_end;
}    

############################################################

=head2 _compare_Transcripts()

 Title: _compare_Transcripts()
 Usage: compares the exons of two transcripts according to overlap and returns 1 if they have at least
        one exon overlap, and 0 otherwise

=cut

sub _compare_Transcripts {        
  my ($self,$transcript1,$transcript2) = @_;
  my @exons1   = @{$transcript1->get_all_Exons};
  my @exons2   = @{$transcript2->get_all_Exons};
  my $overlaps = 0;
  
  foreach my $exon1 (@exons1){
    foreach my $exon2 (@exons2){
      if ( ($exon1->overlaps($exon2)) && ($exon1->strand == $exon2->strand) ){
	return 1;
      }
    }
  }
  return 0;
}    

#########################################################################

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

# version copied from BlatToGenes, needs work on the reverse strand

=head2 check_splice_sites

We want introns of the form:
    
    ...###GT...AG###...   ...###AT...AC###...   ...###GC...AG###...
    
=cut

sub check_splice_sites{
  my ($self,$transcript,$strand) = @_;
 
  my $verbose = 1;
  
  my @exons = sort { $a->start <=> $b->end } @{$transcript->get_all_Exons};
  
  # no need to check single-exon ones
  my $introns  = scalar(@exons) - 1 ; 
  if ( $introns <= 0 ){
    return 1;
  }

  if ($verbose){
   print STDERR "checking splice sites in transcript:\n";
   Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($transcript);
  #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_TranscriptEvidence($transcript);
  }
  
  my $correct  = 0;
  my $wrong    = 0;
  my $other    = 0;
  
  # all exons in the transcripts are in the same seqname coordinate system:
  my $slice = $transcript->start_Exon->contig;
  
  if ($strand == 1 ){
    
  INTRON:
      for (my $i=0; $i<$#exons; $i++ ){
	  my $upstream_exon   = $exons[$i];
	  my $downstream_exon = $exons[$i+1];
	  
	  #my $upstream_start   = $upstream_exon->end + 1;
	  #my $upstream_end     = $upstream_exon->end + 2;
	  #my $downstream_start = $downstream_exon->start - 2;
	  #my $downstream_end   = $downstream_exon->start - 1;
	  
	  #my $up_start   = $slice->chr_start + $upstream_start - 1;
	  #my $up_end     = $slice->chr_start + $upstream_end - 1;
	  #my $down_start = $slice->chr_start + $downstream_start - 1;
	  #my $down_end   = $slice->chr_start + $downstream_end - 1;
	  #my $upstream_site   = 
	  #  $self->get_chr_subseq($slice->chr_name, $up_start, $up_end, $strand );
	  #my $downstream_site = 
	  #  $self->get_chr_subseq($slice->chr_name, $down_start, $down_end, $strand );
	  

	  my $upstream_site;
	  my $downstream_site;
	  eval{
	    $upstream_site = 
		$slice->subseq( ($upstream_exon->end     + 1), ($upstream_exon->end     + 2 ) );
            $downstream_site = 
		$slice->subseq( ($downstream_exon->start - 2), ($downstream_exon->start - 1 ) );

	    #print STDERR "forward; upstream: $upstream_site, downstream:$downstream_site\n" if $verbose;
	  };
	  unless ( $upstream_site && $downstream_site ){
	      print STDERR "problems retrieving sequence for splice sites\n$@";
	      next INTRON;
	  }
	  
	  print STDERR "upstream $upstream_site, downstream: $downstream_site\n" if $verbose;
	  ## good pairs of upstream-downstream intron sites:
	  ## ..###GT...AG###...   ...###AT...AC###...   ...###GC...AG###.
	  # AT-AC is the U12 spliceosome
	  if (  ($upstream_site eq 'GT' && $downstream_site eq 'AG') ||
		($upstream_site eq 'AT' && $downstream_site eq 'AC') ||
		($upstream_site eq 'GC' && $downstream_site eq 'AG') ){
	      $correct++;
	  }
	  else{
	      $other++;
	  }
      } # end of INTRON
  }
  elsif ( $strand == -1 ){
      
    #  example:
    #                                     --------CT...AC-------... 
    #  transcript in reverse strand -> 3' --######GA...TG#####--... 5' 
    #
    #  the slice has been reverse complemented so that the exon looks like:
    #  transcript now in the forward-> 5' --#####CA...TC######--...3' 
    #                                    --------GT...AG-------... 
    
    # we see CA-TC in the slice and we need 
    # to find the complementary sequence to look for a good site

    #my $chr_length = $slice->adaptor->db->get_ChromosomeAdaptor->fetch_by_chr_name( $slice->chr_name )->length;
    
  INTRON:
    for (my $i=0; $i<$#exons; $i++ ){
      my $upstream_exon   = $exons[$i];
      my $downstream_exon = $exons[$i+1];
      
      #my $upstream_start   = $upstream_exon->end + 1;
      #my $upstream_end     = $upstream_exon->end + 2;
      #my $downstream_start = $downstream_exon->start - 2;
      #my $downstream_end   = $downstream_exon->start - 1;
      
     	  my $upstream_site;
	  my $downstream_site;
      eval{
	$upstream_site = 
	  $slice->subseq( ($upstream_exon->end     + 1), ($upstream_exon->end     + 2 ) );
	$downstream_site = 
	  $slice->subseq( ($downstream_exon->start - 2), ($downstream_exon->start - 1 ) );
      };
      
      #my $up_start   = $chr_length  - ($slice->chr_start + $upstream_start - 1)   + 1;
      #my $up_end     = $chr_length  - ($slice->chr_start + $upstream_end - 1  )   + 1;
      #my $down_start = $chr_length  - ($slice->chr_start + $downstream_start - 1) + 1;
      #my $down_end   = $chr_length  - ($slice->chr_start + $downstream_end - 1)   + 1;
      #my $upstream_site   = 
#	$self->get_chr_subseq($slice->chr_name, $up_start, $up_end, $strand );
#      my $downstream_site = 
#	$self->get_chr_subseq($slice->chr_name, $down_start, $down_end, $strand );
#            
      unless ( $upstream_site && $downstream_site ){
	print STDERR "problems retrieving sequence for splice sites\n$@";
	next INTRON;
      }
      #$upstream_site   =~ tr/ACGTacgt/TGCAtgca/;
      #$downstream_site =~ tr/ACGTacgt/TGCAtgca/;
      
      print STDERR "reverse; upstream $upstream_site, downstream: $downstream_site\n" if $verbose;
      if (  ($upstream_site eq 'GT' && $downstream_site eq 'AG') ||
	    ($upstream_site eq 'AT' && $downstream_site eq 'AC') ||
	    ($upstream_site eq 'GC' && $downstream_site eq 'AG') ){
	  $correct++;
      }
      else{
	    $other++;
      }
      
    } # end of INTRON
  }
  unless ( $introns == $other + $correct ){
    print STDERR "STRANGE: introns:  $introns, correct: $correct, wrong: $wrong, other: $other\n";
    return 0;
  }
  if ( $other ){
    print STDERR "rejecting for having non-canonical splice-sites\n" if $verbose;
    return  0;
  }
  else{
    #print STDERR "accepting\n" if $verbose;
    return 1;
  }
}
############################################################

=head2 get_chr_subseq

It return a piece of chromosome sequence specified
by start, end and strand. Its purpose is to
check the splice site sequences and it needs to know
where the dumps of the chromosomes are (fasta files),
which reads from the variable EST_GENOMIC in EST_GeneBuilder_conf.pm
strand is not used.

=cut

sub get_chr_subseq{
  my ( $self, $chr_name, $start, $end, $strand ) = @_;

  my $chr_file = $EST_GENOMIC."/".$chr_name.".fa";
  my $command = "chr_subseq $chr_file $start $end |";
  #print STDERR "command: $command\n";
  open( SEQ, $command ) || $self->throw("Error running chr_subseq within ExonerateToGenes");
  my $seq = uc <SEQ>;
  chomp $seq;
  close( SEQ );

  return $seq;
}

############################################################

=head2 _put_Transcript

 Title   : _put_Transcript
 Usage   : $self->add_Transcript
 Function: method to add transcripts into the array $self->{'_transcripts'} 
 Returns : nothing
 Args    : Bio::EnsEMBL::Transcript

=cut

sub _put_Transcript {
  my ($self,$transcript) = @_;
  $self->throw("No transcript input") unless defined($transcript);
  $self->throw("Input must be Bio::EnsEMBL::Transcript") unless $transcript->isa("Bio::EnsEMBL::Transcript");
  if ( !defined( $self->{'_transcripts'} ) ){
    @{ $self->{'_transcripts'} } = ();
  }
  push( @{ $self->{'_transcripts'} }, $transcript );
}

############################################################

=head2 _get_all_Transcripts

 Title   : _get_all_Transcripts
 Usage   : my @transcripts = $self->_get_all_Transcripts;
 Function: method to get all the transcripts stored in @{ $self->{'_transcripts'} } 
 Example : 
 Returns : Bio::EnsEMBL::Transcript
 Args    : nothing

=cut

sub _get_all_Transcripts {
  my ($self) = @_;
  if ( !defined( $self->{'_transcripts'} ) ) {
    @{ $self->{'_transcripts'} } = ();
    print STDERR "The transcript array you're trying to get is empty\n";
  }
  my @trans = @{ $self->{'_transcripts'} };
  return @trans;
}

############################################################

sub add_runnable{
  my ($self, $value, $strand) = @_;

  if (!defined($self->{'_forward_runnables'})) {
    $self->{'_forward_runnables'} = [];
  }
  if (!defined($self->{'_reverse_runnables'})) {
    $self->{'_reverse_runnables'} = [];
  } 
  if (defined($value)) {
    
    if ($value->isa("Bio::EnsEMBL::Pipeline::RunnableI")) {
      
      if( $strand == -1 ){
	push(@{$self->{'_reverse_runnables'}},$value);
      }
      elsif( $strand == 1){
	push(@{$self->{'_forward_runnables'}},$value);
      }
      else{
	$self->throw( "Cannot add a runnable with strand = $strand" );
      }
    
    } 
    else {
      $self->throw("[$value] is not a Bio::EnsEMBL::Pipeline::RunnableI");
    }
  }
}

############################################################

sub each_runnable{
  my ($self, $strand) = @_;
  
  if (!defined($self->{'_forward_runnables'})) {
    $self->{'_forward_runnables'} = [];
  }
  
  if (!defined($self->{'_reverse_runnables'})) {
    $self->{'_reverse_runnables'} = [];
  } 

  if( $strand == -1 ){
    return @{$self->{'_reverse_runnables'}};
  }
  elsif ($strand == 1){
    return @{$self->{'_forward_runnables'}};
  }
  else{
    $self->throw( "there are no runnables with strand = $strand" );
  }
  
}

############################################################

# run genomewise 

sub run {
  my ($self) = @_;
  my $strand;
    
  # sort out analysis here or we will get into trouble with duplicate analyses
  my $analysis = $self->analysis;
  
  unless ( $analysis ){
      $self->throw("You need an analysis to run this");
  }

  print STDERR "Analysis is: ".$self->analysis->dbID."\n";
  ############################################################
  #### forward strand 
  ############################################################
  $strand = 1;
  my $tcount=0;

  my @f_transcripts1 = $self->_forward_transcripts;
  my @f_transcripts2  = $self->_process_Transcripts(\@f_transcripts1,$strand);
  my @forward_transcripts;
  

  ############################################################
  # run genomewise?
  if ( $USE_GENOMEWISE ){
    foreach my $tran (@f_transcripts2){
      my $runnable = 
	new Bio::EnsEMBL::Pipeline::Runnable::MiniGenomewise(
							     -genomic  => $self->query,
							     -analysis => $self->analysis,
							     -smell    => $GENOMEWISE_SMELL,
							    );
      $self->add_runnable($runnable,$strand);
      $runnable->add_Transcript($tran);
    }
    
  RUN1:
    foreach my $gw_runnable( $self->each_runnable($strand) ){
      $tcount++;
      eval{
	$gw_runnable->run;
      };
      if ($@){
	print STDERR "Your runnable $gw_runnable didn't succeed!\n";
	print STDERR $@;
	next RUN1;
      }
      foreach my $transcript ( $gw_runnable->output ){
	#$self->_set_evidence( $transcript );
	push (@forward_transcripts, $transcript);
      }
    }
    print STDERR $tcount." transcripts run in genomewise in the forward strand\n";

  }
  ############################################################
  # else use a simpler method
  else{
    foreach my $tran (@f_transcripts2){
      my $new_tran = Bio::EnsEMBL::Pipeline::Tools::TranslationUtils->compute_translation( $tran );
      if ( $new_tran ){
	push ( @forward_transcripts, $new_tran );
      }
    }
  }
  
  my @checked_forward_transcripts   = $self->_check_Translations(\@forward_transcripts,$strand);
  
  # cluster them into genes
  my @forward_genes                 = $self->_cluster_into_Genes(@checked_forward_transcripts);
  
  # take the best ten transcripts from each gene
  my @selected_forward_transcripts  = $self->_select_best_transcripts(@forward_genes);

  # cluster again into genes
  my @selected_forward_genes        = $self->_cluster_into_Genes( @selected_forward_transcripts );

  # important: make shared exons unique
  my @ready_forward_genes           = $self->_make_shared_exons_unique( @selected_forward_genes );

  # map them to raw contig coordinates
  my @forward_remapped              = $self->remap_genes(\@ready_forward_genes, $strand);

  $self->output(@forward_remapped);


  ############################################################
  #### reverse strand 
  ############################################################
  $strand = -1;
  my $tcount2=0;
  
  my @reverse_transcripts;
  my @r_transcripts1 = $self->_reverse_transcripts;
  my @r_transcripts2 = $self->_process_Transcripts(\@r_transcripts1,$strand);  
  
  ############################################################
  # use genomewise?
  if ( $USE_GENOMEWISE ){
    foreach my $tran (@r_transcripts2) {
      my $runnable = 
	new Bio::EnsEMBL::Pipeline::Runnable::MiniGenomewise(
							     -genomic  => $self->revcomp_query,
							     -analysis => $self->analysis,
							     -smell    => $GENOMEWISE_SMELL,
							    );
	    $self->add_runnable($runnable, $strand);
	    $runnable->add_Transcript($tran);
    }
    
  RUN2:
    foreach my $gw_runnable( $self->each_runnable($strand)) {
      $tcount2++;
      eval{
	$gw_runnable->run;
      };
      if ($@){
	print STDERR "Your runnable $gw_runnable didn't succeed!\n";
	print STDERR $@;
	next RUN2;
      }
      foreach my $transcript ( $gw_runnable->output ){
	#$self->_set_evidence( $transcript );
	push (@reverse_transcripts, $transcript);
      }
    }
    print STDERR $tcount2." transcripts run in genomewise in the reverse strand\n";
    
  }
  ############################################################
  # else use a simpler method
  else{
    foreach my $tran (@r_transcripts2){
      my $new_tran = Bio::EnsEMBL::Pipeline::Tools::TranslationUtils->compute_translation( $tran );
      if ( $new_tran ){
	push ( @reverse_transcripts, $new_tran );
      }
    }
  }
  
  my @checked_reverse_transcripts = $self->_check_Translations(\@reverse_transcripts,$strand);
  
  # cluster them into genes
  my @reverse_genes               = $self->_cluster_into_Genes(@checked_reverse_transcripts);
  
  # take the best ten transcripts from each gene
  my @selected_reverse_transcripts = $self->_select_best_transcripts(@reverse_genes);
  
  # cluster them again into genes
  my @selected_reverse_genes       = $self->_cluster_into_Genes( @selected_reverse_transcripts );

  # important: make shared exons unique
  my @ready_reverse_genes          = $self->_make_shared_exons_unique( @selected_reverse_genes );
  
  # map them to raw contig coordinates
  my @reverse_remapped             = $self->remap_genes(\@ready_reverse_genes, $strand);
  
  $self->output(@reverse_remapped);

  foreach my $gene ( $self->output ){
    print STDERR "gene with ".scalar(@{$gene->get_all_Transcripts})." transcripts\n";
  }

}

############################################################

sub _select_best_transcripts{
  my ( $self, @genes ) = @_;
  my @selected_transcripts;
 GENE:
  foreach my $gene ( @genes ){
    my $count = 0;
    my @trans = sort { my $result = ( $self->_mean_coverage($b) <=> $self->_mean_coverage($a) );
		       if ( $result){
			 return $result;
		       }
		       else{
			 return ( scalar(@{$b->get_all_Exons}) <=> scalar(@{$a->get_all_Exons}) );
		       }
		     }  @{$gene->get_all_Transcripts};
    
    print STDERR " $gene has ".scalar(@trans)." transcripts\n";
  TRAN:
    foreach my $tran ( @trans ){
      $count++;
      print STDERR "count: $count\n";
      next GENE if $count > $MAX_TRANSCRIPTS_PER_GENE;
      push ( @selected_transcripts, $tran );
      print STDERR "transcript accepted\n";
    }
  }
  return @selected_transcripts;
}

############################################################

sub _mean_coverage{
  my ($self,$tran) = @_;
  my %score;
 EXON:
  foreach my $exon (@{$tran->get_all_Exons} ){
  EVI:
    foreach my $evi ( @{$exon->get_all_supporting_features} ){
      $score{ $evi->hseqname} = $evi->score;
    }
  }
  my $mean_coverage;
  my @ids =  keys %score;
  foreach my $id ( @ids  ){
    $mean_coverage += $score{$id};
  }
  #print STDERR "Mean coverage = ".($mean_coverage/scalar(@ids))."\n";
  return $mean_coverage/scalar(@ids);
}

###################################################################

sub genetype {
  my ($self, $genetype) = @_;

  if(defined $genetype){
    $self->{_genetype} = $genetype;
  }

  return $self->{_genetype};
}

############################################################

sub analysis {
  my ($self, $analysis) = @_;

  if(defined $analysis){
    $self->throw("$analysis is not a Bio::EnsEMBL::Analysis") unless $analysis->isa("Bio::EnsEMBL::Analysis");
    $self->{_analysis} = $analysis;
  }
  return $self->{_analysis};
}

############################################################
# this method set the slice in the exons for each transcript
# and check the translation, and adds the next codon at the end
# codon if it is inside the transcript and it is a stop codon: taa/tag/tga 

sub _check_Translations {
  my ($self,$transcripts,$strand) = @_;
  
  my @trans = @$transcripts;
  my @good_transcripts;

  my $slice = $self->query;
  # are we working on the reverse strand?
  if( $strand == -1 ){
    $slice = $self->revcomp_query;
  }
  
 TRANSCRIPT:
  foreach my $transcript (@trans) {

    # sort the exons 
    $transcript->sort;
    my @exons = @{$transcript->get_all_Exons};

    # at this point, if necessary, only accepts transcripts with more than one exon
    # although we have checked this already, genomewise sometimes bridges 
    # over introns making one exon out of two
    if ( $REJECT_SINGLE_EXON_TRANSCRIPTS && scalar(@exons) == 1 ){
      print STDERR "Rejected a single-exon transcript\n";
      next TRANSCRIPT;
    }
    
  EXON:
    foreach my $exon(@exons){
      $exon->contig($slice);
      $exon->attach_seq($slice);
      
      # if strand = -1 we have inverted the contig, thus
      $exon->strand(1);
 
      # when the gene gets stored, the strand is flipped automatically
      print STDERR "Exon ".$exon->start."-".$exon->end." phase: ".$exon->phase." end phase: ".$exon->end_phase."\n";
    }

    my $translation = $transcript->translation;
    
    # store only genes that translate ( to check it, we get the Bio::Seq )
    my $sequence = $transcript->translate;
    
    unless ( $sequence ){
      print STDERR "TRANSCRIPT WITHOUT A TRANSLATION!!\n";
    }
    if ( $sequence ){
      # length of the translation
      my $length   = $sequence->length;
      
      # total length of the exons in the transcript
      my $t_length = $transcript-> length;
      print STDERR "Translation length : $length\n";
      print STDERR "Exons length       : $t_length\n";
      
      # 5' UTR is usually shorter than 3' UTR, the latter can be very long compared
      # with the translation length ( 5000 vs. 500 ) see e.g. gene SCL ( aka TAL1)
      my $five_prime  = $transcript->five_prime_utr  or print STDERR "No five prime UTR";
      my $three_prime = $transcript->three_prime_utr or print STDERR "No three prime UTR";
      
      # UTRs above are Bio::Seq
      if ( $five_prime ){
	my $length5    = $five_prime->length;
	print STDERR "5prime UTR length  : $length5\n";
      }
      if ( $three_prime ){
	my $length3    = $three_prime->length;
	print STDERR "3prime UTR length  : $length3\n";
      }
      
      # only store the genes whose translation has no stop codons
      my $peptide = $sequence->seq;
      #print STDERR "peptide: $peptide\n";
      if ( $peptide =~ /\*/ ){
	print STDERR "TRANSLATION HAS STOP CODONS!! - skipping\n";
	next TRANSCRIPT
      }
      
      ############################################################
      # put possible stop at the end:
      eval{
	$transcript = Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->set_stop_codon( $transcript );
      };
      if($@){
	print STDERR "there was a problem with the trancript: [$@]\n"
      }
      else{
	push(@good_transcripts, $transcript);
      }
    }
  } # end of TRANSCRIPT
  
  return @good_transcripts;
  
}


###################################################################

=head2 remap_genes

    Usage   :   $self->remap_genes(@genes)
    Function:   Remaps predicted genes into genomic coordinates
    Returns :   array of Bio::EnsEMBL::Gene
    Args    :   arrayref of Bio::EnsEMBL::Gene, INT (strand)

=cut

sub remap_genes {
    my ($self, $genes, $strand) = @_;
    #  my $slice = $self->vcontig;
    #  if ( $strand == -1 ){
    #    $slice = $slice->invert;
    #  }
    
    my @newf;
    my $trancount=1;
    foreach my $gene (@$genes) {
	my @trans = @{$gene->get_all_Transcripts};
	my $newgene;
	# convert to raw contig coords
	eval {
	    # transforming gene to raw contig coordinates.
	    $newgene = $gene->transform;
	};
	if ($@) {
	    print STDERR "Couldn't reverse map gene [$@]\n";
	    foreach my $t ( @{$gene->get_all_Transcripts} ){
		$self->_print_Transcript($t);
	    }
	    next;
	}
	$newgene->type($gene->type);
	$newgene->analysis($gene->analysis);
	push(@newf,$newgene);
	
	
    }
    return @newf;
}

############################################################

sub output{
  my ($self,@genes) = @_;
  
  if (!defined($self->{_output})) {
    $self->{_output} = [];
  }
  
  if(@genes){
    push(@{$self->{_output}},@genes);
  }
  
   return @{$self->{_output}};
}

#########################################################################
  
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
    print $exon->start."-".$exon->end."[".$exon->phase.",".$exon->end_phase."] ";
  }
  print STDERR "\n";
  print STDERR "translation start exon: ".
    $transcript->translation->start_Exon->start."-".$transcript->translation->start_Exon->end.
      " start: ".$transcript->translation->start."\n";
  print STDERR "translation end exon: ".
    $transcript->translation->end_Exon->start."-".$transcript->translation->end_Exon->end.
      " end: ".$transcript->translation->end."\n";
}


############################################################
#
# METHODS FOR CLUSTERING THE TRANSCRIPTS INTO GENES
#
############################################################
#
# similar but not the same as the ones in GeneBuilder
#

=head2 _cluster_into_Genes

    Example :   my @genes = $self->cluster_into_Genes(@transcripts);
Description :   it clusters transcripts into genes according to exon overlap.
                It will take care of difficult cases like transcripts within introns.
                It also unify exons that are shared among transcripts.
    Returns :   a beautiful list of geen objects
    Args    :   a list of transcript objects

=cut

sub _cluster_into_Genes{
  my ($self, @transcripts_unsorted) = @_;
  
  my $num_trans = scalar(@transcripts_unsorted);
  print STDERR "clustering $num_trans transcripts into genes\n";

  #foreach my $t ( @transcripts_unsorted ){
  #  Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($t);
  #}
  # flusold genes
  #$self->flush_Genes;
  
  my @transcripts = sort { my $result = ( $self->transcript_low($a) <=> $self->transcript_low($b) );
			   if ($result){
			     return $result;
			   }
			   else{
			     return ( $self->transcript_high($b) <=> $self->transcript_high($a) );
			   }
			 } @transcripts_unsorted;
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
  
  print STDERR scalar(@clusters)." clusters created, making genes...\n";
  my @genes;
  foreach my $cluster(@clusters){
    my $count = 0;
    my $gene = new Bio::EnsEMBL::Gene;
    my $genetype = $self->genetype;
    my $analysis = $self->analysis;
    $gene->type($genetype);
    $gene->analysis($analysis);
    
    foreach my $transcript( @{$cluster} ){
      $gene->add_Transcript($transcript);
    }
    push ( @genes, $gene);
  }
  
  # test
  #foreach my $gene (@genes){
  #  foreach my $tran ( @{$gene->get_all_Transcripts} ){
  #    foreach my $exon ( @{$tran->get_all_Exons} ){
  #	foreach my $evi ( @{$exon->get_all_supporting_features} ){
  #	  print STDERR "analysis= ".$evi->analysis->logic_name."\n";
  #	}
  #    }
  #  }
  #}

  return @genes;
}

############################################################

sub _make_shared_exons_unique{
  my ( $self, @genes ) = @_;
  my @pruned_genes;
  foreach my $gene ( @genes ){
 
    # make different exon objects that are shared between transcripts 
    # ( regarding attributes: start, end, etc )
    # into unique exon objects 
    my $new_gene = $self->prune_Exons($gene);
    push( @pruned_genes, $new_gene );
  }
  return @pruned_genes;
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

sub transcript_high{
  my ($self,$tran) = @_;
  my $high;
  if ( $tran->start_Exon->strand == 1){
    $high = $tran->end_Exon->end;
  }
  else{
    $high = $tran->start_Exon->end;
  }
  return $high;
}

############################################################

sub transcript_low{
  my ($self,$tran) = @_;
  my $low;
  if ( $tran->start_Exon->strand == 1){
    $low = $tran->start_Exon->start;
  }
  else{
    $low = $tran->end_Exon->start;
  }
  return $low;
}

############################################################

sub by_transcript_high {
  my $alow;
  my $blow;

  my $ahigh;
  my $bhigh;
  
  # alow and ahigh are the left most and right most coordinates for transcript $a 
  if ($a->start_Exon->strand == 1) {
    $alow  = $a->start_Exon->start;
    $ahigh = $a->end_Exon->end;
  } 
  else {
    $alow  = $a->end_Exon->start;
    $ahigh = $a->start_Exon->end;
  }

  # blow and bhigh are the left most and right most coordinates for transcript $b 
  if ($b->start_Exon->strand == 1) {
    $blow  = $b->start_Exon->start;
    $bhigh = $b->end_Exon->end;
  } 
  else {
    $blow  = $b->end_Exon->start;
    $bhigh = $b->start_Exon->end;
  }

  # return the ascending comparison of the right-most coordinates if they're different
  if ($ahigh != $bhigh) {
    return $ahigh <=> $bhigh;
  } 
  # if they'r equal, return the ascending comparison of the left most coordinate
  else {
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
    UNI:
      foreach my $uni (@unique_Exons) {
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
      } 
      else {
	push(@newexons,$exon);
	push(@unique_Exons, $exon);
      }
    }          
    $tran->flush_Exons;
    foreach my $exon (@newexons) {
      $tran->add_Exon($exon);
    }
  }
  # we return the same gene object, but with modified exons in the transcripts
  return $gene;
}

############################################################


1;


