#
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
use Bio::EnsEMBL::Pipeline::Runnable::ClusterMerge;
use Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster;

# config file; parameters searched for here if not passed in as @args
use Bio::EnsEMBL::Pipeline::ESTConf qw (
					EST_INPUTID_REGEX
					EST_REFDBHOST
					EST_REFDBUSER
					EST_REFDBNAME
					EST_REFDBPASS
					EST_E2G_DBNAME
					EST_E2G_DBHOST
					EST_E2G_DBUSER
					EST_E2G_DBPASS     
					EST_GENEBUILDER_INPUT_GENETYPE
					EST_EVIDENCE_TAG
					EST_MIN_EVIDENCE_SIMILARITY
					EST_MAX_EVIDENCE_DISCONTINUITY
					EST_GENOMEWISE_GENETYPE
					USE_cDNA_DB
					cDNA_DBNAME
					cDNA_DBHOST
					cDNA_DBUSER
					cDNA_DBPASS
					cDNA_GENETYPE
				       );

# use new Adaptor to get some extra info from the ESTs
#use Bio::EnsEMBL::Pipeline::DBSQL::ESTFeatureAdaptor;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);

    ## db input_id mandatory and read in by BlastableDB
    #if (!defined $self->seqfetcher) {
    #  my $seqfetcher = $self->make_seqfetcher();
    #  $self->seqfetcher($seqfetcher);
    #}
    # db needs a reference dna database
    my $refdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
						   -host             => $EST_REFDBHOST,
						   -user             => $EST_REFDBUSER,
						   -dbname           => $EST_REFDBNAME,
						   -pass             => $EST_REFDBPASS,
						  );
    
    my $est_e2g_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
							-host             => $EST_E2G_DBHOST,
							-user             => $EST_E2G_DBUSER,
							-dbname           => $EST_E2G_DBNAME,
							-pass             => $EST_E2G_DBPASS,
						       );
    
    $est_e2g_db->dnadb($refdb);
    $self->est_e2g_db($est_e2g_db);
    
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

    $self->genetype($EST_GENOMEWISE_GENETYPE);



    return $self; 
}

############################################################

sub est_e2g_db{
    my ($self, $est_e2g_db) = @_;
    if ($est_e2g_db){
	$self->{_est_e2g_db} = $est_e2g_db;
    }
    return $self->{_est_e2g_db};
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

    Title   :   write_output
    Usage   :   $self->write_output
    Function:   Writes output data to db
    Returns :   array of exons (with start and end)
    Args    :   none

=cut

sub write_output {
  my ($self) = @_;
    
  my $gene_adaptor = $self->db->get_GeneAdaptor;
    
 GENE: 
  foreach my $gene ($self->output) {	
    eval {
      $gene_adaptor->store($gene);
      print STDERR "wrote gene " . $gene->dbID . "\n";
      
      my @transcripts = @{ $gene->get_all_Transcripts};
      foreach my $tran (@transcripts){
	Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Evidence($tran);
      }
      
      
    }; 
    if( $@ ) {
      print STDERR "UNABLE TO WRITE GENE\n\n$@\n\nSkipping this gene\n";
    }
  }
}

############################################################

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data from the database
    Returns :   nothing
    Args    :    string: chr1.1-10000

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

    my $slice = $self->db->get_SliceAdaptor->fetch_by_chr_start_end($chrname,$chrstart,$chrend);    
    $slice->chr_name($chrname);
    $self->query($slice);
    
    ############################################################
    # forward strand
    ############################################################

    $strand = 1;
    print STDERR "\n****** forward strand ******\n\n";

    # get genes
    my $genes  = $slice->get_all_Genes_by_type($genetype);
    
    print STDERR "Number of genes from ests  = " . scalar(@$genes) . "\n";
    
    my $cdna_slice;
    if ( $USE_cDNA_DB ){
	my $cdna_db = $self->cdna_db;
	
	$cdna_slice = $cdna_db->get_SliceAdaptor->fetch_by_chr_start_end($chrname,$chrstart,$chrend);
	my $cdna_genes  = $cdna_slice->get_all_Genes_by_type($cDNA_GENETYPE);
	print STDERR "Number of genes from cdnas = " . scalar(@$cdna_genes) . "\n";
	push (@$genes, @$cdna_genes);
	
    }

    my @plus_transcripts;
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
		push (@plus_transcripts, $transcript );
	    }
	}
    }
    
    print STDERR "In EST_GeneBuilder.fetch_input(): ".scalar(@plus_transcripts) . " forward strand genes\n";
    print STDERR "($single single exon genes NOT thrown away)\n";

    # process transcripts in the forward strand
    
    if( scalar(@plus_transcripts) ){

      my @transcripts  = $self->_process_Transcripts(\@plus_transcripts,$strand);
      
      # make a genomewise runnable for each cluster of transcripts
      foreach my $tran (@transcripts){
	  
	  # use MiniSeq
	  my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::MiniGenomewise(
									      -genomic  => $slice,
									      -analysis => $self->analysis,
									      );
	  
	  $self->add_runnable($runnable,$strand);
	  # we only have one transcript per runnable
	  $runnable->add_Transcript($tran);
      }
  }
    
    
    # minus strand - flip the vc and hope it copes ...
    # but SLOW - get the same genes twice ...
    
    print STDERR "\n****** reverse strand ******\n\n";

    $strand = -1;
    
    # this will return a slice which corresponds to the reversed complement of $slice:
    my $rev_slice = $slice->invert;
    $self->revcomp_query($rev_slice);
    my $revgenes  = $rev_slice->get_all_Genes_by_type($genetype);
    my @minus_transcripts;
    
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
		push (@minus_transcripts, $transcript);
	    }
	}
    }
    print STDERR "In EST_GeneBuilfer.fetch_input(): ".scalar(@minus_transcripts) . " reverse strand genes\n";
    print STDERR "($single single exon genes NOT thrown away)\n";
    
    if(scalar(@minus_transcripts)){
      
	my @transcripts = $self->_process_Transcripts(\@minus_transcripts,$strand);  
	
	foreach my $tran (@transcripts) {
	    
	    my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::MiniGenomewise(
									    -genomic  => $rev_slice,
										-analysis => $self->analysis,
										);
	    $self->add_runnable($runnable, $strand);
	    $runnable->add_Transcript($tran);
	  }
    }
}

############################################################

=head2 _flush_Transcripts

    Title   :   _flush_Transcripts
    Usage   :   $self->_flush_Transcripts
    Function:   it empties out the array of transcripts $self->{'_transcripts'}
    
=cut  

sub _flush_Transcripts {
  my ($self) = @_;
  $self->{'_transcripts'} = [];
  return;
}

############################################################

=head2 _process_Transcripts

    Title   :   _process_Transcripts
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
  # (GT-AG, AT-AC, GC-AG) or if the only intron they have is non standard.
  my @checked_transcripts = $self->check_splice_sites( \@transcripts , $strand );
  
  if ( scalar(@checked_transcripts) == 0 ){
      print STDERR "No transcripts left from the splice-site check\n";
      return;
  }
  my $merge_object 
      = Bio::EnsEMBL::Pipeline::Runnable::ClusterMerge->new(
							    -transcripts => \@checked_transcripts,
							    );
  
  $merge_object->run;
  my @merged_transcripts = $merge_object->output;
  
  # reject the single exon transcripts
  my @filtered_transcripts = @{$self->_reject_single_exon_Transcripts(@merged_transcripts)};
  print STDERR scalar(@filtered_transcripts)." transcripts left after rejecting single-exon transcripts\n";
  
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
    
    Usage   :   @transcripts = $self->_check_Transcripts(@transcripts)
    Function:   checks transcripts obtained from EST2Genome for consistency among exons
                in strand, hit_name (hid), exon content,
                and also checks that the hits associated to consecutive exons do not have a
                discontinuity (in hit-coordinates) larger than than a certain limit
    Returns :   @Bio::EnsEMBL::Transcript
    Args    :   @Bio::EnsEMBL::Transcript, ref to hash for linking hid to transcript, ref to hash for linking 
                exon to hid, 

=cut

sub _check_Transcripts {
  my ($self, $ref_transcripts, $strand) = @_;

  # the source_tag of the supporting evidence is specified in Bio/EnsEMBL/Pipeline/EST_conf.pl
  my $evidence_tag    = $EST_EVIDENCE_TAG;
  
  # the minimum allowed perc. identity of the evidence per exon 
  my $min_similarity  = $EST_MIN_EVIDENCE_SIMILARITY;

  # the maximum allowed discontinuity in EST hits
  my $max_est_gap     = $EST_MAX_EVIDENCE_DISCONTINUITY;

  print STDERR "EST_GeneBuilder: checking consistency of transcripts...\n";

  my @allexons;       # here we'll put all exons that pass the check
  my @alltranscripts; # here we'll put all the transcripts that pass the check
  my %hid_trans;
  my $exon_adaptor    = $self->db->get_ExonAdaptor;
  my $total_rejected        = 0;

  my $slice;
  if ( $strand == +1 ){
      $slice = $self->query;
  }
  else{
      $slice = $self->revcomp_query;
  }

  print STDERR "transcripts:\n";
  foreach my $t (@$ref_transcripts){
      print $t->dbID."\n";
  }


 TRANSCRIPT: 
  foreach my $transcript (@$ref_transcripts){
      
      # reject the transcripts that fall off the slice at the lower end
      unless ( Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_check_Transcript($transcript,$slice) ){
	  next TRANSCRIPT;
      }
      
      $transcript->sort;
      my $exons = $transcript->get_all_Exons;
      #print STDERR "Transcript with ".scalar(@$exons)." exons\n";
      my $hid;
      my $this_strand;
      my @accepted_exons; # here we hold the good exons in this transcript
      my $rejected = 0;
      my $exon_count = 0;
      my $seqname;
      
      my $previous_exon;
    EXON:
      foreach my $exon (@$exons){
	  $exon_count++;

	  my $hstart;
	  my $hend;
	  #print STDERR " --- Exon $exon_count ---\n";
	  # get the supporting_evidence for each exon
	  my @sf = sort { $a->hstart <=> $b->hstart } @{$exon->get_all_supporting_features};
	  
	  # check that you get suporting_evidence at all
	  if ( scalar( @sf ) == 0 ){
	      $self->warn("exon $exon with no supporting evidence, possible sticky exon, ".
			  "exon_id =".$exon->dbID." transcript_id = ".$transcript->dbID."\n");
	  }
	  
	  ####### check the gap with the evidence of the next exon
	  # if the ESTs are of good quality, this should not reject any
	  if ( $exon_count > 1 ){
	      my $est_gap = 0;
	      
	      my @previous_sf = sort { $a->hstart <=> $b->hstart } @{$previous_exon->get_all_supporting_features};
	      
	      if ( scalar( @previous_sf ) != 0 ){
		  
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
		  # test:
		  #print STDERR "EST evidence with gap: $est_gap\n";
		  #print STDERR "prev_exon: ".$previous_exon->start."-".$previous_exon->end.  "\texon : ".$exon->start."-".$exon->end."\n";
		  #print STDERR "prev_sf  : ".$previous_sf[0]->hstart."-".$previous_sf[$#previous_sf]->hend."\tsf[0]: ".$sf[0]->hstart."-".$sf[0]->hend."\n";
		  
		  # check the evidence gap between both exons, do not reject anything yet
		  if ( $est_gap > $max_est_gap ){
		      print STDERR "EST evidence with gap too large: $est_gap\n";
		      print STDERR "prev_exon: ".$previous_exon->start."-".$previous_exon->end.  "\texon : ".$exon->start."-".$exon->end."\n";
		      print STDERR "prev_sf  : ".$previous_sf[0]->hstart."-".$previous_sf[$#previous_sf]->hend."\tsf[0]: ".$sf[0]->hstart."-".$sf[0]->hend."\n";
		      # print STDERR "EST with too large gaps skipping it\n";
		      # next TRANSCRIPT;
		  }
	      }
	  }
     
	  # reject transcript with too large intron length
	  my $intron_length;
	  
	  # 100000 bases is quite tight, we better keep it low for ESTs
	  my $max_intron_length = 100000;
	  if ($exon_count > 1 ){
	      my ( $s, $e, $intron_length);
	      #if ($strand == 1){
		  $s             = $previous_exon->end;
		  $e             = $exon->start;
		  $intron_length = $e - $s - 1;
	      #}
	      #elsif ( $strand == -1 ){
	#	  $s             = $previous_exon->start;
	#	  $e             = $exon->end;
	#	  $intron_length = $s - $e - 1;
	#      }
	      #print STDERR "strand: $strand\n";
	      #print STDERR " $s  - $e  - 1 = $intron_length\n";
	      
	      if ( $intron_length > $max_intron_length ){
		  print STDERR "Rejecting transcript $transcript for having intron too long: $intron_length\n";
		  next TRANSCRIPT;
	      }
	  }
	  
	  
	  $previous_exon = $exon;



	  
      } # end of EXON


      # if the transcript made it to this point, keep it
      push (@alltranscripts, $transcript);
  }
  return @alltranscripts;
  
}

############################################################

=head2 _cluster_Transcripts

    Title   :   _cluster_Transcripts
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
  
  ## test
  #foreach my $tran (@transcripts){
  #  print STDERR "$tran, start: ".$tran->start_Exon->start."\n";
  #}
  
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

sub check_splice_sites{
    my ($self,$transcripts_ref,$strand) = @_;
    my @checked_transcripts;    

  TRANSCRIPT:
    foreach my $transcript ( @$transcripts_ref ){
	
	# all transcripts are in forward coordinates
	my @exons  = sort{ $a->start <=> $b->start } @{$transcript->get_all_Exons};
	
	#print STDERR "checking splice sites in transcript:\n";
	#Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($transcript);
	#Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_TranscriptEvidence($transcript);
	
	my $introns  = scalar(@exons) - 1 ; 
	if ( $introns <= 0 ){
	    push ( @checked_transcripts, $transcript );
	    next TRANSCRIPT;
	}
	
	my $correct  = 0;
	my $other    = 0;
	
	# in the forward strand, exons are on the original slice
	my $slice = $self->query;
	
      INTRON:
	for (my $i=0; $i<$#exons; $i++ ){
	    my $upstream_exon   = $exons[$i];
	    my $downstream_exon = $exons[$i+1];
	    my $upstream_site;
	    my $downstream_site;
	    if ($strand == 1){
		eval{
		    $upstream_site = 
			$self->query->subseq( ($upstream_exon->end     + 1), ($upstream_exon->end     + 2 ) );
		    $downstream_site = 
			$self->query->subseq( ($downstream_exon->start - 2), ($downstream_exon->start - 1 ) );
		};
		unless ( $upstream_site && $downstream_site ){
		    print STDERR "problems retrieving sequence for splice sites\n$@";
		    next INTRON;
		}
	    }
	    elsif( $strand == -1 ){
		# in the reverse strand, exon coords are forward in the revcomp_slice
		#
		#  example:
		#  exons originaly in - strand:       this is how we read the exons:
		#
		#     ------CT...AC---...            ---> 5...#exon1#CA...TC#exon2#...3'
		#  3' #exon2#GA...TG#exon1#...5'               ------GT...AG-----...
		#
		# we calculate CA..TC in the revcomp_query and 
		# make the complementary to get GT..AG (we do not need to apply the reverse)
		
		eval{
		    $upstream_site = 
			$self->query->subseq( ($upstream_exon->end     + 1), ($upstream_exon->end     + 2 ) );
		    $downstream_site = 
			$self->query->subseq( ($downstream_exon->start - 2), ($downstream_exon->start - 1 ) );
		};
		unless ( $upstream_site && $downstream_site ){
		    print STDERR "problems retrieving sequence for splice sites\n$@";
		    next INTRON;
		}
		$upstream_site   =~ tr/ACGTacgt/TGCAtgca/;
		$downstream_site =~ tr/ACGTacgt/TGCAtgca/;
	    }
	    
	    #print STDERR "upstream $upstream_site, downstream: $downstream_site\n";
	    ## good pairs of upstream-downstream intron sites:
	    ## ..###GT...AG###...   ...###AT...AC###...   ...###GC...AG###.
	    if (  ($upstream_site eq 'GT' && $downstream_site eq 'AG') ||
		  ($upstream_site eq 'AT' && $downstream_site eq 'AC') ||
		  ($upstream_site eq 'GC' && $downstream_site eq 'AG') ){
		$correct++;
	    }
	    else{
		$other++;
	    }
	} # end of INTRON
	
	unless ( $introns == $other + $correct ){
	    print STDERR "STRANGE: introns:  $introns, correct: $correct, other: $other\n";
	}
	
	if ( $other > 1 || $other > $correct ){
	    print STDERR "rejected - splice sites correct = $correct, other = $other";
	  Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_TranscriptEvidence($transcript);
	    next TRANSCRIPT;
	}
	else{
	    push ( @checked_transcripts, $transcript );
	}
	
    } # end of TRANSCRIPT
    
    return @checked_transcripts;
}

############################################################
 
sub _check_splice_Sites{
  my ($self, $ref_transcripts, $strand) = @_;

  if ( scalar( @$ref_transcripts ) == 0 ){
    return (0,0,0,0,0); 
  }

  print STDERR "EST_GeneBuilder: checking splice sites in strand $strand...\n";

  # get the contig being analysed
  my $slice = $self->query;
  
  # for reverse strand,  invert the contig, since the exons were retrieved in revcontig
  my $revslice = $slice->invert;   
  
  # upstream/downstream are with respect to the exon
  my ($upstream_AG, $downstream_GT) = (0,0);  # 99%
  my $downstream_GC                 = 0;       
  my ($upstream_AC, $downstream_AT) = (0,0);  # U12 spliceosome
  my $site_count = 0;
  
 TRANSCRIPT:
  foreach my $transcript (@$ref_transcripts){
    
    # count the number of exons
    my $count = 0; 
    my @exons = @{$transcript->get_all_Exons};
    
  EXON:
    foreach my $exon ( @exons ){

      # take the 2 bases right before the exon and after the exon
      my ($upstream, $downstream);
        
      # forward strand
      if ($strand == 1){
	my $seq = $slice; # a Bio::PrimarySeq object
	
	# catch possible exceptions in gettting the sequence (it might not be there!)
	eval{
	  $upstream   = $seq->subseq( ($exon->start)-2 , ($exon->start)-1 ); 
	};
	if ($@){
	  print STDERR "Unable to get subsequence (".(($exon->start)-2).",".(($exon->start)-1).")\n";
	  print STDERR $@;
	  $upstream = 'NN';
	}
	eval{
	  $downstream = $seq->subseq( ($exon->end)+1   , ($exon->end)+2   ); 
	};
	if ($@){
	  print STDERR "Unable to get subsequence (".(($exon->end)+1).",".(($exon->end)+2).")\n";
	  print STDERR $@;
	  $downstream = 'NN';
	}
	# print-outs to test it
	#      if ( $count ==0 ){
	#	print STDERR "FIRST EXON-->".$downstream;
	#      }
	#      if ( $count != 0 && $count != $#clusters ){
	#	print STDERR $upstream."-->EXON-->".$downstream;
	#      }
	#      if ( $count == $#clusters ){
	#	print STDERR $upstream."-->LAST EXON";
	#      }
	#      print "\n";
	
	# the first and last exon are not checked - potential UTR's
	if ( $count != 0 && $upstream eq 'AG') {       
	  $upstream_AG++;
	}
	elsif ( $count != 0 && $upstream eq 'AC') {       
	  $upstream_AC++;
	}
	if ( $count !=$#exons && $downstream eq 'GT') { 
	  $downstream_GT++;
	}
	elsif ( $count !=$#exons && $downstream eq 'GC') { 
	  $downstream_GC++;
	}
	elsif ( $count !=$#exons && $downstream eq 'AT') { 
	  $downstream_AT++;
	}

	$count++;
      } # end of forward strand
      
      # reverse strand
      if ($strand == -1 ){
	my $seq = $revslice; # a Bio::PrimarySeq object
	
	# catch possible exceptions in gettting the sequence (it might not be there!)
	eval{
	  $upstream = $seq->subseq( ($exon->start)-2 , ($exon->start)-1 ); 
	};
	if ($@){
	  print STDERR "Unable to get subsequence (".(($exon->start)-2).",".(($exon->start)-1).")\n";
	  #print STDERR $@;
	  $upstream ='NN';
	}
	eval{
	  $downstream   = $seq->subseq( ($exon->end)+1   , ($exon->end)+2   ); 
	};
	if ($@){
	  print STDERR "Unable to get subsequence (".(($exon->end)+1).",".(($exon->end)+2).")\n";
	  #print STDERR $@;
	  $downstream = 'NN';
	}
	#  in the reverse strand we're looking at coordinates in the reversed-complement slice:
	#
	#        $slice : --------------------TG----GA------------>  forward strand
	#                                     AC    CT               reverse strand 
	#                           downstream   EXON   upstream
	#
	#
	#                   upstream   EXON   downstream              
	#     $revslice : ----------AG----GT-------------------->    forward strand
	#                                                             reverse strand
	#
	# and take the reverse complement
	( $downstream = reverse( $downstream) ) =~ tr/ACGTacgt/TGCAtgca/;
	( $upstream   = reverse( $upstream )  ) =~ tr/ACGTacgt/TGCAtgca/;
	
	# so the conserved sequence should be AC<-EXON<-CT printed as in the reverse strand
	
	#      if ( $count ==0 ){
	#	print STDERR "LAST EXON<--".$upstream;
	#      }
	#      if ( $count != 0 && $count != $#clusters ){
	#	print STDERR $downstream."<--EXON<--".$upstream;
	#      }
	#      if ( $count == $#clusters ){
	#	print STDERR $downstream."<--FIRST EXON";
	#      }
	#      print "\n";
	
	# the first and last exon are not checked - potential UTR's
	# we count it as if everything was in the forward strand
	if ( $count != $#exons && $downstream eq 'AC') {       
	  $downstream_GT++;
	}
	elsif ( $count != $#exons && $downstream eq 'GC') {       
	  $downstream_GC++;
	}
	elsif ( $count != $#exons && $downstream eq 'AT') {       
	  $downstream_AT++;
	}
	if ( $count != 0 && $upstream eq 'CT') { 
	  $upstream_AG++;
	}
	elsif ( $count != 0 && $upstream eq 'GT') {
	  $upstream_AC++;
	}
	$count++;
      } # end of reverse strand
            
    }   # end of EXON
    
    $site_count += $#exons; # count the number of splice sites
    
  }   # end of TRANSCRIPT
    
  my $perc_upstream_AG   = 100*$upstream_AG/$site_count;
  my $perc_downstream_GT = 100*$downstream_GT/$site_count;
  my $prec_upstream_AC   = 100*$upstream_AC/$site_count;
  my $perc_downstream_AT = 100*$downstream_AT/$site_count;

  # most common spliceosome
  print STDERR "upstream   ( of an exon ) splice-sites AG: ".$upstream_AG. 
    " out of ".($site_count)." splice-sites, percentage: $perc_upstream_AG\n";
  print STDERR "downstream ( of an exon ) splice-sites GT: ".$downstream_GT. 
    " out of ".($site_count)." splice-sites, percentage: $perc_downstream_GT\n\n";

  # U12 spliceosome
  print STDERR "upstream   ( of an exon ) U12-splice-sites AC: ".$upstream_AC. 
    " out of ".($site_count)." splice-sites, percentage: $prec_upstream_AC\n";
  print STDERR "downstream ( of an exon ) U12-splice-sites AT: ".$downstream_AT. 
    " out of ".($site_count)." splice-sites, percentage: $perc_downstream_AT\n\n";

  # eventually, we may modify start/end of exons according to these checks
  return ($site_count, $upstream_AG, $downstream_GT, $upstream_AC, $downstream_AT);    
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

  ############### plus strand ##########
  $strand = 1;
  my $tcount=0;

  my @forward_transcripts;
 
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

    push (@forward_transcripts, $gw_runnable->output);
    # convert_output
    #$self->convert_output($gw_runnable, $strand);
  }
  print STDERR $tcount." transcripts run in genomewise (-smell 8) in the forward strand\n";

  my @checked_forward_transcripts   = $self->_check_Translations(\@forward_transcripts,$strand);
  
  my @forward_genes                 = $self->_cluster_into_Genes(@checked_forward_transcripts);
  
  my @forward_remapped              = $self->remap_genes(\@forward_genes, $strand);

  $self->output(@forward_remapped);


  ############### minus strand ##########
  $strand = -1;
  my $tcount2=0;
 
  my @reverse_transcripts;
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
    
    push (@reverse_transcripts, $gw_runnable->output );
    # convert_output
    #$self->convert_output($gw_runnable, $strand);
  }
  print STDERR $tcount2." transcripts run in genomewise (-smell 8) in the reverse strand\n";
  
  my @checked_reverse_transcripts = $self->_check_Translations(\@reverse_transcripts,$strand);
  
  my @reverse_genes               = $self->_cluster_into_Genes(@checked_reverse_transcripts);
  
  my @reverse_remapped            = $self->remap_genes(\@reverse_genes, $strand);

  $self->output(@reverse_remapped);

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

############################################
### convert genomewise output into genes ###
############################################

#sub convert_output {
#  my ($self, $gwr, $strand) = @_;

#  my @genes    = $self->make_genes($gwr, $strand);

#  my @remapped = $self->remap_genes(\@genes, $strand);

#  # store genes
#  $self->output(@remapped);
#}

############################################################
# this method set the slice in the exons for each transcript
# and check the translation

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

    # at this point, only accepts transcripts with more than one exon
    # although we have checked this already, genomewise sometimes bridges 
    # over introns making one exon out of two
    if ( scalar(@exons) == 1 ){
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
	print STDERR "TRANSLATION HAS STOP CODONS!!\n";
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

    Title   :   remap_genes
    Usage   :   $self->remap_genes(@genes)
    Function:   Remaps predicted genes into genomic coordinates
    Returns :   array of Bio::EnsEMBL::Gene
    Args    :   Bio::EnsEMBL::Virtual::Contig, array of Bio::EnsEMBL::Gene

=cut

sub remap_genes {
  my ($self, $genes, $strand) = @_;#
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
      
      # need to explicitly add back genetype and analysis.
      $newgene->type($gene->type);
      $newgene->analysis($gene->analysis);

      push(@newf,$newgene);
      
    };
    if ($@) {
      print STDERR "Couldn't reverse map gene [$@]\n";
      foreach my $t ( @{$gene->get_all_Transcripts} ){
	$self->_print_Transcript($t);
      }
    }
  }
  return @newf
}

=head2 vcontig

 Title   : vcontig
 Usage   : $obj->vcontig($newval)
 Function: inherited from RunnableDB, it sets/gets the virtual contig on which the analysis is carried out 
 Returns : value of vcontig
 Args    : newvalue (optional)


=cut

#####

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
  if (($end2 > $start1 && $start2 < $end1) ) {
	
	#  we have an overlap so we now need to return 
	#  two numbers reflecting how accurate the span 
	#  is. 
	#  0,0 means an exact match with the exon
	# a positive number means an over match to the exon
	# a negative number means not all the exon bases were matched

	my $left  = ($start2 - $start1);
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
    
   if(@genes){
     push(@{$self->{'_output'}},@genes);
   }

   return @{$self->{'_output'}};
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

  
  # flusold genes
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
    my $genetype = $self->genetype;
    my $analysis = $self->analysis;
    $gene->type($genetype);
    $gene->analysis($analysis);
    
    # sort them, get the longest CDS + UTR first (like in prune_Transcripts() )
    foreach my $transcript( @{$cluster} ){
      $gene->add_Transcript($transcript);
    }
    
    # prune out duplicate exons
    my $new_gene = $self->prune_Exons($gene);
    push( @genes, $new_gene );
   }
  
  #$self->final_genes(@genes);
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

sub transcript_high{
  my ($self,$tran) = @_;
  my $high;
  $tran->sort;
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
  $tran->sort;
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
    $tran->flush_Exons;
    foreach my $exon (@newexons) {
      $tran->add_Exon($exon);
    }
  }
  return $gene;
}

############################################################


1;


