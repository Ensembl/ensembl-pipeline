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

Bio::EnsEMBL::Pipeline::RunnableDB::MapGeneToExpression

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::RunnableDB::MapGeneToExpression->new(
									   -input_id  => $id,
									   );
    $obj->fetch_input;
    $obj->run;
    my %expression_map = %{ $obj->output };
    where @{ $expression_map{$transcript_id} } is an array of ests mapped to this transcript
    ests are here Bio::EnsEMBL::Transcript objects   
    
    one can do:
    my @transcript_ids = keys %expression_map;
    foreach my $transcript_id ( @transcripts_ids ){
	@mapped_ests = $expression_map{ $transcript_id };
    }

    $obj->write_output;


=head1 DESCRIPTION

Class to map genes read from an ensembl database to expression vocabulary via ESTs. 
ESTs are also read from an ensembl database. In principle, the typical situation
is to use ensembl genes and ests mapped to the genome.

=head1 CONTACT

eae@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::RunnableDB::MapGeneToExpression;

use diagnostics;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster;
use Bio::EnsEMBL::Pipeline::GeneComparison::GeneCluster;
use Bio::EnsEMBL::Pipeline::GeneComparison::GeneComparison;
use Bio::EnsEMBL::Pipeline::DBSQL::ExpressionAdaptor;



use Bio::EnsEMBL::Pipeline::ESTConf qw(
				       EST_INPUTID_REGEX
				       EST_REFDBHOST
				       EST_REFDBUSER
				       EST_REFDBNAME
				       EST_REFDBPASS
				       EST_E2G_DBHOST
				       EST_E2G_DBUSER
				       EST_E2G_DBNAME
				       EST_E2G_DBPASS
				       EST_TARGET_DBNAME
				       EST_TARGET_DBHOST
				       EST_TARGET_DBUSER
				       EST_TARGET_DBPASS      
				       EST_TARGET_GENETYPE
				       EST_GENEBUILDER_INPUT_GENETYPE
				       EST_EXPRESSION_DBHOST
				       EST_EXPRESSION_DBNAME
				       EST_EXPRESSION_DBUSER
				       EST_EXPRESSION_DBPASS
				      );

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

######################################################################

sub new{
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  
   
  # where the dna is
  my $refdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
						  -host             => $EST_REFDBHOST,
						  -user             => $EST_REFDBUSER,
						  -dbname           => $EST_REFDBNAME,
						);
  
  # where the genes are
  my $ensembl_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
						      '-host'   => $EST_TARGET_DBHOST,
						      '-user'   => $EST_TARGET_DBUSER,
						      '-dbname' => $EST_TARGET_DBNAME,
						      '-pass'   => $EST_TARGET_DBPASS,
						      '-dnadb'  => $refdb,
						     );
  

  # where the ests are (we actually want exonerate_e2g transcripts )
  unless( $self->db){
    my $est_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
						    '-host'   => $EST_E2G_DBHOST,
						    '-user'   => $EST_E2G_DBUSER,
						    '-dbname' => $EST_E2G_DBNAME,
						    '-dnadb'  => $refdb,
						   ); 
    $self->db($est_db);
  }
  $self->est_db( $self->db);
  $self->est_db->dnadb($refdb);
  $self->ensembl_db( $ensembl_db );
  $self->dna_db( $refdb );
  
  # database where the expression vocabularies are.
  # this is also where we are going to store the results
  my $expression_adaptor = Bio::EnsEMBL::Pipeline::DBSQL::ExpressionAdaptor->new(
									  '-host'   => $EST_EXPRESSION_DBHOST,
									  '-user'   => $EST_EXPRESSION_DBUSER,
									  '-dbname' => $EST_EXPRESSION_DBNAME,
									  #'-pass'   => $EST_EXPRESSION_DBPASS,
									 );
  
  $self->expression_adaptor($expression_adaptor);
  
  return $self;
  
}

#########################################################################
#
# GET/SET METHODS 
#
#########################################################################


sub ensembl_db{
  my ( $self, $db ) = @_;
  if ( $db ){
    $db->isa("Bio::EnsEMBL::DBSQL::DBAdaptor") || $self->throw("Input [$db] is not a Bio::EnsEMBL::DBSQL::DBAdaptor");
    $self->{'_ensembl_db'} = $db;
  }
  return $self->{'_ensembl_db'};
}

############################################################

sub expression_adaptor{
  my ( $self, $db ) = @_;
  if ( $db ){
    $db->isa("Bio::EnsEMBL::DBSQL::DBConnection") || $self->throw("Input [$db] is not a Bio::EnsEMBL::DBSQL::DBConnection");
    $self->{_expression_adaptor} = $db;
  }
  return $self->{_expression_adaptor};
}

############################################################

sub est_db{
  my ( $self, $db ) = @_;
  if ( $db ){
    $db->isa("Bio::EnsEMBL::DBSQL::DBAdaptor") || $self->throw("Input [$db] is not a Bio::EnsEMBL::DBSQL::DBAdaptor");
    $self->{'_est_db'} = $db;
  }
  return $self->{'_est_db'};
}

############################################################

sub dna_db{
  my ( $self, $db ) = @_;
  if ( $db ){
    $db->isa("Bio::EnsEMBL::DBSQL::DBAdaptor") || $self->throw("Input [$db] is not a Bio::EnsEMBL::DBSQL::DBAdaptor");
    $self->{'_dna_db'} = $db;
  }
  return $self->{'_dna_db'};
}

#############################################################

sub ensembl_slice{
  my ($self,$slice) = @_;
  if ( $slice ){
    $self->{'_ensembl_slice'} = $slice;
  }
  return $self->{'_ensembl_slice'};
}

#############################################################

sub est_slice{
  my ($self,$slice) = @_;
  if ( $slice ){
    $self->{'_est_slice'} = $slice;
  }
  return $self->{'_est_slice'};
}

############################################################

sub ensembl_genes{
  my ( $self, @genes ) = @_;
  unless( $self->{_ensembl_genes} ){
    $self->{_ensembl_genes} =[];
  }
  if ( @genes ){
    $genes[0]->isa("Bio::EnsEMBL::Gene") || $self->throw("$genes[0] is not a Bio::EnsEMBL::Gene");
    push ( @{ $self->{_ensembl_genes} }, @genes );
  }
  return @{ $self->{_ensembl_genes} };
}

#############################################################

sub ests{
  my ( $self, @genes ) = @_;

  unless( $self->{_ests} ){
    $self->{_ests} =[];
  }
  if ( @genes ){
    $genes[0]->isa("Bio::EnsEMBL::Gene") || $self->throw("$genes[0] is not a Bio::EnsEMBL::Gene");
    push ( @{ $self->{_ests} }, @genes );
  }
  return @{ $self->{_ests} };
}

#########################################################################


sub gene_Clusters {
  my ($self, @clusters) = @_;
  if (@clusters){
    push ( @{$self->{'_gene_clusters'} }, @clusters);
  }
  return @{ $self->{'_gene_clusters'} };
}

#########################################################################

# this holds the gene_id for each transcript

sub _gene_ID{
  my($self,$transcript_id,$gene_id) = @_;
  unless ( $transcript_id && $gene_id ){
    $self->warn("Need two parameters, transcript_id: $transcript_id, gene_id: $gene_id ");
    unless (  $self->{_gene_id}{$transcript_id} ){
      $self->{_gene_id}{$transcript_id} ='none';
    }
    if ($gene_id){
      $self->{_gene_id}{$transcript_id} = $gene_id;
    }
    return $self->{_gene_id}{$transcript_id};
  }
}  

############################################################
#
# FETCH INPUT
#
############################################################

sub fetch_input {
  my( $self) = @_;
  
  # get genomic region 
  my $input_id    = $self->input_id;
  unless ($input_id =~ /$EST_INPUTID_REGEX/ ){
    $self->throw("input $input_id not compatible with EST_INPUTID_REGEX $EST_INPUTID_REGEX");
    }
  my $chrname  = $1;
  my $chrstart = $2;
  my $chrend   = $3;
  
  print STDERR "Chromosome id = $chrname , range $chrstart $chrend\n";

  my $ensembl_sa = $self->ensembl_db->get_SliceAdaptor();
  my $est_sa     = $self->est_db->get_SliceAdaptor();

  my $ensembl_slice  = $ensembl_sa->fetch_by_chr_start_end($chrname,$chrstart,$chrend);
  my $est_slice      = $est_sa->fetch_by_chr_start_end($chrname,$chrstart,$chrend);

  $self->ensembl_slice( $ensembl_slice );
  $self->est_slice( $est_slice );

  # get ests (mapped with Filter_ESTs_and_E2G )
  print STDERR "getting genes of type $EST_GENEBUILDER_INPUT_GENETYPE\n";
  $self->ests(@{ $self->est_slice->get_all_Genes_by_type( $EST_GENEBUILDER_INPUT_GENETYPE, 'evidence' ) });
  print STDERR "got ".scalar( $self->ests )." ests\n";



  # get ensembl genes (from GeneBuilder)
  $self->ensembl_genes(@{ $self->ensembl_slice->get_all_Genes_by_type( $EST_TARGET_GENETYPE, 'evidence' ) });

}
  
############################################################
#
# RUN METHOD
#
############################################################

sub run{
  my ($self,@args) = @_;

  my @genes = $self->ensembl_genes;
  my @est   = $self->ests;

  # first cluster genes by locus
  # calculate on each cluster

  # for each gene in the cluster
  # for each transcript
  # calcultate the ests that map to this transcript

  # if there no genes, we finish a earlier
  unless ( $self->ensembl_genes ){
    print STDERR "no genes found in this region, leaving...\n";
    exit(0);
  }
  print STDERR scalar( $self->ensembl_genes )." genes retrieved\n";
  unless ( $self->ests ){
    print STDERR "No ests in this region, leaving...\n";
    exit(0);
  }
   
  # cluster the genes:
  my @clusters = $self->cluster_Genes( $self->ensembl_genes, $self->ests );
  
 CLUSTER:
  foreach my $cluster ( @clusters ){
      
      # get genes of each type
      my @genes = $cluster->get_Genes_of_Type( $EST_TARGET_GENETYPE );
      my @ests  = $cluster->get_Genes_of_Type( $EST_GENEBUILDER_INPUT_GENETYPE );
      
      # if we have genes of either type, let's try to match them
      if ( @genes && @ests ){
	  print STDERR "Trying to match ".scalar(@genes)." ensembl genes and ".scalar(@ests)." ests\n"; 
	  
	  ############################################################
	  #
	  #
	  # See if you can map genes to clone libraries instead of ESTs
	  #
	  #
	  ############################################################
	  
	  my @est_transcripts;
	  foreach my $est ( @ests ){
	      my @est_trans = @{$est->get_all_Transcripts};
	      push ( @est_transcripts, $self->in_SANBI( @est_trans ));
	  }
	  
	  foreach my $gene ( @genes ){
	      push ( @ensembl_transcripts,  @{$gene->get_all_Transcripts} );
	  }
	  
	  my $matcher = 
	    Bio::EnsEMBL::Pipeline::GeneComparison::GenericTranscriptMatcher->new(
										  -reference_set => \@ensembl_transcripts,
										  -match_set => \@est_transcripts,
										  );
	  
	  $matcher->run;
	  
	  my $matching_map = $matcher->output;
	  $self->output( $matching_map );
      }
      
      
      
      # else we could have only ensembl genes
      elsif(  @genes && !@ests ){
	  # we have nothing to modify them, hence we accept them...
	  print STDERR "Skipping cluster with no ests\n";
	  next CLUSTER;
      }
      # else we could have only ests
      elsif( !@genes && @ests ){
	  print STDERR "Cluster with no genes\n";
	  next CLUSTER;
      }
      # else we could have nothing !!?
      elsif( !@genes && !@ests ){
	  print STDERR "empty cluster, you must be kidding!\n";
	  next CLUSTER;
      }
  } # end of CLUSTER
  
  # before returning, check that we have written anything
  unless( $self->output ){
      print STDERR "No matches found, I'm so sorry\n";
      exit(0);
  }
  return;
}

############################################################
#
# METHODS CALLED FROM RUN METHOD... DOING ALL THE MAGIC
#
############################################################

# this method cluster genes only according to genomic extent
# covered by the genes. The proper clustering of transcripts
# to give rise to genes occurs in _cluster_into_Genes()

sub cluster_Genes{
  my ($self, @genes) = @_;

  # first sort the genes by the left-most position coordinate ####
  my %start_table;
  my $i=0;
  foreach my $gene (@genes){
    $start_table{$i} = $self->_get_start_of_Gene( $gene );
    $i++;
  }
  my @sorted_genes=();
  foreach my $k ( sort { $start_table{$a} <=> $start_table{$b} } keys %start_table ){
    push (@sorted_genes, $genes[$k]);
  }
  
  # we can start clustering
  print "Clustering ".scalar( @sorted_genes )." genes...\n";
  
  # create a new cluster 
  my $cluster = Bio::EnsEMBL::Pipeline::GeneComparison::GeneCluster->new();
  my $cluster_count = 1;
  my @clusters;
  
  # before putting any genes, we must declare the types
  my $ensembl    = [$EST_TARGET_GENETYPE];
  my $est        = [$EST_GENEBUILDER_INPUT_GENETYPE];
  $cluster->gene_Types($ensembl,$est);

  # put the first gene into these cluster
  $cluster->put_Genes( $sorted_genes[0] );
  push (@clusters, $cluster);
  
  # loop over the rest of the genes
 LOOP:
  for (my $c=1; $c<=$#sorted_genes; $c++){
    my $found=0;
    
    # treat the clusters as ranges, so we only need to check if ranges overlap
    # for the moment this is enough
    my $gene_start = $self->_get_start_of_Gene( $sorted_genes[$c] );
    my $gene_end   = $self->_get_end_of_Gene(   $sorted_genes[$c] );
    
    # we need to do this each time, so that start/end get updated
    my $cluster_start = $cluster->start;
    my $cluster_end   = $cluster->end;

    if ( !( $gene_end < $cluster_start || $gene_start > $cluster_end ) ){
      $cluster->put_Genes( $sorted_genes[$c] );
    }
    else{
      # else, create a new cluster
      $cluster = new Bio::EnsEMBL::Pipeline::GeneComparison::GeneCluster; 
      $cluster->gene_Types($ensembl,$est);
      $cluster->put_Genes( $sorted_genes[$c] );
      $cluster_count++;
      push( @clusters, $cluster );
    }
  }

  print STDERR "returning ".scalar(@clusters)." clusters\n";
  return @clusters;
}			


#########################################################################

# this gives the left-most exon coordinate in a gene

sub _get_start_of_Gene{  
  my ($self,$gene) = @_;
  my $start;
  foreach my $tran ( @{$gene->get_all_Transcripts} ){
    foreach my $exon ( @{$tran->get_all_Exons} ){
      unless ($start){
	$start = $exon->start;
      }
      if ( $exon->start < $start ){
	$start = $exon->start;
      }
    }
  }
  return $start;
}


#########################################################################

# this gives the right-most exon coordinate in a gene

sub _get_end_of_Gene{  
  my ($self,$gene) = @_;
  my $end;
  foreach my $tran ( @{$gene->get_all_Transcripts} ){
    foreach my $exon ( @{$tran->get_all_Exons} ){
      unless ($end){
	$end = $exon->end;
      }
      if ( $exon->end > $end ){
	$end = $exon->end;
      }
    }
  }
  return $end;
}

#########################################################################

sub _in_SANBI{
  my ($self,@ests) = @_;
  
  # @ests are transcript objects
  my %id_to_transcript;

  my @est_ids;
 EST:
  foreach my $est ( @ests ){
    my $est_id = $self->_find_est_id($est);
    unless ($est_id){
      #print STDERR "No accession found for ".$est->dbID."\n";
      next EST;
    }
    #print STDERR "est: $est, est_id: $est_id\n";
    if ( $est_id =~/(\S+)\.(\d+)/ ){
      $est_id = $1;
    }
    push( @est_ids, $est_id );
    $id_to_transcript{$est_id} = $est;
  }
  my $expression_adaptor = $self->expression_adaptor;
  my @pairs = $expression_adaptor->get_libraryId_by_estarray( @est_ids );
  
  my @found_ests;
  foreach my $pair ( @pairs ){
    if ( $$pair[1] ){
      push ( @found_ests, $id_to_transcript{ $$pair[0] } );
    }
  }
  return @found_ests;
}


#########################################################################

sub expression_Map{
  my ($self,$transcript,$est) = @_;
  if ( $transcript ){
    my $transcript_id;
      if ($transcript->stable_id){
	$transcript_id = $transcript->stable_id;
      }
      elsif( $transcript->dbID ){
	$transcript_id = $transcript->dbID;
      }
    unless ( $self->{_est_map}{$transcript_id} ){
      $self->{_est_map}->{$transcript_id} = [];
    }
    if ($est){
      push ( @{  $self->{_est_map}->{$transcript_id} }, $est );
    }
  }
  return $self->{_est_map};
}

#########################################################################

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

# method to calculate the exonic length of a transcript which is inside a gene

sub _transcript_exonic_length{
  my ($self,$tran) = @_;
  my $exonic_length = 0;
  foreach my $exon (@{$tran->get_all_Exons}){
    $exonic_length += ($exon->end - $exon->start + 1);
  }
  return $exonic_length;
}

#########################################################################
# method to calculate the length of a transcript in genomic extent, 

sub _transcript_length{
    my ($self,$tran) = @_;
    my @exons= @{$tran->get_all_Exons};
    my $genomic_extent = 0;
    if ( $exons[0]->strand == -1 ){
      @exons = sort{ $b->start <=> $a->start } @exons;
      $genomic_extent = $exons[0]->end - $exons[$#exons]->start + 1;
    }
    elsif( $exons[0]->strand == 1 ){
      @exons = sort{ $a->start <=> $b->start } @exons;
      $genomic_extent = $exons[$#exons]->end - $exons[0]->start + 1;
    }
    return $genomic_extent;
}

#########################################################################
#
# METHODS INVOLVED IN WRITTING THE RESULTS
#
#########################################################################


#########################################################################

=head2 output

the output is a list of Bio::EnsEMBL::Pipeline::GeneComparison::ObjectMap objects

=cut

sub output{
  my ($self, $map)= @_;
  unless ($self->{_output}){
      $self->{_output} = [];
  }
  if ($map){
      push(@{$self->{_output}}, $map );
  }
  return @{$self->{_output}};
}

#########################################################################
# get the est id throught the supporting evidence of the transcript

sub _find_est_id{
  my ($self, $est) = @_;
  my %is_evidence;
  foreach my $exon (@{$est->get_all_Exons}){
    foreach my $evidence ( @{$exon->get_all_supporting_features} ){
      $is_evidence{ $evidence->hseqname } = 1;
    }
  }
  my @evidence = keys %is_evidence;
  unless ( $evidence[0] ){
    print STDERR "No evidence for ".$est->dbID.", hmm... possible sticky single exon gene\n";
  }
  return $evidence[0];
}
 
#########################################################################

sub write_output {
    my ($self) = @_;
    my $expression_adaptor = $self->expression_adaptor;
    
    my @maps = $self->output;
    
    foreach my $map ( @maps ){
	my @list1 = $map->list1;
	foreach my $transcript ( @list1 ){
	    
	    my $t_id;
	    if ($transcript->stable_id){
		$t_id = $transcript->stable_id;
	    }
	    elsif($transcript->dbID){
		$t_id = $transcript->dbID;
	    }
	    
	    my @est_matches = $map->partners($transcript);
	    my @est_ids;
	    foreach my $est ( @est_matches ){
		my $est_id = $self->_find_est_id($est);
		if ( $est_id){
		    my $est_id_no_version;
		    if ( $est_id =~/(\S+)\.(\d+)/){
			$est_id_no_version = $1;
		    }
		    else{
			$est_id_no_version = $est_id;
		    }
		    push (@est_ids, $est_id_no_version);
		}
	    }
	    print STDERR "Storing pairs $transcript_id, @est_ids\n";
	    $expression_adaptor->store_ensembl_link($transcript_id,\@est_ids);
	}
    }
}

########################################################################


1;
