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

Bio::EnsEMBL::Pipeline::RunnableDB::Combine_GeneBuilder_ESTGenes.pm

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::RunnableDB::GeneCombiner->new(
								    -input_id  => $id
								   );
    $obj->fetch_input
    $obj->run

    my @genes = $obj->output;


=head1 DESCRIPTION

It combines the genes produced by GeneBuilder ( from protein information )
with the genes produced by EST_GeneBuilder (from EST and cDNA information)

=head1 CONTACT

eae@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::RunnableDB::GeneCombiner;

use diagnostics;
use vars qw(@ISA);
use strict;

use Bio::Range;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::GeneComparison::GeneCluster;
use Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster;
use Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptComparator;
use Bio::EnsEMBL::Pipeline::GeneCombinerConf;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Bio::EnsEMBL::Pipeline::Tools::ExonUtils;

# config file; parameters searched for here if not passed in as @args

use Bio::EnsEMBL::Pipeline::GeneCombinerConf qw(
						ESTGENE_DBHOST
						ESTGENE_DBUSER
						ESTGENE_DBNAME
						ESTGENE_DBPASS 
						ESTGENE_TYPE
						ESTGENE_MAX_INTRON_LENGTH
						
						ENSEMBL_DBHOST
						ENSEMBL_DBUSER
						ENSEMBL_DBNAME
						ENSEMBL_DBPASS
						ENSEMBL_TYPE
						
						REF_DBHOST
						REF_DBUSER
						REF_DBNAME
						REF_DBPASS
						
						FINAL_DBHOST
						FINAL_DBNAME
						FINAL_DBUSER
						FINAL_DBPASS
						FINAL_TYPE
					       
						GENECOMBINER_INPUTID_REGEX
					       
					       );




@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

######################################################################

sub new{
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  
   
  my $refdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
						  -host             => $REF_DBHOST,
						  -user             => $REF_DBUSER,
						  -dbname           => $REF_DBNAME,
						);
  
  my $ensembl_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
						      '-host'   => $ENSEMBL_DBHOST,
						      '-user'   => $ENSEMBL_DBUSER,
						      '-dbname' => $ENSEMBL_DBNAME,
						      '-dnadb' => $refdb,
						     );
  

  my $estgene_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
						      '-host'   => $ESTGENE_DBHOST,
						      '-user'   => $ESTGENE_DBUSER,
						      '-dbname' => $ESTGENE_DBNAME,
						      '-dnadb' => $refdb,
						     ); 

  
  # dbobj is read by the parent class RunnableDB and it holds the FINAL_DB database
  
  unless( $self->db ){


      my $final_db = new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'   => $FINAL_DBHOST,
							'-user'   => $FINAL_DBUSER,
							'-dbname' => $FINAL_DBNAME,
							'-pass'   => $FINAL_DBPASS,
							'-dnadb'  => $refdb,
						   ); 
      $self->db($final_db);
    }
  
  $self->final_db( $self->db ); 
  $self->final_db->dnadb($refdb);


  # needs to read from two databases and write into another one (possibly a third?)
  
  $self->ensembl_db( $ensembl_db );
  $self->estgene_db( $estgene_db );
  
  return $self;
  
}

#########################################################################
#
# GET/SET METHODS 
#
#########################################################################

sub final_db{
  my ( $self, $db ) = @_;
  if ( $db ){
    $db->isa("Bio::EnsEMBL::DBSQL::DBAdaptor") || $self->throw("Input [$db] is not a Bio::EnsEMBL::DBSQL::DBAdaptor");
    $self->{_final_db} = $db;
  }
  return $self->{_final_db};
}


#############################################################


sub ensembl_db{
  my ( $self, $db ) = @_;
  if ( $db ){
    $db->isa("Bio::EnsEMBL::DBSQL::DBAdaptor") || $self->throw("Input [$db] is not a Bio::EnsEMBL::DBSQL::DBAdaptor");
    $self->{'_ensembl_db'} = $db;
  }
  return $self->{'_ensembl_db'};
}

############################################################

sub estgene_db{
  my ( $self, $db ) = @_;
  if ( $db ){
    $db->isa("Bio::EnsEMBL::DBSQL::DBAdaptor") || $self->throw("Input [$db] is not a Bio::EnsEMBL::DBSQL::DBAdaptor");
    $self->{'_estgene_db'} = $db;
  }
  return $self->{'_estgene_db'};
}


#############################################################

sub ensembl_vc{
  my ($self,$vc) = @_;
  if ( $vc ){
    $self->{'_ensembl_vc'} = $vc;
  }
  return $self->{'_ensembl_vc'};
}

#############################################################

sub estgene_vc{
  my ($self,$vc) = @_;
  if ( $vc ){
    $self->{'_estgene_vc'} = $vc;
  }
  return $self->{'_estgene_vc'};
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

sub estgenes{
  my ( $self, @genes ) = @_;
  unless( $self->{_estgenes} ){
    $self->{_estgenes} =[];
  }
  if ( @genes ){
    $genes[0]->isa("Bio::EnsEMBL::Gene") || $self->throw("$genes[0] is not a Bio::EnsEMBL::Gene");
    push ( @{ $self->{_estgenes} }, @genes );
  }
  return @{ $self->{_estgenes} };
}

#########################################################################

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

############################################################
#
# FETCH INPUT
#
############################################################

sub fetch_input {
  my( $self) = @_;
  
  # get genomic region 
  my $chrid    = $self->input_id;
  print STDERR "input_id: $chrid\n";
   
  my $chrname;
  my $chrstart;
  my $chrend;
  if ( $chrid =~/$GENECOMBINER_INPUTID_REGEX/ ){
    $chrname  = $1;
    $chrstart = $2;
    $chrend   = $3;
  }
  else{
    $self->throw("Not a valid input_id... $chrid");
  }
  print STDERR "Chromosome id = $chrname , range $chrstart $chrend\n";

  my $ensembl_gpa = $self->ensembl_db->get_SliceAdaptor();
  my $estgene_gpa = $self->estgene_db->get_SliceAdaptor();

  my $ensembl_vc  = $ensembl_gpa->fetch_by_chr_start_end($chrname,$chrstart,$chrend);
  my $estgene_vc  = $estgene_gpa->fetch_by_chr_start_end($chrname,$chrstart,$chrend);
  
  $self->ensembl_vc( $ensembl_vc );
  $self->estgene_vc( $estgene_vc );

  print STDERR $self->ensembl_vc."\t".$self->estgene_vc."\n";

  
}
  
############################################################
#
# RUN METHOD
#
############################################################

sub run{
  my ($self,@args) = @_;

  # get estgenes ( from EST_GeneBuilder )
  print STDERR "getting genes of type $ESTGENE_TYPE\n";
  my @est_genes = @{ $self->estgene_vc->get_all_Genes_by_type( $ESTGENE_TYPE, 'evidence' ) };
  
  # get ensembl genes (from GeneBuilder)
  print STDERR "getting genes of type $ENSEMBL_TYPE\n";
  my @ensembl_genes = @{ $self->ensembl_vc->get_all_Genes_by_type( $ENSEMBL_TYPE, 'evidence' ) };

  # if there are no genes, we finish earlier
  unless ( @est_genes ){
    unless (@ensembl_genes){
	print STDERR "no genes found, leaving...\n";
	  exit(0);
    }
    print STDERR "No estgenes found, writing ensembl genes as they are\n";
    my @transcripts;
    foreach my $gene (@ensembl_genes){
      push(@transcripts, @{$gene->get_all_Transcripts});
    }
    my @newgenes = $self->_make_Genes(\@transcripts);
      
    my @remapped = $self->_remap_Genes(\@newgenes);
    
    $self->output(@remapped);
    return;
  }
  
  # need to CLONE all genes, as their transcripts may share exon objects
  # which makes it dangerous when modifying exon coordinates

  my @cloned_ensembl_genes;
  foreach my $gene (@ensembl_genes){
      my $newgene = $self->_clone_Gene($gene);
      push( @cloned_ensembl_genes, $newgene);
  }
  $self->ensembl_genes( @cloned_ensembl_genes );
  
  my @cloned_est_genes;
  foreach my $gene ( @est_genes){
      my $newgene = $self->_clone_Gene($gene);
      push( @cloned_est_genes, $newgene );
  }
  $self->estgenes( @cloned_est_genes );

  # store the original transcripts, just for the numbers
  my @original_ens_transcripts;
  foreach my $gene ( $self->ensembl_genes ){
    my @transcripts =@{ $gene->get_all_Transcripts};
  TRAN1:
    foreach my $transcript (@transcripts){
      $transcript->type($gene->type);
      $self->_transcript_Type($transcript,$ENSEMBL_TYPE);
      push (  @original_ens_transcripts, $transcript );
    }
  }
  my @original_est_transcripts;
  foreach my $gene ( $self->estgenes ){
    my @transcripts = @{$gene->get_all_Transcripts};
  TRAN2:
    foreach my $transcript (@transcripts){
      $transcript->type($gene->type);
      $self->_transcript_Type($transcript,$ESTGENE_TYPE);
      push (  @original_est_transcripts, $transcript );
    }
  }

  # cluster estgenes and ensembl genes
  my @genes             = ( $self->ensembl_genes, $self->estgenes );
  my @gene_clusters     = $self->cluster_Genes( @genes );
  my @unclustered_genes = $self->unclustered_Genes;
  
  # the accepted set of transcripts
  my @transcripts;
 
  ############################################################
  # at this stage we could separate the clusters into three sorts:

  # 1) only ensembl genes --> we leave them as they are
  # 2) ensembl+est genes  -->  first: extend UTRs
  #                           second: include alternative forms (do checks) 
  # 3) only est genes     --> if only 1 in the cluster, take it if good coverage 
  #                       --> if more: take a set that share exons among themselves (exact matches)
  #                                    not necessarily all the same exon, it should be the
  #                                    largest connected graph where vertices are trnascripts and
  #                                    edges is the relation 'share an exon'
  ############################################################

 CLUSTER:
  foreach my $cluster ( @gene_clusters ){
    
    # get genes of each type
    my @ens_genes = $cluster->get_Genes_of_Type( $ENSEMBL_TYPE );
    my @est_genes = $cluster->get_Genes_of_Type( $ESTGENE_TYPE );
    
    ############################################################
    # if we have genes of either type, let's try to match them
    ############################################################
    if ( @ens_genes && @est_genes ){
      print STDERR "Matching ".scalar(@ens_genes)." ensembl genes and ".scalar(@est_genes)." est genes\n"; 
      my @ens_transcripts;
      my @est_transcripts;
      
      print STDERR "=== ensembl genes ===\n";
      foreach my $gene ( @ens_genes ){	
      TRAN1:
	foreach my $transcript ( @{$gene->get_all_Transcripts} ){
	  $transcript->type($gene->type);
	  Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Evidence($transcript);
	  unless (Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_check_Transcript($transcript, $self->ensembl_vc)
		  && Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_check_Translation($transcript) ){
	    print STDERR "skipping this transcript\n";
	    next TRAN1;
	  }
	  push ( @ens_transcripts, $transcript );
	}
      }
      
      print STDERR "=== est genes ===\n";
      foreach my $gene ( @est_genes ){
      TRAN2:
	foreach my $transcript ( @{$gene->get_all_Transcripts} ){
	  $transcript->type($gene->type);
	  Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Evidence($transcript);
	  # reject those with too long introns
	  if ( $self->_too_long_intron_size( $transcript ) ){
	    next TRAN2;
	  }
	  unless (Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_check_Transcript($transcript, $self->estgene_vc)
		  && Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_check_Translation($transcript) ){
	    print STDERR "skipping this transcript\n";
	    next TRAN2;
	  }
	  push ( @est_transcripts, $transcript );
	}
      }
      
      my ( $new_ens, $accepted_est ) = $self->_pair_Transcripts( \@ens_transcripts, \@est_transcripts );
      if ( $new_ens ){
	push ( @transcripts, @$new_ens );
      }
      if ( $accepted_est ){
	push ( @transcripts, @$accepted_est );
      }
    }
    
    ############################################################
    # else we could have only ensembl genes
    ############################################################
    elsif(  @ens_genes && !@est_genes ){
      # we have nothing to modify them, hence we accept them...
      my @ens_transcripts;
      
      # but before we check, just in case, you know
      print STDERR "=== ensembl genes only ===\n";
      foreach my $gene ( @ens_genes ){
      TRAN3:
	foreach my $transcript ( @{$gene->get_all_Transcripts} ){
	  $transcript->type($gene->type);
	  Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Evidence($transcript);
	  unless (Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_check_Transcript($transcript, $self->ensembl_vc)
		  && Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_check_Translation($transcript) ){
	    print STDERR "skipping this transcript\n";
	    next TRAN3;
	  }
	  push ( @ens_transcripts, $transcript );
	}
      }
      print STDERR "Accepting ".scalar(@ens_transcripts)." ens-transcripts\n";
      push ( @transcripts, @ens_transcripts );
    }
    
    ############################################################
    # else we could have only est genes
    ############################################################
    elsif( !@ens_genes && @est_genes ){
      print STDERR "Checking ".scalar(@est_genes)." est genes\n";
      my @est_transcripts;
      
      print STDERR "=== est genes only ===\n";
      foreach my $gene ( @est_genes ){
      TRAN4:
	foreach my $transcript ( @{$gene->get_all_Transcripts} ){
	  $transcript->type($gene->type);
	  Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Evidence($transcript);
	  # reject those with too long introns
	  if ( $self->_too_long_intron_size( $transcript ) ){
	    next TRAN4;
	  }
	  unless (Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_check_Transcript($transcript,$self->estgene_vc)
		  && Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_check_Translation($transcript) ){
	    print STDERR "skipping this transcript\n";
	    next TRAN4;
	  }
	  push ( @est_transcripts, $transcript );
	}
      }
      
      # we have to check the est genes to see whether they are ok
      # they must have some common properties in order
      # to form a proper set of alternative forms
      my @accepted_trans = $self->_check_est_Cluster( @est_transcripts );
      push ( @transcripts, @accepted_trans );
    }
    
    # else we could have nothing !!?
    elsif( !@ens_genes && !@est_genes ){
      print STDERR "empty cluster, you must be kidding!\n";
    }
  } # end of CLUSTER
  
  # make the genes 
  my @newgenes = $self->_make_Genes(\@transcripts);

  print STDERR "================ GENES MADE =====================\n";
  print STDERR scalar($self->ensembl_genes)." ensembl genes (".scalar(@original_ens_transcripts)." transcripts) and ".
      scalar($self->estgenes)." estgenes (".scalar(@original_est_transcripts)." transcripts)\n";    
  print STDERR "have produced ".scalar(@transcripts)." transcripts, which are clustered into ".scalar(@newgenes)." genes:\n";
  my $count = 0;
  foreach my $gene (@newgenes){
      $count++;
      print STDERR "gene $count: ".$gene->type."\n";
      foreach my $transcript (@{$gene->get_all_Transcripts}){
	Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Evidence($transcript);
      }
  }

  print STDERR scalar(@newgenes)." genes 2b remapped\n";
  # remap them to raw contig coordinates
  my @remapped = $self->_remap_Genes(\@newgenes);
  print STDERR "==================== REMAPPED GENES =====================\n";
  foreach my $gene ( @remapped ){
      #print STDERR "type    : ".$gene->type."\n";
      #print STDERR "analysis: ".$gene->analysis->dbID." ".$gene->analysis->logic_name."\n";
      
      # test:
      foreach my $t ( @{$gene->get_all_Transcripts} ){
	Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Evidence($t);
      }
  }
  
  print STDERR scalar(@remapped)." genes remapped\n";
  
  # store the genes
  $self->output(@remapped);
  print STDERR scalar($self->output)." stored, to be written on the db\n";


  return @remapped;

}


############################################################
#
# METHODS CALLED FROM RUN METHOD... DOING ALL THE MAGIC
#
############################################################




# this method cluster genes only according to genomic extent
# covered by the genes. The proper clustering of transcripts
# to give rise to genes occurs in _cluster_into_Genes()

sub cluster_Genes {
  my ($self) = @_;
  my @genes = ( $self->ensembl_genes, $self->estgenes );
 
  my @forward_genes;
  my @reverse_genes;
 
 GENE:
  foreach my $gene ( @genes ){
    foreach my $transcript (@{$gene->get_all_Transcripts}){
      if ( $transcript->start_Exon->strand == 1 ){
	push( @forward_genes, $gene );
      }
      else{
	push( @reverse_genes, $gene );
      }
      next GENE;
    }
  }
  my @clusters;
  if ( @forward_genes ){
    my @forward_clusters = $self->_cluster_Genes_by_genomic_range( @forward_genes );
    if ( @forward_clusters){
      push( @clusters, @forward_clusters);
    }
  }
  if ( @reverse_genes ){
    my @reverse_clusters = $self->_cluster_Genes_by_genomic_range( @reverse_genes );
    if ( @reverse_clusters){
      push( @clusters, @reverse_clusters);
    }
  }
  return @clusters;
}
############################################################

sub _cluster_Genes_by_genomic_range{
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
  my $ensembl    = [$ENSEMBL_TYPE];
  my $genomewise = [$ESTGENE_TYPE];
  $cluster->gene_Types($ensembl,$genomewise);

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
      $cluster->gene_Types($ensembl,$genomewise);
      $cluster->put_Genes( $sorted_genes[$c] );
      $cluster_count++;
      push( @clusters, $cluster );
    }
  }

  print STDERR "returning ".scalar(@clusters)." clusters\n";
  return @clusters;
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

sub unclustered_Genes{
  my ($self,@genes ) = @_;
  unless( $self->{_unclustered_genes}){
    $self->{_unclustered_genes} = [];
  }
  if (@genes){
    push( @{ $self->{_unclustered_genes} }, @genes );
  }
  return @{ $self->{_unclustered_genes} };
}

#########################################################################

# this gives the left-most exon coordinate in a gene

sub _get_start_of_Gene{  
  my ($self,$gene) = @_;
  my $start;
  foreach my $tran ( @{$gene->get_all_Transcripts}){
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
  foreach my $tran ( @{$gene->get_all_Transcripts}){
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

sub _pair_Transcripts{
  my ($self,$ens_transcripts,$est_transcripts) = @_;
  my @potential_isoforms;

  my @ens_transcripts = @$ens_transcripts;
  my @est_transcripts = @$est_transcripts;
  my %used_est_transcript;
  my @accepted_isoforms;
  
  # sort the transcripts by their genomic and exonic length
  @ens_transcripts = sort {  my $result = ( $self->_transcript_length($b) <=>
					    $self->_transcript_length($a) );
			     unless ($result){
			       return ( $self->_transcript_exonic_length($b) <=>
					$self->_transcript_exonic_length($a) );
			     }
			     return $result;
			   } @ens_transcripts;
  
  
  @est_transcripts = sort {  my $result = ( $self->_transcript_length($b) <=>
					    $self->_transcript_length($a) );
			     unless ($result){
			       return ( $self->_transcript_exonic_length($b) <=>
					$self->_transcript_exonic_length($a) );
			     }
			     return $result;
			   } @est_transcripts;
  

  ###### LOOK FOR ESTGENES THAT EXTEND UTRs

  # matrix holding the number and length of exon overlap for each pair of transcripts
  my $overlap_number_matrix;
  my $overlap_length_matrix;
  
  my %merge_list;
  my %selected; # keeps track of the est-transcripts that have been used
  
  # we base everything on the ensembl transcripts
 ENS_TRANSCRIPT:
  foreach my $ens_tran ( @ens_transcripts ){
    
    # we store them according to whether there is a est/cdna matching both UTRs or just one
    my @list_both_ends;
    my @list_one_end;
    my @lists;

  EST_TRANSCRIPT:
    foreach my $est_tran ( @est_transcripts ){
      
      # check with which est_transcripts it can merge
      my ($merge,$overlaps) = $self->_test_for_Merge( $ens_tran, $est_tran );
      if ($merge == 1 && $overlaps > 0 ){
	
	print STDERR "Can merge:\n";
	Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Evidence($ens_tran);
	Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Evidence($est_tran);
	
	# we prefer those with both UTR ends matched
	my ($match_5prime, $match_3prime) = $self->_check_UTRMatches( $ens_tran,$est_tran);
	
	# calculate then how much they overlap
	my ($overlap_number,$overlap_length, $exact_matches) = _compare_Transcripts( $ens_tran, $est_tran );
	
	if ( $match_5prime == $match_3prime && $match_5prime == 1){
	  push( @list_both_ends, [ $overlap_number, $overlap_length, $exact_matches, $est_tran ]);
	}
	elsif( $match_5prime == 1 || $match_3prime == 1 ){
	  push( @list_one_end, [ $overlap_number, $overlap_length, $exact_matches, $est_tran ]);
	}
      }
    }
    print STDERR "possible merges: both_ends: ".scalar(@list_both_ends)."\n";
    foreach my $element ( @list_both_ends ){
	print STDERR "dbID: ".$element->[3]->dbID."\n";
    }
    print STDERR "possible merges: one_end: ".scalar(@list_one_end)."\n";
    foreach my $element ( @list_one_end ){
	print STDERR "dbID: ".$element->[3]->dbID."\n";
    }
    push( @lists, \@list_both_ends, \@list_one_end);
    
    # take the est_transcript which overlap the most
    my @sorted_lists;
    foreach my $list ( @lists ){
      my @sorted_pairs = sort { my $result = ( $$b[2] <=> $$a[2] );
				if ($result){
				  return $result;
				}
				else{
				  $result = ( $$b[0] <=> $$a[0] );
				  if ($result){
				    return $result;
				  }
				  else{
				    return ( $$b[1] <=> $$a[1] );
				  }
				}
			      } @{ $list };
      push( @sorted_lists, @sorted_pairs );
    }
    # test:
    print STDERR "Candidates:\n";
    foreach my $pair ( @sorted_lists ){
      Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript( $pair->[3] );
    }
    
    # try to merge it, if succeeded, mark the chosen one and leave the rest
  MODIFY:
    for (my $i = 0; $i< scalar(@sorted_lists); $i++){
	print STDERR "trying to modify ".$ens_tran->dbID."\n";
	my $est_tran = $sorted_lists[$i][3];
	unless( $used_est_transcript{ $est_tran } && $used_est_transcript{ $est_tran } == 1){
	    my $modified;
	    ( $ens_tran, $modified ) = $self->_extend_UTRs($ens_tran,$est_tran);
	    if ($modified == 1 ){
		$used_est_transcript{$est_tran} = 1;
		print STDERR "Merged\n";
		print STDERR "transfering supporting evidence\n";
		$self->_transfer_transcript_supporting_evidence($est_tran,$ens_tran);

                # we only allow one cdna to modify each transcript
		last MODIFY;
	    }
	}
    }
  }  # end of ENS_TRANSCRIPT
  
  ###### LOOK FOR POSSIBLE ALTERNATIVE TRANSCRIPTS
  # those which have not been used for extending UTRs are candidate isoforms
  #
  # This analysis relies on the fact that the set of cdna-genes being used
  # is non redundant, i.e. it forms a putative set of isoforms
  my @candidates;
  foreach my $est_tran ( @est_transcripts ){
    unless (  $used_est_transcript{ $est_tran } && $used_est_transcript{$est_tran} == 1 ){
      push ( @candidates, $est_tran );
    }
  }
  print STDERR "\n".scalar(@candidates)." est transcripts left to be tested as isoforms\n\n";
  if (@candidates){
      
      # at the moment we just check whether there is an ens_transcript
      # that shares one exon, part of the protein product and
      # there is one intron matching one exon or vice versa.
    
  CANDIDATE:
    foreach my $est_tran ( @candidates ){
      my $one_exon_shared = 0;
      my $similar_protein = 0;
      my $exon_in_intron  = 0;
      my $exact_merge     = 0;
      my $fuzzy_merge     = 0;
      my $overlaps = 0;

      # check first it does not merge with anything in the cluster
      foreach my $ens_tran ( @ens_transcripts ){
	my $comparator1 = Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptComparator->new(
											    -comparison_level         => 3,
											    -splice_mismatch          => 1,
											   );
	my $comparator2 = Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptComparator->new(
											    -comparison_level         => 2,
											   );
	
	($exact_merge,$overlaps) = $comparator2->compare($est_tran, $ens_tran);
	($fuzzy_merge,$overlaps) = $comparator1->compare($est_tran, $ens_tran);
	if ($exact_merge == 1){
	  print STDERR "est transcript ".$est_tran->dbID." is redundant, it merges exactly with an ensembl transcript, skipping it\n";
	  next CANDIDATE;
	}
	if ( $fuzzy_merge == 1){
	      print STDERR "est transcript ".$est_tran->dbID." is redundant, it merges 'fuzzily' with an ensembl transcript, skipping it\n";
	      next CANDIDATE;
	    }
	  }
	  
      foreach my $ens_tran ( @ens_transcripts ){
	$one_exon_shared = $self->_check_exact_exon_Match(  $est_tran, $ens_tran);
	$similar_protein = $self->_check_protein_Match(     $est_tran, $ens_tran);
	$exon_in_intron  = $self->_check_exon_Skipping(     $est_tran, $ens_tran);
	
	print "Comparing\n";
	Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($est_tran);
	Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($ens_tran);
	print STDERR "shared_exon    : $one_exon_shared\n";
	print STDERR "similar protein: $similar_protein\n";
	print STDERR "exon_in_intron : $exon_in_intron\n";
	
	print STDERR "not using the protein info\n";
	if ( $one_exon_shared == 1 
	     #&& $similar_protein == 1
	     && $exon_in_intron  == 1
	     #&& $exact_merge     == 0 
	   ){
	    print STDERR "BINGO, and ISOFORM found !!\n";
	    
	    # if they don't have the same start/end translation, try to lock phases:
	    if ( !( $est_tran->translation->start == $ens_tran->translation->start) ||  
		 !( $est_tran->translation->end   == $ens_tran->translation->end  )   ){ 
		print STDERR "Transcripts don't have the same start/end translation. We are not trying to lock the start/end of translation:\n";
		#Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($est_tran);
		#$est_tran = $self->_lock_Phases($est_tran, $ens_tran) ;
	    }
	    else{
		print STDERR "EST-transcript has the same translation start and end:\n";
		print STDERR "EST translation start: ".$est_tran->translation->start." end: ".$est_tran->translation->end."\n";
		print STDERR "ensembl transl  start: ".$ens_tran->translation->start." end: ".$ens_tran->translation->end."\n";
	    }
	  push ( @accepted_isoforms, $est_tran );
	  
	  next CANDIDATE; 
	  ############################################################
	  # NOTE: we accept the est gene if it does not merge with any
	  # of the ensembl genes and if it is a valid isoform
	  # with respect to just ONE ensembl gene.
	  # We could try to extend that to every ensembl gene
	  # in the cluster.
	  # This would make it harder to lock phases
	  ############################################################
	}
      }
    }  # end of CANDIDATE 
  }    # end of 'if (@candidates)'
  else{
      # nothing to do then
  }
  
  # return the ens_transcripts (modified or not)
  return ( \@ens_transcripts, \@accepted_isoforms );
  
}
#########################################################################

sub _lock_Phases{
  my ($self,$est_tran,$ens_tran) = @_;
  
  # keep the original est_tran translation
  my $original_translation = $est_tran->translation;
  my $original_transcript  = $self->_clone_Transcript($est_tran);

  # we go first for the easy-peasy case: when the est gene
  # overlaps the start exon and end exon of the ensembl translation
  my $ens_start_exon = $ens_tran->translation->start_Exon;
  my $ens_end_exon   = $ens_tran->translation->end_Exon;
  my $est_start_exon;
  my $est_end_exon;
  my $est_start_translation;
  my $est_end_translation;
  
  foreach my $exon (@{$est_tran->get_all_Exons}){
      if ( $exon->overlaps( $ens_start_exon ) ){
	  $est_start_exon = $exon;
	  
	  if ( $exon->strand == 1 && $ens_start_exon->strand == 1 ){
	      $est_start_translation = 
		  $ens_start_exon->start - $exon->start + $ens_tran->translation->start;
	  }
	  elsif( $exon->strand == -1 && $ens_start_exon->strand == -1 ){
	      $est_start_translation = 
		  $exon->end - $ens_start_exon->end + $ens_tran->translation->start;
	  }
      }
      
      if ( $exon->overlaps( $ens_end_exon ) ){
	  $est_end_exon = $exon;
	  
	  if ( $exon->strand == 1 && $ens_end_exon->strand == 1 ){
	      $est_end_translation = 
		  $ens_end_exon->start - $exon->start + $ens_tran->translation->end;
	  }
	  elsif( $exon->strand == -1 && $ens_end_exon->strand == -1 ){
	      $est_end_translation =
		  $exon->end - $ens_end_exon->end + $ens_tran->translation->end;
	  }
      }
  }
  
  # check that we get something useful:
  unless ( $est_start_exon && $est_end_exon ){
    print STDERR "could not lock phases. Translation start/end exon not overlapping\n";
    return $est_tran;
  }
  unless ( $est_start_translation > 0 && $est_start_translation <= $est_start_exon->length ){
    print STDERR "could not lock phases. Bad start translation: $est_start_translation\n";
    return $est_tran;
  }
  unless ( $est_end_translation > 0 && $est_end_translation <= $est_end_exon->length ){
    print STDERR "could not lock phases. Bad end translation: $est_end_translation\n";
    return $est_tran;
  }

  # recompute phases:
  $est_tran->sort;
  my $seen_start = 0;
  my $seen_end   = 0;
  my $previous_exon;
  foreach my $exon (@{$est_tran->get_all_Exons}){
    if ( $exon == $est_start_exon && $exon == $est_end_exon ){
	$exon->phase(-1);
	$exon->end_phase(-1);
	$seen_start = 1;
	$seen_end   = 1;
    }
    elsif ( $exon == $est_start_exon ){
	$exon->phase(-1);
	$exon->end_phase( ( $exon->length - $est_start_translation )%3 );
	$seen_start = 1;
    }
    elsif( $exon == $est_end_exon ){
      $exon->end_phase(-1);
      $exon->phase( $previous_exon->end_phase );
      $seen_end = 1;
    }
    elsif( $seen_start == 1 && $seen_end == 0 ){
      $exon->phase( $previous_exon->end_phase );
      $exon->end_phase( ( $exon->phase + $exon->length )%3 );
    }
    elsif( ($seen_start == 0 && $seen_end == 0 ) || ( $seen_start == 1 && $seen_end == 1) ){
      $exon->phase(-1);
      $exon->end_phase(-1);
    }
    $previous_exon = $exon;
  }
  
  my $new_est_translation = Bio::EnsEMBL::Translation->new();
  $new_est_translation->start_Exon( $est_start_exon );
  $new_est_translation->end_Exon($est_end_exon);
  $new_est_translation->start( $est_start_translation );
  $new_est_translation->end( $est_end_translation );
  
  $est_tran->translation( $new_est_translation );
  unless ( Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_check_Translation( $est_tran )){
      print STDERR "Oh, shit we have to put back the original phases!!!!!!!!!!\n";
      print STDERR "Getting back to the original translation:\n";
    Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($original_transcript);
      return $original_transcript;
  }
  
  print STDERR "Succeeded in setting phases Returning:\n";
 Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($est_tran);
  return $est_tran;
}
    
  




############################################################
# this function takes est_transcripts that have been clustered together
# but not with any ensembl transcript and tries to figure out whether they
# make an acceptable set of alt-forms

sub _check_est_Cluster{
  my ($self,@est_transcripts) = @_;
  my %color;

  #if ( scalar(@est_transcripts) == 1 ){
  print STDERR "cluster with ".scalar(@est_transcripts)." transcripts\n";
  #}

  # adjacency lists:
  my %adj;
  
  for(my $i=0;$i<scalar(@est_transcripts);$i++){
    for(my $j=0;$j<scalar(@est_transcripts);$j++){
      
      next if $j==$i;
      print STDERR "Comparing transcripts:\n";
      Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($est_transcripts[$i]);
      Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($est_transcripts[$j]);
      if ( $self->_check_exact_exon_Match( $est_transcripts[$i], $est_transcripts[$j]) &&
	   $self->_check_protein_Match(    $est_transcripts[$i], $est_transcripts[$j])    ){
	push ( @{ $adj{$est_transcripts[$i]} } , $est_transcripts[$j] );
      }
    }
  }
  
  print STDERR "adjacency lists:\n";
  foreach my $tran (@est_transcripts){
    print STDERR $tran->dbID." -> ";
    foreach my $link ( @{ $adj{ $tran } } ){
      print STDERR $link.",";
    }
    print STDERR "\n";
  }
  
  foreach my $tran ( @est_transcripts ){
    $color{$tran} = "white";
  }
  
  my @potential_genes;
  
  # find the connected components doing a depth-first search
  foreach my $tran ( @est_transcripts ){
    if ( $color{$tran} eq 'white' ){
      my @potential_gene;
      $self->_visit( $tran, \%color, \%adj, \@potential_gene);
      push ( @potential_genes, \@potential_gene );
    }
  }
  print STDERR scalar(@potential_genes)." potential genes created\n";
  
  # take only the sets with more than one transcript?:
  #@potential_genes = sort { scalar( @{ $b } ) <=> scalar( @{ $a } ) } @potential_genes;
  my @accepted;
  my $first_gene = shift @potential_genes;
  while ( $first_gene && scalar( @{ $first_gene } ) > 0 ){
    if ( scalar(@{ $first_gene } ) == 1 ){
      my $ok = $self->_check_Completeness( $first_gene );
    }
    push( @accepted, @{ $first_gene } );
    $first_gene = shift @potential_genes;
  }
  return @accepted;
}

#########################################################################

sub _visit{
  my ($self, $node, $color, $adj, $potential_gene) = @_;
  
  # node is a transcript object;
  $color->{ $node } = 'gray';

  foreach my $trans ( @{ $adj->{$node} } ){
    if ( $color->{ $trans } eq 'white' ){
      $self->_visit( $trans, $color, $adj, $potential_gene );
    }
  }
  unless ( $color->{$node} eq 'black'){
    push( @{ $potential_gene }, $node);
  }
  $color->{ $node } = 'black';    
  return;
}

#########################################################################

sub _check_Completeness{
  my ($self, $array_ref ) = @_;
  my @transcripts = @{ $array_ref };
  my $transcript = $transcripts[0];
  
  #take the 3' exon:
  $transcript->sort;
  my @exons = @{$transcript->get_all_Exons};
  my $three_prime = $exons[$#exons];
  my $seq;
  eval{
    $seq = $three_prime->seq;
  };
  if ($seq){
    if ( $seq =~/AATAAA/ || $seq =~/ATTAAA/ ){
      print STDERR "polyA signal included\n";
      return 1;
    }
    else{
      print STDERR "polyA signal NOT included\n";
    }
  }
  else{
    print STDERR "cannot retrieve exon sequence\n";
  }
  return 0;
}
    
#########################################################################

sub _clone_Gene{
  my ($self,$gene) = @_;
  
  my $newgene = new Bio::EnsEMBL::Gene;
  if ($gene->type){
      $newgene->type( $gene->type);
  }
  if ( defined $gene->dbID ){
      $newgene->dbID($gene->dbID);
  }
  if ( defined $gene->analysis ){
      $newgene->analysis($gene->analysis);
  }
  foreach my $transcript (@{$gene->get_all_Transcripts}){
    my $newtranscript = $self->_clone_Transcript($transcript);
    $newgene->add_Transcript($newtranscript);
  }
  return $newgene;
}

#########################################################################

sub _clone_Transcript{
  my ($self,$transcript) = @_;
  
  #print STDERR "Cloning:\n";
  #$self->_print_Transcript($transcript);
  my $newtranscript  = new Bio::EnsEMBL::Transcript;
  my $newtranslation = new Bio::EnsEMBL::Translation;
  
  my $translation_start_exon = $transcript->translation->start_Exon;
  my $translation_end_exon   = $transcript->translation->end_Exon; 
  
  foreach my $exon ( @{$transcript->get_all_Exons} ){
    my $newexon = $self->_clone_Exon($exon);
    if ($exon == $translation_start_exon){
      $newtranslation->start_Exon($newexon);
      $newtranslation->start($transcript->translation->start);
    }
    if ($exon == $translation_end_exon){
      $newtranslation->end_Exon($newexon);
      $newtranslation->end($transcript->translation->end);
    }
    $newtranscript->add_Exon($newexon);
  }
  #$newtranscript->sort;
  $newtranscript->dbID($transcript->dbID);
  if (defined $transcript->type ){
      $newtranscript->type($transcript->type);
  }
  $self->_transcript_Type($newtranscript,$self->_transcript_Type($transcript));
  $newtranscript->translation($newtranslation);
 
  return $newtranscript;
}

#########################################################################

sub _clone_Exon{
  my ($self,$exon) = @_;
  my $newexon = new Bio::EnsEMBL::Exon;
  $newexon->start      ($exon->start);
  $newexon->end        ($exon->end);
  $newexon->phase      ($exon->phase);
  $newexon->end_phase  ($exon->end_phase);
  $newexon->strand     ($exon->strand);
  $newexon->dbID       ($exon->dbID);
  $newexon->contig     ($exon->contig);
  $newexon->sticky_rank($exon->sticky_rank);
  $newexon->seqname    ($exon->seqname);
  $newexon->attach_seq ($self->ensembl_vc);

  foreach my $evidence ( @{$exon->get_all_supporting_features} ){
    $newexon->add_supporting_features( $evidence );
  }
  return $newexon;
}

#########################################################################

sub _transcript_Type{
  my($self,$transcript,$type) = @_;
  unless (  $self->{type}{$transcript} ){
    $self->{type}{$transcript} ='none';
  }
  if ($type){
    $self->{type}{$transcript} = $type;
  }
  return $self->{type}{$transcript};
}

#########################################################################

# having a set of est_genes only, we have no reference transcript (ensembl one),
# so to determine whether two transcripts are two alternative forms of the same gene
# we check whether they share at least an exon

sub _check_exact_exon_Match{
 my ($self, $tran1, $tran2 ) = @_;
 my @exons1 = @{$tran1->get_all_Exons};
 my @exons2 = @{$tran2->get_all_Exons};
 my $exact_match = 0;
 
 # how many exact matches we need (maybe 1 is enough)
 foreach my $exon1 (@exons1){
     foreach my $exon2 (@exons2){
	 return 1 if  ( $exon1->start == $exon2->start && $exon2->end == $exon2->end );
     }
 }
 return 0;
}

#########################################################################
#
# having a set of est_genes only, we have no reference transcript (ensembl one),
# so to determine whether two transcripts are two alternative forms of the same gene
# we check whether they have a similar protein product

sub _check_protein_Match{
 my ($self, $tran1, $tran2 ) = @_;

 my $seq1;
 my $seq2;

 my $compatible_proteins = 0;
 eval{
     $seq1 = $tran1->translate;
     $seq2 = $tran2->translate;
 };
 if ( $seq1 && $seq2 ){
     
   if ( $seq1 =~/\*/ || $seq2 =~/\*/ ){ 
     print STDERR "On of the peptides has a stop codon\n";
     return 0;
   }
   if ( $seq1 eq $seq2 ){
     print STDERR "Identical translation\n";
     $compatible_proteins = 1;
   }
   elsif( $seq1 =~/$seq2/ || $seq2 =~/$seq1/ ){
     $compatible_proteins = 1;
   }
 }
 if ( $compatible_proteins != 0 ){
     return 1;
 }
 else{
     return 0;
 }
}
#########################################################################
#
# this function is another check used to find out if two given transcripts
# are two alternative variants. This function returns 1 if there is at least
# one 'skipped exon' between the two

# Note that overlaps exon-intron is not checked to be an intron 
# containing the complete exon. This allows also one partial exon skipping
# where the different in the two transcripts will be in the size of some exons

sub _check_exon_Skipping{
  my ($self,$tran1,$tran2) = @_;

  my @exons1 = @{$tran1->get_all_Exons};
  my @exons2 = @{$tran2->get_all_Exons};	
  
  # order is important here
  @exons1 = sort { $a->start <=> $b->start } @exons1;
  @exons2 = sort { $a->start <=> $b->start } @exons2;

  # we actually compare ranges:
  # we create two lists per transcript: exon_list and intron_list
  my @exon_list1;
  my @intron_list1;
  for(my $i=0; $i<scalar(@exons1);$i++){
    my $exon_range = Bio::Range->new();
    $exon_range->start($exons1[$i]->start);
    $exon_range->end($exons1[$i]->end);
    push( @exon_list1, $exon_range);
    if ( $i+1 < scalar(@exons1)){
      my $intron_range = Bio::Range->new();
      my $intron_start = $exons1[$i]->end + 1;
      my $intron_end   = $exons1[$i+1]->start - 1;
      $intron_range->start($intron_start);
      $intron_range->end($intron_end);
      push( @intron_list1, $intron_range);
    }
  }
  my @exon_list2;
  my @intron_list2;
  for(my $i=0; $i<scalar(@exons2);$i++){
    my $exon_range = Bio::Range->new();
    $exon_range->start($exons2[$i]->start);
    $exon_range->end($exons2[$i]->end);
    push( @exon_list2, $exon_range);
    if ( $i+1 < scalar(@exons2)){
      my $intron_range = Bio::Range->new();
      my $intron_start = $exons2[$i]->end + 1;
      my $intron_end   = $exons2[$i+1]->start - 1;
      $intron_range->start($intron_start);
      $intron_range->end($intron_end);
      push( @intron_list2, $intron_range);
    }
  }
  # probably, the list with more introns will be more likely to have one exon
  # in a big intron of the other transcript, it won't be too many comparisons if we
  # don't optimize that
  
  for (my $i=0; $i<scalar(@exon_list1);$i++){
    for ( my $j=0; $j<scalar(@intron_list2);$j++){
      
      # if one exon overlaps one intron
      if ( $exon_list1[$i]->overlaps( $intron_list2[$j] ) ){
	
	# either the exon falls completely in the intron:
	if ( $exon_list1[$i]->start >= $intron_list2[$j]->start &&
	     $exon_list1[$i]->end   <= $intron_list2[$j]->end ){
	  return 1;
	}
	# or the exon overlaps the 2 flanking exons
	elsif ( $exon_list1[$i]->overlaps( $exon_list2[$j]) &&
		defined( $exon_list2[$j+1] )                &&
		$exon_list1[$i]->overlaps( $exon_list2[$j+1] ) ){
	  return 1;
	}
	# else, we could allow and exon to partially overlap one intron, but how much overlap?
	# (the previous OR the next exon overlaps one of the flanking exons )
	#if ( ( $i>0            && $exon_list1[$i-1]->overlaps( $exon_list2[$j] ) ) ||
	#     ( $i<$#exon_list1 && 
	#       $j<$#exon_list2 &&
	#       $exon_list1[$i+1]->overlaps( $exon_list2[$j+1] ) )
	#   ){
	#  return 1;
	#}
      }
    }
  }
  
  for (my $i=0; $i<scalar(@exon_list2);$i++){
    for ( my $j=0; $j<scalar(@intron_list1);$j++){
      
      # if one exon overlaps one intron
      if ( $exon_list2[$i]->overlaps( $intron_list1[$j] ) ){
		
	# either the exon falls completely in the intron:
	if ( $exon_list2[$i]->start >= $intron_list1[$j]->start &&
	     $exon_list2[$i]->end   <= $intron_list1[$j]->end ){
	  return 1;
	}
	# or the exon overlaps the 2 flanking exons
	elsif ( $exon_list2[$i]->overlaps( $exon_list1[$j]) &&
		defined( $exon_list1[$j+1] )                &&
		$exon_list2[$i]->overlaps( $exon_list1[$j+1] ) ){
	  return 1;
	}
	# else, we could allow and exon to partially overlap one intron, but how much overlap?
	# and the previous OR the next exon overlaps one of the flanking exons
	#if ( ( $i>0            && $exon_list2[$i-1]->overlaps( $exon_list1[$j] ) ) ||
	#     ( $i<$#exon_list2 &&
	#       $j<$#exon_list1 &&
	#       $exon_list2[$i+1]->overlaps( $exon_list1[$j+1] ) )
	#   ){
	#  return 1;
	#}
      }
    }
  }
  
    return 0;
}

#########################################################################



# this function checks whether two transcripts merge
# according to consecutive exon overlap
# it only considers 1-to-1 matches, so things like
#                        ____     ____        
#              exons1 --|____|---|____|------ etc... $j
#                        ____________  
#              exons2 --|____________|------ etc...  $k
#
# are considered a mismatch

sub _test_for_Merge{
  my ($self,$tran1,$tran2) = @_;
  my @exons1 = @{$tran1->get_all_Exons};
  my @exons2 = @{$tran2->get_all_Exons};	
 
  my $foundlink = 0; # flag that gets set when starting to link exons
  my $start     = 0; # start looking at the first one
  my $overlaps  = 0; # independently if they merge or not, we compute the number of exon overlaps
  my $merge     = 0; # =1 if they merge

  my $one2one_overlap = 0;
  my $one2two_overlap = 0;
  my $two2one_overlap = 0;
 EXON1:
  for (my $j=0; $j<=$#exons1; $j++) {
    
  EXON2:
    for (my $k=$start; $k<=$#exons2; $k++){
    #print STDERR "comparing ".($j+1)." and ".($k+1)."\n";
	    
      # if exon 1 is not the first, check first whether it matches the previous exon2 as well, i.e.
      #                        ____     ____        
      #              exons1 --|____|---|____|------ etc... $j
      #                        ____________  
      #              exons2 --|____________|------ etc...  $k
      #
      if ($foundlink == 1 && $j != 0){
	if ( $k != 0 && $exons1[$j]->overlaps($exons2[$k-1]) ){
	  #print STDERR ($j+1)." <--> ".($k)."\n";
	  $overlaps++;
	  $two2one_overlap++;
	  next EXON1;
	}
      }
      
      # if texons1[$j] and exons2[$k] overlap go to the next exon1 and  next $exon2
      if ( $exons1[$j]->overlaps($exons2[$k]) ){
	#print STDERR ($j+1)." <--> ".($k+1)."\n";
        $overlaps++;
	
        # in order to merge the link always start at the first exon of one of the transcripts
        if ( $j == 0 || $k == 0 ){
          $foundlink = 1;
        }
      }          
      else {  
	# if you haven't found an overlap yet, look at the next exon 
	if ( $foundlink == 0 ){
	  next EXON2;
	}
	# leave if we stop finding links between exons before the end of transcripts
	if ( $foundlink == 1 ){
	  $merge = 0;
	  last EXON1;
	}
      }
      
      # if foundlink = 1 and we get to the end of either transcript, we merge them!
      if ( $foundlink == 1 && ( $j == $#exons1 || $k == $#exons2 ) ){
	
	# and we can leave
        $merge = 1;
	last EXON1;
      }
      # if foundlink = 1 but we're not yet at the end, go to the next exon 
      if ( $foundlink == 1 ){
	
	# but first check whether in exons2 there are further exons overlapping exon1, i.e.
        #                       ____________        
	#             exons1 --|____________|------ etc...
	#                       ____     ___  
	#             exons2 --|____|---|___|------ etc...
	# 
	my $addition = 0;
	while ( $k+1+$addition < scalar(@exons2) && $exons1[$j]->overlaps($exons2[$k+1+$addition]) ){
	  #print STDERR ($j+1)." <--> ".($k+2+$addition)."\n";
	  $one2two_overlap++;
	  $overlaps++;
          $addition++;
	}      
	$start = $k+1+$addition;
	next EXON1;
      }    
      
    } # end of EXON2 
    
    # if you haven't found any match for this exon1, start again from the first exon2:
    if ($foundlink == 0){
      $start = 0;
    }
 
  }   # end of EXON1      

  # we only make them merge if $merge = 1 and the 2-to-1 and 1-to-2 overlaps are zero;
  if ( $merge == 1 && $one2two_overlap == 0 && $two2one_overlap == 0 ){
    return ( 1, $overlaps );
  }
  else{
    return ( 0, $overlaps);
  }
}
  
#########################################################################
# this function checks whether two transcripts merge
# with exact exon matches, except for
# possible mismatches in the extremal exons

sub test_for_semiexact_Merge{
  my ($self,$est_tran,$ens_tran) = @_;
  
  my @exons1 = @{$est_tran->get_all_Exons};
  my @exons2 = @{$ens_tran->get_all_Exons};	
  
  @exons1 = sort {$a->start <=> $b->start} @exons1;
  @exons2 = sort {$a->start <=> $b->start} @exons2;

  my $foundlink = 0; # flag that gets set when starting to link exons
  my $start     = 0; # start looking at the first one
  my $merge     = 0; # =1 if they merge
  
 EXON1:
  for (my $j=0; $j<=$#exons1; $j++) {
      
    EXON2:
      for (my $k=$start; $k<=$#exons2; $k++){
	  #print STDERR "comparing j = $j : ".$exons1[$j]->start."-".$exons1[$j]->end.
	  #    " and k = $k : ".$exons2[$k]->start."-".$exons2[$k]->end."\n";
	  
	  # we allow some mismatches at the extremities
	  #                        ____     ____     ___   
	  #              exons1   |____|---|____|---|___|  $j
	  #                         ___     ____     ____  
	  #              exons2    |___|---|____|---|____|  $k
	  
	  # if there is no overlap, go to the next EXON2
	  if ( $foundlink == 0 && !($exons1[$j]->overlaps($exons2[$k]) ) ){
	      #print STDERR "foundlink = 0 and no overlap --> go to next EXON2\n";
	      next EXON2;
	  }
	  # if there is no overlap and we had found a link, there is no merge
	  if ( $foundlink == 1 && !($exons1[$j]->overlaps($exons2[$k]) ) ){
	      #print STDERR "foundlink = 1 and no overlap --> leaving\n";
	      $merge = 0;
	      last EXON1;
	  }	
	  
	  # the first exon can have a mismatch in the start
	  if ( ($k == 0 || $j == 0) && $exons1[$j]->end == $exons2[$k]->end ){
	      
	      # but if it is also the last exon
	      if ( ( ( $k == 0 && $k == $#exons2 )   || 
		     ( $j == 0 && $j == $#exons1 ) ) ){
		  
		  # we force it to match the start
		  if ( $exons1[$j]->start == $exons2[$k]->start ){
		      $foundlink  = 1;
		      $merge      = 1;
		      #print STDERR "merged single exon transcript\n";
		      last EXON1;
		  }
		  # we call it a non-merge
		  else{
		      $foundlink = 0;
		      $merge     = 0;
		      #print STDERR "non-merged single exon transcript\n";
		      last EXON1;
		  }
	      }
	      else{
		  #else, we have a link
		  $foundlink = 1;
		  $start = $k+1;
		  #print STDERR "found a link\n";
		  next EXON1;
	      }
	  }
	  # the last one can have a mismatch on the end
	  elsif ( ( $k == $#exons2 || $j == $#exons1 ) &&
		  ( $foundlink == 1 )                  &&
		  ( $exons1[$j]->start == $exons2[$k]->start ) 
		  ){
	      #print STDERR "link completed, merged transcripts\n";
	      $merge = 1;
	      last EXON1;
	  }
	  # the middle one must have exact matches
	  elsif ( ($k != 0 && $k != $#exons2) && 
		  ($j != 0 && $j != $#exons1) &&
		  ( $foundlink == 1)          &&
		  ( $exons1[$j]->start == $exons2[$k]->start ) &&
		  ( $exons1[$j]->end   == $exons2[$k]->end   )
		  ){
	      $start = $k+1;
	      #print STDERR "continue link\n";
	      next EXON1;
	  }

      } # end of EXON2 
    
      if ($foundlink == 0){
	  $start = 0;
      }
      
  }   # end of EXON1      
  
  return $merge;
}

#########################################################################
# this function checks whether two transcripts merge
# with fuzzy exon matches: there is consecutive exon overlap 
# but there are mismatches of 2 base allowed at the edges of any exon pair
#
# Why 2 bases: 2 bases is perhaps not meaningful enough to be considered
# a biological difference, and it is possibly an artifact of any of the
# analysis previously run: genomewise, est2genome,... it is more likely to
# happen 

sub _test_for_fuzzy_semiexact_Merge{
  my ($self,$est_tran,$ens_tran) = @_;
  
  my @exons1 = @{$est_tran->get_all_Exons};
  my @exons2 = @{$ens_tran->get_all_Exons};	
  
  @exons1 = sort {$a->start <=> $b->start} @exons1;
  @exons2 = sort {$a->start <=> $b->start} @exons2;

  my $foundlink = 0; # flag that gets set when starting to link exons
  my $start     = 0; # start looking at the first one
  my $merge     = 0; # =1 if they merge
  
 EXON1:
  for (my $j=0; $j<=$#exons1; $j++) {
      
    EXON2:
      for (my $k=$start; $k<=$#exons2; $k++){
	#print STDERR "comparing j = $j : ".$exons1[$j]->start."-".$exons1[$j]->end.
	#  " and k = $k : ".$exons2[$k]->start."-".$exons2[$k]->end."\n";
	
	  # we allow some mismatches at the extremities
	  #                        ____     ____     ___   
	  #              exons1   |____|---|____|---|___|  $j
	  #                         ___     ____     ____  
	  #              exons2    |___|---|____|---|____|  $k
	  
	  # if there is no overlap, go to the next EXON2
	  if ( $foundlink == 0 && !($exons1[$j]->overlaps($exons2[$k]) ) ){
	    #print STDERR "foundlink = 0 and no overlap --> go to next EXON2\n";
	    next EXON2;
	  }
	  # if there is no overlap and we had found a link, there is no merge
	  if ( $foundlink == 1 && !($exons1[$j]->overlaps($exons2[$k]) ) ){
	      #print STDERR "foundlink = 1 and no overlap --> leaving\n";
	      $merge = 0;
	      last EXON1;
	  }	
	  
	  # the first exon can have a mismatch ( any number of bases) in the start
	  # and a 2base mismatch at the end
	  if ( ($k == 0 || $j == 0) && abs($exons1[$j]->end - $exons2[$k]->end)<3 ){
	      
	      # but if it is also the last exon
	      if ( ( ( $k == 0 && $k == $#exons2 )   || 
		     ( $j == 0 && $j == $#exons1 ) ) ){
		  
		  # we force it to match the start (with a mismatch of 2bases allowed)
		  if ( abs($exons1[$j]->start - $exons2[$k]->start)< 3 ){
		      $foundlink  = 1;
		      $merge      = 1;
		      #print STDERR "merged single exon transcript\n";
		      last EXON1;
		  }
		  # we call it a non-merge
		  else{
		      $foundlink = 0;
		      $merge     = 0;
		      #print STDERR "non-merged single exon transcript\n";
		      last EXON1;
		  }
	      }
	      else{
		  #else, we have a link
		  $foundlink = 1;
		  $start = $k+1;
		  #print STDERR "found a link\n";
		  next EXON1;
	      }
	  }
	  # the last one can have any mismatch on the end
	  # but must have a match at the start (wiht 2bases mismatch allowed)
	  elsif ( ( $k == $#exons2 || $j == $#exons1 ) &&
		  ( $foundlink == 1 )                  &&
		  ( abs($exons1[$j]->start - $exons2[$k]->start)<3 ) 
		  ){
	      #print STDERR "link completed, merged transcripts\n";
	      $merge = 1;
	      last EXON1;
	  }
	# the middle one must have exact matches
	# (up to a 2base mismatch)
	elsif ( ($k != 0 && $k != $#exons2) && 
		($j != 0 && $j != $#exons1) &&
		( $foundlink == 1)          &&
		abs( $exons1[$j]->start - $exons2[$k]->start )<3 &&
		abs( $exons1[$j]->end   - $exons2[$k]->end   )<3
		  ){
	      $start = $k+1;
	      #print STDERR "continue link\n";
	      next EXON1;
	  }

      } # end of EXON2 
    
      if ($foundlink == 0){
	  $start = 0;
      }
      
  }   # end of EXON1      
  
  return $merge;
}

#########################################################################
# this gets the left most start coordinate for the transcript, regardless of the strand

sub _get_start_of_Transcript{
  my ($self,$transcript) = @_;
  my @exons = @{$transcript->get_all_Exons};
  @exons    = sort { $a->start <=> $b->start } @exons;
  my $start = $exons[0]->start;
  
  return $start;
}

#########################################################################
   
# this compares both transcripts and calculate the number of overlapping exons,
# the length of the overlap, and the number of exact overlapping exons

sub _compare_Transcripts {         
  my ($tran1, $tran2) = @_;
  my @exons1   = @{$tran1->get_all_Exons};
  my @exons2   = @{$tran2->get_all_Exons};

  my $overlaps = 0;
  my $overlap_length = 0;
  my $exact = 0;

  foreach my $exon1 (@exons1){
      foreach my $exon2 (@exons2){
	  if ( ($exon1->overlaps($exon2)) && ($exon1->strand == $exon2->strand) ){
	      $overlaps++;
	      
	      # exact?
	      if ( $exon1->start == $exon2->start && $exon1->end == $exon2->end ){
		  $exact++;
	      }


	      # calculate the extent of the overlap
	      if ( $exon1->start > $exon2->start && $exon1->start <= $exon2->end ){
		  if ( $exon1->end < $exon2->end ){
		      $overlap_length += ( $exon1->end - $exon1->start + 1);
		  }
		  elsif ( $exon1->end >= $exon2->end ){
		      $overlap_length += ( $exon2->end - $exon1->start + 1);
		  }
	      }
	      elsif( $exon1->start <= $exon2->start && $exon2->start <= $exon1->end ){
		  if ( $exon1->end < $exon2->end ){
		      $overlap_length += ( $exon1->end - $exon2->start + 1);
		  }
		  elsif ( $exon1->end >= $exon2->end ){
		      $overlap_length += ( $exon2->end - $exon2->start + 1);
		  }
	      }
	  }
      }
  }
  
  return ($overlaps,$overlap_length,$exact);
}    

#########################################################################
#
# this method just checks which UTR ends match between the two trancripts:
# ensembl genes are supposed to have UTRs most of them, so we probably need
# to check for the overlap with the first and the last exons of the transcript

sub _check_UTRMatches{
  # the order here is important
  my ($self,$ens_tran,$est_tran) = @_;
  #my $ens_translation = $ens_tran->translation;
  my @ens_exons       = @{$ens_tran->get_all_Exons};
  my $strand          = $ens_exons[0]->strand;

  # we only look at the start and end of the translations in the ensembl gene
  #my $ens_t_start_exon = $ens_translation->start_Exon;
  #my $ens_t_end_exon   = $ens_translation->end_Exon;
  my @est_exons = @{$est_tran->get_all_Exons};
  
  my $ens_start_exon = $ens_tran->start_Exon;
  my $ens_end_exon   = $ens_tran->end_Exon;

  # sorting paranoia
  if ( $strand == 1 ){
    @est_exons = sort { $a->start <=> $b->start } @est_exons;
    #@ens_exons = sort { $a->start <=> $b->start } @ens_exons;
  }
  if ( $strand == -1 ){
    @est_exons = sort { $b->start <=> $b->start } @est_exons;
    #@ens_exons = sort { $b->start <=> $b->start } @ens_exons;
  }
  my $match_3prime = 0;
  my $match_5prime = 0;

  my $ens_t_start_exon = $ens_start_exon;
  my $ens_t_end_exon   = $ens_end_exon;
  
 EST_EXON:
  for ( my $i=0; $i<= $#est_exons; $i++ ){
    
    # forward:
    if ( $strand == 1 ){
      
      # check first the 5' UTR
      if ( $est_exons[$i]->end == $ens_t_start_exon->end &&
	   $est_exons[$i]->start <= $ens_t_start_exon->start ){
	$match_5prime = 1;
      }
      # 3' UTR:
      if ( $est_exons[$i]->start == $ens_t_end_exon->start &&
	   $est_exons[$i]->end >= $ens_t_end_exon->end ){
	$match_3prime = 1;
      }
    }
    # reverse:
    elsif( $strand == -1){
      # 5' UTR:
      if ( $est_exons[$i]->start == $ens_t_start_exon->start &&
	   $est_exons[$i]->end >= $ens_t_start_exon->end ){
	$match_5prime = 1;
      }
      # 3' UTR
      if ( $est_exons[$i]->end == $ens_t_end_exon->end &&
	   $est_exons[$i]->start <= $ens_t_end_exon->start ){
	$match_3prime = 1;
      }
    }
  }  # end of EST_EXON
  
  return ( $match_5prime, $match_3prime);
}
#########################################################################
#
# we try to extend UTRs a bit more when possible
#

sub _extend_UTRs{
  # the order here is important
  my ($self,$ens_tran,$est_tran) = @_;

  # clone the translation
  my $ens_translation = new Bio::EnsEMBL::Translation;
  $ens_translation->start($ens_tran->translation->start);
  $ens_translation->end($ens_tran->translation->end);
  $ens_translation->start_Exon($ens_tran->translation->start_Exon);
  $ens_translation->end_Exon($ens_tran->translation->end_Exon);

  $ens_tran->sort;
  $est_tran->sort;
  
  my @ens_exons        = @{$ens_tran->get_all_Exons};
  my $strand           = $ens_exons[0]->strand;

  my $ens_t_start_exon = $ens_tran->start_Exon;
  my $ens_t_end_exon   = $ens_tran->end_Exon;

  my $ens_translation_start_exon = $ens_translation->start_Exon;
  my $ens_translation_end_exon   = $ens_translation->end_Exon;

  my @est_exons = @{$est_tran->get_all_Exons};
  
  # sorting paranoia
  if ( $strand == 1 ){
    @est_exons = sort { $a->start <=> $b->start } @est_exons;
    @ens_exons = sort { $a->start <=> $b->start } @ens_exons;
  }
  if ( $strand == -1 ){
    @est_exons = sort { $b->start <=> $b->start } @est_exons;
    @ens_exons = sort { $b->start <=> $b->start } @ens_exons;
  }

  my $modified = 0;
 EST_EXON:
  for ( my $i=0; $i<= $#est_exons; $i++ ){
    
  FORWARD:
    if ( $strand == 1 ){
      
      # check first the 5' UTR
      # if one est_exon has a coinciding end with the ens_exon where the translation starts
      if ( $est_exons[$i]->end == $ens_t_start_exon->end &&
	   $est_exons[$i]->start < $ens_t_start_exon->start ){
	
	# check that this est_exon does not start on top of another previous ens_exon
       	# we accept this situation:               Not this one:    
	#      __      __                                 __      __
	#     |__|----|__|--- ...   ens_tran             |__|----|__|--- ...
	#            ____                                   ________
	#           |____|--- ...   est_tran               |________|---  ...
	my $ok = 1;
	#foreach my $prev_exon ( @ens_exons ){
	#  if ( $prev_exon eq $ens_t_start_exon ){
	#    last;
	#  }
	#  if( $prev_exon->end >= $est_exons[$i]->start ){
	#    $ok = 0;
	#  }
	#}
	if ( $ok == 0 ){
	  # the est_exon starts on top of another exon, what should we do?
	}
	if ( $ok ){	  
	  # if this is the start of translation as well, need to modify translation
	  if ( $ens_t_start_exon == $ens_translation_start_exon ){
	    my $tstart = $ens_translation->start;
	    $tstart += ( $ens_t_start_exon->start - $est_exons[$i]->start );
	    $ens_translation->start($tstart);
	    #print STDERR "setting start of translation $ens_translation to $tstart\n";
	    $ens_t_start_exon->phase(-1);
	    $ens_translation->start_Exon($ens_t_start_exon);
	  }
	  
	  # modify the start coordinate
	  $ens_t_start_exon->start( $est_exons[$i]->start );

	  # add the est evidence to be able to see why we modified this exon
	  foreach my $evidence ( @{$est_exons[$i]->get_all_supporting_features} ){
	    $ens_t_start_exon->add_supporting_features( $evidence );
	  }
	  
	  # we need to add any possible extra UTR est_exons
	  $ens_tran = $self->_add_5prime_exons($ens_tran, $est_tran, $est_exons[$i], $i);
	  $modified = 1;
	}
      }
      
      # now check the 3' UTR
      if ( $est_exons[$i]->start == $ens_t_end_exon->start &&
	   $est_exons[$i]->end > $ens_t_end_exon->end ){
	
	# check that this est_exon does not end on top of another of the following ens_exon's
       	my $ok = 1;
	# we sort them in the opposite direction
	@ens_exons = sort {$b->start <=> $a->start} @ens_exons;
	#foreach my $next_exon ( @ens_exons ){
	#  if ( $next_exon eq $ens_t_end_exon ){
	#    last;
	#  }
	#  if( $next_exon->start <= $est_exons[$i]->end ){
	#    $ok = 0;
	#  }
	#}
	if ( $ok == 0 ){
	  # the est_exon end on top of another following exon, what should we do?
	}
	if ( $ok ){
	  if ( $ens_t_end_exon == $ens_translation_end_exon ){
	    $ens_t_end_exon->end_phase(-1);
	    $ens_translation->end_Exon( $ens_t_end_exon );
	    
	    # since the start coordinate does not change as we are in the forward strand
	    # we don't need to change the translation start/end
	  }
	  # modify the end coordinate
	  $ens_t_end_exon->end( $est_exons[$i]->end );

	  # add the est evidence to be able to see why we modified this exon
	  foreach my $evidence (@{ $est_exons[$i]->get_all_supporting_features} ){
	    $ens_t_end_exon->add_supporting_features( $evidence );
	  }
	  
	  # we need to add any possible extra UTR est_exons
	  $ens_tran = $self->_add_3prime_exons($ens_tran, $est_tran, $est_exons[$i], $i);
	  $modified = 1;
	}
      }
      
    } # end of strand == 1

  REVERSE:
    if ( $strand == -1 ){
      
      # check first the 5' UTR
      # if one est_exon has a coinciding end with the ens_exon where the translation starts
	if ( $est_exons[$i]->start == $ens_t_start_exon->start &&
	   $est_exons[$i]->end > $ens_t_start_exon->end ){
	    
        # check that this est_exon does not end on top an ens_exon on the right of $ens_t_start_exon, i.e.
       	# we accept a situation like:               but not like:
	#           __      __                             __      __
	#  ...  ---|__|----|__|   ens_tran         ... ---|__|----|__|
	#           ____                                   ________
	#       ---|____|         est_tran         ... ---|________|
	my $ok = 1;
	#foreach my $prev_exon ( @ens_exons ){
	#  if ( $prev_exon eq $ens_t_start_exon ){
	#    last;
	#  }
	#  if( $prev_exon->start <= $est_exons[$i]->end ){
	#    $ok = 0;
	#  }
	#}
	if ( $ok == 0 ){
	  # the est_exon ends on top of another ens_exon, what should we do?
	}
	if ( $ok ){
	  if ( $ens_t_start_exon == $ens_translation_start_exon ){
	    my $tstart = $ens_translation->start;
	    $tstart += ( $est_exons[$i]->end - $ens_t_start_exon->end );
	    $ens_translation->start($tstart);
	    $ens_t_start_exon->phase(-1);
	    $ens_translation->start_Exon($ens_t_start_exon);
	  }
	  # let's merge them
	  $ens_t_start_exon->end( $est_exons[$i]->end );
	  
	  # add the est evidence to be able to see why we modified this exon
	  foreach my $evidence ( @{$est_exons[$i]->get_all_supporting_features} ){
	    $ens_t_start_exon->add_supporting_features( $evidence );
	  }
	  
	  # we need to add any possible extra UTR est_exons
	  # but need to check possible incompatibilities with  previous ens_exons
	  $ens_tran = $self->_add_5prime_exons($ens_tran, $est_tran, $est_exons[$i], $i);
	  $modified = 1;
	}
      }
      
      # now check the 3' UTR
	if ( $est_exons[$i]->end == $ens_t_end_exon->end &&
	   $est_exons[$i]->start < $ens_t_end_exon->start ){
	  
	# check that this est_exon does not start on top of another ens_exon on the left of $ens_t_end_exon
       	my $ok = 1;
	# we sort them in the opposite direction
	@ens_exons = sort {$a->start <=> $b->start} @ens_exons;
	#foreach my $next_exon ( @ens_exons ){
	#  if ( $next_exon eq $ens_t_end_exon ){
	#    last;
	#  }
	#  if( $next_exon->end >= $est_exons[$i]->start ){
	#    $ok = 0;
	#  }
	#}
	if ( $ok == 0 ){
	  # the est_exon end on top of another following exon, what should we do?
	}
	if ( $ok ){
	  if ( $ens_t_end_exon == $ens_translation_end_exon ){
	    $ens_translation->end_Exon( $ens_t_end_exon );
	    $ens_t_end_exon->end_phase(-1);
	    # no need to change ens of translation as this is counted from the end of the
	    # translation->end_exon
	  }
	  # let's merge them
	  $ens_t_end_exon->start( $est_exons[$i]->start );
	  
	  # add the est evidence to be able to see why we modified this exon
	  foreach my $evidence ( @{$est_exons[$i]->get_all_supporting_features} ){
	    $ens_t_end_exon->add_supporting_features( $evidence );
	  }
	  
	  # we need to add any possible extra UTR est_exons
	  $ens_tran = $self->_add_3prime_exons($ens_tran, $est_tran, $est_exons[$i], $i);
	  $modified = 1;
	}
      }
      
    } # end of strand == -1
    
  }   # end of EST_EXON
  
  if ( $modified ){
    $self->_transcript_Type($ens_tran,'modified');
  }

  $ens_tran->translation($ens_translation);
  return ($ens_tran,$modified);
  
}


#########################################################################

sub _add_5prime_exons{ 

# the same idea as in Combine_Genewises_and_E2Gs.pm, we include the extra exons in the
# 5prime UTR region, previous to the one we've used to modify the
# exon. Maybe at the start of the translation in the ensembl transcript, 
# but also possibly just at the end of the transcript, as the ensembl transcript could have
# UTRs originally.
# 
# However, we also need to check that the exons we are going to add do not clash
# with other exons that may be already there. Thus
#
# we accept this situation:               Not this one:    
#              __                                  __      __
#             |__|---  ...   ens_tran             |__|----|__|--- ...
#     __     ____                                ___     ____
#    |__|---|____|---  ...   est_tran           |___|---|____|--- ...
#
# the reason so far is because I still haven't thought what to do with it
# ( we also accept when est_tran does not add any extra exons and ens has any exons before that one
#   but this of course does not change anything in the transcript )

  my ($self, $ens_tran, $est_tran, $est_exon_UTR, $est_exon_position ) = @_;
  my @est_exons = @{$est_tran->get_all_Exons};
  my $strand = $est_exons[0]->strand;
  
  # we have est_exons to add if $est_exon_position > 0
  if ( $est_exon_position > 0 ){

    # check whether there is any ens_exon previous to the est_exon_UTR 
    # (second case in the picture above)
    my $overlap = 0;
    if ( $strand == 1 ){
      my $start_range = $ens_tran->start_Exon->start;
      my $end_range   = $est_exon_UTR->start;
      foreach my $ens_exon ( @{$ens_tran->get_all_Exons} ){
	if ( $ens_exon->start >= $start_range && $ens_exon->end < $end_range ){
	  $overlap = 1;
	}
      }
    }
    # the same check but in the reverse strand
    if ( $strand == -1 ){
      my $start_range = $est_exon_UTR->end;
      my $end_range   = $ens_tran->start_Exon->end;
      foreach my $ens_exon ( @{$ens_tran->get_all_Exons} ){
	if ( $ens_exon->start > $start_range && $ens_exon->end <= $end_range ){
	  $overlap = 1;
	}
      }
    }
    if ( $overlap ){
      # forget it for the time being
    }
    else{
      # add the extra exons
      my $count = 0;
      while ( $count < $est_exon_position ){
	my $new_exon = $self->_clone_Exon( $est_exons[ $count ] );
	$new_exon->phase(-1);
	$new_exon->end_phase(-1);
	$ens_tran->add_Exon($new_exon);
	$ens_tran->sort;
	$count++;
      }
    }
  }
  return $ens_tran;
}
#########################################################################

sub _add_3prime_exons{ 

  # same as with the 5prime exons but left-right mirrored
  
  my ($self, $ens_tran, $est_tran, $est_exon_UTR, $est_exon_position ) = @_;
  my @est_exons = @{$est_tran->get_all_Exons};
  my $strand = $est_exons[0]->strand;
  my $overlap = 0;

  # we have est_exons to add if $est_exon_position < $#est_exons
  if ( $est_exon_position < $#est_exons ){

    # check whether there is any ens_exon previous to the est_exon_UTR 
    # (second case in the picture above)
    if ( $strand == 1){
      my $start_range = $est_exon_UTR->end;
      my $end_range   = $ens_tran->end_Exon->end;
      foreach my $ens_exon ( @{$ens_tran->get_all_Exons} ){
	if ( $ens_exon->start > $start_range && $ens_exon->end <= $end_range ){
	  $overlap = 1;
	}
      }
    }
    # the same check but in the reverse strand
    if ( $strand == -1 ){
      my $start_range = $ens_tran->end_Exon->start;
      my $end_range   = $est_exon_UTR->start;
      foreach my $ens_exon ( @{$ens_tran->get_all_Exons} ){
	if ( $ens_exon->start >= $start_range && $ens_exon->end < $end_range ){
	  $overlap = 1;
	}
      }
    }

    if ( $overlap ){
      # forget it for the time being
    }
    else{
      # add the extra exons
      my $count = $#est_exons;
      while ( $count > $est_exon_position ){
	my $new_exon = $self->_clone_Exon( $est_exons[ $count ] );
	$new_exon->phase(-1);
	$new_exon->end_phase(-1);
	$ens_tran->add_Exon($new_exon);
	$ens_tran->sort;
	$count--;
      }
    }
  }
  return $ens_tran;
}

#########################################################################
#
# METHODS INVOLVED IN WRITTING THE RESULTS
#
#########################################################################


# so far we just create one gene per transcript
# here we should include the funky clustering algorithm to do things properly

sub _make_Genes{
  my ( $self, $transcripts ) = @_;
  my @transcripts = @{ $transcripts };

  my $genetype = $FINAL_TYPE;
  my $analysis = $self->_analysis;

  my @selected_transcripts;
  my $count  = 0;
  
 TRANSCRIPT:
  foreach my $transcript (@transcripts) {
      $count++;
      unless (Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_check_Transcript($transcript,$self->ensembl_vc)
	      && Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_check_Translation($transcript) ){
	  print STDERR "skipping this transcript\n";
	  next TRANSCRIPT;
    }
    push(@selected_transcripts,$transcript);
  } 
  my @genes = $self->_cluster_into_Genes(@transcripts);
  
  foreach my $gene ( @genes ){
      $gene->type($genetype);
      $gene->analysis($analysis);
      foreach my $t (@{$gene->get_all_Transcripts}){
	  unless (  Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_check_Translation($t) ){
	      $self->warn("transcript with stop codons!");
	  }
      }
  }

  return @genes;
}

###################################################################

sub _remap_Genes {
  my ($self,$genes) = @_;
  my @genes = @$genes;
  my @new_genes;

  my $final_db  = $self->final_db; 
  my $final_gpa = $self->final_db->get_SliceAdaptor();
  my $chrid     = $self->input_id;
  if ( !( $chrid =~ s/\.(.*)-(.*)// ) ){
    $self->throw("Not a valid input_id... $chrid");
  }
  $chrid       =~ s/\.(.*)-(.*)//;
  my $chrstart = $1;
  my $chrend   = $2;
  my $final_vc = $final_gpa->fetch_by_chr_start_end($chrid,$chrstart,$chrend);
  my $genetype = $FINAL_TYPE;
  
 GENE:  
  foreach my $gene (@genes) {
    
    $gene->analysis($self->_analysis);
    $gene->type($genetype);
    my @trans = @{$gene->get_all_Transcripts};
    my $new_gene;
    # convert to raw contig coords
    eval {
      # transforming gene to raw contig coordinates.
      print STDERR "****************about to transform***********************\n";
      foreach my $transcript ( @{$gene->get_all_Transcripts} ){
	Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Evidence($transcript);
      }

      $new_gene = $gene->transform;
    };
    if ($@) {
      print STDERR "Couldn't reverse map gene [$@]\n";
      foreach my $t ( @{$gene->get_all_Transcripts} ){
	$self->_print_Transcript($t);
      }
      next;
    }
    $new_gene->type($gene->type);
    $new_gene->analysis($gene->analysis);  
    print STDERR "****************transformed gene***********************\n";
    foreach my $transcript ( @{$new_gene->get_all_Transcripts} ){
      Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Evidence($transcript);
    }
    push( @new_genes, $new_gene);
    
  }
  return @new_genes;
}





#########################################################################

sub _cluster_into_Genes{
  my ($self,@transcripts) = @_;
  
  @transcripts = sort by_transcript_high @transcripts;
  my @clusters;

  # clusters transcripts by exon overlap 
  foreach my $tran (@transcripts) {
    
    # store all clusters to which this transcript can match
    my @matching_clusters;
    my ($trans_start, $trans_end, $trans_strand) = get_transcript_start_end_strand($tran);
    
  CLUSTER: 
    foreach my $cluster (@clusters) {
      
      if (!($trans_start > $cluster->end || $trans_end < $cluster->start) &&
	  $trans_strand == $cluster->strand) {
	
	foreach my $cluster_transcript (@{$cluster->get_Transcripts}) {
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
    }
    
    if (scalar(@matching_clusters) == 0) {
      my $newcluster = new Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster;
      $newcluster->put_Transcripts($tran);
      push(@clusters,$newcluster);
    } 
    elsif (scalar(@matching_clusters) == 1) {
      $matching_clusters[0]->put_Transcripts($tran);
    } 
    else {
      # Merge the matching clusters into a single cluster
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
      @clusters =  @new_clusters;
    }
  }
  
  # safety and sanity checks
  $self->_check_Clusters(scalar(@transcripts), \@clusters);

  # @clusters is an array of arrayrefs, and each of those arrayrefs contain
  # the set of transcripts that have been clustered together

  # make and store genes
  my @genes;
  foreach my $cluster (@clusters){
    my $gene = new Bio::EnsEMBL::Gene;
    foreach my $transcript ( @{$cluster->get_Transcripts} ){
      $gene->add_Transcript($transcript);
    }
    
    # prune out duplicate exons
    $self->prune_Exons($gene);
    
    push( @genes, $gene );
  }
  
  return @genes;
}

#########################################################################


sub _check_Clusters{
    my ($self, $num_transcripts, $clusters) = @_;
    
    #Safety checks
    my $ntrans = 0;
    my %trans_check_hash;
    foreach my $cluster (@$clusters) {
	$ntrans += scalar(@{$cluster->get_Transcripts});
	foreach my $trans (@{$cluster->get_Transcripts}) {
	    if (defined($trans_check_hash{$trans})) {
		$self->throw("Transcript " . $trans->dbID . " added twice to clusters\n");
	    }
	    $trans_check_hash{$trans} = 1;
	}
	if (!scalar(@{$cluster->get_Transcripts})) {
	  $self->throw("Empty cluster");
	}
      }
    if ($ntrans != $num_transcripts) {
      $self->throw("Not all transcripts have been added into clusters $ntrans and " . $num_transcripts. " \n");
    } 
    return;
}

#########################################################################


sub by_transcript_high {
  my $alow;
  my $blow;
  my $ahigh;
  my $bhigh;

  if ($a->start_Exon->strand == 1) {
    $alow = $a->start_Exon->start;
    $ahigh = $a->end_Exon->end;
  } 
  else {
    $alow = $a->end_Exon->start;
    $ahigh = $a->start_Exon->end;
  }

  if ($b->start_Exon->strand == 1) {
    $blow = $b->start_Exon->start;
    $bhigh = $b->end_Exon->end;
  } 
  else {
    $blow = $b->end_Exon->start;
    $bhigh = $b->start_Exon->end;
  }

  if ($ahigh != $bhigh) {
    return $ahigh <=> $bhigh;
  } 
  else {
    return $alow <=> $blow;
  }
}

#########################################################################

sub get_transcript_start_end_strand {
  my ($transcript) = @_;
  my $start;
  my $end;
  
  my $start_exon = $transcript->start_Exon;
  my $end_exon = $transcript->end_Exon;
  
  if ($start_exon->strand == 1) {
    $start = $start_exon->start;
    $end   = $end_exon->end;
  } 
  else {
    $end   = $start_exon->end;
    $start = $end_exon->start;
  }
  return ($start, $end, $start_exon->strand);
}

#########################################################################

# when two exon-objects represent the same physical exon, we want to store just
# one of the objects, i.e. the transcript-objects will share the exon-object as well

sub prune_Exons {
  my ($self,$gene) = @_;
  my @unique_Exons; 
  
  # keep track of all unique exons found so far to avoid making duplicates
  # need to be very careful about translation->start_exon and translation->end_exon
  foreach my $tran (@{$gene->get_all_Transcripts}) {
    my @newexons;
    foreach my $exon (@{$tran->get_all_Exons}) {
      my $found;
    UNI:
      foreach my $uni (@unique_Exons) {
	if ($uni->start     == $exon->start  &&
	    $uni->end       == $exon->end    &&
	    $uni->strand    == $exon->strand &&
	    $uni->phase     == $exon->phase  &&
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
 } # end of TRANSCRIPT
  return;
}


#########################################################################

sub _analysis {
  my ($self, $analysis) = @_;

  if(defined $analysis){
    $self->throw("$analysis is not a Bio::EnsEMBL::Analysis") unless $analysis->isa("Bio::EnsEMBL::Analysis");
    $self->{'_analysis'} = $analysis;
  }

  return $self->{'_analysis'};
}

#########################################################################

# writes data into the db specified in ...


sub write_output {
  my ($self, @genes) = @_;
  unless (@genes){
    @genes = $self->output;
  }
  #print STDERR "about to write ".scalar(@genes)." into the db\n";
  
  # dbobj holds a reference to FINAL_DB
  my $gene_adaptor = $self->db->get_GeneAdaptor;
  
 GENE: 
  foreach my $gene (@genes) {	
    
    unless ( $gene->analysis ){
      $gene->analysis( $self->_analysis );
    }
    unless ( $gene->type ){
      $gene->type( $FINAL_TYPE );
    }
    
    print STDERR "about to write a gene\n";
    print STDERR "type    : ".$gene->type."\n";
    print STDERR "analysis: ".$gene->analysis->dbID." ".$gene->analysis->logic_name."\n";
    #foreach my $t ( @{$gene->get_all_Transcripts} ){
    #  Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Evidence($t);
    #}
    
    eval {
	$gene_adaptor->store($gene);
	print STDERR "wrote gene dbID " . $gene->dbID . "\n";
	foreach my $transcript ( @{$gene->get_all_Transcripts} ){
	Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Evidence($transcript);
      }
      
    }; 
    if( $@ ) {
      print STDERR "UNABLE TO WRITE GENE\n\n$@\n\nSkipping this gene\n";
      foreach my $transcript ( @{$gene->get_all_Transcripts} ){
	Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Evidence($transcript);
      }
     } 
}
}


########################################################################

=head2 _too_long_intron_size

Return: returns 1 if it finds an intron of size >  $ESTGENE_MAX_INTRON_LENGTH in the passed in transcript

=cut


sub _too_long_intron_size{
  my ($self,$transcript ) = @_;
  
  my $prev_exon;
  my $exon_count = 0;
  foreach my $exon ( @{$transcript->get_all_Exons} ){
    $exon_count++;

    if ( $exon_count > 1 ){
      
      my $intron_size = 0;
      if ($exon->strand == 1) {
	$intron_size = abs($exon->start - $prev_exon->end   - 1);
      } 
      else {
	$intron_size = abs($exon->end   - $prev_exon->start - 1);
      }
      
      if ($intron_size > $ESTGENE_MAX_INTRON_LENGTH ){ 
	print STDERR "intron size too long: $intron_size for transcript $transcript\n";
	return 1;
      }
    }
    $prev_exon = $exon;
  }
  return 0;
}


############################################################

sub _transfer_transcript_supporting_evidence{
    my ($self,$trans_source,$trans_target) = @_;
    
    my @exons_source  = @{$trans_source->get_all_Exons};
    my @exons_target  = @{$trans_target->get_all_Exons};
    
    # most of the time the numbers of exons doesn't vary
    if ( scalar( @exons_source ) == scalar ( @exons_target ) ){
	#print STDERR "passing evi info between 2 transcripts with same number of exons\n";
	while ( scalar ( @exons_source ) > 0 ){
	    my $exon_in  = shift( @exons_source  );
	    my $exon_out = shift( @exons_target );  
	    
	    # check just in case
	    if ( $exon_in->overlaps( $exon_out ) ){
	      Bio::EnsEMBL::Pipeline::Tools::ExonUtils->_transfer_supporting_evidence( $exon_in,$exon_out);
	    }
	    else{
		$self->warn("Trying to pass evidence between exons that do not overlap, this won't work!");
	    }
	}
    }
    else{
	# if not the same number of exons, we cannot know how the split happened
	print STDERR "passing evi info between 2 transcripts with different number of exons\n";
	foreach my $exon_in ( @exons_source ){
	    foreach my $exon_out( @exons_target ){
		if ( $exon_out->overlaps($exon_in) ){
		  Bio::EnsEMBL::Pipeline::Tools::ExonUtils->_transfer_supporting_evidence( $exon_in,$exon_out);
		}
	    }
	}
    }
}



############################################################

=head2 _transfer_supporting_evidence

 Title   : _transfer_supporting_evidence
 Usage   : $self->transfer_supporting_evidence($source_exon, $target_exon)
 Function: Transfers supporting evidence from source_exon to target_exon, 
           after checking the coordinates are sane and that the evidence is not already in place.
 Returns : nothing, but $target_exon has additional supporting evidence

=cut

sub _transfer_supporting_evidence{
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



1;
