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

use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::GeneComparison::GeneCluster;
use Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster;
use Bio::EnsEMBL::Pipeline::GeneCombinerConf;

# config file; parameters searched for here if not passed in as @args

use Bio::EnsEMBL::Pipeline::GeneCombinerConf qw(
						ESTGENE_DBHOST
						ESTGENE_DBUSER
						ESTGENE_DBNAME
						ESTGENE_DBPASS 
						ESTGENE_TYPE
						
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
					       );




@ISA = qw(Bio::Root::RootI Bio::EnsEMBL::Pipeline::RunnableDB);

######################################################################

sub new{
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  my $refdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
						  -host             => $REF_DBHOST,
						  -user             => $REF_DBUSER,
						  -dbname           => $REF_DBNAME,
						  -pass             => $REF_DBPASS,
						);
  
  my $ensembl_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
						      '-host'   => $ENSEMBL_DBHOST,
						      '-user'   => $ENSEMBL_DBUSER,
						      #'-pass'   => $ENSEMBL_DBPASS,
						      '-dbname' => $ENSEMBL_DBNAME,
						     );
  

  my $estgene_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
						      '-host'   => $ESTGENE_DBHOST,
						      '-user'   => $ESTGENE_REFDBUSER,
						      '-dbname' => $ESTGENE_DBNAME,
						      '-dnadb'  => $refdb,
						     ); 

  my $final_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
						      '-host'   => $FINAL_DBHOST,
						      '-user'   => $FINAL_REFDBUSER,
						      '-dbname' => $FINAL_DBNAME,
						      '-dnadb'  => $refdb,
						     ); 


  # needs to read from two databases and write into another one (possibly a third?)
  
  $self->ensembl_db( $ensembl_db );
  $self->estgene_db( $estgene_db );
  $self->final_db( $final_db );

  return $self;
  
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

sub final_db{
  my ( $self, $db ) = @_;
  if ( $db ){
    $db->isa("Bio::EnsEMBL::DBSQL::DBAdaptor") || $self->throw("Input [$db] is not a Bio::EnsEMBL::DBSQL::DBAdaptor");
    $self->{_final_db} = $db;
  }
  return $self->{_final_db};
}


#############################################################

sub fetch_input {
  my( $self) = @_;
  
  # get genomic region 
  my $chrid    = $self->input_id;
  if ( !( $ $chrid =~ s/\.(.*)-(.*)// ) ){
    $self->throw("Not a valid input_id... $chrid");
  }
  $chrid       =~ s/\.(.*)-(.*)//;
  my $chrstart = $1;
  my $chrend   = $2;
  print STDERR "Chromosome id = $chrid , range $chrstart $chrend\n";

  my $ensembl_gpa = $self->ensembl_db->get_StaticGoldenPathAdaptor();
  my $estgene_gpa = $self->estgene_db->get_StaticGoldenPathAdaptor();

  my $ensembl_vc  = $ensembl_gpa->fetch_VirtualContig_by_chr_start_end($chrid,$chrstart,$chrend);
  my $estgene_vc  = $estgene_gpa->fetch_VirtualContig_by_chr_start_end($chrid,$chrstart,$chrend);

  $self->ensembl_vc( $ensembl_vc );
  $self->estgene_vc( $estgene_vc );

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
  if ( @genes ){
    $genes[0]->isa("Bio::EnsEMBL::Gene") || $self->throw("$genes[0] is not a Bio::EnsEMBL::Gene");
    push ( @{ $self->{'_ensembl_genes'} }, @genes );
  }
  return @{ $self->{'_ensembl_genes'} };
}

#############################################################

sub estgenes{
  my ( $self, @genes ) = @_;
  if ( @genes ){
    $genes[0]->isa("Bio::EnsEMBL::Gene") || $self->throw("$genes[0] is not a Bio::EnsEMBL::Gene");
    push ( @{ $self->{'_estgenes'} }, @genes );
  }
  return @{ $self->{'_estgenes'} };
}
  
############################################################
#
# RUN METHOD
#
############################################################

sub run{
  my ($self,@args) = @_;

  # get ensembl genes (from GeneBuilder)
  $self->ensembl_genes( $self->ensembl_vc->get_Genes_by_Type( $ENSEMBL_TYPE, 'evidence' ) );
  
  # get estgenes ( from EST_GeneBuilder )
  $self->estgenes( $self->estgene_vc->get_Genes_by_Type( $ESTGENE_TYPE, 'evidence' ) ); 
  
  # cluster estgenes and ensembl genes
  my @genes             = ( $self->ensembl_genes, $self->estgenes );
  my @gene_clusters     = $self->cluster_Genes( @genes );
  my @unclustered_genes = $self->unclustered_Genes;
  
  # the accepted set of transcripts
  my @transcripts;
 
  # at this stage we could separate the clusters into three sorts:

  # 1) only ensembl genes --> we leave them as they are
  # 2) ensembl+est genes  -->  first: extend UTRs
  #                           second: include alternative forms (do checks) 
  # 3) only est genes     --> if only 1 in the cluster, take it if good coverage 
  #                       --> if more: take a set that share exons among themselves (exact matches)
  #                                    not necessarily all the same exon, it should be the
  #                                    largest connected graph where vertices are trnascripts and
  #                                    edges is the relation 'share an exon'
 
 CLUSTER:
  foreach my $cluster ( @gene_clusters ){
    
    # get genes of each type
    my @ens_genes = $cluster->get_Genes_of_Type( $GB_FINAL_GENETYPE );
    my @est_genes = $cluster->get_Genes_of_Type( $EST_GENOMEWISE_GENETYPE );
    
    # if we have genes of either type, let's try to match them
    if ( @ens_genes && @est_genes ){
      my @ens_transcripts;
      my @est_transcripts;
      foreach my $gene ( @ens_genes ){
	push ( @ens_transcripts, $gene->each_Transcript );
      }
      foreach my $gene ( @est_genes ){
	push ( @est_transcripts, $gene->each_Transcript );
      }
      my ( $new_ens, $accepted_est ) = $self->_pair_Transcripts( \@ens_transcripts, \@est_transcripts );
      if ( $new_ens ){
	push ( @transcripts, @$new_ens );
      }
      if ( $accepted_est ){
	push ( @transcripts, @$accepted_est );
      }
    }
    
    # else we could have only ensembl genes
    elsif(  @ens_genes && !@est_genes ){
      # we have nothing to modify them, hence we accept them...
      foreach my $gene ( @ens_genes ){
	push ( @transcripts, $gene->each_Transcript );
      }
    }

    # else we could have only est genes
    elsif( !@ens_genes && @est_genes ){
      my @est_transcripts;
      foreach my $gene ( @est_genes ){
	push ( @est_transcripts, $gene->each_Transcript );
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
  
  # make the genes ( this goes through a clustering as in GeneBuilder )
  my @newgenes = $self->_make_Genes(\@transcripts);

  
  my @remapped  = $self->_remap_Genes(\@new_genes);

  # store genes
  $self->output(@remapped);

}


############################################################
#
# METHODS CALLED FROM RUN METHOD... DOING ALL THE MAGIC
#
############################################################

# this method cluster genes according to genomic extent
# covered by the genes. The proper clustering of transcripts
# to give rise to genes occurs in _make_Genes()

sub cluster_Genes{
  my ($self) = @_;

  my @genes = ( $self->ensembl_genes, $self->estgenes );
  
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
  
  # put the first gene into these cluster
  $cluster->put_Genes( $sorted_genes[0] );
  push (@clusters, $cluster);
  
  # loop over the rest of the genes
 LOOP:
  for (my $c=1; $c<=$#sorted_genes; $c++){
    my $found=0;
    
    # treat the clusters as ranges, so we only need to check if ranges overlap
    # for the moment this is enough
    my $gene_start = $self->get_start_of_Gene( $sorted_genes[$c] );
    my $gene_end   = $self->get_end_of_Gene( $sorted_genes[$c] );
    
    # we need to do this each time, so that start/end get updated
    my $cluster_start = $cluster->start;
    my $cluster_end   = $cluster->end;

    if ( !( $gene_end < $cluster_start || $gene_start > $cluster_end ) ){
      $cluster->put_Genes( $sorted_genes[$c] );
    }
    else{
      # else, create a new cluster
      $cluster = new Bio::EnsEMBL::Pipeline::GeneComparison::GeneCluster; 
      $cluster->put_Genes( $sorted_genes[$c] );
      $cluster_count++;
      push( @clusters, $cluster );
    }
  }

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
  if (@genes){
    push( @{ $self->{'_unclustered_genes'} }, @genes );
  }
  return @{ $self->{'_unclustered_genes'} };
}

#########################################################################

# this gives the left-most exon coordinate in a gene

sub _get_start_of_Gene{  
  my ($self,$gene) = @_;
  my @exons        = $gene->get_all_Exons;
  @exons           = sort{ $a->start <=> $b->start } @exons;
  my $start        = $exons[0]->start;
  return $start;
}

#########################################################################

# this gives the right-most exon coordinate in a gene

sub _get_end_of_Gene{  
  my ($self,$gene) = @_;
  my @exons        = $gene->get_all_Exons;
  @exons           = sort{ $b->end <=> $a->end } @exons;
  my $end          = $exons[0]->end;
  return $end;
}

#########################################################################

sub _pair_Transcripts{
  my ($self,$ens_transcripts,$est_transcripts) = @_;
  my @potential_isoforms;

  my @ens_transcripts = @$ens_transcripts;
  my @est_transcripts = @$est_transcripts;
  my %used_est_transcript;
  my @accepted_isoforms;
  
  # sort the transcripts by the number of exons in descending order
  # for equal number of exons, order according to start (left most) coordinate
  @ens_transcripts = sort { my $result = ( scalar( $b->get_all_Exons ) <=> scalar( $a->get_all_Exons ) );
			    if ($result){ 
			      return $result;
			    }
			    else{
			      return ( $self->_get_start_of_Transcript($a) 
				       <=>
				       $self->_get_start_of_Transcript($b) 
				     );
			    }
			  } @ens_transcripts;
  
  
  @est_transcripts = sort { my $result = ( scalar( $b->get_all_Exons ) <=> scalar( $a->get_all_Exons ) );
			    if ($result){
			      return $result;
			    }
			    else{
			      return ( $self->_get_start_of_Transcript($a) 
				       <=>
				       $self->_get_start_of_Transcript($b) 
				     );
			    }
			  } @est_transcripts;
  

  ###### LOOK FOR ESTs THAT EXTEND UTRs

  # matrix holding the number and length of exon overlap for each pair of transcripts
  my $overlap_number_matrix;
  my $overlap_length_matrix;
  
  my %merge_list;
  my %selected; # keeps track of the est-transcripts that have been used

  # we base everything on the ensembl transcripts
 ENS_TRANSCRIPT:
  foreach my $ens_tran ( @ens_transcripts ){
      my @list;
      
    EST_TRANSCRIPT:
      foreach my $est_tran ( @est_transcripts ){
	  
	  # check with which est_transcripts it can merge
	  my ($merge,$overlaps) = $self->_test_for_Merge( $ens_tran, $est_tran );
	  if ($merge == 1 && $overlaps > 0 ){
	      
	      # calculate then how much they overlap
	      my ($overlap_number,$overlap_length, $exact_matches) = _compare_Transcripts( $ens_tran, $est_tran );
	      
	      push( @list, [ $overlap_number, $overlap_length, $exact_matches, $est_tran ]);
	  }
      }
      # take the est_transcript which overlap the most
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
			    } @list;
      
      # try to merge it, if succeeded, mark the chosen one and leave the rest
      
    MODIFY:
      for (my $i = 0; $i< scalar(@list); $i++){
	  my $est_tran = $$list[$i][3];
	  my $modified;
	  ( $ens_tran, $modified ) = $self->_extend_UTRs($ens_tran,$est_tran);
	  if ($modified == 1 ){
	      $used_est_transcript{$est_tran} = 1;
	      
	      # we only allow one cdna to modify each transcript
	      last MODIFY;
	  }
      }
      
  }   # end of ENS_TRANSCRIPT
  
  ###### LOOK FOR POSSIBLE ALTERNATIVE TRANSCRIPTS
  # those which haven't been used for extending UTRs are candidate isoforms
  #
  my @candidates;
  foreach my $est_tran ( @est_transcripts ){
      unless (  $used_est_transcript{ $est_tran } && $used_est_transcript{$est_tran} == 1 ){
	  push ( @candidates, $est_tran);
      }
  }
  if (@candidates){
 
      # at the moment we just check whether there is an ens_transcript
      # that shares one exon, part of the protein product and
      # there is one intron matching one exon or vice versa.

      foreach my $est_tran ( @candidates ){
	  my $one_exon_shared = 0;
	  my $similar_protein = 0;
	  my $exon_in_intron  = 0;
	  
	  foreach my $ens_tran ( @ens_transcripts ){
	      $one_exon_shared = $self->_check_exact_exon_Match( $est_tran, $ens_tran);
	      $similar_protein = $self->_check_protein_Match(    $est_tran, $ens_tran);
	      $exon_in_intron  = $self->_check_alternative_Form( $est_tran, $ens_tran);
	      
	      if ( $one_exon_shared == 1 && 
		   $similar_protein == 1 &&
		   $exon_in_intron  == 1 ){
		  push ( @accepted_isoforms, $est_tran );
	      }
	  }
      }      
  }  # end of 'if (@candidates)'
  else{
      # nothing to do then
  }
  
  # return the ens_transcripts (modified or not)
  return ( @ens_transcripts, @accepted_isoforms );
  
}
#########################################################################

# this function takes est_transcripts that have been clustered together
# but not with any ensembl transcript and tries to figure out whether they
# make an acceptable set of alt-forms

sub _check_est_Cluster{
  my ($self,@est_transcripts) = @_;
  my %color;

  # adjacency lists:
  my %adj;
  
  for(my $i=1;$i<=scalar(@est_transcripts);$i++){
    for(my $j=1;$j<=scalar(@est_transripts);$j++){
      
      next if $j==$i;
      if ( $self->_check_exact_exon_Match( $est_transcripts[$i], $est_transcripts[$j]) &&
	   $self->_check_protein_Match(    $est_transcripts[$i], $est_transcripts[$j])    ){
	  push ( @{ $adj{$est_transcripts[$i]} } , $est_transcripts[$j] );
      }
    }
  }
  
  foreach my $tran ( @est_transcripts ){
    $color{$trans} = "white";
  }
  
  my @potential_genes;
  
  # find the connected components doing a depth-first search
  foreach my $tran ( @est_transcripts ){
    if ( $color{$trans} eq 'white' ){
      my $potential_gene = [];
      $self->_visit( $tran, \%color, \%adj, $potential_gene);
      push ( @potential_genes, $potential_gene );
    }
  }
  print STDERR scalar(@potential_genes)." genes created\n";
  
  # take only the sets with more than one transcript:
  @potential_genes = sort { scalar( @{ $b } ) <=> scalar( @{ $a } ) } @potential_genes;
  my @accepted;
  my $first_gene = shift @potential_genes;
  while ( scalar( @{ $first_gene } ) > 0 ){
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
    push( @{ $potential_gene }, $node );
  }
  $color->{ $node } = 'black';    
  #$self->_print_Transcript($node);
  return;
}

#########################################################################
  
sub _print_Transcript{
  my ($self,$transcript) = @_;
  print STDERR "transcript: $transcript\n";
  foreach my $exon ( $transcript->get_all_Exons ){
    print $exon->start."-".$exon->end." ";
  }
  print STDERR "\n";
}


#########################################################################

# having a set of est_genes only, we have no reference transcript (ensembl one),
# so to determine whether two transcripts are two alternative forms of the same gene
# we check whether they share at least an exon

sub _check_exact_exon_Match{
 my ($self, $tran1, $tran2 ) = @_;
 my @exons1 = $tran1->get_all_Exons;
 my @exons2 = $tran2->get_all_Exons;
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

 my $compatible_protein = 0;
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

sub _check_alternative_Form{
    my ($self,$tran1,$tran2) = @_;
    my @exons1 = $tran1->get_all_Exons;
    my @exons2 = $tran2->get_all_Exons;	
    
    # we actually compare ranges:
    my @exon_list1;
    my @intron_list1;
    for(my $i=0; $i<scalar(@exons1);$i++){
	my $exon_range = new Bio::Range();
	$exon_range->start($exons1[$i]->start);
	$exon_range->end($exons1[$i]->end);
	push( @exon_list1, $exon_range);
	if ( $i+1 < scalar(@exons1)){
	    my $intron_range = new Bio::Range();
	    $intron_range->start($exons1[$i]->end + 1);
	    $intron_range->end($exons1[$i+1]->start - 1);
	    push( @intron_list1, $intron_range);
	}
    }
    my @exon_list2;
    my @intron_list2;
    for(my $i=0; $i<scalar(@exons2);$i++){
	my $exon_range = new Bio::Range();
	$exon_range->start($exons2[$i]->start);
	$exon_range->end($exons2[$i]->end);
	push( @exon_list2, $exon_range);
	if ( $i+1 < scalar(@exons2)){
	    my $intron_range = new Bio::Range();
	    $intron_range->start($exons2[$i]->end + 1);
	    $intron_range->end($exons2[$i+1]->start - 1);
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
		
		# and the previous or the next exon overlaps one of the flanking exons
		if ( ( $i>0            && 
		       $exon_list1[$i-1]->overlaps( $exon_list2[$j] ) ) 
		     ||
		     ( $i<$#exon_list1 &&
		       $j<$#exon_list2 &&
		       $exon_list1[$i+1]->overlaps( $exon_list2[$j+1] ) )
		     ){
		    return 1;
		}
	    }
	}
    }
	
    for (my $i=0; $i<scalar(@exon_list2);$i++){
	for ( my $j=0; $j<scalar(@intron_list1);$j++){
	    
	    # if one exon overlaps one intron
	    if ( $exon_list2[$i]->overlaps( $intron_list1[$j] ) ){
		
		# and the previous or the next exon overlaps one of the flanking exons
		if ( ( $i>0            && 
		       $exon_list2[$i-1]->overlaps( $exon_list1[$j] ) ) 
		     ||
		     ( $i<$#exon_list2 &&
		       $j<$#exon_list1 &&
		       $exon_list2[$i+1]->overlaps( $exon_list1[$j+1] ) )
		     ){
		    return 1;
		}
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
  my @exons1 = $tran1->get_all_Exons;
  my @exons2 = $tran2->get_all_Exons;	
 
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

# this gets the left most start coordinate for the transcript, regardless of the strand

sub _get_start_of_Transcript{
  my ($self,$transcript) = @_;
  my @exons = $transcript->get_all_Exons;
  @exons    = sort { $a->start <=> $b->start } @exons;
  my $start = $exons[0]->start;
  
  return $start;
}

#########################################################################
   
# this compares both transcripts and calculate the number of overlapping exons,
# the length of the overlap, and the number of exact overlapping exons

sub _compare_Transcripts {         
  my ($tran1, $tran2) = @_;
  my @exons1   = $tran1->get_all_Exons;
  my @exons2   = $tran2->get_all_Exons;

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

sub _extend_UTRs{
  # the order here is important
  my ($self,$ens_tran,$est_tran) = @_;
  my $ens_translation = $ens_tran->translation;
  my @ens_exons       = $ens_tran->get_all_Exons;
  my $strand          = $ens_exons[0]->strand;

  # we only look at the start and end of the translations in the ensembl gene
  my $ens_t_start_exon = $ens_translation->start_exon;
  my $ens_t_end_exon   = $ens_translation->end_exon;
  
  my @est_exons = $est_tran->get_all_Exons;
  
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
  for ( my $i=0; $i< $#est_exons; $i++ ){
    
  FORWARD:
    if ( $strand == 1 ){
      
      # check first the 5' UTR
      # if one est_exon has a coinciding end with the ens_exon where the translation starts
      if ( $est_exons[$i]->end == $ens_t_start_exon->end &&
	   $est_exons[$i]->start <= $ens_t_start_exon->start ){

	# check that this est_exon does not start on top of another previous ens_exon
       	# we accept this situation:               Not this one:    
	#      __      __                                 __      __
	#     |__|----|__|--- ...   ens_tran             |__|----|__|--- ...
	#            ____                                   ________
	#           |____|--- ...   est_tran               |________|---  ...
	my $ok = 1;
	foreach my $prev_exon ( @ens_exons ){
	  if ( $prev_exon eq $ens_t_start_exon ){
	    last;
	  }
	  if( $prev_exon->end >= $est_exons[$i]->start ){
	    $ok = 0;
	  }
	}
	if ( $ok == 0 ){
	  # the est_exon starts on top of another exon, what should we do?
	}
	if ( $ok ){
	  # let's merge them
	  my $tstart = $ens_translation->start;
	  $tstart += ( $ens_t_start_exon->start - $est_exons[$i]->start );
	  $ens_t_start_exon->start( $est_exons[$i]->start );
	  $ens_tran->translation($tstart);
	  
	  # add the est evidence to be able to see why we modified this exon
	  foreach my $evidence ( $est_exons[$i]->each_Supporting_Feature ){
	    $ens_t_start_exon->add_Supporting_Feature( $evidence );
	  }

	  # we need to add any possible extra UTR est_exons
	  $ens_tran = $self->_add_5prime_exons($ens_tran, $est_tran, $est_exons[$i], $i);
	  $modified = 1;
	}
      }

      
      # now check the 3' UTR
      if ( $est_exons[$i]->start == $ens_t_end_exon->start &&
	   $est_exons[$i]->end >= $ens_t_end_exon->end ){
	
	# check that this est_exon does not end on top of another of the following ens_exon's
       	my $ok = 1;
	# we sort them in the opposite direction
	@ens_exons = sort {$b->start <=> $a->start} @ens_exons;
	foreach my $next_exon ( @ens_exons ){
	  if ( $next_exon eq $ens_t_end_exon ){
	    last;
	  }
	  if( $next_exon->start <= $est_exons[$i]->end ){
	    $ok = 0;
	  }
	}
	if ( $ok == 0 ){
	  # the est_exon end on top of another following exon, what should we do?
	}
	if ( $ok ){
	  # let's merge them
	  $ens_t_end_exon->end( $est_exons[$i]->end );
	  $ens_translation->end_exon( $ens_t_end_exon );
	  # since the start coordinate does not change as we are in the forward strand
	  # we don't need to change the translation start/end
	  
	  # add the est evidence to be able to see why we modified this exon
	  foreach my $evidence ( $est_exons[$i]->each_Supporting_Feature ){
	    $ens_t_end_exon->add_Supporting_Feature( $evidence );
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
	   $est_exons[$i]->end >= $ens_t_start_exon->end ){
	
	# check that this est_exon does not end on top an ens_exon on the right of $ens_t_start_exon, i.e.
       	# we accept a situation like:               but not like:
	#           __      __                             __      __
	#  ...  ---|__|----|__|   ens_tran         ... ---|__|----|__|
	#           ____                                   ________
	#       ---|____|         est_tran         ... ---|________|
	my $ok = 1;
	foreach my $prev_exon ( @ens_exons ){
	  if ( $prev_exon eq $ens_t_start_exon ){
	    last;
	  }
	  if( $prev_exon->start <= $est_exons[$i]->end ){
	    $ok = 0;
	  }
	}
	if ( $ok == 0 ){
	  # the est_exon ends on top of another ens_exon, what should we do?
	}
	if ( $ok ){
	  # let's merge them
	  my $tstart = $ens_translation->start;
	  $tstart += ( $est_exons[$i]->end - $ens_t_start_exon->end );
	  $ens_t_start_exon->end( $est_exons[$i]->end );
	  $ens_tran->translation($tstart);
	  
	  # add the est evidence to be able to see why we modified this exon
	  foreach my $evidence ( $est_exons[$i]->each_Supporting_Feature ){
	    $ens_t_start_exon->add_Supporting_Feature( $evidence );
	  }
	  
	  # we need to add any possible extra UTR est_exons
	  # but need to check possible incompatibilities with  previous ens_exons
	  $ens_tran = $self->_add_5prime_exons($ens_tran, $est_tran, $est_exons[$i], $i);
	  $modified = 1;
	}
      }
            
      # now check the 3' UTR
      if ( $est_exons[$i]->end == $ens_t_end_exon->end &&
	   $est_exons[$i]->start <= $ens_t_end_exon->start ){
	
	# check that this est_exon does not start on top of another ens_exon on the left of $ens_t_end_exon
       	my $ok = 1;
	# we sort them in the opposite direction
	@ens_exons = sort {$a->start <=> $b->start} @ens_exons;
	foreach my $next_exon ( @ens_exons ){
	  if ( $next_exon eq $ens_t_end_exon ){
	    last;
	  }
	  if( $next_exon->end >= $est_exons[$i]->start ){
	    $ok = 0;
	  }
	}
	if ( $ok == 0 ){
	  # the est_exon end on top of another following exon, what should we do?
	}
	if ( $ok ){
	  # let's merge them
	  $ens_t_end_exon->start( $est_exons[$i]->start );
	  $ens_translation->end_exon( $ens_t_end_exon );
	  # no need to change ens of translation as this is counted from the end of the
	  # translation->end_exon

	  # add the est evidence to be able to see why we modified this exon
	  foreach my $evidence ( $est_exons[$i]->each_Supporting_Feature ){
	    $ens_t_end_exon->add_Supporting_Feature( $evidence );
	  }
	  
	  # we need to add any possible extra UTR est_exons
	  $ens_tran = $self->_add_3prime_exons($ens_tran, $est_tran, $est_exons[$i], $i);
	  $modified = 1;
	}
      }
	
    } # end of strand == -1
    
    
  }   # end of EST_EXON

  return ($ens_tran,$modified);

}
	    
#########################################################################

sub _add_5prime_exons{ 

# the same idea as in Combine_Genewises_and_E2Gs.pm, we include the extra exons in the
# 5prime UTR region, previous to the one we've used to modify the
# exon at the start of the translation in the ensembl transcript. 
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
  my @est_exons = $est_tran->get_all_Exons;
  my $strand = $est_exons[0]->strand;
  
  # we have est_exons to add if $est_exon_position > 0
  if ( $est_exon_position > 0 ){

    # check whether there is any ens_exon previous to the est_exon_UTR 
    # (second case in the picture above)
    my $overlap = 0;
    if ( $strand == 1 ){
      my $start_range = $ens_tran->start_exon->start;
      my $end_range   = $est_exon_UTR->start;
      foreach my $ens_exon ( $ens_tran->get_all_exons ){
	if ( $ens_exon->start >= $start_range && $ens_exon->end < $end_range ){
	  $overlap = 1;
	}
      }
    }
    # the same check but in the reverse strand
    if ( $strand == -1 ){
      my $start_range = $est_exon_UTR->end;
      my $end_range   = $ens_tran->start_exon->end;
      foreach my $ens_exon ( $ens_tran->get_all_exons ){
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
  my @est_exons = $est_tran->get_all_Exons;
  my $strand = $est_exons[0]->strand;
  my $overlap = 0;

  # we have est_exons to add if $est_exon_position < $#est_exons
  if ( $est_exon_position < $#est_exons ){

    # check whether there is any ens_exon previous to the est_exon_UTR 
    # (second case in the picture above)
    if ( $strand == 1){
      my $start_range = $est_exon_UTR->end;
      my $end_range   = $ens_tran->end_exon->end;
      foreach my $ens_exon ( $ens_tran->get_all_exons ){
	if ( $ens_exon->start > $start_range && $ens_exon->end <= $end_range ){
	  $overlap = 1;
	}
      }
    }
    # the same check but in the reverse strand
    if ( $strand == -1 ){
      my $start_range = $ens_tran->end_exon->start;
      my $end_range   = $est_exon_UTR->start;
      foreach my $ens_exon ( $ens_tran->get_all_exons ){
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
	$ens_tran->add_Exon($new_exon);
	$ens_tran->sort;
	$count--;
      }
    }
  }
  return $ens_tran;
}

#########################################################################

sub _clone_Exon{
  my ($self,$old_exon) = @_;
  my $new_exon = new Bio::EnsEMBL::Exon;
  $new_exon->start( $old_exon->start );
  $new_exon->end( $old_exon->end );
  $new_exon->strand( $old_exon->strand );
  $new_exon->phase( $old_exon->phase ); # this should be -1 as is in UTR
  $new_exon->contig_id( $old_exon->contig_id );
  $new_exon->attach_seq($self->vc);
  
  # transfer supporting evidence
  foreach my $evidence ( $old_exon->each_Supporting_Feature ){
    $new_exon->add_Supporting_Feature( $evidence );
  }
  return $new_exon;
}


#########################################################################
#
# GET/SET METHODS 
#
#########################################################################

# ... put them here ...



#########################################################################
#
# METHODS INVOLVED IN WRITTING THE RESULTS
#
#########################################################################


# so far we just create one gene per transcript
# here we should include the funky clustering algorithm to do things properly

sub _make_Genes{
  my ( $self, @transcripts ) = @_;

  my $genetype = $FINAL_TYPE;

  # sort out analysis here or we will get into trouble with duplicate analyses
  my $ana_Adaptor = $self->dbobj->get_AnalysisAdaptor;
  my @analyses = $ana_Adaptor->fetch_by_logic_name($genetype);
  my $analysis_obj;

  if(scalar(@analyses) > 1){
    $self->throw("panic! > 1 analysis for $genetype\n");
  }
  elsif(scalar(@analyses) == 1){
    $analysis_obj = $analyses[0];
  }
  else{
    # make a new analysis object
    $analysis_obj = new Bio::EnsEMBL::Analysis( -db              => 'NULL',
					        -db_version      => 1,
					        -program         => $genetype,
					        -program_version => 1,
					        -gff_source      => $genetype,
					        -gff_feature     => 'gene',
					        -logic_name      => $genetype,
					        -module          => 'GeneCombiner',
					      );
  }

  $self->_analysis($analysis_obj);
  my @genes;
  my $contig = $self->vc;
  my $count  = 0;

 TRANSCRIPT:
  foreach my $transcript (@transcripts) {
    $count++;
    
    # create a new gene object out of this transcript
    my $gene   = new Bio::EnsEMBL::Gene;
    $gene->type($genetype);
    $gene->temporary_id($contig->id . ".$genetype.$count");
    $gene->analysis($self->_analysis);

    # add transcript to gene
    $transcript->temporary_id($contig->id . ".$genetype.$count");
    $gene->add_Transcript($transcript);

    # sort the exons 
    $transcript->sort;
    my @exons = $transcript->get_all_Exons;

    my $translation = $transcript->translation;
    $translation->temporary_id($contig->id . ".$genetype.$count");
  
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
      my $five_prime  = $transcript->five_prime_utr  or warn "No five prime UTR";
      my $three_prime = $transcript->three_prime_utr or warn "No three prime UTR";

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
      if ( $peptide =~ /\*/ ){
	print STDERR "TRANSLATION HAS STOP CODONS!!\n";
      }
      else{
	push(@genes,$gene);
      }
      print STDERR "\n";
    }
  
  } # end of TRANSCRIPT

  my @clusters = $self->_cluster_into_Genes(@genes);

  my @newgenes;

  make genes here....



  return @newgenes;
}

###################################################################

sub _remap_Genes {
  my ($self,$genes) = @_;
  my @genes = @$genes;
  my @new_genes;
  my $contig = $self->vc;
  
GENE:  
  foreach my $gene (@genes) {
    
 #   eval {
#      my $genetype = $gene->type;
#      my $new_gene = $contig->convert_Gene_to_raw_contig($gene);
      
#      # need to explicitly add back genetype and analysis.
#      $new_gene->type($genetype);
#      $new_gene->analysis($gene->analysis);
      
#      # DO WE REALLY NEED TO DO THIS???
#      # sort out supporting feature coordinates
#      foreach my $tran ($new_gene->each_Transcript) {
#	foreach my $exon($tran->get_all_Exons) {
#	  foreach my $sf ($exon->each_Supporting_Feature) {
	    
#	    # this should be sorted out by the remapping to rawcontig ... strand is fine
#	    if ($sf->start > $sf->end) {
#	      my $tmp = $sf->start;
#	      $sf->start($sf->end);
#	      $sf->end($tmp);
#	    }
#	  }
#	}
#      }

#      # is this a special case single coding exon gene with UTRS?
#      if($tran->translation->start_exon() eq $tran->translation->end_exon() 
#	 && $gene->type eq 'est_combined'){
#	#	print STDERR "single coding exon, with UTRs\n";
	
#	# problems come about when we switch from + strand on FPC contig to - strand on raw contig.
#	my $fpc_strand;
	
#	foreach my $exon($tran->get_all_Exons) {
#	  if ($exon eq $tran->translation->start_exon()) {
#	    $fpc_strand = $exon->strand;
#	    last;
#	  }
#	}
	
#	foreach my $tran ($new_gene->each_Transcript) {
#	  foreach my $exon($tran->get_all_Exons) {
	    
#	    # oh dear oh dear oh dear
#	    # this is still giving some problems
#	    print STDERR "oh dear, strand problems with single-conding-exon gene\n";
#	    if ($exon eq $tran->translation->start_exon()) {
#	      if($fpc_strand == 1 && $exon->strand == -1){
#		print STDERR "fpc strand 1, raw strand -1 - flipping translation start/end\n";
#		$exon->end($exon->end - ($tran->translation->start -1));
#		$exon->start($exon->end - ($tran->translation->end -1));
#	      }
#	    }
#	  }
#	}
	
#      } # end special case single coding exon
      
#      # final exon coord sanity check
#      foreach my $exon($new_gene->get_all_Exons){
#	# make sure we deal with stickies!
#	if($exon->isa("Bio::EnsEMBL::StickyExon")){
#	  foreach my $ce($exon->each_component_Exon){
#	    # exon start and end must both be within the raw contig!!!
#	    if($ce->start < 1){
#	      $self->throw("can't set exon->start < 1 (" . $ce->start . ") - discarding gene\n");
#	    }
	    
#	    if($ce->end > $ce->contig->primary_seq->length){
#	      $self->throw("exon extends beyond end of contig - discarding gene\n");
#	    }
#	  }
#	}
#	else{
#	  # regular exon
#	  # exon start and end must both be within the raw contig!!!
#	  if($exon->start < 1){
#	    $self->throw("can't set exon->start < 1 (" . $exon->start . ") - discarding gene\n");
#	  }
	  
#	  if($exon->end > $exon->contig->primary_seq->length){
#	    $self->throw("exon extends beyond end of contig - discarding gene\n");
#	  }
#	}
#      }
#      # if we get to here, the gene is fine, so push it onto the array to be returned
#      push( @new_genes, $new_gene );
      
#    };
    
    # did we throw exceptions?
    if ($@) {
      print STDERR "Couldn't reverse map gene:  [$@]\n";
    }
  }
  
  return @new_genes;
}





#########################################################################

sub _cluster_into_Genes{
  my ($self,$genes2cluster) = @_;
  my $num_old_genes = 0;
  my @genes2cluster = @$gene2cluster;
  my @transcripts_unsorted;
  
  foreach my $gene ( @$gene2cluster ){
    $num_old_genes++;
    foreach my $tran ($gene->each_Transcript) {
      push(@transcripts_unsorted, $tran);
    }
  }
  
  my @transcripts = sort by_transcript_high @transcripts_unsorted;
  my @clusters;

  # clusters transcripts by exon overlap 
  foreach my $tran (@transcripts) {
    
    # store all clusters to which this transcript can match
    my @matching_clusters;
  CLUSTER: 
    foreach my $cluster (@clusters) {
      foreach my $cluster_transcript (@$cluster) {
        foreach my $exon1 ($tran->get_all_Exons) {
	  
          foreach my $cluster_exon ($cluster_transcript->get_all_Exons) {
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
      push (@{ $matching_clusters[0] }, $tran);
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
  $self->check_Clusters(scalar(@transcripts), $num_old_genes, \@clusters);

  # @clusters is an array of arrayrefs, and each of those arrayrefs contain
  # the set of transcripts that have been clustered together

  # check the version by Steve where he uses TranscriptCluster objects, etc...

  return @clusters;
}

sub check_Clusters{
  my ($self, $num_transcripts, $num_old_genes, $clusters) = @_;
  
  #Safety checks
  my $ntrans = 0;
  my %trans_check_hash;
  foreach my $cluster (@$clusters) {
    $ntrans += scalar(@$cluster);
    foreach my $trans (@$cluster) {
      if (defined($trans_check_hash{"$trans"})) {
        $self->throw("Transcript " . $trans->dbID . " added twice to clusters\n");
      }
      $trans_check_hash{"$trans"} = 1;
    }
    if (!scalar(@$cluster)) {
      $self->throw("Empty cluster");
    }
  }
  if ($ntrans != $num_transcripts) {
    $self->throw("Not all transcripts have been added into clusters $ntrans and " . $num_transcripts. " \n");
  } 
  #end safety checks
  
  if (scalar(@$clusters) < $num_old_genes) {
    $self->warn("Reclustering reduced number of genes from " . 
		$num_old_genes . " to " . scalar(@$clusters). "\n");
  } elsif (scalar(@$clusters) > $num_old_genes) {
    $self->warn("Reclustering increased number of genes from " . 
		$num_old_genes . " to " . scalar(@$clusters). "\n");
  }

}

sub by_transcript_high {
  my $alow;
  my $blow;
  my $ahigh;
  my $bhigh;

  if ($a->start_exon->strand == 1) {
    $alow = $a->start_exon->start;
    $ahigh = $a->end_exon->end;
  } 
  else {
    $alow = $a->end_exon->start;
    $ahigh = $a->start_exon->end;
  }

  if ($b->start_exon->strand == 1) {
    $blow = $b->start_exon->start;
    $bhigh = $b->end_exon->end;
  } 
  else {
    $blow = $b->end_exon->start;
    $bhigh = $b->start_exon->end;
  }

  if ($ahigh != $bhigh) {
    return $ahigh <=> $bhigh;
  } 
  else {
    return $alow <=> $blow;
  }
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
  my($self) = @_;
  
  my $gene_adaptor = $self->dbobj->get_GeneAdaptor;
  
 GENE: 
  foreach my $gene ($self->output) {	
    # do a per gene eval...
    eval {
      $gene_adaptor->store($gene);
      print STDERR "wrote gene dbID " . $gene->dbID . "\n";
    }; 
    if( $@ ) {
      print STDERR "UNABLE TO WRITE GENE\n\n$@\n\nSkipping this gene\n";
    }
  }
}

########################################################################

1;
