=head1 NAME

GeneCluster

=head1 SYNOPSIS


=head1 DESCRIPTION

This object holds one or more genes which has been clustered according to 
comparison criteria external to this class (for instance, in the 
methods compare and _compare_Genes methods of the class GeneComparison).
Each GeneCluster object holds the IDs of the genes clustered and the beginning and end coordinates
of each one (taken from the start and end coordinates of the first and last exon in the correspondig
get_all_Exons array)

=head1 CONTACT

eae@sanger.ac.uk

=cut

# Let the code begin ...

package Bio::EnsEMBL::Pipeline::GeneComparison::GeneCluster;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster;

@ISA = qw(Bio::EnsEMBL::Root);

=head1 METHODS

=cut

#########################################################################


=head2 new()

new() initializes the attributes:

$self->{'_benchmark_types'}
$self->{'_prediction_types'}
$self->{'_benchmark_genes'}
$self->{'_prediction_genes'}

=cut

sub new {
  my ($class,$whatever)=@_;

  if (ref($class)){
    $class = ref($class);
  }
  my $self = {};
  bless($self,$class);

if ($whatever){
    $self->throw( "Can't pass an object to new() method. Use put_Genes() to include Bio::EnsEMBL::Gene in cluster");
  }

  $self->{'_ys_gene'}= {};
  %{$self->{'_statistics'}}=();
  
  return $self;
}

#########################################################################

=head2 put_Genes()

  function to include one or more genes in the cluster.
  Useful when creating a cluster. It takes as argument an array of genes, it returns nothing.

=cut

sub put_Genes {
  my ($self, @new_genes)= @_;
  if ( !defined( $self->{'_benchmark_types'} ) || !defined(  $self->{'_prediction_types'} ) ){
    $self->throw( "Cluster lacks references to gene-types, unable to put the gene");
  }

 GENE:
  foreach my $gene (@new_genes){
    foreach my $type ( @{ $self->{'_benchmark_types'} } ){
      if ($gene->type eq $type){
	push ( @{ $self->{'_benchmark_genes'} }, $gene );
	next GENE; 
      }
    }
    foreach my $type ( @{ $self->{'_prediction_types'} } ){
      if ($gene->type eq $type){
	push ( @{ $self->{'_prediction_genes'} }, $gene );
	next GENE;
      }
    }
  }
}

#########################################################################

=head2 get_Genes()

  it returns the array of genes in the GeneCluster object

=cut

sub get_Genes {
  my $self = shift @_;
  my @genes;
  if ( !defined( $self->{'_benchmark_genes'} ) && !defined( $self->{'_prediction_genes'} ) ){
    $self->warn("The gene array you try to retrieve is empty");
    @genes = ();
  }
  if ( $self->{'_benchmark_genes'} ){
    push( @genes, @{ $self->{'_benchmark_genes'} } );
  }
  if ( $self->{'_prediction_genes'} ){
    push( @genes, @{ $self->{'_prediction_genes'} } );
  }
  return @genes;
}

############################################################

sub strand{
  my $self = shift;
  my @genes = $self->get_Genes;
  unless (@genes){
    $self->warn("cannot retrieve the strand in a cluster with no genes");
  }
  my $strand;
  foreach my $gene (@genes){
    foreach my $transcript (@{$gene->get_all_Transcript}){
      unless (defined($strand)){
	$strand = $transcript->start_Exon->strand;
	next;
      }
      if ( $transcript->start_Exon->strand != $strand ){
	$self->throw("You have a cluster with genes on opposite strands");
      }
    }
  }
  return $strand;
}

#########################################################################

=head2 get_separated_Genes()

  Handy method to get the genes in the genes in the cluster separated by type.
  It returns two arrayrefs.

=cut


sub get_separated_Genes {
  my ($self) = @_;
  return ( $self->{'_benchmark_genes'}, $self->{'_prediction_genes'} );
}

#########################################################################

=head2 get_Gene_Count()

  it returns the number of genes in the GeneCluster object

=cut

sub get_Gene_Count {
  my $self = shift @_;
  my $count =0;
  if ( $self->{'_benchmark_genes'} ){
    $count += scalar( @{ $self->{'_benchmark_genes'}  } ); 
  }
  if ( $self->{'_prediction_genes'} ){
    $count += scalar( @{ $self->{'_prediction_genes'} } );
  }
  #print STDERR "In GeneCluster.get_Gene_Count(), Count = ".$count."\n";
  return $count;
}

#########################################################################

=head2 gene_Types()
  
  It accepts two array references to set the types. One array holds the gene-types for the 
  benchmark genes and the other on for the predicted genes. 
  It can also be used to get the two type-arrays: ($types1, $types2) = $cluster->gene_Types;
  The conventions throughout are (first entry: benchmark, second entry: prediction)

=cut

sub gene_Types {
  my ($self, $benchmark_types, $prediction_types) = @_;
  if ( $benchmark_types && $prediction_types ) {
    $self->{'_benchmark_types'}  = $benchmark_types;
    $self->{'_prediction_types'} = $prediction_types;
  }
  return ($self->{'_benchmark_types'},$self->{'_prediction_types'});
}

#########################################################################

=head2 get_Genes_of_Type()

  We can get the genes in each cluster of one type. 
  We pass one string identifying the genetype.
  The difference with get_Genes_by_Type is that that as an arrayref as argument.

=cut
  
sub get_Genes_of_Type() {
    my ($self,$type) = @_;
    unless ($type){
      $self->throw( "must provide a type");
    }
    my @genes = $self->get_Genes;  # this should give them in order, but we check anyway
    my @selected_genes;
    push ( @selected_genes, grep { $_->type eq $type } @genes );
    return @selected_genes;
  }

#########################################################################

=head2 get_Genes_by_Type()

  We can get the genes in each cluster of a given type. 
  We pass an arrayref containing the types we want to retrieve.

=cut
  
  sub get_Genes_by_Type() {
    my ($self,$types) = @_;
    unless ($types){
      $self->throw( "must provide a type");
    }
    my @genes = $self->get_Genes;  # this should give them in order, but we check anyway
    my @selected_genes;
    foreach my $type ( @{ $types } ){
      push ( @selected_genes, grep { $_->type eq $type } @genes );
    }
    return @selected_genes;
  }

#########################################################################

=head2 pair_Transcripts()

  Title   : pair_Transcripts()
  Usage   : my @array_of_pairs = $gene_cluster->pair_Transcripts
  Function: This method make pairs of transcripts according to maximum reciprocal exon overlap. 
  Example : look for instance in the method find_missing_Exons
  Returns : three arrayrefs =  
            1.- a list of Bio::EnsEMBL::Pipeline::GeneComparison::Transcripts, each holding a pair of transcripts, 
            2.- a list with the unpaired transcripts, and 
            3.- a list those transcripts which have been paired up twice or more
  Args    : nothing

=cut
  
sub pair_Transcripts {
  my ($self) = @_;
  
  # get the genes separated by type list 'benchmark'-like or 'prediction'-like
  my ( $ann_genes, $pred_genes ) = $self->get_separated_Genes;
  my (@ann_transcripts,@pred_transcripts);
  my (@ann_trans,@pred_trans);
  foreach my $gene ( @$pred_genes ){
    # print STDERR "gene ".$gene->type." put in pred_trans array\n";
    push( @pred_trans, @{$gene->get_all_Transcripts} );
  }
  foreach my $gene ( @$ann_genes ){
    # print STDERR "gene ".$gene->type." put in ann_trans array\n";
    push( @ann_trans, @{$gene->get_all_Transcripts} );
  }
  
  # tran1 are predicted genes
  # first sort the transcripts by their start position coordinate
  my %start_table_pred;
  my %start_table_ann;
  my $i=0;
  foreach my $tran ( @pred_trans ) {
    $start_table_pred{$i} = $tran->start_Exon->start;
    $i++;
  }
  my $j=0;
  foreach my $tran ( @ann_trans  ) {
    $start_table_ann{$j} = $tran->start_Exon->start;
    $j++;
  }
  foreach my $pos ( sort { $start_table_pred{$a} <=> $start_table_pred{$b} } keys %start_table_pred ){
    push (@pred_transcripts, $pred_trans[$pos]);
  }
  foreach my $pos ( sort { $start_table_ann{$a}  <=> $start_table_ann{$b}  } keys %start_table_ann ){
    push (@ann_transcripts, $ann_trans[$pos]);
  }

  # pair the transcripts, but first, some variable definition

  my %seen_pred;           # these keep track of those transcript linked and with how much overlap
  my %seen_ann;           # ditto, for @transcripts2
  my @pairs;           # list of (Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster) transcript-pairs being created 
  my @unpaired_ann;        # list of prediction Bio::EnsEMBL::Transcript which are left unpaired
  my @unpaired_pred;        # list of annotation Bio::EnsEMBL::Transcript which are left unpaired
  my @doubled;         # those which have been paired up twice
  my $overlap_matrix;  # matrix holding the number of exon overaps for each pair of transcripts
  my $link;            # matrix with 1 for the pairs linked and undef otherwise 
  my @overlap_pairs;   # each entry holds an array with the overlap and the two transcripts being compared
  my %repeated;        # to keep track of repeated transcripts

  # first calculate all possible overlaps
  foreach my $pred_tran ( @pred_transcripts ){
    foreach my $ann_tran ( @ann_transcripts ){
      my ($overlap_number,$overlap_length) = _compare_Transcripts( $pred_tran, $ann_tran );
      $$overlap_matrix{ $pred_tran }{ $ann_tran } = $overlap_number;
      my @list = ( $$overlap_matrix{ $pred_tran }{ $ann_tran }, $pred_tran, $ann_tran , $overlap_length);
      push ( @overlap_pairs, \@list );
      # print STDERR "Overlap( ".$ann_tran->stable_id.",".$pred_tran->stable_id." ) = "
      #.$$overlap_matrix{ $pred_tran }{ $ann_tran }."\n";
    }
  }
  # @overlap_pairs contains an array of lists
  # each list contains: ( number_of_exon_overlaps, prediction, annotation, length_of_overlap )
  
  if ( @overlap_pairs ){
    # sort the list of @overlap_pairs on the overlap
    my @sorted_pairs = sort { my $result = ( $$b[0] <=> $$a[0] );
			      if ($result){
				return $result;
			      }
			      else{
				return ( $$b[3] <=> $$a[3] );
			      }
			    } @overlap_pairs;
    
    #foreach my $pair ( @sorted_pairs ){
    #  print STDERR "pred: $$pair[1], ann: $$pair[2], overlaps: $$pair[0], length: $$pair[3]\n";
    #}


    # take the first pair of the list
    my $first = shift @sorted_pairs;
    my ($max_overlap,$pred_tran,$ann_tran,$max_overlap_length) =  @$first;
    $seen_pred{ $pred_tran } = $max_overlap;
    $seen_ann{ $ann_tran }   = $max_overlap;
    $$link{ $pred_tran }{ $ann_tran } = 1;
    # print STDERR "putting together $pred_tran and $ann_tran\n";
    
    # scan through each overlap
  PAIR:
    foreach my $list ( @sorted_pairs ){
      # each list contains @$list = ( overlap, transcript1, transcript2 , overlap_length)
      
      # first of all, if the overlap is zero, ditch it
      if ( $$list[0] == 0 ){
	next PAIR;
      }
      
      # if we've seen both transcripts already reject them
      if ( $$link{ $$list[1] }{ $$list[2] } && defined( $seen_pred{ $$list[1] } ) && defined( $seen_ann{ $$list[2] } ) ){
	next PAIR;
      }
      
      # if the same score...
      if ( $$list[0] == $max_overlap ) {
	
	# if we've seen both transcripts already, check they have the highest score
	#if ( defined( $seen1{ $$list[1] } ) && defined( $seen2{ $$list[2] } ) ){
	#  if ( $$list[0] == $seen1{ $$list[1] } && $$list[0] == $seen2{ $$list[2] } ){
	#    $$link{ $$list[1] }{ $$list[2] } = 1;
	#  }
	#  next PAIR;
	#}
	
	# if the pair is entirely new, we accept it
	if ( !defined( $seen_pred{ $$list[1] } ) && !defined( $seen_ann{ $$list[2] } ) ){
	  $$link{ $$list[1] }{ $$list[2] } = 1;
	  #print STDERR "putting together $$list[1] and $$list[2]\n";
	  $seen_pred{ $$list[1] } = $$list[0];
	  $seen_ann{ $$list[2] } = $$list[0];
	  next PAIR;
	}
	
	# we accept repeats only if this is their maximum overlap as well
	#if ( !defined( $seen_ann{$$list[2]} ) && defined( $seen_pred{$$list[1]} ) && $$list[0] == $seen_pred{$$list[1]} ){
	#  $$link{ $$list[1] }{ $$list[2] } = 1;
	#  $seen_ann{ $$list[2] } = $$list[0];
	#  if ( !defined( $repeated{ $$list[1] } ) ){
	#    push( @doubled, $$list[1] );
	#    $repeated{ $$list[1] } = 1;
	#  }
	#  next PAIR;
        #}
	#if ( !defined( $seen_pred{$$list[1]} ) && defined( $seen_ann{$$list[2]} ) && ($$list[0] == $seen_ann{$$list[2]}) ){ 
	#  $$link{ $$list[1] }{ $$list[2] } = 1;
	#  $seen_pred{ $$list[1] } = $$list[0];
	#  if ( !defined( $repeated{ $$list[2] } ) ){
	#    push( @doubled, $$list[2] );
	#    $repeated{ $$list[2] } = 1;
	#  }
	#  next PAIR;
	#}
      }
      
      # if the score is lower, only accept if the pair is completely new
      if ( $$list[0] < $max_overlap ){
	if ( !defined( $seen_pred{ $$list[1] } ) && !defined( $seen_ann{ $$list[2] } ) ){
	  $$link{ $$list[1] }{ $$list[2] } = 1;
	  $seen_pred{ $$list[1] } = $$list[0];
	  $seen_ann{  $$list[2] } = $$list[0];
	  $max_overlap = $$list[0];
	  #print STDERR "putting together $$list[1] and $$list[2]\n";
	  next PAIR;
	}
      }
    }
    
    # create a new cluster for each pair linked
    foreach my $pred_tran ( @pred_transcripts ){
      foreach my $ann_tran ( @ann_transcripts ){
	if ( $$link{ $pred_tran }{ $ann_tran } && $$link{ $pred_tran }{ $ann_tran } == 1 ){
	  my $pair = Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster->new();
	  $pair->put_Transcripts( $pred_tran, $ann_tran );
	  push( @pairs, $pair );
	}
      }
    }
    
    # finally, check for the unseen ones
    foreach my $pred_tran ( @pred_transcripts ){
      if ( !defined( $seen_pred{ $pred_tran } ) ){
	push( @unpaired_pred, $pred_tran );
      }
    }
    foreach my $ann_tran ( @ann_transcripts ){
      if ( !defined( $seen_ann{ $ann_tran } ) ){
	push( @unpaired_ann, $ann_tran );
      }
    }
    #print STDERR scalar(@pairs)." transcript pairs created\n";
    #my $count2=1;
    #foreach my $pair ( @pairs ){
    #  print STDERR "pair $count2:\n".$pair->to_String;
    #  $count2++;
    #}
    #$count2=1;
    #print STDERR scalar(@unpaired)." unpaired transcripts\n";
    #foreach my $unpaired ( @unpaired ){
    #  print STDERR "unpaired $count2: ".$unpaired->stable_id."\n";
    #}
    # return the pairs, the unpaired transcripts, and those transcript which have been taken twice or more
    
  }
  else{
    print STDERR "no pairs could be created from comparing ".scalar(@pred_transcripts). " predicted transcripts and "
      .scalar(@ann_transcripts)." annotated transcripts\n";
    @pairs = ();
    @unpaired_ann  = ();
    @unpaired_pred = ();
    @doubled  = ();
  }
  return (\@pairs,\@unpaired_ann,\@unpaired_pred);
}

#########################################################################
   

=head2 _compare_Transcripts()

 Title: _compare_Transcripts()
 Usage: this internal function compares the exons of two transcripts according to overlap
        and returns the number of overlaps
=cut

sub _compare_Transcripts {         
  my ($tran1, $tran2) = @_;
  my @exons1   = @{ $tran1->get_all_Exons };
  my @exons2   = @{ $tran2->get_all_Exons };
  my $overlaps = 0;
  my $overlap_length = 0;
  foreach my $exon1 (@exons1){
    foreach my $exon2 (@exons2){
      if ( $exon1->overlaps($exon2) && ( $exon1->strand == $exon2->strand ) ){
	$overlaps++;
	
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
  
  return ($overlaps,$overlap_length);
			
} 
	


#########################################################################

=head2 get_first_Gene()

  it returns the first gene in the cluster, which usually would be the benchmark gene 

=cut

sub get_first_Gene {
  my $self = shift @_;
  return @{$self->{'_benchmark_genes'}}[0];
}

#########################################################################

=head2 to_String()

  it returns a string containing the information about the genes contained in the
  GeneCluster object

=cut

sub to_String {
  my $self = shift @_;
  my $data='';
  foreach my $gene ( $self->get_Genes ){
    my @exons = @{ $gene->get_all_Exons };
     
    $data .= sprintf "Id: %-16s"      , $gene->stable_id;
    $data .= sprintf "Contig: %-20s"  , $exons[0]->contig->id;
    $data .= sprintf "Exons: %-3d"    , scalar(@exons);
    $data .= sprintf "Start: %-9d"    , $self->_get_start($gene);
    $data .= sprintf "End: %-9d"      , $self->_get_end  ($gene);
    $data .= sprintf "Strand: %-2d\n" , $exons[0]->strand;
  }
  return $data;
}

#########################################################################

=head2 _get_start()

 function to get the start position of a gene - it reads the gene object and it returns
 the start position of the first exon

=cut

sub _get_start {
  my ($self,$gene) = @_;
  my @exons = @{ $gene->get_all_Exons };
  my $st;
  
  if ($exons[0]->strand == 1) {
    @exons = sort {$a->start <=> $b->start} @exons;
    $st = $exons[0]->start;
  } else {
    @exons = sort {$b->start <=> $a->start} @exons; # they're read in opposite direction (from right to left)
    $st = $exons[0]->end;                           # the start is the end coordinate of the right-most exon
  }                                                 # which is here the first of the list of sorted @exons
  return $st;
}

#########################################################################

=head2 _get_end()

 function to get the end position of a gene - it reads the gene object and it returns
 the end position of the last exon

=cut

sub _get_end {
  my ($self,$gene) = @_;
  my @exons = @{ $gene->get_all_Exons };
  my $end;
  
  if ($exons[0]->strand == 1) {
    @exons = sort {$a->start <=> $b->start} @exons;
    $end = $exons[$#exons]->end;
  } else {
    @exons = sort {$b->start <=> $a->start} @exons; # they're read in opposite direction (from right to left)
    $end = $exons[$#exons]->start;                  # the end is the start coordinate of the left-most exon  
  }                                                 # which is here the last of the list @exons
  return $end;
}

#########################################################################

=head2 statistics()

 returns a hash containing the statistics of each gene in this cluster

=cut

sub statistics {
  my ($self,%stats) = @_;
  if (%stats){
    %{$self->{'_statistics'}}=%stats;
  }
  return  %{$self->{'_statistics'}};
}


#########################################################################

=head2 _translateable_exon_length()

 internal function that returns the length of the translateable exons 

=cut

sub _translateable_exon_length {
	my ($trans)= @_;
	my @exons = $trans->translateable_exons;
    my $length = 0;
    foreach my $ex (@exons) {
        $length += $ex->length;
    }
    return $length;

}

# method to get the start of the cluster, which we take to be the left_most exon_coordinate 
# i.e. the start coordinate of the first exon ordered as { $a->start <=> $b->start }, regardless of the strand

sub start{
  my ($self) = @_;
  my @genes = $self->get_Genes;
  my $start;
  foreach my $gene ( @genes ) {
    my @exons = @{ $gene->get_all_Exons};
    @exons = sort { $a->start <=> $b->start } @exons;
    my $this_start = $exons[0]->start;
    unless ( $start ){
      $start = $this_start;
    }
    if ( $this_start < $start ){
      $start = $this_start;
    }
  }
  return $start;
}
      
# method to get the end of the cluster, which we take to be the right_most exon_coordinate
# this being the end coordinate of the first exon ordered as { $b->end <=> $a->end }, regardless of the strand

sub end{
  my ($self) = @_;
  my @genes = $self->get_Genes;
  my $end;
  foreach my $gene ( @genes ) {
    my @exons = @{$gene->get_all_Exons};
    @exons = sort { $b->end <=> $a->end } @exons;
    
    # this is the largest end of all exons
    my $this_end = $exons[0]->end;
    unless ( $end ){
      $end = $this_end;
    }
    if ( $this_end > $end ){
      $end = $this_end;
    }
  }
  return $end;
}
      



1;
