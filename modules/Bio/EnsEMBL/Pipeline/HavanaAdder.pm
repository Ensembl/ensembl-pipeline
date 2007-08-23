#
#
# BioPerl module for GeneBuilder
#
# Cared for by EnsEMBL <ensembl-dev@ebi.ac.uk>
#
# You may distribute this module under the same terms as perl itself
#
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

This module reads your favourite annotations (genewise, combined_genes, est2genome, genomewise, ... )
on the one hand, and ab initio predictions plus features on the other hand. Ab initio predictions and features
are passed to Bio::EnsEMBL::Pipeline::Runnable::PredictionGeneBuilder which generates putative transcripts
from supported prediction exons (see documentation in that module for details).
The product of Bio::EnsEMBL::Pipeline::Runnable::PredictionGeneBuilder is combined with
all the other annotations and redundant transcripts are eliminated in the method prune_Transcripts().
The resulting transcripts are combined into genes. For more details, follow the list of methods called
by build_Genes() method and the description in each one.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX


The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::HavanaAdder;

use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::Attribute;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster;

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::General     qw (
							       GB_INPUTID_REGEX
							      );

use Bio::EnsEMBL::Pipeline::Config::HavanaAdder qw (
                                                    GB_ENSEMBL_INPUT_GENETYPE
                                                    GB_HAVANA_INPUT_GENETYPE
                                                    GB_HAVANA_INPUT_TRANSCRIPTTYPES
                                                    GB_MERGED_TRANSCRIPT_TYPE
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

    $self->throw("Must input a slice to HavanAdder") unless defined($slice);
    $self->{_final_genes} = [];
    $self->{_gene_types}  = [];

    $self->query($slice);
    $self->gene_types($GB_ENSEMBL_INPUT_GENETYPE);
    $self->gene_types($GB_HAVANA_INPUT_GENETYPE);
  
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
  my @all_transcripts = $self->combined_Transcripts;
  
  # do a preliminary clustering
  my @preliminary_genes = $self->cluster_into_Genes(@all_transcripts);

  # merge redundant ensembl transcripts which match a havana one
  $self->_merge_redundant_transcripts(\@preliminary_genes);

  # make shared exons unique objects
  my @genes =  $self->_make_shared_exons_unique( @preliminary_genes );
  
  print STDERR scalar(@genes)." genes built\n";
  
  $self->final_genes( @genes );
}


sub _merge_redundant_transcripts{
  my ($self, $genes) = @_;

 GENE:
  foreach my $gene(@$genes){
    my @transcripts = @{$gene->get_all_Transcripts};
    my @havana;
    my @ensembl;

    
    # are there any havana transcripts?
    foreach my $transcript(@transcripts){
      my $havana_t = 0;

      foreach my $htranscript (@{$GB_HAVANA_INPUT_TRANSCRIPTTYPES}){
        if($transcript->biotype eq  $htranscript){
        #  print "IM HERE\n";
          push(@havana, $transcript);
          $havana_t = 1;
        }
      }
      if($havana_t == 0){
       # print "IM with ENSEMBL\n";
        push(@ensembl, $transcript);
      }
    }
    if (!scalar(@havana)){
#      print "No havana transcripts to consider\n";
      next GENE;
    }
    

# coul dthis lead to trying to remove a transcript >once?

    # compare each havana transcript to each ensembl one
    foreach my $ht(@havana){
      print "Deleting havana transcript supporting features\n";
      $ht->flush_supporting_features;
      foreach my $et(@ensembl){
        my $delete_t = 0;
        if ($delete_t = $self->are_matched_pair($ht, $et)){
          my @t_pair = ($ht, $et);
          $self->set_transcript_relation($delete_t, @t_pair);
          $ht->biotype($GB_MERGED_TRANSCRIPT_TYPE);
          $et->biotype($GB_MERGED_TRANSCRIPT_TYPE);
          $self->_remove_transcript_from_gene($gene, $delete_t);          
        }
      }
    }
  }

}

sub are_matched_pair {
  my($self, $havana, $ensembl) = @_;


  # Fetch all exons in each transcript 
  my @hexons = @{$havana->get_all_Exons};
  my @eexons = @{$ensembl->get_all_Exons};

  my @thexons = @{$havana->get_all_translateable_Exons};
  my @teexons = @{$ensembl->get_all_translateable_Exons};


#  print "_________________\n";
  # Fetch only coding exons
  #$self->clear_coding_exons_cache;
  #my @thexons = @{get_coding_exons_for_transcript($havana)};
  #my @teexons = @{get_coding_exons_for_transcript($ensembl)};

  
  # Check that the number of exons is the same in both transcripts
  return 0 unless scalar(@hexons) == scalar(@eexons);


  #print "____________________________________\n";
  #print "HAVANA ID: ",$havana->dbID, " ENSEMBL: ",$ensembl->dbID,"\n";

 # Check return 3 different possible values:
 # 0 means keep both transcript
 # return ($ensembl) means keep havana transcript and remove ensembl 
 # return ($havana) means keep ensembl transcript and remove hanana

  # double check translation coords
  #print "HAVANA TRANS START: ",$havana->translation->genomic_start()," END: ",$havana->translation->genomic_end,"\n";
  #print "ENSEMBL TRANS START: ",$ensembl->translation->genomic_start," END: ",$ensembl->translation->genomic_end,"\n";

  return 0 unless($havana->translation->genomic_start == $ensembl->translation->genomic_start);
  return 0 unless($havana->translation->genomic_end   == $ensembl->translation->genomic_end);



  # special case for single exon genes
  if(scalar(@hexons ==1)){
    print "SINGLE EXONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
 
    if ($hexons[0]->start     == $eexons[0]->start &&
        $hexons[0]->end       == $eexons[0]->end &&
        $hexons[0]->strand    == $eexons[0]->strand
        ){
      # Both are exactly the same so we delete Ensembl one
      return $ensembl;
      
    }elsif($hexons[0]->start     <= $eexons[0]->start &&
           $hexons[0]->end       >= $eexons[0]->end &&
           $hexons[0]->strand    == $eexons[0]->strand &&
           $eexons[0]->start     == $teexons[0]->coding_region_start($ensembl) &&
           $eexons[0]->end       == $teexons[0]->coding_region_end($ensembl)
           ){
      # Ensembl gene don't have UTR and Havana has then delete Ensembl one
      return $ensembl;
      
    }elsif((($hexons[0]->start    != $eexons[0]->start ||
            $hexons[0]->end       != $eexons[0]->end) &&
            $hexons[0]->strand    == $eexons[0]->strand) &&
            ($eexons[0]->start    != $teexons[0]->coding_region_start($ensembl) ||
            $eexons[0]->end       != $teexons[0]->coding_region_end($ensembl))
           ){
      # Both contain UTR keep ENSEMBL
      return $havana;
      
    }else{
      print "Keep single both\n";
      return 0;
      
    }
  }
  # if is a multi exons transcript
  else{
    # First we check the internal structure of the transcript where everything has to be exactly equal
    print "CHECKING INTERNAL EXONS \n";
    for(my $i=1; $i<=($#hexons-1); $i++){
      return 0 unless ($hexons[$i]->start     == $eexons[$i]->start &&
                       $hexons[$i]->end       == $eexons[$i]->end &&
                       $hexons[$i]->strand    == $eexons[$i]->strand &&
                       $hexons[$i]->phase     == $eexons[$i]->phase &&
                       $hexons[$i]->end_phase == $eexons[$i]->end_phase);
    }
    print "INTERNAL EXONS ARE OK \n";
    # Then check the first an last exon to check if they are the same. If just start and end of UTR are different keep ensembl one

    # CASE 1: Both coding and UTR are the same, keep Havana and delete Ensembl
    if ($hexons[0]->start     == $eexons[0]->start &&
        $hexons[0]->end       == $eexons[0]->end &&
        $hexons[0]->strand    == $eexons[0]->strand &&
        $hexons[-1]->start    == $eexons[-1]->start &&
        $hexons[-1]->end      == $eexons[-1]->end &&
        $hexons[-1]->strand   == $eexons[-1]->strand 
        ){
      print "MULTIEXON DELETE ENSEMBL\n";
      return $ensembl;
      
    }elsif (#CASE 2": HAVANA HAS UTR AND ENSEMBL DOESNT, KEEP HAVANA
            $hexons[0]->strand == 1 &&
            $hexons[0]->end       == $eexons[0]->end &&
            $hexons[0]->strand    == $eexons[0]->strand &&
            $hexons[-1]->start    == $eexons[-1]->start &&
            $hexons[-1]->strand   == $eexons[-1]->strand &&
            $eexons[0]->start     == $teexons[0]->coding_region_start($ensembl) &&
            $eexons[-1]->end      == $teexons[-1]->coding_region_end($ensembl) &&
            ($hexons[-1]->end     != $eexons[-1]->end ||
             $hexons[0]->start    != $eexons[0]->start)  
            ){
      print "MULTIEXON DELETE ENSEMBL\n";      
      return $ensembl;
      
    }elsif (# CASE 3: BOTH ENSEMBL AND HAVANA HAVE UTR KEEP ENSEMBL
            $hexons[0]->strand == 1 &&
            $hexons[0]->end       == $eexons[0]->end &&
            $hexons[0]->strand    == $eexons[0]->strand &&
            $hexons[-1]->start    == $eexons[-1]->start &&
            $hexons[-1]->strand   == $eexons[-1]->strand &&
            ($eexons[0]->start    != $teexons[0]->coding_region_start($ensembl) ||
            $eexons[-1]->end      != $teexons[-1]->coding_region_end($ensembl)) &&
            ($hexons[-1]->end     != $eexons[-1]->end ||
            $hexons[0]->start     != $eexons[0]->start)
            ){
      print "MULTIEXON DELETE HAVANA\n";      
      return $havana;
      
    }elsif (# CASE 4: Same as case 2 but in reverse strand
            $hexons[0]->strand == -1 &&
            $hexons[0]->start     == $eexons[0]->start &&
            $hexons[0]->strand    == $eexons[0]->strand &&
            $hexons[-1]->end      == $eexons[-1]->end &&
            $hexons[-1]->strand   == $eexons[-1]->strand &&
            $eexons[-1]->start    == $teexons[-1]->coding_region_start($ensembl) &&
            $eexons[0]->end       == $teexons[0]->coding_region_end($ensembl) &&
            ($hexons[0]->end      != $eexons[0]->end ||
             $hexons[-1]->start   != $eexons[-1]->start)
            
            ){
      print "MULTIEXON DELETE ENSEMBL\n";      
      return $ensembl;
      
    }elsif (# CASE 5: Same as case 3 but in reverse strand
            $hexons[0]->strand == -1 &&
            $hexons[0]->start     == $eexons[0]->start &&
            $hexons[0]->strand    == $eexons[0]->strand &&
            $hexons[-1]->end      == $eexons[-1]->end &&
            $hexons[-1]->strand   == $eexons[-1]->strand &&
            ($eexons[-1]->start   != $teexons[-1]->coding_region_start($ensembl) ||
            $eexons[0]->end       != $teexons[0]->coding_region_end($ensembl)) &&
            ($hexons[0]->end      != $eexons[0]->end ||
             $hexons[-1]->start   != $eexons[-1]->start)

            ){
      print "MULTIEXON DELETE HAVANA\n";      
      return $havana;
      
    }else{

      print "Keep MULTIEXON BOTH\n";
      return 0;
      
    }
    
  }
  
  print " WEIRD CASE NOT WE DID NOT THOUGHT ABOUT, CHECK RULES!\n";
  return 0;
  
}
sub set_transcript_relation {
  
  my($self, $delete_t, @t_pair) = @_;
  
 # print " To delete: ",$delete_t," havana ", $t_pair[0],"\n";
  
  # If transcript to delete is Havana we create an xref for the entry say that the transcript is CDS equal to ENSEMBL
  if ($delete_t == $t_pair[0]){
    # transfer OTT ID and/or ENST
    foreach my $entry(@{ $t_pair[0]->get_all_DBEntries}){
      if ($entry->dbname eq 'Vega_transcript'){
        my $newentry = new Bio::EnsEMBL::DBEntry
            (
              -primary_id => $entry->primary_id,
              -display_id => $entry->display_id,
              -priority => 1,
              -xref_priority => 0,
              -version => 1,
              -release => 1,
              -dbname => 'shares_CDS_with'
              );
        
        $newentry->status("XREF");
        
        $t_pair[1]->add_DBEntry($newentry);
      }
    }
    
    # We add a transcript attribute to the ensembl transcript with the start and end coords of the Havana transcript that we will delete
    my $attrib_value = $t_pair[0]->slice->coord_system_name.":".$t_pair[0]->slice->coord_system->version.":".$t_pair[0]->slice->seq_region_name.":".
                       $t_pair[0]->start.":".$t_pair[0]->end.":1";
   # print "ATTRIB VALUE:---------- ",$attrib_value,"\n";
    my $attribute = Bio::EnsEMBL::Attribute->new
       (-CODE => 'TranscriptEdge',
        -NAME => 'Transcript Edge',
        -DESCRIPTION => '',
        -VALUE => $attrib_value);

   $t_pair[1]->add_Attributes($attribute);

    #$t_pair[1]->get_all_supporting_features;
    # When we delete a Havana transcript we want to transfer the exon supporting features to the transcript we keep    
    my @delete_e = @{$delete_t->get_all_Exons};
    my @exons    = @{$t_pair[1]->get_all_Exons};
    
    if (scalar(@delete_e) == scalar(@exons)){
      # Single exon genes
      if (scalar(@delete_e) == 1){
        $self->transfer_supporting_evidence($delete_e[0], $exons[0]); 
      }else{
        #  print "IM GETTING HERE\n";
        my $e;
        for ($e = 0, $e<scalar(@delete_e), $e++){
          if($delete_e[$e] && $exons[$e]){
            $self->transfer_supporting_evidence($delete_e[$e], $exons[$e]); 
          }
        }
      }
    }
  }else{
    # If the transcript to delete is ENSEMBL we add an xref entry say that both transcripts are exact matches (including UTR)
    foreach my $entry(@{ $t_pair[0]->get_all_DBEntries}){
      if ($entry->dbname eq 'Vega_transcript'){
        my $enstentry = new Bio::EnsEMBL::DBEntry
            ( 
              -primary_id => $entry->primary_id,
              -display_id => $entry->display_id,
              -version => 1,
              -release => 1,
              -priority => 1,
              -xref_priority => 0,
              -dbname => 'shares_CDS_and_UTR_with'
              );
        
        $enstentry->status("XREF");
        
        $t_pair[0]->add_DBEntry($enstentry);
      }
    } 
  # Transfer the supporting features both for transccript and exon of the transcript to delete to the transcript we keep
    $self->transfer_supporting_features($delete_t,$t_pair[0]);
  }
}

sub transfer_supporting_features{
  my ($self, $delete_t, $transcript) = @_;
  
  print "TRANSCRIPT IS:::: ", $transcript,"\n";
  
  my @exon_features;

  # Delete all the supporting features for the Havana Transcript 
  #$transcript->flush_supporting_features;
  
  my @delete_tsf = @{ $delete_t->get_all_supporting_features };
  #my @transcript_sf = @{ $transcript->get_all_supporting_features };
  
  # print "NUMBER OF TRANSCRIPT SF: ",scalar(@transcript_sf),"\n";
  # print " AND DELETE TSF: ", scalar(@delete_tsf),"\n";
 DTSF: foreach my $dtsf (@delete_tsf){
   next DTSF unless $dtsf->isa("Bio::EnsEMBL::FeaturePair");
   #TSF: foreach my $tsf (@transcript_sf){
   #  next TSF unless $tsf->isa("Bio::EnsEMBL::FeaturePair");
   
   # Check that the supporting feature to transfer does not already exist
   #if($dtsf->start    == $tsf->start &&
   #   $dtsf->end      == $tsf->end &&
   #   $dtsf->strand   == $tsf->strand &&
   #   $dtsf->hseqname eq $tsf->hseqname &&
   #   $dtsf->hstart   == $tsf->hstart &&
   #   $dtsf->hend     == $tsf->hend){
   
   #  print STDERR "feature already exists\n";
   #  next DTSF;
   # Now check that the supporting feature is not longer (to avoid write_gene check to fail)
   #}elsif($dtsf->start   <= $tsf->start &&
   #       $dtsf->end     >= $tsf->end &&
   #       $dtsf->strand  == $tsf->strand){
   
   #print STDERR "feature already\n";
   #  next DTSF;
   #}
   #}
   
   # Transfer supporting feature if we reach this point
   #  print "TRANSFErring SUpporting Feature\n";
   $transcript->add_supporting_features($dtsf);
 }
  
  my @delete_e = @{$delete_t->get_all_Exons};
  my @exons    = @{$transcript->get_all_Exons};
  
  if (scalar(@delete_e) == scalar(@exons)){
    # Single exon genes
    if (scalar(@delete_e) == 1){
     $self->transfer_supporting_evidence($delete_e[0], $exons[0]); 
    }else{
      #  print "IM GETTING HERE\n";
      my $e;
      for ($e = 0, $e<scalar(@delete_e), $e++){
        if($delete_e[$e] && $exons[$e]){
          $self->transfer_supporting_evidence($delete_e[$e], $exons[$e]); 
        }
      }
    }
  }
  
#  print "NUMBER AFT ADDITTION: ", scalar(@{ $transcript->get_all_supporting_features }),"\n";
}

sub _remove_transcript_from_gene {
  my ($self, $gene, $trans_to_del)  = @_;
  
  my @newtrans;
  foreach my $trans (@{$gene->get_all_Transcripts}) {
    if ($trans != $trans_to_del) {
      push @newtrans,$trans;
    }
  }

# The naughty bit!
  $gene->{_transcript_array} = [];

  foreach my $trans (@newtrans) {
    $gene->add_Transcript($trans);
  }

  return scalar(@newtrans);
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


=head2 get_Genes

 Description: retrieves ensembl and havana gene annotations with supporting evidence. 
 ReturnType : none, but $self->combined_Transcripts is filled
 Args       : none

=cut

sub get_Genes {
  my ($self) = @_;
  my @transcripts;
 
  my $ensemblslice = $self->fetch_sequence($self->input_id, $self->ensembl_db);
  my $havanaslice = $self->fetch_sequence($self->input_id, $self->havana_db);
  print STDERR "Fetching ensembl genes\n";  
  my @genes = @{$ensemblslice->get_all_Genes_by_type($GB_ENSEMBL_INPUT_GENETYPE)};
  print STDERR "Retrieved ".scalar(@genes)." genes of type ".$GB_ENSEMBL_INPUT_GENETYPE."\n";
  print STDERR "Fetching havana genes\n";  
  my @hgenes = @{$havanaslice->get_all_Genes_by_type($GB_HAVANA_INPUT_GENETYPE)};
    print STDERR "Retrieved ".scalar(@hgenes)." genes of type ".$GB_HAVANA_INPUT_GENETYPE."\n";

  push(@genes, @hgenes);

  foreach my $gene(@genes){
  TRANSCRIPT:
    foreach my $tran (@{$gene->get_all_Transcripts}) {
      push(@transcripts, $tran);
    }
  }


  print STDERR "Finished fetching genes\n";
  $self->combined_Transcripts(@transcripts);
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
  
  my @forward_clusters;
  my @reverse_clusters;
  
  if ( @forward_transcripts ){
    @forward_clusters = $self->_cluster_Transcripts_by_genomic_range( @forward_transcripts );
  }
  if ( @reverse_transcripts ){
    @reverse_clusters = $self->_cluster_Transcripts_by_genomic_range( @reverse_transcripts );
  }
  my @clusters;
  if ( @forward_clusters ){
    push( @clusters, @forward_clusters);
  }
  if ( @reverse_clusters ){
    push( @clusters, @reverse_clusters);
  }
  return @clusters;
}

############################################################

=head2 _cluster_Transcripts_by_genomic_range

 Description : It clusters transcripts according to genomic overlap
  Args       : Array of Bio::EnsEMBL::Transcript
  Return     : Array of Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster

=cut

sub _cluster_Transcripts_by_genomic_range{
  my ($self,@mytranscripts) = @_;
  # first sort the transcripts

  my @transcripts = sort { $a->start <=> $b->start ? $a->start <=> $b->start : $b->end <=> $a->end } @mytranscripts;


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
    #print STDERR "\nIn cluster ".($count+1)."\n";
    #print STDERR "start: $cluster_starts[$count] end: $cluster_ends[$count]\n";
    #print STDERR "comparing:\n";
    #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript( $transcripts[$c] );
    
    if ( !( $transcripts[$c]->end < $cluster_starts[$count] ||
	    $transcripts[$c]->start > $cluster_ends[$count] ) ){
      $cluster->put_Transcripts( $transcripts[$c] );
      
      # re-adjust size of cluster
      if ($transcripts[$c]->start < $cluster_starts[$count]) {
	$cluster_starts[$count] = $transcripts[$c]->start;
      }
      if ( $transcripts[$c]->end > $cluster_ends[$count]) {
	$cluster_ends[$count] =  $transcripts[$c]->end;
      }
    }
    else{
      # else, create a new cluster with this feature
      $count++;
      $cluster = Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptCluster->new();
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

=head2 cluster_into_Genes

    Example :   my @genes = $self->cluster_into_Genes(@transcripts);
Description :   it clusters transcripts into genes according to exon overlap.
                It will take care of difficult cases like transcripts within introns.
                It also unify exons that are shared among transcripts.
    Returns :   a beautiful list of geen objects
    Args    :   a list of transcript objects

=cut


sub no_cluster_into_Genes{
  my ($self, @transcripts_unsorted) = @_;
  
  my $num_trans = scalar(@transcripts_unsorted);

  # First clean the coding exon cache in case it has any exons stored from previous called to the cluster_into_Genes function.
  #$self->clear_coding_exons_cache;

  my @transcripts = sort { $a->start <=> $b->start ? $a->start <=> $b->start  : $b->end <=> $a->end } @transcripts_unsorted;
  my @clusters;

  # clusters transcripts by whether or not any coding exon overlaps with a coding exon in 
  # another transcript (came from original prune in GeneBuilder)
  foreach my $tran (@transcripts) {

    my @matching_clusters;
  CLUSTER: 
    foreach my $cluster (@clusters) {
      foreach my $cluster_transcript (@$cluster) {
        if ($tran->end  >= $cluster_transcript->start &&
            $tran->start <= $cluster_transcript->end) {
          
          my $exons1 = $tran->get_all_Exons;
          my $cluster_exons = $cluster_transcript->get_all_Exons;
          
          foreach my $exon1 (@{$exons1}) {
            foreach my $cluster_exon (@{$cluster_exons}) {
              
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
  #print STDERR scalar(@clusters)." created, turning them into genes...\n";
  my @genes;
  foreach my $cluster(@clusters){
    my $count = 0;
    my $gene = new Bio::EnsEMBL::Gene;
    foreach my $transcript (@$cluster){
      $gene->add_Transcript($transcript);
    }
    push( @genes, $gene );
  }
  return @genes;
}

############################################################

=head2 re_cluster_into_Genes

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

  # First clean the coding exon cache in case it has any exons stored from previous called to the cluster_into_Genes function.
  $self->clear_coding_exons_cache;

  my @transcripts = sort { $a->coding_region_start <=> $b->coding_region_start ? $a->coding_region_start <=> $b->coding_region_start  : $b->coding_region_end <=> $a->coding_region_end } @transcripts_unsorted;
  my @clusters;

  # clusters transcripts by whether or not any coding exon overlaps with a coding exon in 
  # another transcript (came from original prune in GeneBuilder)
  foreach my $tran (@transcripts) {
  # First clean the coding exon cache in case it has any exons stored from previous called to the cluster_into_Genes function.
 # $self->clear_coding_exons_cache;

    my @matching_clusters;
  CLUSTER: 
    foreach my $cluster (@clusters) {
      
     # $self->clear_coding_exons_cache;

      foreach my $cluster_transcript (@$cluster) {
        if ($tran->coding_region_end  >= $cluster_transcript->coding_region_start &&
            $tran->coding_region_start <= $cluster_transcript->coding_region_end) {
          
          # foreach my $exon1 (@{$tran->get_all_Exons}) {
          # foreach my $cluster_exon (@{$cluster_transcript->get_all_Exons}) {
          my $exons1 = get_coding_exons_for_transcript($tran);
          my $cluster_exons = get_coding_exons_for_transcript($cluster_transcript);

          foreach my $exon1 (@{$exons1}) {
            foreach my $cluster_exon (@{$cluster_exons}) {
              
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
  #print STDERR scalar(@clusters)." created, turning them into genes...\n";
  my @genes;
  foreach my $cluster(@clusters){
    my $count = 0;
    my $gene = new Bio::EnsEMBL::Gene;
    foreach my $transcript (@$cluster){
      $gene->add_Transcript($transcript);
    }
    push( @genes, $gene );
  }
  return @genes;
}

############################################################

=head2 get_coding_exons_for_transcript

    Example :    my $exons1 = $self->get_coding_exons_for_gene($tran);
Description :   It returns the coding exons of a transcript and stores 
                them in a hash to safe computer time                
    Returns :   An ArrayRef than contain Exon objects.
    Args    :   a transcript object

=cut

{
  my %coding_exon_cache;

  sub clear_coding_exons_cache {
    %coding_exon_cache = ();
  }


sub get_coding_exons_for_transcript {
    my ($trans) = @_;

    if (exists($coding_exon_cache{$trans})) {
      return $coding_exon_cache{$trans};
    } else {
      my %coding_hash;
      
      next if (!$trans->translation);
      foreach my $exon (@{$trans->get_all_translateable_Exons}) {
        $coding_hash{$exon} = $exon;
      }

      my @coding = sort { $a->start <=> $b->start } values %coding_hash;
      #my @coding = values %coding_hash;

      $coding_exon_cache{$trans} = \@coding;
      return $coding_exon_cache{$trans};
    }
  }
}

############################################################

sub check_Clusters{
  my ($self, $num_transcripts, $clusters) = @_;
  #Safety checks
  my $ntrans = 0;

  my $cluster_num = 0;

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
  #$tran->sort;
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
  #$tran->sort;
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

=head2 prune_features

 Description: prunes out duplicated features
 Returntype : array of Bio::EnsEMBL::SeqFeature
 Args       : array of Bio::EnsEMBL::SeqFeature

=cut
    
sub prune_features {
    my ($self,$feature_hash)  = @_;
    my @pruned;
 
  ID:
    foreach my $id (keys %{ $feature_hash }) {
	my @features = @{$feature_hash->{$id}};
	@features = sort {$a->start <=> $b->start} @features;
	
	unless ( @features ){
	    print STDERR "No features here for id: $id\n";
	    next ID;
	}
	while ( @features && !defined $features[0] ){
	    #print STDERR "jumping an undefined feature\n";
	    shift @features;
	}
	
	my $prev = -1;
	
      FEATURE: 
	foreach  my $f (@features) {
	    if ($prev != -1 && $f->hseqname eq $prev->hseqname &&
		$f->start   == $prev->start &&
		$f->end     == $prev->end   &&
		$f->hstart  == $prev->hstart &&
		$f->hend    == $prev->hend   &&
		$f->strand  == $prev->strand &&
		$f->hstrand == $prev->hstrand) 
	    {
		#keep the one with highest score
		if ( $f->score > $prev->score ){
		    $prev->score( $f->score );
		}
		#print STDERR "pruning duplicated feature\n";
		#print STDERR "previous: ".$prev->gffstring."\n";
		#print STDERR "thisone : ".$f->gffstring."\n";
		next FEATURE;
	    } 
	    else {
		push(@pruned,$f);
		$prev = $f;
	    }
	}
    }
    return @pruned;
}

############################################################


############################################################
#
# GETSET METHODS
#
############################################################

# get/set method holding a reference to the db with genewise and combined genes
# this reference is set in Bio::EnsEMBL::Pipeline::RunnableDB::Gene_Builder

sub ensembl_db{
 my ($self,$ensembl_db) = @_;
 if ( $ensembl_db ){
   $self->{_ensembl_db} = $ensembl_db;
 }
 
 return $self->{_ensembl_db};
}

sub havana_db{
 my ($self,$havana_db) = @_;
 if ( $havana_db ){
   $self->{_havana_db} = $havana_db;
 }
 
 return $self->{_havana_db};
}

############################################################

sub combined_Transcripts {
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
  my ($self, @genes) = @_;
  
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

  if(!$self->{_predictions}){
    $self->{_predictions} = [];
  }
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

sub query {
  my ($self,$slice) = @_;
  
  if (defined($slice)) {
    $self->{_query} = $slice;
  }
  return $self->{_query};
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
  
  #print "NOW IM HERE\n";
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
    #print STDERR "from ".$source_exon->dbID." to ".$target_exon->dbID."\n";
    #$self->print_FeaturePair($feat);
    # I may need to add a paranoid check to see that no exons longer than the current one are transferred 
    $target_exon->add_supporting_features($feat);
    $unique_evidence{ $feat } = 1;
    $hold_evidence{ $feat->hseqname }{ $feat->start }{ $feat->end }{ $feat->hstart }{ $feat->hend } = 1;
  }
}


#fetches sequence from appropriate database

sub fetch_sequence{
  my ($self, $name, $db) = @_;

  my $sa = $db->get_SliceAdaptor; 

  my $slice = $sa->fetch_by_name($name);

  return $slice;
}

1;
