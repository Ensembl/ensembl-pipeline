##
#
# Cared for by Ensembl  <ensembl-dev@ebi.ac.uk>
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::FPC_BlastMiniGenewise

=head1 SYNOPSIS

my $obj = Bio::EnsEMBL::Pipeline::RunnableDB::MiniGenewise->new(
					     -db        => $db,
					     -input_id  => $id,
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

package Bio::EnsEMBL::Pipeline::RunnableDB::FPC_BlastMiniGenewise;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DnaPepAlignFeature;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Pipeline::Runnable::Protein::Seg;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Bio::EnsEMBL::Pipeline::Tools::GeneUtils;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Databases qw (
							     GB_GW_DBNAME
							     GB_GW_DBHOST
							     GB_GW_DBUSER
							     GB_GW_DBPASS
							     GB_GW_DBPORT
							    );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Similarity qw (
							      GB_SIMILARITY_DATABASES
							      GB_SIMILARITY_INPUTID_REGEX
							      GB_SIMILARITY_MULTI_EXON_COVERAGE
							      GB_SIMILARITY_SINGLE_EXON_COVERAGE
							      GB_SIMILARITY_MAX_INTRON
							      GB_SIMILARITY_MIN_SPLIT_COVERAGE
							      GB_SIMILARITY_GENETYPE
							      GB_SIMILARITY_MAX_LOW_COMPLEXITY 
							      GB_SIMILARITY_MASKING
							      GB_SIMILARITY_SOFTMASK
							      GB_SIMILARITY_GENETYPEMASKED
							      GB_SIMILARITY_POST_GENEMASK
							      GB_SIMILARITY_POST_EXONMASK
							      GB_SIMILARITY_BLAST_FILTER
							      
							    );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Targetted  qw (
							     GB_TARGETTED_GW_GENETYPE
							    );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::General    qw (
							     GB_INPUTID_REGEX
							    );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Scripts    qw (
							     GB_KILL_LIST
							    );

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB );

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);    
  
  my ($genewise_db) = $self->_rearrange([qw(GENEWISE_DB)], @args);
  # the following will create a genewise db from config varibles is none was given
  $self->genewise_db($genewise_db);

  # make sure at least one protein source database has been defined
  &throw("No protein source databases defined in Config::GeneBuild::Similarity::GB_SIMILARITY_DATABASES\n") 
      unless scalar(@{$GB_SIMILARITY_DATABASES});

  # make all seqfetchers
  foreach my $db(@{$GB_SIMILARITY_DATABASES}){
    my $type = $db->{"type"};
    my $seqfetcher =  $self->make_seqfetcher($db->{index}, $db->{seqfetcher});  
    $self->add_seqfetcher_by_type($type, $seqfetcher);
  }
  
  return $self; 
}


=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output
    Function:   Writes output data to db
    Returns :   array of exons (with start and end)
    Args    :   none

=cut

sub write_output {
  my($self,@features) = @_;
  
  my $gene_adaptor = $self->genewise_db->get_GeneAdaptor;
  my @genes = $self->output;
  
  my $unwritten_genes = 0;
  GENE: foreach my $gene (@genes) {	
    # do a per gene eval...
    eval {
      $gene_adaptor->store($gene);
      print STDERR "Wrote gene " . $gene->dbID . "\n";
    }; 
    if( $@ ) {
      &warning("UNABLE TO WRITE GENE\n\n$@\n\nSkipping this gene");
      $unwritten_genes++;
    }	
  }
  
  if ($unwritten_genes) {
    &throw("Failed to write " . $unwritten_genes . "out of ". scalar(@genes) . "\n");
  }
  
}

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data from the database
    Returns :   nothing
    Args    :   none

=cut

sub fetch_input {
  my($self) = @_;
    
  print STDERR "Fetching input id : " . $self->input_id. " \n";  
  &throw("No input id") unless defined($self->input_id);
    
  my ($main_bit_of_input_id, $extra_iid_one, $extra_iid_two) = 
      ($self->input_id =~ /$GB_SIMILARITY_INPUTID_REGEX/);
  
  &throw("Input id '" . $self->input_id . "' could not be parsed with " . 
         $GB_SIMILARITY_INPUTID_REGEX) if not $main_bit_of_input_id;
  
  $self->input_id($main_bit_of_input_id);
  # fetching the repeat-masked sequence from the genewise db works
  # because the Bio::EnsEMBL::DBSQL::DBAdaptor looks to the contained 
  # dnadb when asked for a RepeatFeatureAdaptor.
  $self->fetch_sequence($GB_SIMILARITY_MASKING, $self->genewise_db);
  
  my ($single_pid_db, $single_pid, $id_pool_bins, $id_pool_index);
  
  if ($extra_iid_one and $extra_iid_two) {

    if ($extra_iid_one =~ /^\d+$/ and 
        $extra_iid_two =~ /^\d+$/) {
      # assume to be a pair of number for splitting protein set. See later
      
      ($id_pool_bins, $id_pool_index) = ($extra_iid_one, $extra_iid_two);
      
      if ($id_pool_index > $id_pool_bins or
          $id_pool_index < 1) {
        
        &warning("Could not get sensible values for id_pool_bins ".
                 "('$id_pool_bins') and id_pool_index ".
                 "('$id_pool_index'); doing all proteins in region"); 
        ($id_pool_bins, $id_pool_index) = (0,0);
      }
    } else {
      # assume to be a database type name (from config) and protein id
      ($single_pid_db, $single_pid) = ($extra_iid_one, $extra_iid_two);
    }
  }
  
  if ($single_pid_db and $single_pid) {
    # force a run on this protein, regardless of whether it is masked or 
    #killed;  primarily useful for testing and clean-up after a main run
    
    my ($database) = grep { $_->{'type'} eq $single_pid_db } 
    @{$GB_SIMILARITY_DATABASES};
    
    if (not $database) {
      &throw("Your input id refers to a database ".
             "($single_pid_db) that isnt present in config\n");
    }
    
    printf(STDERR "Restricting to $single_pid from $single_pid_db based ".
           "on input id\n");
    
    my $seqfetcher =  $self->get_seqfetcher_by_type($database->{'type'});
    my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise
        ('-genomic'  => $self->query,
         '-ids'      => [$single_pid],
         '-seqfetcher' => $seqfetcher);	
    $self->runnable($runnable);
  } else {	
    # Features will be masked before using them as seeds to BlastMiniGenewise. 
    
    my ($ex_msk_reg_ref) = $self->mask_gene_region_lists($self->query);
    my @exonmask_regions = @$ex_msk_reg_ref;
    
    my %kill_list = %{$self->fill_kill_list};

    DATABASE: foreach my $database(@{$GB_SIMILARITY_DATABASES}){
      my (%features);
      
      my $pafa = $self->db->get_ProteinAlignFeatureAdaptor();
      
      my @all_feats = @{$pafa->fetch_all_by_Slice_and_score($self->query, 
                                                            $database->{'threshold'}, 
                                                            $database->{'type'})};
                 
      printf(STDERR "Fetched %d features for %s with score ".
             "above %d from %s\@%s", 
             scalar(@all_feats),
             $database->{'type'}, 
             $database->{'threshold'}, 
             $self->db->dbname, 
             $self->db->host);
      
      foreach my $f (@all_feats) {
        my $name = $f->hseqname;
        if ($name =~ /(\S+)\.\d+/) { 
          $f->hseqname($1);
        }
        
        push @{$features{$f->hseqname}}, $f;
      }
      
      printf(STDERR " (feats come from %d proteins)\n", scalar(keys %features));
      
      # flag IDs that have a feature that overlaps with a mask feature
      
      my @ids_to_ignore;
      SEQID: foreach my $sid (keys %features) {
        my $ex_idx = 0;
        #print STDERR "Looking at $sid\n";
        FEAT: foreach my $f (sort {$a->start <=> $b->start} 
                             @{$features{$sid}}) {
          #printf STDERR "Feature: %d %d\n", $f->start, $f->end;
          for( ; $ex_idx < @exonmask_regions; ) {
            my $mask_exon = $exonmask_regions[$ex_idx];
            
            # printf STDERR " Mask exon %d %d\n", $mask_exon->{'start'}, 
            $mask_exon->{'end'};
            if ($mask_exon->{'start'} > $f->end) {
              # no exons will overlap this feature
              next FEAT;
            } elsif ( $mask_exon->{'end'} >= $f->start) {
              # overlap
              push @ids_to_ignore, $f->hseqname;
              # printf STDERR "Ignoring %s\n", $f->hseqname;
              next SEQID;
            }  else {
              $ex_idx++;
            }
          }
        }
      }
      
      # remove those IDs that are either masked targetted genes or killed; 
	    
      foreach my $dud_id (@ids_to_ignore, keys %kill_list) {
        if (exists $features{$dud_id}) {
          delete $features{$dud_id};
        }
      }
	    	    
      printf (STDERR "There are %d prots left after removal of masked/killed proteins\n", scalar(keys %features));

      if ($id_pool_bins and $id_pool_index) {
        my @local_ids = sort keys %features;
		
        my (@restricted_list);
        for (my $i = $id_pool_index - 1; $i < @local_ids; $i += $id_pool_bins) {
          push @restricted_list, $local_ids[$i];
        }
        
        %features = map { $_ => $features{$_} } @restricted_list;
        
        printf(STDERR "Restricting to %d prots based on input id : @restricted_list\n", scalar(@restricted_list));
      }
      
      my @ids = sort keys %features;
      
      if ($GB_SIMILARITY_BLAST_FILTER) {
        my @sortedids = $self->sort_hids_by_coverage($database,\%features);
        my @newids = $self->prune_features($self->query,\@sortedids,\%features);
		
        @ids = sort @newids;
        
        printf (STDERR "There are %d prots left after Anopheles-style filter\n", scalar(@ids)); 
      }

      
      if (@ids) {
        # only make a runnable if there are any ids to do
        my $seqfetcher =  $self->get_seqfetcher_by_type($database->{'type'});
        my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise
        ('-genomic'  => $self->query,
         '-ids'      => \@ids,
         '-seqfetcher' => $seqfetcher);
        
        $self->runnable($runnable);	
        # at present, we'll only ever have one ...
      }
    }
  }
}    
  

=head2 run

    Title   :   run
    Usage   :   $self->run
    Function:   calls the run method on each runnable, and then calls convert_output
    Returns :   nothing, but $self->output contains results
    Args    :   none

=cut

sub run {
  my ($self) = @_;
  
  $self->genewise_db->disconnect_when_inactive(1);
  foreach my $runnable ($self->runnable) {
    $runnable->run;
  }
  $self->genewise_db->disconnect_when_inactive(0);
  
  $self->convert_output;
}

=head2 mask_gene_region_lists

    Title   :   mask_gene_region_list
    Usage   :   @list = $self->mask_region_list
    Function:   Gets, records and returns a list of start/ends that 
                correspond to regions already covered by exons
    Returns :   a list
    Args    :   none

=cut

sub mask_gene_region_lists {
  my ($self, $slice) = @_;
  
  if (not defined($self->{'_mask_region_lists'})) { 
    my (@mask_gene_types);
    
    if (@{$GB_SIMILARITY_GENETYPEMASKED}) {
      @mask_gene_types = @{$GB_SIMILARITY_GENETYPEMASKED};
    } else {
      @mask_gene_types = ($GB_TARGETTED_GW_GENETYPE);
    }
    
    my (@mask_gene_regions, @mask_exon_regions);
    
    foreach my $type (@mask_gene_types) {
      #print STDERR "Fetching gene type : $type\n";
      
      foreach my $mask_genes (@{$slice->get_all_Genes_by_type($type)}) {
        my @mask_exons = grep { $_->seqname eq $slice->id } (sort {$a->start <=> $b->start} @{$mask_genes->get_all_Exons});
        push @mask_gene_regions, { start => $mask_exons[0]->start, 
                                   end   => $mask_exons[-1]->end };
        
        foreach my $mask_exon (@mask_exons) {
          push @mask_exon_regions, { start => $mask_exon->start,
                                     end   => $mask_exon->end };
        }
      }
      #printf STDERR "  Initial mask gene list %d long\n", scalar(@mask_gene_regions);
      #printf STDERR "  Initial mask exon list %d long\n", scalar(@mask_exon_regions);
    }
    # make the mask list non-redundant. Much faster when checking against features
    my (@nr_mask_exon_regions, @nr_mask_gene_regions);
    
    foreach my $mask_exon_reg (sort {$a->{'start'} <=> $b->{'start'}} @mask_exon_regions) {
      if (@nr_mask_exon_regions and $nr_mask_exon_regions[-1]->{'end'} > $mask_exon_reg->{'start'}) {
        if ($mask_exon_reg->{'end'} > $nr_mask_exon_regions[-1]->{'end'}) {
          $nr_mask_exon_regions[-1]->{'end'} = $mask_exon_reg->{'end'};
        }
      } else {
        push @nr_mask_exon_regions, $mask_exon_reg;		
      }
    }
    foreach my $mask_gene_reg (sort {$a->{'start'} <=> $b->{'start'}} @mask_gene_regions) {
      if (@nr_mask_gene_regions and $nr_mask_gene_regions[-1]->{'end'} > $mask_gene_reg->{'start'}) {
        if ($mask_gene_reg->{'end'} > $nr_mask_gene_regions[-1]->{'end'}) {
          $nr_mask_gene_regions[-1]->{'end'} = $mask_gene_reg->{'end'};
        }
      } else {
        push @nr_mask_gene_regions, $mask_gene_reg;		
      }
    }
    
    $self->{'_mask_region_lists'} = [\@nr_mask_exon_regions, \@nr_mask_gene_regions];
  }
  
  return @{$self->{'_mask_region_lists'}};
}
    
    


=head2 convert_output

  Title   :   convert_output
  Usage   :   $self->convert_output
  Function:   converts output from each runnable into gene predictions
  Returns :   nothing, but $self->output contains results
  Args    :   none

=cut

sub convert_output {
  my ($self) =@_;
  my $trancount = 1;
  my $genetype = $GB_SIMILARITY_GENETYPE;

  foreach my $runnable ($self->runnable) {
    &throw("I don't know what to do with $runnable") 
        unless $runnable->isa("Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise");
										 
    if(!defined($genetype) || $genetype eq ''){
      $genetype = 'similarity_genewise';
      &warning("Setting genetype to $genetype\n");
    }
    
    my $anaAdaptor = $self->db->get_AnalysisAdaptor;
    # print STDERR $anaAdaptor . "\n";

    my $analysis_obj = $anaAdaptor->fetch_by_logic_name($genetype);

    if ( !defined $analysis_obj ){
      # make a new analysis object
      $analysis_obj = new Bio::EnsEMBL::Analysis
	(-db              => 'NULL',
	 -db_version      => 1,
	 -program         => $genetype,
	 -program_version => 1,
	 -gff_source      => $genetype,
	 -gff_feature     => 'gene',
	 -logic_name      => $genetype,
	 -module          => 'FPC_BlastMiniGenewise',
      );
    }

    my @results = $runnable->output;

    # filter out masked genes here if appropriate
    if ($GB_SIMILARITY_POST_GENEMASK or $GB_SIMILARITY_POST_EXONMASK) {
      my (@mask_regions, @filtered_results);
      
      my ($exonmask_regions, $genemask_regions) = $self->mask_gene_region_lists;
      
      if ($GB_SIMILARITY_POST_GENEMASK) {
        @mask_regions = @$genemask_regions;
        
      } else {
        @mask_regions = @$exonmask_regions;
      }

      GENE: foreach my $gene (@results) {
        my $keep_gene = 1;
        my $mask_reg_idx = 0;
        
        my @exons = sort {$a->start <=> $b->start} ($gene->sub_SeqFeature);
        my @test_regions;
        
        if ($GB_SIMILARITY_POST_GENEMASK) {
          @test_regions = ({start => $exons[0]->start, end => $exons[-1]->end});
        }
        else {
          @test_regions = map { { start => $_->start, end => $_->end } } @exons;
        }
        
        FEAT: foreach my $f (@test_regions) {
          for( ; $mask_reg_idx < @mask_regions; ) {
            my $mask_reg = $mask_regions[$mask_reg_idx];
            
            if ($mask_reg->{'start'} > $f->{'end'}) {
              # no mask regions will overlap this feature
              next FEAT;
		    }
            elsif ( $mask_reg->{'end'} >= $f->{'start'}) {
              # overlap			
              $keep_gene = 0;
              last FEAT;
            }			
            else {
              $mask_reg_idx++;
            }
		}
        }
        
        if ($keep_gene) {
          push @filtered_results, $gene;
        }
      }
      
      @results = @filtered_results;
    }

    my $genes = $self->make_genes($genetype, $analysis_obj, \@results);
    $self->output(@$genes);
  }
}


=head2 make_genes

  Title   :   make_genes
  Usage   :   $self->make_genes
  Function:   makes Bio::EnsEMBL::Genes out of the output from runnables
  Returns :   array of Bio::EnsEMBL::Gene  
  Args    :   $genetype: string
              $analysis_obj: Bio::EnsEMBL::Analysis
              $results: reference to an array of FeaturePairs

=cut

sub make_genes {
  my ($self, $genetype, $analysis_obj, $results) = @_;
  my @genes;

  &throw("[$analysis_obj] is not a Bio::EnsEMBL::Analysis\n") 
      unless defined($analysis_obj) && $analysis_obj->isa("Bio::EnsEMBL::Analysis");

  my $count = 0;
  foreach my $tmpf (@{$results}) {
    $count++;

    my ($first_exon) = $tmpf->sub_SeqFeature;
    my ($first_supp_feat) = $first_exon->sub_SeqFeature;
    my $prot_id = $first_supp_feat ? $first_supp_feat->hseqname : "Unknown protein";

    my $unchecked_transcript = 
        Bio::EnsEMBL::Pipeline::Tools::GeneUtils->SeqFeature_to_Transcript($tmpf, 
                                                                           $self->query, 
                                                                           $analysis_obj, 
                                                                           $self->genewise_db, 
                                                                           0);
    if (not defined $unchecked_transcript) {
      printf(STDERR "   Transcript %d (%s) : REJECTED (could not make from SeqFeatures)\n", $count, $prot_id);
      next;
    }
        
    # validate transcript
    my @seqfetchers = $self->each_seqfetcher;
    
    my $valid_transcripts = 
        Bio::EnsEMBL::Pipeline::Tools::GeneUtils->validate_Transcript($unchecked_transcript, 
                                                                      $self->query, 
                                                                      $GB_SIMILARITY_MULTI_EXON_COVERAGE, 
                                                                      $GB_SIMILARITY_SINGLE_EXON_COVERAGE, 
                                                                      $GB_SIMILARITY_MAX_INTRON, 
                                                                      $GB_SIMILARITY_MIN_SPLIT_COVERAGE, 
                                                                      \@seqfetchers, 
                                                                      $GB_SIMILARITY_MAX_LOW_COMPLEXITY);

    if (not defined $valid_transcripts) {
      printf (STDERR "   Transcript %d (%s) : REJECTED (validation failed)\n", $count, $prot_id);
      next;
    }
    else {
      printf (STDERR "   Transcript %d (%s) : ACCEPTED\n", $count, $prot_id);
    }
      
    # make genes from valid transcripts
    foreach my $checked_transcript(@$valid_transcripts){
      # add start codon if appropriate
      $checked_transcript = Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->set_start_codon($checked_transcript);

      # add terminal stop codon if appropriate
      $checked_transcript = Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->set_stop_codon($checked_transcript);

      # flip reverse strand supporting features. This used to be done 
      # during remapping, whch is not necessary any more

      foreach my $exon(@{$checked_transcript->get_all_Exons}) {
        foreach my $sf(@{$exon->get_all_supporting_features}) {
          # this should be sorted out by the remapping to rawcontig ... strand is fine
          if ($sf->start > $sf->end) {
            my $tmp = $sf->start;
            $sf->start($sf->end);
            $sf->end($tmp);
          }
        }
      }

      my $gene = new Bio::EnsEMBL::Gene;
      $gene->type($genetype);
      $gene->analysis($analysis_obj);
      $gene->add_Transcript($checked_transcript);

      
      push (@genes, $gene);
    }
  }
  
  return \@genes;

}



=head2 output

 Title   : output
 Usage   :
 Function: get/set for output array
 Example :
 Returns : array of Bio::EnsEMBL::Gene
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
   #print STDERR "have ".@{$self->{'_output'}}." as output\n";
   return @{$self->{'_output'}};
}

=head2 make_seqfetcher

 Title   : make_seqfetcher
 Usage   :
 Function: makes a Bio::EnsEMBL::SeqFetcher to be used for fetching protein sequences. If 
           $index is specified, then a Getseqs fetcher is made, otherwise it throws
 Example :
 Returns : Bio::EnsEMBL::SeqFetcher
 Args    : string representing path to sequence index


=cut

sub make_seqfetcher {
  my ( $self, $index, $seqfetcher_class ) = @_;
  my $seqfetcher;

  (my $class = $seqfetcher_class) =~ s/::/\//g;
  require "$class.pm";

  if(defined $index && $index ne ''){
    my @db = ( $index );
    
    # make sure that your class is compatible with the index type
    $seqfetcher = "$seqfetcher_class"->new('-db' => \@db, );
    
  }
  else{
    &throw("Can't make seqfetcher\n");
  }

  return $seqfetcher;

}

=head2 each_seqfetcher

  Title   :   each_seqfetcher
  Usage   :   my @seqfetchers = $self->each_seqfetcher
  Function:   Returns an array of Bio::DB::RandomAccessI representing the various sequence indices 
              listed in Config::GeneBuild::Similarity::GB_SIMILARITY_DATABASES
  Returns :   Array of Bio::DB::RandomAccessI
  Args    :   none

=cut

sub each_seqfetcher {
  my ($self) = @_;
  
  if(!defined($self->{'_seqfetchers'})) {
     $self->{'_seqfetchers'} = {};
   }
    
  # NB array of seqfetchers
   return values(%{$self->{'_seqfetchers'}});
}

=head2 each_seqfetcher_by_type

  Title   :   each_seqfetcher_by_type
  Usage   :   my %seqfetchers_by_type = $self->each_seqfetcher_by_type
  Function:   Returns a hash of Bio::DB::RandomAccessI representing the various sequence indices 
              listed in Config::GeneBuild::Similarity::GB_SIMILARITY_DATABASES keyed by type listed therein.
  Returns :   Hash of all seqfetchers linking db_type to Bio::DB::RandomAccessI
  Args    :   none

=cut

sub each_seqfetcher_by_type {
  my ($self) = @_;
  
  if(!defined($self->{'_seqfetchers'})) {
     $self->{'_seqfetchers'} = {};
   }
    
  # NB hash of seqfetchers
   return %{$self->{'_seqfetchers'}};
}

=head2 add_seqfetcher_by_type

  Title   :   add_seqfetcher_by_type
  Usage   :   $self->add_seqfetcher_by_type('swall', $seqfetcher)
  Function:   Adds a Bio::DB::RandomAccessI into $self->{'_seqfetchers'} keyed by type
  Returns :   Nothing
  Args    :   $type - string representing db type
              $seqfetcher - Bio::DB::RandomAccesI

=cut

sub add_seqfetcher_by_type{
  my ($self, $type, $seqfetcher) = @_;

  &throw("no type specified\n") unless defined ($type); 
  &throw("no suitable seqfetcher specified: [$seqfetcher]\n") 
    unless defined ($seqfetcher) && $seqfetcher->isa("Bio::DB::RandomAccessI"); 
  $self->{'_seqfetchers'}{$type} = $seqfetcher;
}


=head2 get_seqfetcher_by_type

  Title   :   get_seqfetcher_by_type
  Usage   :   my $seqfetcher = $self->get_seqfetcher_by_type('swall')
  Function:   Fetches the seqfetcher associated with a particular db type as specified in Config::GeneBuild::Similarity::GB_SIMILARITY_DATABASES
  Returns :   Bio::DB::RandomAccessI
  Args    :   $type - string representing db type

=cut

sub get_seqfetcher_by_type{
  my ($self, $type) = @_;
  my %seqfetchers = $self->each_seqfetcher_by_type;
  foreach my $dbtype(keys %seqfetchers){
    if ($dbtype eq $type){
      return $seqfetchers{$dbtype};
    }
  }
}


=head2 fill_kill_list

 Title   : fill_kill_list
 Usage   : 
 Function: 
           
 Returns : 
 Args    : 

=cut

sub fill_kill_list {
  my ($self) = @_;
  my %kill_list;
  
  if (defined($GB_KILL_LIST) && $GB_KILL_LIST ne '') {
    open (KILL_LIST, "< $GB_KILL_LIST") or die "can't open $GB_KILL_LIST";
    while (<KILL_LIST>) {
      
      chomp;
      my @list = split;
      next unless scalar(@list); 	# blank or empty line
      $kill_list{$list[0]} = 1;
    }
    
    close KILL_LIST or die "error closing $GB_KILL_LIST\n";
  }
  return \%kill_list;
}


=head2 sort_hids_by_coverage

 Title   : sort_hids_by_coverage
 Usage   : my @sorted_ids = $self->sort_hids_by_coverage($database,\%hids)
 Function: This is supposed to order each hid in function of there coverage with the sequence they
    have been build with
 Example :
 Returns : an array of ordered hids
 Args    : Database description(where to fetch the feature from). A hash refererence containing hid
    and its features.
    
=cut

sub sort_hids_by_coverage{
  my ($self,$database,$hash_ref) = @_;
  my @id =  keys %{$hash_ref};
  my $seqfetcher =  $self->get_seqfetcher_by_type($database->{'type'});
  my $forward_start;
  my $forward_end;
  my $reverse_start;
  my $reverse_end;
  my $f_matches = 0;
  my $r_matches = 0;
  my $protname;
  my $seq;
  my $plength;
  my $features;
  my $best_cov;
  my %idsreturned;
  my @sorted;

  IDS: foreach my $id (@id) {
       
    $protname = undef;
    $seq = undef;
    $plength = 0;
    $forward_start = 0;
    $forward_end = 0;
    $reverse_start = 0;
    $reverse_end = 0;
    $best_cov = 0;
    $features = $hash_ref->{$id};
    
    FEAT: foreach my $feat(@$features) {
      if ($feat->strand == 1) {
        if (!defined($protname)){
          $protname = $feat->hseqname;
        }
        if($protname ne $feat->hseqname){
          &warning("$protname ne " . $feat->hseqname . "\n");
        }
        
        if((!$forward_start) || $forward_start > $feat->hstart){
          $forward_start = $feat->hstart;
        }
        
        if((!$forward_end) || $forward_end < $feat->hend){
          $forward_end= $feat->hend;
        }
      }
      
      if ($feat->strand == -1) {
        if (!defined($protname)){
          $protname = $feat->hseqname;
        }
        if($protname ne $feat->hseqname){
          &warning("$protname ne " . $feat->hseqname . "\n");
        }
        
        if((!$reverse_start) || $reverse_start > $feat->hstart){
          $reverse_start = $feat->hstart;
        }
        
        if((!$reverse_end) || $reverse_end < $feat->hend){
          $reverse_end= $feat->hend;
        }
      }      
    }
    
    $f_matches = ($forward_end - $forward_start + 1);
    $r_matches = ($reverse_end - $reverse_start + 1);
    
    if($self->get_length_by_id($protname)) {
      $plength = $self->get_length_by_id($protname);
    }
    
    else {
      
      SEQFETCHER:
      foreach my $seqfetcher( $self->each_seqfetcher){
        eval{
          $seq = $seqfetcher->get_Seq_by_acc($protname);
        };
	       
        if ($@) {
          &warning("FPC_BMG:Error fetching sequence for [$protname] - trying next seqfetcher:[$@]\n");
        }
        
        if (defined $seq) {
          last SEQFETCHER;
        }
        
      }
	   
      if(!defined $seq){
        &warning("FPC_BMG:No sequence fetched for [$protname] - can't check coverage - set coverage to 1 (entry will be considered as having low coverage...sorry\n");
      }
      
      eval {
        $plength = $seq->length;
      };
      
      if($plength) {
        $self->add_length_by_id($protname,$plength);
      }
    }	 
    
    my $f_coverage;
    my $r_coverage;
    
    #check the low complexity
    my $valid = 0;
    if ($seq) {	   
      eval {
        # Ugh! 
        my $analysis = Bio::EnsEMBL::Analysis->new(
                                                   -db           => 'low_complexity',
                                                   -program      => '/usr/local/ensembl/bin/seg',
                                                   -program_file => '/usr/local/ensembl/bin/seg',
                                                   -gff_source   => 'Seg',
                                                   -gff_feature  => 'annot',
                                                   -module       => 'Seg',
                                                   -logic_name   => 'Seg'
                                                   
                                                   );
        
        my $seg = new  Bio::EnsEMBL::Pipeline::Runnable::Protein::Seg(    
                                                                          -query    => $seq,
                                                                          -analysis => $analysis,
                                                                          );
        
        $seg->run;
        
        
        if($seg->get_low_complexity_length > $GB_SIMILARITY_MAX_LOW_COMPLEXITY){
          &warning("discarding feature too much low complexity\n");
          $valid = 0;
        }
        $valid = 1;	 	 
      };
      
      if($@){
        print STDERR "problem running seg: \n[$@]\n";
        $valid = 1;		# let transcript through
      }
    }
    
    if(!defined($plength) || $plength == 0 || $valid == 0){
      &warning("no sensible length for $protname - can't get coverage - or too much low complexity\n");
      
      #Can't check the length of the protein, set the coverage to 1 for both strand
      $f_matches = 1;
      $r_matches = 1;
      $plength = 100;
    }
    
    $f_coverage = $f_matches/$plength;
    $r_coverage = $r_matches/$plength;
    $f_coverage *= 100;
    $r_coverage *= 100;
    
    if ($f_coverage > $r_coverage) {
      $best_cov = $f_coverage;
    }
    else {
      $best_cov = $r_coverage;
    }
    
    $idsreturned{$protname} = $best_cov;
  }
  
  @sorted = sort { $idsreturned{$b} <=> $idsreturned{$a} }  keys %idsreturned;
  
  return @sorted;
}


=head2 prune_features

 Title   : prune_features
 Usage   : my @pruned_ids = $self->prune_features($genseq,\@sortedids,\%idhash);
 Function: This method as for goal to limit the number of hids sent to BlastMiniGenewise.
 Example :
 Returns : An array of pruned hids
 Args    : genomic sequence, an array of sorted hids, an array reference of all the hids and 
    their features

=cut

sub prune_features{
  my ($self,$genseq,$sortedids_array_ref,$hash_ref) = @_;
  my $forwardcountstring = '0' x $genseq->length;
  my $reversecountstring = '0' x $genseq->length;
  my %ids;
  my %hidcount;
  my %finalids;
  
  my @ident = @$sortedids_array_ref;
  
  IDS: foreach my $id (@ident) {
    my $features = $hash_ref->{$id};
    
    my @exons;
    
    next IDS unless (ref($features) eq "ARRAY");
    next IDS unless (scalar(@$features) >= 1);
    
    FEAT: foreach my $feat(@$features) {
      my $length = $feat->end - $feat->start + 1;
      
      if($feat->strand == 1) {
        my $count = 0;
        
        my $str = substr($forwardcountstring,$feat->start,$length);
        
        foreach my $byte (split //, $str) {
          $count = $count + $byte;
        }
        
        $hidcount{$id}->{'lengthforward'} = $hidcount{$id}->{'lengthforward'} + $length;
        $hidcount{$id}->{'coverageforward'} = $hidcount{$id}->{'coverageforward'} + $count;
        
        if ($count > 10) {
          next FEAT;
        }
        else {
          substr($forwardcountstring,$feat->start,$length) = '1' x $length; 
          $ids{$feat->hseqname} = 1;
        }
      }
      elsif ($feat->strand == -1) {
        my $count = 0;
        my $length = $feat->end - $feat->start + 1;
        
        my $str = substr($reversecountstring,$feat->start,$length);
        
        foreach my $byte (split //, $str) {
          $count = $count + $byte;
        }
        $hidcount{$id}->{'lengthreverse'} = $hidcount{$id}->{'lengthreverse'} + $length;
        $hidcount{$id}->{'coveragereverse'} = $hidcount{$id}->{'coveragereverse'} + $count;
        
        if ($count > 10) {
          next FEAT;
        } 
        else {
          substr($reversecountstring,$feat->start,$length) = '1' x $length; 		    
          $ids{$feat->hseqname} = 1;		    
        }
      }
    }
  }
  
  foreach my $hidk (keys %ids) {
    my $percforward;
    my $percreverse;
    
    if ($hidcount{$hidk}->{'lengthforward'} > 0) {
      $percforward = $hidcount{$hidk}->{'coverageforward'} * 100 / $hidcount{$hidk}->{'lengthforward'};
    }
    
    if ($hidcount{$hidk}->{'lengthreverse'} > 0) {
      $percreverse = $hidcount{$hidk}->{'coveragereverse'} * 100 / $hidcount{$hidk}->{'lengthreverse'};
    }
    
    if (($percforward <= 90)&&($percreverse <= 90)) {
      $finalids{$hidk} = 1;
    }
  }
  
  my @returnids = keys %finalids;
  return @returnids;
}


=head2 add_length_by_id

 Title   : add_length_by_id
 Usage   : $self->add_length_by_id("QE345",345)
 Function: Caches the length of a protein
 Example :
 Returns : nothing
 Args    : protein id and its length


=cut

sub add_length_by_id{
  my ($self,$id,$length) = @_;

  &throw("no id specified\n") unless defined ($id); 
  &throw("no length specified\n") unless defined ($length);
  
  if (!defined $self->{'_idlength'}{$id}) {
    $self->{'_idlength'}{$id} = $length;
  }
}

=head2 get_length_by_id

 Title   : get_length_by_id
 Usage   : $self->get_length_by_id("QGY7980")
 Function: retrieves the length of a given protein
 Example :
 Returns : Length of a protein
 Args    : Protein id


=cut

sub get_length_by_id{
  my ($self,$id) = @_;
  my %idlength = $self->{'_idlength'};
  if ($idlength{$id}) {
    return $idlength{$id};
  }
}

=head2 genewise_db

 Title   : genewise_db
 Usage   : retrieves /sets the db pointed at by host GB_GW_HOST -
           where preliminary genes are being written
 Function:

 Returns :
 Args    :

=cut

sub genewise_db {
  my( $self, $genewise_db ) = @_;
  
  if ($genewise_db){
    $genewise_db->isa("Bio::EnsEMBL::DBSQL::DBAdaptor")
        || &throw("Input [$genewise_db] isn't a ".
                  "Bio::EnsEMBL::DBSQL::DBAdaptor");
    $self->{_genewise_db} = $genewise_db;
  }
  
  if(!$self->{_genewise_db}){
    $GB_GW_DBHOST = $self->db->host     if (!$GB_GW_DBHOST);
    $GB_GW_DBUSER = $self->db->username if (!$GB_GW_DBUSER);
    $GB_GW_DBPASS = $self->db->password if (!$GB_GW_DBPASS);
    $GB_GW_DBNAME = $self->db->dbname   if (!$GB_GW_DBNAME);
    $GB_GW_DBPORT = $self->db->dbname   if (!$GB_GW_DBPORT);
        
    $self->{_genewise_db} = new Bio::EnsEMBL::DBSQL::DBAdaptor
        (
         '-host'   => $GB_GW_DBHOST,
         '-user'   => $GB_GW_DBUSER,
         '-pass'   => $GB_GW_DBPASS,
         '-port'   => $GB_GW_DBPORT,
         '-dbname' => $GB_GW_DBNAME,
         '-dnadb'  => $self->db,
         );
  }
  return $self->{_genewise_db};
}




1;



