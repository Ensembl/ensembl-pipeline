#
# Written by Eduardo Eyras
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code
#
=head1 NAME

Bio::EnsEMBL::Pipeline::GeneComparison::ComparativeTools

=head1 SYNOPSIS

=head1 DESCRIPTION

Class containing some useful methods to exploit the info in compara.
It contains methods to be used for checking genes for orthology on syntenic slices, 
checking for synteny breaking, etc... more information can be founf in the POD.

=head1 CONTACT

ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::Pipeline::GeneComparison::ComparativeTools;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Root;

@ISA = qw(Bio::EnsEMBL::Root);

############################################################

=head2 get_all_syntenic_slices


=cut

=head2 test_for_synteny_breaking

   Args       : $focus_slice     -> slice of which we want to retrieve the syntenic slices
                $focus_db        -> db used by compara for the dna matches
                $focus_species   -> the species for the transcript $transcript
                $compara_db      -> compara db holding the dna-dna matches between
                                    the focus and the target species
                $target_db       -> db used by compara for the dna matches
                $target_species  -> the species against which we want to compare $transcript
  Description : this method reads a slice and one of the target species 
                specified in compara and retrieves the syntenic slice in that species
                if any. It gets all the dna align features and thread them
                into a bigger slice. The syntenic slice is made to be at least
                as big as the original slice.
  Returntype  : an arrayref of slices (the syntenic slices)
  
=cut

sub get_all_syntenic_slices{
  my ($self, $focus_slice, $focus_db, $focus_species, $compara_db, $target_db, $target_species) = @_;
  
  ############################################################
  # store the size of the focus slice:
  my $focus_length = $focus_slice->chr_end - $focus_slice->chr_start + 1 ;
  print STDERR "search syntenic slice of ".$focus_slice->chr_name."."
    .$focus_slice->chr_start."-".$focus_slice->chr_end."\n";
  
  # get the adaptor
  my $adaptor = $compara_db->get_DnaAlignFeatureAdaptor;
  
  # it returns an array reference of Bio::EnsEMBL::DnaDnaAlignFeature objects
  #my @features = @{$adaptor->fetch_all_by_Slice($slice, $focus_species, $focus_db->assembly_type)};
  
  my @features = @{$adaptor->fetch_all_by_species_region($focus_species,
							 $focus_db->assembly_type,
							 $target_species,
							 $target_db->assembly_type,
							 $focus_slice->chr_name,
							 $focus_slice->chr_start, 
							 $focus_slice->chr_end
							)};
	
  
  ############################################################
  # chain features into longer ones:
  print STDERR scalar(@features)." syntenic features found\n";
  my %chr_features;
  foreach my $feature_pair ( @features ){
    #print STDERR $feature_pair->hseqname."\t".
    #$feature_pair->hstart."\t".
    #$feature_pair->hend."\t".
    #$feature_pair->hstrand."\n";
    push( @{ $chr_features{$feature_pair->hseqname} }, $feature_pair );
  }
  
  my @slices;
  foreach my $chr ( keys %chr_features ){
    my @chr_features = sort { $a->hstart <=> $b->hstart } @{ $chr_features{$chr} };
    my $start = $chr_features[0]->hstart;
    my $end   = $chr_features[0]->hend;
    foreach (my $i=1; $i<scalar(@chr_features); $i++ ){
      if ( $chr_features[$i]->start - $chr_features[$i-1] - 1 < $focus_length ){
	$end = $chr_features[$i]->hend;
      }
      else{
	############################################################
	# if the slice is not as big as the original, we make it as big
	# it's ok if they overlap
	my $length = $end - $start + 1;
	if ( $focus_length > $length ){
	  $start -= ( $focus_length - $length )/2;
	  $end   += ( $focus_length - $length )/2;
	}
	
	my $target_slice = $target_db->get_SliceAdaptor->fetch_by_chr_start_end($chr,$start,$end);
	push(@slices, $target_slice);
	$start = $chr_features[$i]->hstart;
	$end   = $chr_features[$i]->hend;
      }
    }
    my $length = $end - $start + 1;
    if ( $focus_length > $length ){
      $start -= int( ( $focus_length - $length )/2 );
      $end   += int( ( $focus_length - $length )/2 );
    }
    my $target_slice = $target_db->get_SliceAdaptor->fetch_by_chr_start_end($chr,$start,$end);
    push(@slices, $target_slice);
  }
  if (@slices){
    print STDERR "Produced slices:\n";
    foreach my $slice (@slices){
      print STDERR $slice->chr_name.".".$slice->chr_start."-".$slice->chr_end."\n";
    }
  }
  return \@slices;
}

############################################################

=head2 test_for_synteny_breaking

   Args       : $transcript      -> a transcript to test
                $db              -> a db with an assembly to which the transcript is associated
                $focus_db        -> db used by compara for the dna matches
                $focus_species   -> the species for the transcript $transcript
                $compara_db      -> compara db holding the dna-dna matches between
                                    the focus and the target species
                $target_db       -> db used by compara for the dna matches
                $target_species  -> the species against which we want to compare $transcript
                $threshold       -> it will return matches with coverage above this threshold
  Description : method to check for orthology for the flaking transcripts of a given transcript
                in a region syntenic to the piece containing these transcripts.
                This is the logic:
  
              left-flaking transcript   transcript         right-flaking transcript
                      _____              ______             ______
                 ____|_____|____________|______|___________|______|___   slice
                        |                   |                  |
                        |                   |                  |
                 ______\|/_________________\|/________________\|/_____   syntenic slice
            exonerate   ?                   ?                  ? 
             match?
                                                                         synteny breaking?
                       Yes                 Yes                Yes              No
                       Yes                 No                 Yes              Yes
                       No                  Yes                Yes              No
                       No                  No                 Yes              Yes
                       Yes                 Yes                No               No
                       Yes                 No                 No               Yes
                       No                  No                 No               Yes
                       No                  Yes                No               No


  Returntype  : a boolean, whether ot not synteny is broken

=cut

sub test_for_synteny_breaking{
  my ($self, $transcript, $db, $focus_db, $focus_species, $compara_db, $target_db, $target_species, $threshold ) = @_;

  ############################################################
  # transcript is a transcript
  my @exons = sort { $a->start <=> $b->start } @{$transcript->get_all_Exons};
  my $low    = $exons[0]->start;
  my $high   = $exons[$#exons]->end;
  my $strand = $exons[0]->strand;

  ############################################################
  # get a slice for the transcript with extra sequence on both sides
  my $focus_slice;
  if ($transcript->dbID){
    $focus_slice = $db->get_SliceAdaptor->fetch_by_transcript_id( $transcript->dbID, 500000 );
  }

  ############################################################
  # no dbID, need to get a slice from the slice where transcript sits:
  if (  !$transcript->dbID || !$focus_slice  ){
    my $chr_name  = $exons[0]->contig->chr_name;
    my $chr_start = $exons[0]->contig->chr_start;
    my $chr_end   = $exons[0]->contig->chr_end;
    
    my $start = $chr_start + $low  - 1;
    my $end   = $chr_start + $high - 1;
    
    ############################################################
    # get slice from the same db where our transcript is
    $focus_slice = 
      $db->get_SliceAdaptor->fetch_by_chr_start_end( $chr_name, $start - 500000, $end + 500000);
  }

  my ($left, $right) = $self->get_flanking_genes($transcript,$focus_slice);
  my @left  = @$left;
  my @right = @$right;
  
  print STDERR scalar(@left)." genes found on the left\n";
  print STDERR scalar(@right)." genes found on the right\n";  
  my $left_gene  = pop @left;
  my $right_gene = shift @right;

  ################################################################
  # find the slice(s) syntenic to this focus slice where the genes are  
  my @target_slices;
  my $target_slices = 
    $self->get_all_syntenic_slices( $focus_slice, $focus_db, $focus_species, $compara_db, $target_db, $target_species);
  if ( $target_slices && @{$target_slices} ){
    @target_slices = @{$target_slices};
  }
  else{
    print STDERR "Could not find any syntenic region\n";
    return 1; #synteny is broken
  }
  
  ############################################################
  # run exonerate for each gene: left-flanking, transcript and right-flanking gene
  # take the orthologue with best coverage:
  my $synteny_breaking = 0;
  foreach my $target_slice (@target_slices){
    print STDERR "running exonerate with ".$target_slice->chr_name.".".$target_slice->chr_start."-".$target_slice->chr_end."\n";
    
    my @left_orthology_genes;
    my @right_orthology_genes;
    my @orthology_transcripts;
    
    ############################################################
    if ( $left_gene ){    
      my @sorted_left  = sort{ $self->_coverage($b) <=> $self->_coverage($a) } $self->align_with_exonerate($left_gene, $target_slice );
    
      if (  $sorted_left[0] ){
	print STDERR "left ortholog - coverage: ".$self->_coverage( $sorted_left[0] )."\n";
	Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript( $sorted_left[0],1 );
      }
      @left_orthology_genes  =  @{$self->filter_by_coverage( \@sorted_left, 80 )};
    }
    
    ############################################################
    if( $right_gene ){
      my @sorted_right = sort{ $self->_coverage($b) <=> $self->_coverage($a) } $self->align_with_exonerate($right_gene, $target_slice );
      
      if (  $sorted_right[0] ){
	print STDERR "right ortholog - coverage: ".$self->_coverage( $sorted_right[0] )."\n";
	Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript( $sorted_right[0],1 );
      }
      
      @right_orthology_genes = @{$self->filter_by_coverage( \@sorted_right, 80 )};
    }

    ############################################################
    if( $transcript ){
      my @sorted       = sort{ $self->_coverage($b) <=> $self->_coverage($a) } $self->align_with_exonerate($transcript, $target_slice );
      
      if ( $sorted[0] ){
	print STDERR "\"transcript\" ortholog - coverage: ".$self->_coverage( $sorted[0] )." \n";
	Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript( $sorted[0],1 );
      }
      @orthology_transcripts =  @{$self->filter_by_coverage( \@sorted, 80 )};
    }
    
    ############################################################
    if ( ( $left_orthology_genes[0] && !$orthology_transcripts[0] && $right_orthology_genes[0] )
	 ||
	 (  $left_orthology_genes[0] && !$orthology_transcripts[0] )
	 ||
	 (  !$orthology_transcripts[0] && $right_orthology_genes[0] )
       ){
      print STDERR "There is synteny breaking\n";
      $synteny_breaking = 1;
    }
    
    if (  $left_orthology_genes[0] && $orthology_transcripts[0] && $right_orthology_genes[0] ){
      print STDERR "all genes have orthologs\n";
      if ( $left_orthology_genes[0]->end < $orthology_transcripts[0]->start 
	   &&
	   $right_orthology_genes[0]->start > $orthology_transcripts[0]->end 
	 ){
	print STDERR "synteny is preserved\n";
	$synteny_breaking = 0;
      }
      else{
	print STDERR "synteny is reversed\n";
	$synteny_breaking = 0;
      }
    }
    
    if (  !$left_orthology_genes[0] && !$orthology_transcripts[0] && !$right_orthology_genes[0] ){
      print STDERR "there is no synteny at all\n";
      $synteny_breaking = 1;
    }
  }
  return $synteny_breaking;
}
############################################################

sub _coverage{
    my ($self,$tran) = @_;
    my @exons = @{$tran->get_all_Exons};
    my @evi = @{$exons[0]->get_all_supporting_features};
    return $evi[0]->score;
}

############################################################

sub _percent_id{
    my ($self,$tran) = @_;
    my @exons = @{$tran->get_all_Exons};
    my @evi = @{$exons[0]->get_all_supporting_features};
    return $evi[0]->percent_id;
}

############################################################

=head2 get_flanking_genes

  Args       :  $transcript      -> a transcript to test
                $slice           -> a slice where we want to try to align the transcript
  Description:  get the genes on either side of a transcript
                (not necessarily on the same strand)
                given the slice and given the transcript which must be on that slice coordinates.
                The convetion here is to see the dna horizontally, forward stand in the upper-side and
                reverse strand underneath. The Coordinate system runs from left to right, so left
                of the transcript means 'smaller coordinates' and right of the gene means 'larger coordinates'
  Returntype :  2 arrayrefs containing the left-flanking and right flanking-genes respectively

=cut

sub get_flanking_genes{
  my ($self, $transcript, $focus_slice ) = @_;
  
  ############################################################
  my @exons = sort { $a->start <=> $b->start } @{$transcript->get_all_Exons};
  my $low    = $exons[0]->start;
  my $high   = $exons[$#exons]->end;
  my $strand = $exons[0]->strand;
  #print STDERR "Finding flanking genes in slice: ".$focus_slice->chr_name.".".$focus_slice->chr_start."-".$focus_slice->chr_end."\n";
  #print STDERR "transcript: $low |---| $high  strand:$strand\n";

  ############################################################
  # find the N=2 flanking genes
  # first, get all the genes in the focus_slice
  my @genes = @{$focus_slice->get_all_Genes};
  my @left;
  my @right;
  my %start;

  ############################################################
  # the retrieved genes are in the slice coordinates!
  foreach my $gene ( @genes ){
    foreach my $trans ( @{$gene->get_all_Transcripts} ){
      
      if ( $trans->type ){
	my $type = $trans->type;
	next if ( $type eq 'Putative' || $type eq 'Novel_Transcript' || $type eq 'Pseudogene' );
      }
      my @exons = sort { $a->start <=> $b->end } @{$trans->get_all_Exons};
      #print STDERR "candidate: ".
      #($exons[0]->start + $focus_slice->chr_start - 1).
      #" |---| ".
      #($exons[$#exons]->end + $focus_slice->chr_start - 1).
      #" strand:".$exons[0]->strand."\n";
      
      ############################################################
      # do not restrict the strand
      $start{$trans} = $exons[0]->start + $focus_slice->chr_start - 1;
      if ( $exons[$#exons]->end + $focus_slice->chr_start - 1< $low ){
	push( @left, $trans );
	#print STDERR "left\n";
      }
      elsif ( $exons[0]->start + $focus_slice->chr_start - 1 > $high ){
	push( @right, $trans );
	#print STDERR "right\n";
      }
      else{
	#print STDERR ($exons[0]->start + $focus_slice->chr_start - 1)." > $high\n";
	#print STDERR "skipped\n";
	next;
      }
    }
  }
  print STDERR "left: ".scalar(@left)." right: ".scalar(@right)."\n";
  
  # sort them:
  @left = map { $_->[1] } sort { $a->[0] <=> $b->[0] } map { [ $start{$_}, $_ ] } @left;
  @right= map { $_->[1] } sort { $a->[0] <=> $b->[0] } map { [ $start{$_}, $_ ] } @right;
  
  return (\@left, \@right);
}


############################################################

=head2 test_for_orthology

  Args       :  $transcript      -> a transcript to test
                $db              -> a db with an assembly to which the transcript is associated
                $focus_db        -> db used by compara for the dna matches
                $focus_species   -> the species for the transcript $transcript
                $compara_db      -> compara db holding the dna-dna matches between
                                    the focus and the target species
                $target_db       -> db used by compara for the dna matches
                $target_species  -> the species against which we want to compare $transcript
                $threshold       -> it will return matches with coverage above this threshold
  Description: method to try to align a transcript on a syntenic region for a give target species.
               It has many arguments as I have tried to make it independent of the
               way to get the focus/target databases from compara
  Returntype : an arrayref with the orthologs in chr coordinates
               If there are no orthologs it returns undef.
               If there is no syntenic region it will also return undef.
=cut
  
  
sub test_for_orthology{
  my ($self, $transcript, $db, $focus_db, $focus_species, $compara_db, $target_db, $target_species, $threshold ) = @_;
  
  ############################################################
  my @exons = sort { $a->start <=> $b->start } @{$transcript->get_all_Exons};
  my $low    = $exons[0]->start;
  my $high   = $exons[$#exons]->end;
  my $strand = $exons[0]->strand;
  

  ############################################################
  # get a slice for the transcript with extra sequence on both sides

  my $focus_slice;
  if ($transcript->dbID){
    $focus_slice = $db->get_SliceAdaptor->fetch_by_transcript_id( $transcript->dbID, 1000 );
  }

  ############################################################
  # if no dbID, need to get a slice from the slice where transcript sits:
  if (  !$transcript->dbID || !$focus_slice  ){
    my $chr_name  = $exons[0]->contig->chr_name;
    my $chr_start = $exons[0]->contig->chr_start;
    my $chr_end   = $exons[0]->contig->chr_end;
    
    my $start = $chr_start + $low  - 1;

    my $end   = $chr_start + $high - 1;
    
    ############################################################
    # get slice from the same db where our transcript is
    $focus_slice = 
      $db->get_SliceAdaptor->fetch_by_chr_start_end( $chr_name, $start - 1000, $end + 1000);
  }
  
  ################################################################
  # find the slice(s) syntenic to this focus slice where the gene is  
  my @target_slices;
  my $target_slices = 
    $self->get_all_syntenic_slices( $focus_slice, $focus_db, $focus_species, $compara_db, $target_db, $target_species);
  if ( $target_slices ){
    @target_slices = @{$target_slices};
  }
  else{
    print STDERR "Could not find any syntenic region\n";
    return undef;
  }
  ############################################################
  # do some checks on the slices:
  #print STDERR scalar(@target_slices)." target slices found\n";
  #foreach my $slice ( @target_slices ){
  #  print STDERR $slice->chr_name.".".$slice->chr_start."-".$slice->chr_end."\n";
  #}
  
  ############################################################
  # run genewise for this gene
  my @orthologues;
  foreach my $target_slice ( @target_slices ){
    push( @orthologues, $self->align_with_exonerate( $transcript, $target_slice ) );
  }
  print STDERR "Found ".scalar(@orthologues)." orthologues\n";
  
  if ( @orthologues ){
    foreach my $ortho ( @orthologues ){
      my @exons    = sort { $a->start <=> $b->end } @{$ortho->get_all_Exons};
      my $start    = $exons[0]->contig->chr_start + $exons[0]->start - 1;
      my $end      = $exons[0]->contig->chr_start + $exons[$#exons]->end - 1;
      my $strand   = $exons[0]->strand;
      my $seqname  = $exons[0]->seqname;
      $seqname     =~ s/\.\d+-\d+$//;
      my $coverage = $self->_coverage(   $ortho );
      my $perc_id  = $self->_percent_id( $ortho );
      my $id;
      if ( $transcript->dbID || $transcript->stable_id ){
	$id = $transcript->stable_id  || $transcript->dbID;
      }
      else{
	$id = "no id";
      }
      
      print STDERR "$focus_species $id $target_species coverage:$coverage percent_id:$perc_id extent:$seqname.$start-$end strand:$strand\n";
      Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript( $ortho, 1);
    }
    
    if ($threshold){
      my @selected;
      @orthologues = sort { $self->_coverage($b) <=> $self->_coverage($a) } @orthologues;
      
      #print STDERR "best match: coverage = ".$self->_coverage($orthologues[0])."\n";
      #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript( $orthologues[0], 1 );
      
      @selected = @{$self->filter_by_coverage( \@orthologues, $threshold )};
      print STDERR scalar(@selected)." matches with >= $threshold % coverage\n";
      
      return \@selected;
    }
    else{
      return \@orthologues;
    }
  }
}

############################################################
# filter according to coverage

sub filter_by_coverage{
  my ($self, $orthologues, $coverage_cutoff ) = @_;
  my @selected;
  foreach my $ortholog ( @$orthologues ){
    my $coverage = $self->_coverage($ortholog);
    if ( $coverage >= $coverage_cutoff ){
      push( @selected, $ortholog );
    }
  }
  return \@selected;
}
  

############################################################


=head2 align_with_exonerate

  Args       :  $transcript      -> a transcript to test
                $slice           -> a slice where we want to try to align the transcript
  Description:  method to try to align a transcript on a slice with exonerte.
                There is no restriction on the origins of the slice, as only its sequece is used.
                This method is used by 'test_for_orthology' to try to
                align a transcript on a slice which is syntenic to the slice on which the transcript
                is located.
  Returntype :  a list of transcripts ( all the result from the exonerate run)
                If there is no mathces it will return an empty list
=cut


sub align_with_exonerate{
  my ($self, $transcript, $slice ) = @_;
  
  ############################################################
  # create database
  my $database = "/tmp/db_seqs.$$";
  open( DB_SEQ,">$database") || $self->throw("Could not open $database $!");
  
  my $seqout = Bio::SeqIO->new('-format' => 'Fasta',
			       '-fh'     => \*DB_SEQ);
  
  $seqout->write_seq($slice);
  close( DB_SEQ );
  
  ############################################################
  # create runnable
  my $options = " --forcegtag FALSE ";
  my $runnable = Bio::EnsEMBL::Pipeline::Runnable::NewExonerate->new(
								     -database    => $database,
								     -query_seqs  => [$transcript->seq],
								     -query_type => 'dna',
								     -target_type=> 'dna',
								     -exonerate   => 'exonerate-0.6.7',
								     -options     => $options,
								    );
  
  $runnable->run;
  my @transcripts = $runnable->output;

  ### put the results in the slice ###
  foreach my $t ( @transcripts ){
    foreach my $exon ( @{$t->get_all_Exons} ){
      $exon->contig($slice);
    }
  }
  return @transcripts;
}

############################################################
# compare upstream regions for two genes
# this requires two genes or two transcripts in slice coordinates
sub check_upstream_region_of_genes{
  my ($self, $gene1, $db1, $gene2, $db2, $length) = @_;
  unless($length){
    $length = 1000;
  }
  my @exons1  = sort { $a->start <=> $b->end } @{$gene1->get_all_Exons};
  my $strand1 = $exons1[0]->strand;
  
  my @exons2  = sort { $a->start <=> $b->end } @{$gene2->get_all_Exons};
  my $strand2 = $exons2[0]->strand;
  
  ############################################################
  # get upstream region for 1
  my $start1;
  my $upstream1;
  if ($strand1 == 1 ){
    my $chr_name  =  $exons1[0]->contig->chr_name;
    my $new_start =  $exons1[0]->contig->chr_start + $exons1[0]->start - 1 - $length;
    my $new_end   =  $exons1[0]->contig->chr_start + $exons1[0]->start - 1;
    my $seq = $db1->get_SliceAdaptor->fetch_by_chr_start_end( $chr_name, $new_start, $new_end )->seq;    
    
    my $display_id = "gene1 strand:$strand1 chr".$chr_name.
      " start:".($new_start)." end: ". ($new_end);
    $upstream1 = Bio::Seq->new( -seq      => $seq,
				-moltype  => 'dna',
				-alphabet => 'dna',
				-id       => $display_id );
  }
  else{
    my $chr_name  =  $exons1[0]->contig->chr_name;
    my $new_start =  $exons1[$#exons1]->contig->chr_start + $exons1[$#exons1]->end;
    my $new_end   =  $exons1[$#exons1]->contig->chr_start + $exons1[$#exons1]->end + $length - 1;
    my $slice= $db1->get_SliceAdaptor->fetch_by_chr_start_end( $chr_name, $new_start, $new_end )->invert;
    my $seq = $slice->seq;
    
    my $display_id = "gene1 strand:$strand1 chr".$slice->chr_name.
      " start:".($slice->chr_start)." end: ".($slice->chr_end);
    $upstream1 = Bio::Seq->new( -seq      => $seq,
				-moltype  => 'dna',
				-alphabet => 'dna',
				-id       => $display_id );
  }
  
  
  ############################################################
  # get upstream region for 2
  my $start2;
  my $upstream2;
  if ($strand2 == 1 ){
    #$start2 = $exons2[0]->start;
    #my $seq = $exons2[0]->contig->subseq( $start2 - $length, $start2 - 1, $strand2);    
    my $chr_name  =  $exons2[0]->contig->chr_name;
    my $new_start =  $exons2[0]->contig->chr_start + $exons2[0]->start - 1 - $length;
    my $new_end   =  $exons2[0]->contig->chr_start + $exons2[0]->start - 1;
    my $seq = $db2->get_SliceAdaptor->fetch_by_chr_start_end( $chr_name, $new_start, $new_end )->seq;
    
    my $display_id = "gene2 strand:$strand2 chr".$chr_name.
      " start:".($new_start)." end: ". ($new_end);
    $upstream2 = Bio::Seq->new( -seq      => $seq,
				-moltype  => 'dna',
				-alphabet => 'dna',
				-id       => $display_id );
  }
  else{
    my $chr_name  =  $exons2[0]->contig->chr_name;
    my $new_start =  $exons2[$#exons1]->contig->chr_start + $exons2[$#exons1]->end;
    my $new_end   =  $exons2[$#exons1]->contig->chr_start + $exons2[$#exons1]->end + $length - 1;
    my $slice= $db2->get_SliceAdaptor->fetch_by_chr_start_end( $chr_name, $new_start, $new_end )->invert;
    my $seq = $slice->seq;
    
    my $display_id = "gene2 strand:$strand2 chr".$slice->chr_name.
      " start:".($slice->chr_start)." end: ".($slice->chr_end);
    $upstream2 = Bio::Seq->new( -seq      => $seq,
			       -moltype  => 'dna',
			       -alphabet => 'dna',
			       -id => $display_id );
  }


  ############################################################
  # run exonerate with these two sequences
  my @output = $self->run_exonerate( $upstream1, $upstream2 );

  # output is a set of transcripts
  my $coverage;
  foreach my $output ( @output ){
    foreach my $exon ( @{$output->get_all_Exons} ){
      foreach my $evi ( @{ $exon->get_all_supporting_features} ){
	unless( $coverage){
	  $coverage = $evi->coverage;
	}
	if ( $coverage && $evi->coverage > $coverage ){
	  $coverage = $evi->coverage;
	}
      }
    }
  }
}



1;
