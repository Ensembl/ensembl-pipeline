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

use Bio::EnsEMBL::Pipeline::Runnable::NewExonerate;
use Bio::SeqIO;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Bio::EnsEMBL::Pipeline::Runnable::Blast;
use Bio::EnsEMBL::Root;

use strict;


use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Root);


############################################################

=head2 get_all_syntenic_slices

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
  my $focus_length = $focus_slice->chr_end - $focus_slice->chr_start + 1;
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
							 $focus_slice->chr_end,
							 'WGA'
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
      
      ############################################################
      # if the distance between features is smaller than the original length, we bridge them over
      if ( $chr_features[$i]->start - $chr_features[$i-1] - 1 < $focus_length ){
	$end = $chr_features[$i]->hend;
      }
      ############################################################
      # else we create a slice with the current (start,end) feature
      else{

	############################################################
	# if the slice is not as big as the original, we make it as big it's ok if they overlap
	my $length = $end - $start + 1;
	if ( $focus_length > $length ){
	  $start -= ( $focus_length - $length )/2;
	  $end   += ( $focus_length - $length )/2;
	}

        my $target_slice = $target_db->get_SliceAdaptor->fetch_by_chr_start_end($chr,$start,$end);
	push(@slices, $target_slice);

        # update to the latest feature
	$start = $chr_features[$i]->hstart;
	$end   = $chr_features[$i]->hend;
      }

      if ( $i == scalar(@chr_features) - 1 ){
         my $length = $end - $start + 1;
        if ( $focus_length > $length ){
          $start -= int( ( $focus_length - $length )/2 );
          $end   += int( ( $focus_length - $length )/2 );
        }
        my $target_slice = $target_db->get_SliceAdaptor->fetch_by_chr_start_end($chr,$start,$end);
	push(@slices, $target_slice);
      }
    }
    
    #my $length = $end - $start + 1;
    #if ( $focus_length > $length ){
    #  $start -= int( ( $focus_length - $length )/2 );
    #  $end   += int( ( $focus_length - $length )/2 );
    #}
    #my $target_slice = $target_db->get_SliceAdaptor->fetch_by_chr_start_end($chr,$start,$end);
    #push(@slices, $target_slice);
  }

  if (@slices){
    print STDERR "Produced slices:\n";
    foreach my $slice (@slices){
      print STDERR $slice->chr_name.".".$slice->chr_start."-".$slice->chr_end."\n";
    }
  }
  #print STDERR "returning ".scalar(@slices)."\n";
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

  unless ($threshold){
    $threshold = 40;
  }

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
  my $synteny_is_broken = 1;

  foreach my $target_slice (@target_slices){
    print STDERR "running exonerate with ".
      $target_slice->chr_name.".".$target_slice->chr_start."-".$target_slice->chr_end."\n";
    
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
      @left_orthology_genes  =  @{$self->filter_by_coverage( \@sorted_left, $threshold )};
    }
    
    ############################################################
    if( $right_gene ){
      my @sorted_right = sort{ $self->_coverage($b) <=> $self->_coverage($a) } $self->align_with_exonerate($right_gene, $target_slice );
      
      if (  $sorted_right[0] ){
	print STDERR "right ortholog - coverage: ".$self->_coverage( $sorted_right[0] )."\n";
	Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript( $sorted_right[0],1 );
      }
      
      @right_orthology_genes = @{$self->filter_by_coverage( \@sorted_right, $threshold )};
    }

    ############################################################
    if( $transcript ){
      my @sorted       = sort{ $self->_coverage($b) <=> $self->_coverage($a) } $self->align_with_exonerate($transcript, $target_slice );
      
      if ( $sorted[0] ){
	print STDERR "\"transcript\" ortholog - coverage: ".$self->_coverage( $sorted[0] )." \n";
	Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript( $sorted[0],1 );
      }
      @orthology_transcripts =  @{$self->filter_by_coverage( \@sorted, $threshold )};
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
    elsif (  !$left_orthology_genes[0] && !$orthology_transcripts[0] && !$right_orthology_genes[0] ){
      print STDERR "there is no synteny at all\n";
      $synteny_breaking = 1;
    }
    elsif (  $left_orthology_genes[0] && $orthology_transcripts[0] && $right_orthology_genes[0] 
	     ||
	     !$left_orthology_genes[0] && $orthology_transcripts[0] && $right_orthology_genes[0] 
	     ||
	     $left_orthology_genes[0] && $orthology_transcripts[0] && !$right_orthology_genes[0] 
	     ||
	     !$left_orthology_genes[0] && $orthology_transcripts[0] && !$right_orthology_genes[0] 
	  ){

      # one occurrence of synteny conservation is enough to flag it as preserved
      $synteny_is_broken = 0;
      #if ( $left_orthology_genes[0]->end < $orthology_transcripts[0]->start 
#	   &&
#	   $right_orthology_genes[0]->start > $orthology_transcripts[0]->end 
#	 ){
#	print STDERR "synteny is preserved\n";
#	$synteny_breaking = 0;
#      }
#      else{
#	print STDERR "synteny is reversed\n";
#	$synteny_breaking = 0;
#      }
    }

  }
  return $synteny_is_broken;
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
  my ($self, $transcript, $db, $focus_db, $focus_species, $compara_db, $target_db, $target_species, $threshold , $gene_id) = @_;
  
  # $gene_id is a hashref with transcript objects as keys and gene ids as values, handy for printing reports
  #unless( $threshold ){
  #  $threshold = 60;
  #}

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
  if ( $target_slices && @{$target_slices} ){
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
    #push( @orthologues, $self->align_with_exonerate( $transcript, $target_slice ) );
    return $self->align_with_tblastx( $transcript, $target_slice );
  }
  exit(0);
  print STDERR "Found ".scalar(@orthologues)." orthologues\n";
  
  if ( @orthologues ){
    foreach my $ortho ( @orthologues ){
      
      my $coverage = $self->_coverage(   $ortho );

      if ( $threshold && $coverage > $threshold ){
	my @exons    = sort { $a->start <=> $b->end } @{$ortho->get_all_Exons};
	my $start    = $exons[0]->contig->chr_start + $exons[0]->start - 1;
	my $end      = $exons[0]->contig->chr_start + $exons[$#exons]->end - 1;
	my $strand   = $exons[0]->strand;
	my $seqname  = $exons[0]->seqname;
	$seqname     =~ s/\.\d+-\d+$//;
	my $perc_id  = $self->_percent_id( $ortho );
	my $id;
	if ( $transcript->dbID || $transcript->stable_id ){
	  $id = $transcript->stable_id  || $transcript->dbID;
	}
	else{
	  $id = "no id";
	}
	
	my $g_id = $gene_id->{$transcript};
	
	print STDERR "$focus_species $g_id $id $target_species coverage:$coverage percent_id:$perc_id extent:$seqname.$start-$end strand:$strand\n";
	Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_SimpleTranscript( $ortho, 1);
      }
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
  
  my $id;
  if ( $transcript->dbID ){
    $id = $transcript->stable_id || $transcript->dbID;
  }
  else{
    $id = "no id";
  }
   
  my $seq = $transcript->seq;
  unless ( $seq->display_id ){
    $seq->display_id($id);
  }
  #print STDERR "transcript seq siplay_id : ".$seq->display_id."\n";

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
  #my $options = " --forcegtag FALSE ";
  my $options = ""; 
  my $runnable = Bio::EnsEMBL::Pipeline::Runnable::NewExonerate->new(
								     -database    => $database,
								     -query_seqs  => [$seq],
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

sub align_with_tblastx{
  my ($self, $transcript, $slice) =@_;
  
  my $id;
  if ( $transcript->dbID ){
    $id = $transcript->stable_id || $transcript->dbID;
  }
  else{
    $id = "no id";
  }
    
  my $seq = $transcript->seq;
  my $transcript_length = $seq->length;
  unless ( $seq->display_id ){
    $seq->display_id($id);
  }

  ############################################################
  # create database
  my $file = 'seq_'.$$.'.fa';
  my $database = "/tmp/".$file;
  open( DB_SEQ,">$database") || die("Could not open $database $!");
  
  my $seqout = Bio::SeqIO->new('-format' => 'Fasta',
			       '-fh'     => \*DB_SEQ);
  
  $seqout->write_seq($seq);
  close( DB_SEQ );
  
  system("pressdb $database");#/tmp/db_seqs.$$");


   ############################################################
  my $blast =  Bio::EnsEMBL::Pipeline::Runnable::Blast->new ('-query'     => $slice,
							     '-program'   => 'wutblastx',
							     '-database'  => $database,
							     -threshold_type => "PVALUE",
							     '-threshold' => 1e-10,
							     #'-filter'    => $filter,
							     '-options'   => 'V=1000000'
							    );

    
  $blast->add_regex($file,'(\S+)');
  $blast->run();
    
  my @featurepairs = $blast->output();
  my @pos_strand = grep { $_->strand == 1} @featurepairs;  
  my @neg_strand = grep { $_->strand == -1} @featurepairs;  

  foreach my $fp (sort{ $a->hstart <=> $b->hstart} @pos_strand) {
   print $fp->gffstring . "\n";
  }
  foreach my $fp (sort{ $a->hstart <=> $b->hstart} @neg_strand) {
   print $fp->gffstring . "\n";
  }

#  ############################################################
#  # make the blast features into likely transcript structures:
  
#  # cluster the features by genomic:
#  my @pos_clusters = @{$self->cluster_Features_by_genomic(@pos_strand)};
#  my @neg_clusters = @{$self->cluster_Features_by_genomic(@neg_strand)};
    
#  my @pos_features;
#  foreach my $cluster ( @pos_clusters ){
#     my @features = sort { $b->score <=> $a->score } @{$cluster->get_all_sub_SeqFeatures};
#     push ( @pos_features, $features[0] );
#  }

#  my @neg_features;
#  foreach my $cluster ( @neg_clusters ){
#     my @features = sort { $b->score <=> $a->score } @{$cluster->get_all_sub_SeqFeatures};
#     push ( @neg_features, $features[0] );
#  }

  # calculate the coverage
  my ($pos_trans_pos_genom_clusters, $neg_trans_pos_genom_clusters) = $self->cluster_all_Features_by_transcript(@pos_strand);
  my ($pos_trans_neg_genom_clusters, $neg_trans_neg_genom_clusters) = $self->cluster_all_Features_by_transcript(@neg_strand);

  print STDERR "length: $transcript_length\n"; 
  my $length = 0;
  foreach my $cluster ( @$pos_trans_pos_genom_clusters ){
     my ($start,$end);
     foreach my $fp ( @{$cluster} ){
       unless ( $start ){
         $start = $fp->hstart;
       }
       unless ( $end ){
         $end = $fp->hend;
       }
       if ( $fp->hstart < $start ){
         $start = $fp->hstart;
       }
       if ( $fp->hend > $end ){
         $end = $fp->hend;
       }
     }
     $length += ( $end - $start + 1 );
#     print STDERR "cluster: ".$start."-".$end."\t $length\n";
  }
  my $pos_pos_coverage = 100*$length/$transcript_length;
  print STDERR "coverage on forward transcript - positive strand = $pos_pos_coverage\n";

  $length = 0;
  foreach my $cluster ( @$neg_trans_pos_genom_clusters ){
     my ($start,$end);
     foreach my $fp ( @{$cluster} ){
       unless ( $start ){
         $start = $fp->hstart;
       }
       unless ( $end ){
         $end = $fp->hend;
       }
       if ( $fp->hstart < $start ){
         $start = $fp->hstart;
       }
       if ( $fp->hend > $end ){
         $end = $fp->hend;
       }
     }
     $length += ( $end - $start + 1 );
#          print STDERR "cluster: ".$start."-".$end."\t $length\n";
  }
  my $neg_pos_coverage = 100*$length/$transcript_length;
  print STDERR "coverage on reverse transcript - positive strand = $neg_pos_coverage\n";

  $length = 0;
  foreach my $cluster ( @$pos_trans_neg_genom_clusters ){
     my ($start,$end);
     foreach my $fp ( @{$cluster} ){
       unless ( $start ){
         $start = $fp->hstart;
       }
       unless ( $end ){
         $end = $fp->hend;
       }
       if ( $fp->hstart < $start ){
         $start = $fp->hstart;
       }
       if ( $fp->hend > $end ){
         $end = $fp->hend;
       }
     }
     $length += ( $end - $start + 1 );
#     print STDERR "cluster: ".$start."-".$end."\t $length\n";
  }
  my $pos_neg_coverage = 100*$length/$transcript_length;
  print STDERR "coverage on forward transcript - negative strand = $pos_neg_coverage\n";

  $length = 0;
  foreach my $cluster ( @$neg_trans_neg_genom_clusters ){
     my ($start,$end);
     foreach my $fp ( @{$cluster} ){
       unless ( $start ){
         $start = $fp->hstart;
       }
       unless ( $end ){
         $end = $fp->hend;
       }
       if ( $fp->hstart < $start ){
         $start = $fp->hstart;
       }
       if ( $fp->hend > $end ){
         $end = $fp->hend;
       }
     }
     $length += ( $end - $start + 1 );
 #         print STDERR "cluster: ".$start."-".$end."\t $length\n";

  }
  my $neg_neg_coverage = 100*$length/$transcript_length;
  print STDERR "coverage on reverse transcript - negative strand = $neg_neg_coverage\n";
  
  return @featurepairs;  
}

############################################################

sub cluster_Features_by_genomic{
 my ($self,@feat) = @_;

 my @features = sort{ $a->start <=> $b->start} @feat;
 
 # Create the first cluster
 my $cluster = new Bio::EnsEMBL::SeqFeature;
  
 # main cluster feature - holds all clusters
  my $cluster_list = new Bio::EnsEMBL::SeqFeature; 

 # Start off the cluster with the first exon
 $cluster->add_sub_SeqFeature($features[0]);

 $cluster->strand($features[0]->strand);    
 $cluster_list->add_sub_SeqFeature($cluster);
  
 # Loop over the rest of the features
 my $count = 0;
  
 EXON:
  foreach my $f (@features) {
    if ($count > 0) {
      
      # Add to cluster if overlap AND if strand matches
      if ( !( $f->start > $cluster->end || $f->end < $cluster->start )
           && 
           ( $f->strand == $cluster->strand) ) { 
	      $cluster->add_sub_SeqFeature($f);
      }  
      else {
	# Start a new cluster
	$cluster = new Bio::EnsEMBL::SeqFeature;
	$cluster->add_sub_SeqFeature($f);
	$cluster->strand($f->strand);
		
	# and add it to the main_cluster feature
	$cluster_list->add_sub_SeqFeature($cluster);
      }
    }
    $count++;
  }
  return $cluster_list;
}

############################################################

sub cluster_all_Features_by_transcript{
 my ($self,@feat) = @_;
 
 # separate by hstrands
 my @hpos_features;
 my @hneg_features;
 foreach my $feat ( @feat ){
   if ( $feat->hstrand == 1 ){
     push( @hpos_features, $feat );
   }
   else{
     push( @hneg_features, $feat );
   }
 }
 
 my $pos_cluster_list;
 if ( @hpos_features ){
  $pos_cluster_list = $self->cluster_Features_by_Transcript( @hpos_features );
 }
 my $neg_cluster_list; 
 if ( @hneg_features ){
   $neg_cluster_list = $self->cluster_Features_by_Transcript( @hneg_features );
 }
 return ($pos_cluster_list,$neg_cluster_list);
}

############################################################

sub cluster_Features_by_Transcript{
 my ($self,@feat) = @_;
 
 my @clusters;
 my @cluster_hstarts;
 my @cluster_hends;
 my @features = sort{ $a->hstart <=> $b->hstart} @feat;
 
 # create the first cluster
 my $count = 0;
 my $cluster = [];
  
 # start it off with the first feature
 my $first_feat = shift( @features );
 push (@$cluster, $first_feat);
 $cluster_hstarts[$count] = $first_feat->hstart;
 $cluster_hends[$count]   = $first_feat->hend;
  
 # store the list of clusters
 push(@clusters,$cluster);
  
 # loop over the rest of the features
  
 FEATURE:
  foreach my $f ( @features ){
    #print STDERR "trying to place feature: $f ".$f->start."-".$f->end."\n";    
    # add $f to the current cluster if overlaps and strand matched
    #print STDERR "comparing with cluster $count : "
    #  .$cluster_hstarts[$count]."-".$cluster_hends[$count]."\n";
    
    if (!($f->hend < $cluster_hstarts[$count] || $f->hstart > $cluster_hends[$count])) {      
      push(@$cluster,$f);
      
      # re-adjust size of cluster
      if ($f->hstart < $cluster_hstarts[$count]) {
	$cluster_hstarts[$count] = $f->hstart;
      }
      if ($f->hend  > $cluster_hends[$count]) {
	$cluster_hends[$count] = $f->hend;
      }
      
    }
    else{
      # else, start create a new cluster with this feature
      $count++;
      $cluster = [];
      push (@$cluster, $f);
      $cluster_hstarts[$count] = $f->hstart;
      $cluster_hends[$count]   = $f->hend;
          
      # store it in the list of clusters
      push(@clusters,$cluster);
    }
  }
  return \@clusters;
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

############################################################

sub run_genewise{
  my ($self, $transcript, $slices) =@_;
  my @slices = @$slices;
  
  my $tseq;
  eval{
    $tseq = $transcript->translate();
  };
  foreach my $slice ( @slices){
    
    print STDERR "running genewise on slice ".$slice->name."\n";
    my $genewise = new Bio::EnsEMBL::Pipeline::Runnable::Genewise(  -genomic => $slice,
								    -protein => $tseq,
								    -reverse => 1,
								    -endbias => 1,);
    
    $genewise->run;
    my @features = $genewise->output;
    
    foreach my $feature (@features){
      print STDERR $feature->gffstring."\n";
    }
  }
}

############################################################

sub run_blast{
  my ($self, $gene, $slices) =@_;
  my @slices = @$slices;
  
  my @transcripts = sort { $b->length <=> $a->length} @{$gene->get_all_Transcripts};
  my $tseq = $transcripts[0]->translate();
  
  
  foreach my $slice ( @slices){
    
    my $blastdb = Bio::EnsEMBL::Pipeline::Runnable::BlastDB->new(
								 -sequences => [$slice],
								 -type      => 'DNA',
								);
    
    $blastdb->run;
    my $dbname = $blastdb->dbname;
    
    my $blast   = Bio::EnsEMBL::Pipeline::Runnable::Blast->new (
								'-query'     => $tseq,
								'-program'   => 'wublastp',
								'-database'  => $dbname,
								'-threshold' => 1e-6,
								'-filter'    => 1,
								'-options'   => 'V=1000000',
							       );
    
    $blast->run();
    
    my @featurepairs = $blast->output();
    
    foreach my $fp (@featurepairs) {
      print $fp->gffstring . "\n";
    }
  }
}

############################################################


1;
