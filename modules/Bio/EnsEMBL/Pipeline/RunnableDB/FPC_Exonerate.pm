#
#
# Cared for by Val Curwen (vac@sanger.ac.uk)
#
# Copyright Val Curwen
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::FPC_Exonerate

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::RunnableDB::FPC_Exonerate->new(
					     -dbobj     => $db,
					     -input_id  => $id,
					     -analysis  => $analysis
                                             );
    $obj->fetch_input
    $obj->run

    my @newfeatures = $obj->output;

    $obj->write_output();

=head1 DESCRIPTION

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::RunnableDB::FPC_Exonerate;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::Exonerate;
use Bio::EnsEMBL::Pipeline::SeqFetcher;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Gene;
use Bio::SeqIO;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB );

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
           
    $self->{'_fplist'} = []; #create key to an array of feature pairs

    $self->throw("Analysis object required") unless ($self->analysis);

    return $self;
}

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data for exonerate
    Returns :   nothing
    Args    :   none

=cut

sub fetch_input {
    my ($self) = @_;

    print STDERR "Fetching input \n";
    $self->throw("No input id") unless defined($self->input_id);

    # eg ctg25118
    my $contigid  = $self->input_id;

    $self->dbobj->static_golden_path_type('UCSC');

    my $stadaptor = $self->dbobj->get_StaticGoldenPathAdaptor();

    # returns all the contigs making up the fpc contig munged together as a list of virtual contigs.
    my @contig    = $stadaptor->fetch_VirtualContig_list_sized($contigid,500000,10000,1000000,100);

    foreach my $contig (@contig){
	print STDERR "Analysing contig " . $contig->id . "\n";
	foreach my $rc ($contig->_vmap->each_MapContig) {
	    my $strand = "+";
	    if ($rc->orientation == -1) {
		$strand = "-";
	    }
	    print STDERR "orientation: $strand\n";

	}

	my @ests = $self->_get_ests($contig);
	my @genomic = ( $contig->get_repeatmasked_seq() );

	if(scalar(@ests) && scalar(@genomic))
	{
	    my $executable =  $self->analysis->program_file();
	    my $exonerate = new Bio::EnsEMBL::Pipeline::Runnable::Exonerate('-genomic'  => \@genomic,
									    '-est'      => \@ests,
									    '-exonerate' => $executable);
	    $self->runnable($exonerate);
	    $self->{$exonerate} = $contig;
	}
	else { print STDERR "No ests to be analysed for " . $contig->id()  . "\n"; }
    }
}

=head2 output

    Title   :   output
    Usage   :   $self->output()
    Function:   Returns the contents of $self->{_output}, 
                which holds predicted genes.
    Returns :   Array of Bio::EnsEMBL::Gene
    Args    :   None

=head2 run

    Title   :   run
    Usage   :   $self->run()
    Function:   Runs the exonerate analysis, producing 
                Bio::EnsEMBL::Gene predictions
    Returns :   Nothing, but $self{_output} contains the predicted genes.
    Args    :   None

=cut

sub run {
  my ($self) = @_;

  $self->throw("Can't run - no runnable objects") unless defined($self->{_runnables});
  
  foreach my $runnable ($self->runnable) {
        $runnable->run;
  }

  $self->_convert_output();

}

=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output
    Function:   Writes output data to db
    Returns :   array of exons (with start and end)
    Args    :   none

=cut

sub write_output {

  my($self) = @_;
  
  my @features = $self->output();
  my $db = $self->dbobj();
  my $gene_obj = $db->gene_Obj;
  
  my $EXON_ID_SUBSCRIPT       = "EXOE";
  my $TRANSCRIPT_ID_SUBSCRIPT = "EXOT";
  my $GENE_ID_SUBSCRIPT       = "EXOG";
  my $PROTEIN_ID_SUBSCRIPT    = "EXOP";
  
  my $sth = $db->prepare("lock tables gene write, genetype write, exon write, transcript write, exon_transcript write, translation write,dna read,contig read,clone read,feature read,analysis write");

  $sth->execute;
  
  eval {
    (my $gcount = $gene_obj->get_new_GeneID($GENE_ID_SUBSCRIPT))
      =~ s/$GENE_ID_SUBSCRIPT//;
    (my $tcount = $gene_obj->get_new_TranscriptID($TRANSCRIPT_ID_SUBSCRIPT))
      =~ s/$TRANSCRIPT_ID_SUBSCRIPT//;
    (my $pcount = $gene_obj->get_new_TranslationID($PROTEIN_ID_SUBSCRIPT))
      =~ s/$PROTEIN_ID_SUBSCRIPT//;
    (my $ecount = $gene_obj->get_new_ExonID($EXON_ID_SUBSCRIPT))
      =~ s/$EXON_ID_SUBSCRIPT//;
    
    
    foreach my $gene (@features) {
      $gene->id($GENE_ID_SUBSCRIPT . $gcount);
      $gcount++;
      
      # Convert all exon ids and save in a hash
      my %namehash;
      
      foreach my $ex ($gene->each_unique_Exon) {
	$namehash{$ex->id} = $EXON_ID_SUBSCRIPT.$ecount;
	$ex->id($EXON_ID_SUBSCRIPT.$ecount);
	$ecount++;
      }
      
      foreach my $tran ($gene->each_Transcript) {
	$tran->id($TRANSCRIPT_ID_SUBSCRIPT . $tcount);
	$tran->translation->id($PROTEIN_ID_SUBSCRIPT . $pcount);
	
	my $translation = $tran->translation;
	
	$tcount++;
	$pcount++;
	
	foreach my $ex ($tran->each_Exon) {
	  my @sf = $ex->each_Supporting_Feature;

	  if (defined($namehash{$translation->start_exon_id}) && $namehash{$translation->start_exon_id} ne "") {
	    $translation->start_exon_id($namehash{$translation->start_exon_id});
	  }
	  if (defined($namehash{$translation->end_exon_id}) &&$namehash{$translation->end_exon_id} ne "") {
	    $translation->end_exon_id  ($namehash{$translation->end_exon_id});
	  }
	}
	
      }
      
      $gene_obj->write($gene);
    }
  };
  if ($@) {
    $sth = $db->prepare("unlock tables");
    $sth->execute;
    
    $self->throw("Error writing gene for " . $self->input_id . " [$@]\n");
  } else {
    $sth = $db->prepare("unlock tables");
    $sth->execute;
  }
  return 1;
}

=head2 runnable

    Title   :   runnable
    Usage   :   $self->runnable($arg)
    Function:   Sets a runnable for this RunnableDB
    Returns :   Bio::EnsEMBL::Pipeline::RunnableI
    Args    :   Bio::EnsEMBL::Pipeline::RunnableI

=head2 _convert_output

    Title   :   _convert_output
    Usage   :   $self->_convert_output
    Function:   converts exons found by exonerate into Bio::EnsEMBL::Gene 
                ready to be stored in the database
    Returns :   Nothing, but  $self->{_output} contains the Gene objects
    Args    :   None

=cut

# not presently doing anything with gene, intron & splice site features. Quite possibly 
# we ought to ...
sub _convert_output {
  my ($self) =@_;
  my $count=1;

  # make an array of genes for each runnable
  foreach my $runnable ($self->runnable) {
    my $contig = $self->{$runnable};
    my @features   = $runnable->output;
    my %homols = ();
    my @genes = ();
    
    # sort the hits into a hash keyed by hid
    foreach my $f(@features) {
      push (@{$homols{$f->hseqname()}}, $f);
    }
    
    # make one gene per hid, using the exons predicted by exonerate
    foreach my $id (keys %homols) {
      my @exonfeat;
      foreach my $ft( @{$homols{$id}}) {
	if($ft->primary_tag eq 'exon') {
	  push(@exonfeat, $ft);
	}	  
      }
	   
      my $gene   = $self->_make_gene($contig, $count, @exonfeat);
      if(defined($gene)) {
	push (@genes, $gene);
	$count++;
      }
    }

    # map genes back to genomic coordinates
    my @remapped = $self->_remap_genes($contig,@genes);
       
    if (!defined($self->{_output})) {
      $self->{_output} = [];
    }
	
    push(@{$self->{_output}},@remapped);
  }
  
}


=head2 _get_ests

    Title   :   _get_ests
    Usage   :   $self->_get_ests(@features)
    Function:   Screens FeaturePairs in @features for vert EST blast hits, retrieves sequences and
                makes them into an array of Bio::Seq
    Returns :   Array of Bio::EnsEMBL::Seq
    Args    :   None

=cut

sub _get_ests {
  my ($self,$contig) = @_;
  my @mrnafeatures = ();
  my %idhash = ();
  
  foreach my $f ($contig->get_all_SimilarityFeatures()) {
    if (defined($f->analysis)      && defined($f->score) && 
	defined($f->analysis->db)  && $f->analysis->db eq "vert") {
      
      if (!defined($idhash{$f->hseqname})) { 
	push(@mrnafeatures,$f);
	$idhash{$f->hseqname} =1;
      } 
      else {
	# feature pair on this sequence already seen
	print STDERR ("Ignoring feature " . $f->hseqname . "\n");
      }
    }
  }
  
  unless (@mrnafeatures)
    {
      print STDERR ("No EST hits\n");
      return;
    }
  my @seq = $self->_get_Sequences(@mrnafeatures);

  # validate sequences
  my @valid_seq   = $self->_validate_sequence(@seq);
  
  return @valid_seq;
}

=head2 _make_gene

    Title   :   _make_gene
    Usage   :   $self->_make_gene($contig,$count,@exonfeat)
    Function:   Converts exon features predicted by exonerate into a gene
    Returns :   Bio::EnsEMBL::Gene
    Args    :   Bio::EnsEMBL::Virtual::Contig, integer, array of Bio::EnsEMBL::FeaturePair

=cut

sub _make_gene {
  my ($self, $contig, $count, @exonfeat) = @_;
  my $excount = 1;
  my @exons = ();
  my $time   = time; 
  chomp($time); 

  my $gene   = new Bio::EnsEMBL::Gene;
  $gene->type('fpc_exonerate');
  $gene->id($self->input_id . ".exo.$count");
  $gene->created($time);
  $gene->modified($time);
  $gene->version(1);
  
  # make a transcript
  my $tran   = new Bio::EnsEMBL::Transcript;
  $tran->id($self->input_id . ".exo.$count");
  $tran->created($time);
  $tran->modified($time);
  $tran->version(1);
  
  my $transl = new Bio::EnsEMBL::Translation;
  $transl->id($self->input_id . ".exo.$count");
  $transl->version(1);
  
  # add transcript to gene
  $gene->add_Transcript($tran);
  $tran->translation($transl);
  
  foreach my $fp(@exonfeat) {
    $fp->analysis->gff_source('fpc_exonerate');
    
    my $exon = new Bio::EnsEMBL::Exon;
    $exon->id($self->input_id . ".exo.$count.$excount");
    $exon->contig_id($contig->id);
    $exon->created($time);
    $exon->modified($time);
    $exon->version(1);
    
    $exon->start($fp->start);
    $exon->end  ($fp->end);
    
    if($fp->strand == $fp->hstrand) {
      $fp->strand(1);
      $fp->hstrand(1);
    }
    else {
      $fp->strand(-1);
      $fp->hstrand(-1);
    }
    
    $exon->strand($fp->strand);
    $exon->phase(0);
    $exon->attach_seq($contig->primary_seq);
    $exon->add_Supporting_Feature($fp);
  
  
    push(@exons,$exon);
    $excount++;
  }

  if ($#exons < 0) {
    print STDERR "Odd.  No exons found\n";
    return;
  } 
  else {
    # assemble exons
    if ($exons[0]->strand == -1) {
      @exons = sort {$b->start <=> $a->start} @exons;
    } 
    else {
      @exons = sort {$a->start <=> $b->start} @exons;
    }
    
    foreach my $exon(@exons){
      $tran->add_Exon($exon);
    }
    
    $transl->start_exon_id($exons[0]->id);
    $transl->end_exon_id  ($exons[$#exons]->id);
    
    if ($exons[0]->strand == 1) {
      $transl->start($exons[0]->start);
      $transl->end  ($exons[$#exons]->end);
    } else {
      $transl->start($exons[0]->end);
      $transl->end  ($exons[$#exons]->start);
    }
    return $gene;
  }
}

=head2 _remap_genes

    Title   :   _remap_genes
    Usage   :   $self->_remap_genes($contig,@genes)
    Function:   Remaps predicted genes into genomic coordinates
    Returns :   array of Bio::EnsEMBL::Gene
    Args    :   Bio::EnsEMBL::Virtual::Contig, array of Bio::EnsEMBL::Gene

=cut

sub _remap_genes {
  my ($self, $contig, @genes) = @_;
  my @remapped;

  foreach my $gene(@genes) {
    eval {
      # check supporting feature data
      print STDERR "CHECKING SUPPORTING FEATURES BEFORE REMAP\n";
      foreach my $t($gene->each_Transcript){
	foreach my $e($t->each_Exon){
	  foreach my $sf($e->each_Supporting_Feature){
	    $self->_print_FeaturePair($sf);
	  }
	}
      }
	
      my $newgene = $contig->convert_Gene_to_raw_contig($gene);
      $newgene->type('fpc_exonerate');
      push(@remapped,$newgene);
      # check supporting feature data
      print STDERR "CHECKING SUPPORTING FEATURES AFTER REMAP\n";
      foreach my $t($newgene->each_Transcript){
	foreach my $e($t->each_Exon){
	  foreach my $sf($e->each_Supporting_Feature){
	    $self->_print_FeaturePair($sf);
	  }
	}
      }
      
    };

    if ($@) {
      print STDERR "Couldn't reverse map gene " . $gene->id . " [$@]\n";
    }
    
  }

  return @remapped;

}

###############################
# Sequence fetching routines
###############################

# This is shared with various other modules; need a sequence fetching module ...

=head2 _get_Sequences
  
    Title   :   _get_Sequences
    Usage   :   $self->_get_Sequences(@features)
    Function:   Gets a Bio::Seq for each of the hit sequences in @features
    Returns :   Array of Bio::Seq
    Args    :   Array of Bio::EnsEMBL::FeaturePair

=cut
  
sub _get_Sequences {
  my ($self,@pairs) = @_;
  
  my @seq;
  
  foreach my $pair (@pairs) {
    my $id = $pair->hseqname;
    if ($pair->analysis->db eq "vert") {
      
      eval {
	my $seq = $self->_get_Sequence($id);
	push(@seq,$seq);
      };
      if ($@) {
	$self->warn("Couldn't fetch sequence for $id [$@]");
      } 
    }
  }
  return @seq;
}

=head2 _get_Sequence

  Title   : _get_Sequence
  Usage   : my $seq = _get_Sequence($id);
  Function: Fetches sequence that has ID $id. Tries a variety of methods for 
            fetching sequence. If sequence found, it is cached.
  Returns : Bio::PrimarySeq
  Args    : ID string 

=cut

sub _get_Sequence {
  my ($self,$id) = @_;
  if (defined($self->{_seq_cache}{$id})) {
    return $self->{_seq_cache}{$id};
  } 
  
  my $seqfetcher = new Bio::EnsEMBL::Pipeline::SeqFetcher;
  my $seq = $seqfetcher->run_efetch($id);

 if(!defined $seq) {
    $seq = $seqfetcher->run_getz($id, 'embl');
  }
  
  if(!defined $seq) {
    $self->throw("Couldn't find sequence for [$id]");
  }

  $self->{_seq_cache}{$id} = $seq;
  
  return $seq;
  
}

=head2 _validate_sequence

    Title   :   _validate_sequence
    Usage   :   $self->_validate_sequence(@seq)
    Function:   Takes an array of Seq or PrimarySeq objects and 
                returns valid ones, removing invalid characters 
                for nucleic acid sequences and rejecting sequences 
                that are not nucleic acid
    Returns :   Array of Bio::Seq
    Args    :   Array of Bio::Seq

=cut

sub _validate_sequence {
  my ($self, @seq) = @_;
  my @validated;
  foreach my $seq (@seq)
    {
      print STDERR ("$seq is not a Bio::PrimarySeq or Bio::Seq\n") 
	unless ($seq->isa("Bio::PrimarySeq") ||
		$seq->isa("Bio::Seq"));
      my $sequence = $seq->seq;
      if ($sequence !~ /[^acgtn]/i)
        {
	  push (@validated, $seq);
        }
      else 
        {
	  $_ = $sequence;
	  my $len = length ($_);
	  my $invalidCharCount = tr/mrwsykvhdbxMRWSYKVHDBX/n/;
	  #extract invalid characters
	  $sequence =~ s/[ACGTN]//ig;
	  if ($invalidCharCount / $len > 0.05)
            {
	      $self->warn("Ignoring ".$seq->display_id()
			  ." contains more than 5% ($invalidCharCount) "
			  ."odd nucleotide codes ($sequence)\n Type returns "
			  .$seq->moltype().")\n");
            }
	  else
            {
	      $self->warn ("Cleaned up ".$seq->display_id
			   ." for blast : $invalidCharCount invalid chars ($sequence)\n");
	      $seq->seq($_);
	      push (@validated, $seq);
            }
        }
    } 
  return @validated;  
}

=head2 _get_tmp_file

    Title   :   _get_tmp_file
    Usage   :   $self->_get_tmp_file($dir,$stub,$ext)
    Function:   Generates a unique file name in directory $dir with combination 
                of $stub, a random number, and $ext
    Returns :   Array of Bio::EnsEMBL::FeaturePair
    Args    :   None

=cut

sub _get_tmp_file {
  my ($self,$dir,$stub,$ext) = @_;
  
  
  if ($dir !~ /\/$/) {
    $dir = $dir . "/";
  }
  
  # but just ignores the return value!!!
  $self->_check_disk_space($dir);
  
  my $num = int(rand(10000));
  my $file = $dir . $stub . "." . $num . "." . $ext;
  
  # this introduces a race condition; suspect it's the way it's done all over EnsEMBL; could be a problem with big pipeline runs?
  while (-e $file) {
    $num = int(rand(10000));
    $file = $stub . "." . $num . "." . $ext;
  }			
  
  return $file;
}

=head2 _check_disk_space

    Title   :   _check_disk_space
    Usage   :   $self->_check_disk_space($dir, $minimumkb)
    Function:   Checks space in $dir
    Returns :   1 if "enough" space, 0 otherwise
    Args    :   $dir (directory name), $minimumkb (minimumkb needed)

=cut

sub _check_disk_space {
  my ($self,$dir,$minimumkb) = @_;
  
  # what if $minimumkb is unset? check for undef?
  
  $self->throw("No directory entered") unless defined($dir);
  
  open(DF,"df -k $dir |");
  
  my @lines = <DF>;
  $self->throw("Wrong number of lines output from df") unless scalar(@lines) == 2;
  my @f = split(' ',$lines[1]);
  
  my $kbytes = $f[3];
  
  if ($minimumkb && $kbytes > $minimumkb) {
    return 1;
  } else {
    return 0;
  }
}

=head2 _print_FeaturePair

    Title   :   _print_FeaturePair
    Usage   :   $self->_print_FeaturePair($pair)
    Function:   Prints attributes of a Bio::EnsEMBL::FeaturePair
    Returns :   Nothing
    Args    :   A Bio::EnsEMBL::FeaturePair

=cut

sub _print_FeaturePair {
  my ($self,$pair) = @_;
  
  print STDERR $pair->seqname . "\t" . $pair->start . "\t" .
               $pair->end . "\t" . $pair->score . "\t" . 
               $pair->strand . "\t" . $pair->hseqname . "\t" . 
	       $pair->hstart . "\t" . $pair->hend . "\t" . 
	       $pair->hstrand . "\n";
}

1;


