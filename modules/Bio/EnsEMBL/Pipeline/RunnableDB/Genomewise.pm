#
#
# Cared for by EnsEMBL  <ensembl-dev@ebi.ac.uk>
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::Genomewise

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::RunnableDB::Genomewise->new(
					     -dbobj     => $db,
					     -input_id  => $id
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

package Bio::EnsEMBL::Pipeline::RunnableDB::Genomewise;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::RootI;
use Bio::EnsEMBL::Pipeline::RunnableDBI;
use Bio::EnsEMBL::Pipeline::Runnable::Genomewise;
use Bio::EnsEMBL::Pipeline::Runnable::Est2Genome;
use Bio::EnsEMBL::Pipeline::Runnable::Exonerate;
use Bio::EnsEMBL::Pipeline::SeqFetcher::BPIndex;

use Bio::EnsEMBL::Pipeline::GeneConf qw (EXON_ID_SUBSCRIPT
					 TRANSCRIPT_ID_SUBSCRIPT
					 GENE_ID_SUBSCRIPT
					 PROTEIN_ID_SUBSCRIPT
					 );

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);

    # dbobj input_id mandatory and read in by BlastableDB
    if (!defined $self->seqfetcher) {
      my $seqfetcher = new Bio::EnsEMBL::Pipeline::Seqfetcher::BPIndex( '-index'  => '/work2/vac/genomewise/uni.test.inx',
									'-format' => 'Fasta');
      $self->seqfetcher($seqfetcher);
    }
    return $self; 
}


sub write_output {
    my($self,@features) = @_;

    my $db = $self->dbobj;
   
    if( !defined $db ) {
      $self->throw("unable to make write db");
    }
    
    my %contighash;
    my $gene_obj = $db->gene_Obj;

    my @newgenes = $self->output;
    return unless ($#newgenes >= 0);

    # get new ids
    eval {

	my $genecount  = 0;
	my $transcount = 0;
	my $translcount = 0;
	my $exoncount  = 0;

	# get counts of each type of ID we need.

	foreach my $gene ( @newgenes ) {
	    $genecount++;
	    foreach my $trans ( $gene->each_Transcript ) {
		$transcount++;
		$translcount++;
	    }
	    foreach my $exon ( $gene->each_unique_Exon() ) {
		$exoncount++;
	    }
	}

	# get that number of ids. This locks the database

	my @geneids  =  $gene_obj->get_New_external_id('gene',$GENE_ID_SUBSCRIPT,$genecount);
	my @transids =  $gene_obj->get_New_external_id('transcript',$TRANSCRIPT_ID_SUBSCRIPT,$transcount);
	my @translids =  $gene_obj->get_New_external_id('translation',$PROTEIN_ID_SUBSCRIPT,$translcount);
	my @exonsid  =  $gene_obj->get_New_external_id('exon',$EXON_ID_SUBSCRIPT,$exoncount);

	# database locks are over.

	# now assign ids. gene and transcripts are easy. Exons are harder.
	# the code currently assummes that there is one Exon object per unique
	# exon id. This might not always be the case.

	foreach my $gene ( @newgenes ) {
	    $gene->id(shift(@geneids));
	    my %exonhash;
	    foreach my $exon ( $gene->each_unique_Exon() ) {
		my $tempid = $exon->id;
		$exon->id(shift(@exonsid));
		$exonhash{$tempid} = $exon->id;
	    }
	    foreach my $trans ( $gene->each_Transcript ) {
		$trans->id(shift(@transids));
		$trans->translation->id(shift(@translids));
		$trans->translation->start_exon_id($exonhash{$trans->translation->start_exon_id});
		$trans->translation->end_exon_id($exonhash{$trans->translation->end_exon_id});
	    }
	    
	}

	# paranoia!
	if( scalar(@geneids) != 0 || scalar(@exonsid) != 0 || scalar(@transids) != 0 || scalar (@translids) != 0 ) {
	    $self->throw("In id assignment, left with unassigned ids ".scalar(@geneids)." ".scalar(@transids)." ".scalar(@translids)." ".scalar(@exonsid));
	}

    };
    if( $@ ) {
	$self->throw("Exception in getting new ids. Exiting befor write\n\n$@" );
    }


    # this now assummes that we are building on a single VC.



  GENE: foreach my $gene (@newgenes) {	
      # do a per gene eval...
      eval {
	
	  $gene_obj->write($gene);
      }; 
      if( $@ ) {
	  print STDERR "UNABLE TO WRITE GENE\n\n$@\n\nSkipping this gene\n";
      }
	    
  }
   
}

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data from the database
    Returns :   nothing
    Args    :    string: chr1.1-10000

=cut

sub fetch_input {
    my( $self) = @_;
    
    print STDERR "Fetching input: " . $self->input_id. " \n";
    $self->throw("No input id") unless defined($self->input_id);

    my $chrid  = $self->input_id;
       $chrid =~ s/\.(.*)-(.*)//;

    my $chrstart = $1;
    my $chrend   = $2;

    print STDERR "Chromosome id = $chrid , range $chrstart $chrend\n";

    $self->dbobj->static_golden_path_type('UCSC');

    my $stadaptor = $self->dbobj->get_StaticGoldenPathAdaptor();
    my $contig    = $stadaptor->fetch_VirtualContig_by_chr_start_end($chrid,$chrstart,$chrend);

    $contig->_chr_name($chrid);

    foreach my $rc ($contig->_vmap->each_MapContig) {
	my $strand = "+";
	if ($rc->orientation == -1) {
	    $strand = "-";
	}
	
	print STDERR $rc->contig->id . "\tsequence\t" . $rc->contig->id . "\t" . $rc->start . "\t" . $rc->end . "\t100\t" . $strand . "\t0\n";
      }
    
    my $genseq    = $contig->get_repeatmasked_seq;
    
    print STDERR "Length is " . $genseq->length . "\n";
    print STDERR "Fetching features \n";
    print STDERR "contig: " . $contig . " \n";

    my @features  = $contig->get_all_SimilarityFeatures_above_score('unigene.seq',200);
    
    print STDERR "Number of features = " . scalar(@features) . "\n";

    my %idhash;
    
    foreach my $f (@features) {
        print "Feature " . $f . " " . $f->seqname . " " . $f->source_tag . "\n";
      if ($f->isa("Bio::EnsEMBL::FeaturePair") && 
	  defined($f->hseqname) && $f->hseqname =~ /Hs/) {
	# only want human sequences
      $idhash{$f->hseqname} = 1;
    }
  }
    
    my @ids = keys %idhash;

    print STDERR "Feature ids are @ids\n";

    my $seqfetcher = new Bio::EnsEMBL::Pipeline::SeqFetcher;

    foreach my $id(@ids){
      my $seq = $self->seqfetcher->get_Seq_by_acc($id);

      my $e2g = Bio::EnsEMBL::Pipeline::Runnable::Est2Genome->new( -genomic => $contig->primary_seq,
								   -est => $seq
								 );
      $self->add_e2g_Runnable($e2g);
    }
    

    my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::Genomewise();
    $runnable->seq($contig->primary_seq);
    $self->gw_runnable($runnable);

    # at present, we'll only ever have one ...
    $self->vc($contig);
}
  

sub gw_runnable{
  my ($self, $value) = @_;

  if (defined($value)) {
    $self->{'_gw_runnable'} = $value;
  }

  return $self->{'_gw_runnable'};
  
} 

sub add_e2g_Runnable {
  my ($self,$arg) = @_;
  
  if (!defined($self->{'_e2g_runnables'})) {
    $self->{'_e2g_runnables'} = [];
  }
  
  if (defined($arg)) {
    if ($arg->isa("Bio::EnsEMBL::Pipeline::RunnableI")) {
      push(@{$self->{'_e2g_runnables'}},$arg);
    } else {
      $self->throw("[$arg] is not a Bio::EnsEMBL::Pipeline::RunnableI");
    }
  }
}

sub get_e2g_runnables {
  my ($self) = @_;
  
  if (!defined($self->{'_e2g_runnables'})) {
    $self->{'_e2g_runnables'} = [];
  }
  
  return @{$self->{'_e2g_runnables'}};
}


sub run {
  my ($self) = @_;

# run est2genomes  
  foreach my $runnable ($self->get_e2g_runnables) {
    $runnable->run;
    $self->convert_e2g_output($runnable);
  }
  
  # now run genomewise
  my $gw_runnable = $self->gw_runnable;
  my @e2g_transcripts = @{$self->{'_e2g_transcripts'}};
  foreach my $t(@e2g_transcripts){
    $gw_runnable->add_Transcript($t);
  }
  $gw_runnable->run;

  
}

sub convert_e2g_output {
  my ($self, $runnable) = @_;
  my $genetype = 'e2g';
  my $count = 1;
  my $time  = time; chomp($time);
  my @results = $runnable->output;

  my @transcripts;
  
  foreach my $t(@results) {
    my $transcript = new Bio::EnsEMBL::Transcript;
    $transcript->id($self->vc->id . ".$genetype.$count");
    $transcript->version(1);

    my @exons;
    my $excount = 1;

    foreach my $exon_pred($t->sub_SeqFeature){
      # e2g has no concept of phase, but remapping will fail if this is unset
      #      $ex->phase(-1);
      $exon_pred->phase(0);

      # make an exon
      my $exon = new Bio::EnsEMBL::Exon;
      
      $exon->id($self->vc->id . ".$genetype.$count.$excount");
      $exon->contig_id($self->vc->id);
      $exon->created($time);
      $exon->modified($time);
      $exon->version(1);
      
      $exon->start($exon_pred->start);
      $exon->end  ($exon_pred->end);
      $exon->strand($exon_pred->strand);
      
      $exon->phase($exon_pred->{_phase});
      $exon->attach_seq($self->vc->primary_seq);

      # supporting feature data goes here

      push(@exons, $exon);
      $excount++
    }

    if ($#exons < 0) {
      print STDERR "Odd.  No exons found\n";
    } 
    else {
      
      print STDERR "num exons: " . scalar(@exons) . "\n";
      
      if ($exons[0]->strand == -1) {
	@exons = sort {$b->start <=> $a->start} @exons;
      } else {
	@exons = sort {$a->start <=> $b->start} @exons;
      }
      
      foreach my $exon (@exons) {
	$transcript->add_Exon($exon);
      }
      
      push (@transcripts, $transcript);
    }
  }
    
    
  
  if (!defined($self->{'_e2g_transcripts'})) {
    $self->{'_e2g_transcripts'} = [];
  }
  
  print STDERR "e2g transcripts: " . scalar(@transcripts) . "\n";
  
  push(@{$self->{'_e2g_transcripts'}},@transcripts);  
}
  
sub convert_output {
  my ($self) =@_;
  
  my $gwr = $self->gw_runnable;

  my @transcripts = $gwr->output;
 
  my @genes = $self->make_genes;

  # check translations
  foreach my $gene(@genes){
    foreach my $trans ( $gene->each_Transcript ) {
      print STDERR "translation: \n";
      my $seqio = Bio::SeqIO->new(-fh => \*STDERR);
      print STDERR "checktrans: ";
      $seqio->write_seq($trans->translate); 
      print STDERR "\n ";	
    }
  }
  
  my @remapped = $self->remap_genes('genomewise', \@genes);
 # check translations
  foreach my $gene(@remapped){
    foreach my $trans ( $gene->each_Transcript ) {
      print STDERR "translation: \n";
      my $seqio = Bio::SeqIO->new(-fh => \*STDERR);
      print STDERR "checktrans: ";
      $seqio->write_seq($trans->translate); 
      print STDERR "\n ";	
    }
  }
  # store genes
  
  if (!defined($self->{'_output'})) {
    $self->{'_output'} = [];
  }

  print STDERR "remmapped genes: " . scalar(@remapped) . "\n";

  push(@{$self->{'_output'}},@remapped);
}



sub make_genes {

  my ($self) = @_;
  my $count = 0;
  my @genes;
  my $genetype = 'go';
  my $time  = time; chomp($time);

  my @trans = $self->gw_runnable->output;
  print "transcripts: " . scalar(@trans) . "\n";

  foreach my $transcript ($self->gw_runnable->output) {
    $count++;
    my $gene   = new Bio::EnsEMBL::Gene;
    $gene->type($genetype);
    $gene->id($self->vc->id . ".$genetype.$count");
    $gene->version(1);
    
    # add transcript to gene
    $transcript->id($self->vc->id . ".$genetype.$count");
    $gene->add_Transcript($transcript);


    # and store it
    push(@genes,$gene);

    # sort the exons 
    $transcript->sort;
    my $excount = 1;
    my @exons = $transcript->each_Exon;

    foreach my $exon(@exons){
      $exon->id($self->vc->id . ".$genetype.$count.$excount");
      $exon->contig_id($self->vc->id);
      $exon->created($time);
      $exon->modified($time);
      $exon->attach_seq($self->vc->primary_seq);
      $excount++;
      }
    
    print STDERR "exons: " . scalar(@exons) . "\n";

    # sort out translation
    my $translation  = new Bio::EnsEMBL::Translation;    
    $translation->id($self->vc->id . ".$genetype.$count");
    $translation->version(1);    
    
    $translation->start_exon_id($exons[0]->id);

    $translation->end_exon_id  ($exons[$#exons]->id);
    if ($exons[0]->phase == 0) {
      $translation->start(1);
    } elsif ($exons[0]->phase == 1) {
      $translation->start(3);
    } elsif ($exons[0]->phase == 2) {
      $translation->start(2);
    }
    $translation->end  ($exons[$#exons]->end - $exons[$#exons]->start + 1);

    $transcript->translation($translation);
  }
  return @genes;
}

sub remap_genes {
  my ($self,$genetype,$genes) = @_;
  my $contig = $self->vc;
  
  print STDERR "genes before remap: " . scalar(@$genes) . "\n";
  
  my @newf;
  my $trancount=1;
  foreach my $gene (@$genes) {
    eval {
      my $newgene = $contig->convert_Gene_to_raw_contig($gene);
      $newgene->type($genetype);
      foreach my $tran ($newgene->each_Transcript) {
	foreach my $exon($tran->each_Exon) {
	  if ($exon->isa('Bio::EnsEMBL::StickyExon')){
	    # need to deal with component exons
	    foreach my $ce($exon->each_component_Exon){
	      print STDERR $ce->contig_id . "\tgenewise\texon\t" . $ce->start . "\t" . $ce->end . "\t100\t" . $ce->phase . "\n";
	    }
	  }
	  else {
	    print STDERR $exon->contig_id . "\tgenewise\texon\t" . $exon->start . "\t" . $exon->end . "\t100\t" . $exon->phase . "\n";
	  }
	  foreach my $sf($exon->each_Supporting_Feature) {
	    print STDERR "sub_align: " . 
	      $sf->seqname . "\t" .
	      $sf->start . "\t" .
	      $sf->end . "\t" .
	      $sf->strand . "\t" .
	      $sf->hseqname . "\t" .
	      $sf->hstart . "\t" .
	      $sf->hend . "\n";
	  }
	}
      }
      push(@newf,$newgene);
      
    };
    if ($@) {
      
      
      print STDERR "contig: $contig\n";
      foreach my $tran ($gene->each_Transcript) {
	foreach my $exon($tran->each_Exon) {
	  foreach my $sf($exon->each_Supporting_Feature) {
	    print STDERR "hid: " . $sf->hseqname . "\n";
	  }
	}
      }
      
      
      print STDERR "Couldn't reverse map gene " . $gene->id . " [$@]\n";
    }
    
    
  }
  
  return @newf;
  
}


sub output {
    my ($self) = @_;
   
    if (!defined($self->{'_output'})) {
      $self->{'_output'} = [];
    } 
    return @{$self->{'_output'}};
}


=head2 vc

 Title   : vc
 Usage   : $obj->vc($newval)
 Function: 
 Returns : value of vc
 Args    : newvalue (optional)

=cut


1;
