#
#
# Cared for by Michele Clamp  <ensembl-dev@ebi.ac.uk>
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
					     -dbobj     => $db,
					     -input_id  => $id,
					     -golden_path => $gp
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

use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise;
use Bio::EnsEMBL::Pipeline::GeneConf qw (EXON_ID_SUBSCRIPT
					 TRANSCRIPT_ID_SUBSCRIPT
					 GENE_ID_SUBSCRIPT
					 PROTEIN_ID_SUBSCRIPT
					 );
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;

use Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;
use Data::Dumper;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB );

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);    
        
    if(!defined $self->seqfetcher) {
      my $seqfetcher = new Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;
      $self->seqfetcher($seqfetcher);
    }
       
    my ($path) = $self->_rearrange([qw(GOLDEN_PATH)], @args);
    $path = 'UCSC' unless (defined $path && $path ne '');
    $self->dbobj->static_golden_path_type($path);
    return $self; 
}

=head1 RunnableDB implemented methods

=head2 dbobj

    Title   :   dbobj
    Usage   :   $self->dbobj($obj);
    Function:   Gets or sets the value of dbobj
    Returns :   A Bio::EnsEMBL::Pipeline::DB::ObjI compliant object
                (which extends Bio::EnsEMBL::DB::ObjI)
    Args    :   A Bio::EnsEMBL::Pipeline::DB::ObjI compliant object

=cut

=head2 input_id

    Title   :   input_id
    Usage   :   $self->input_id($input_id);
    Function:   Gets or sets the value of input_id
    Returns :   valid input id for this analysis (if set) 
    Args    :   input id for this analysis 

=head2 vc

 Title   : vc
 Usage   : $obj->vc($newval)
 Function: 
 Returns : value of vc
 Args    : newvalue (optional)

=head1 FPC_BlastMiniGenewise implemented methods

=head2 fetch_output

    Title   :   fetch_output
    Usage   :   $self->fetch_output($file_name);
    Function:   Fetchs output data from a frozen perl object
                stored in file $file_name
    Returns :   array of exons (with start and end)
    Args    :   none

=cut

sub fetch_output {
    my($self,$output) = @_;
    
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
	print STDERR $gene->id . "\n";
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
    Args    :   none

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

#    $self->dbobj->static_golden_path_type('UCSC');

    my $stadaptor = $self->dbobj->get_StaticGoldenPathAdaptor();
    my $contig    = $stadaptor->fetch_VirtualContig_by_chr_start_end($chrid,$chrstart,$chrend);

    $contig->_chr_name($chrid);

    foreach my $rc ($contig->_vmap->each_MapContig) {
	my $strand = "+";
	if ($rc->orientation == -1) {
	    $strand = "-";
	}
	
    }

    my $genseq    = $contig->get_repeatmasked_seq;

    print STDERR "Length is " . $genseq->length . "\n";
    print STDERR "Fetching features \n";

    print STDERR "contig: " . $contig . " \n";

    # need to pass in bp value of zero to prevent globbing on StaticContig.
    my @features  = $contig->get_all_SimilarityFeatures_above_score('sptr',200, 0);
    
    # lose version numbers - probably temporary till pfetch indices catch up
    foreach my $f(@features) {
      my $name = $f->hseqname;
      $name =~ /(\S+)\.\d+/;
      $f->hseqname($1);
    }

    print STDERR "Number of features = " . scalar(@features) . "\n";

    my @genes     = $contig->get_Genes_by_Type('TGE_gw');

    print STDERR "Found " . scalar(@genes) . " genewise genes\n";

    my %redids;
    my $trancount = 1;

    foreach my $gene (@genes) {
      print STDERR "Found genewise gene " . $gene->id . "\n";
      foreach my $tran ($gene->each_Transcript) {
	foreach my $exon ($tran->each_Exon) {
	  print STDERR "Exon " . $exon->id . " " . $exon->strand . "\n";

	  if ($exon->seqname eq $contig->id) {
	    
	  FEAT: foreach my $f (@features) {
	      if ($exon->overlaps($f)) {
		$redids{$f->hseqname} = 1;
		print STDERR "ID " . $f->hseqname . " covered by genewise\n";
	      }
	    }
	  }
	}
	$trancount++;
      }
    }

    my %idhash;
    
    foreach my $f (@features) {
        print "Feature " . $f . " " . $f->seqname . " " . $f->source_tag . "\n";
      if ($f->isa("Bio::EnsEMBL::FeaturePair") && 
	  defined($f->hseqname) &&
	    $redids{$f->hseqname} != 1) {
      $idhash{$f->hseqname} = 1;
      
    }
  }
    
    my @ids = keys %idhash;

    print STDERR "Feature ids are @ids\n";

    my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise('-genomic'  => $genseq,
									   '-ids'      => \@ids,
									   '-seqfetcher' => $self->seqfetcher,
									   '-trim'     => 1);
    
    
    $self->runnable($runnable);
    # at present, we'll only ever have one ...
    $self->vc($contig);
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

    # is there ever going to be more than one?
    foreach my $runnable ($self->runnable) {
	$runnable->run;
    }
    
    $self->convert_output;

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
  
  my $count = 1;

  
  my $trancount = 1;
  my $genetype;
  foreach my $runnable ($self->runnable) {
    if ($runnable->isa("Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise")){
      $genetype = "similarity_genewise";
    }
    else{
      $self->throw("I don't know what to do with $runnable");
    }

    my $anaAdaptor = $self->dbobj->get_AnalysisAdaptor;
    my @analyses = $anaAdaptor->fetch_by_logic_name($genetype);
    my $analysis_obj;
    if(scalar(@analyses) > 1){
      $self->throw("panic! > 1 analysis for $genetype\n");
    }
    elsif(scalar(@analyses) == 1){
      $analysis_obj = $analyses[0];
    }
    else{
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
    my @genes = $self->make_genes($count, $genetype, $analysis_obj, \@results);

    my @remapped = $self->remap_genes(@genes);

    $self->output(@remapped);

  }
}


=head2 make_genes

  Title   :   make_genes
  Usage   :   $self->make_genes
  Function:   makes Bio::EnsEMBL::Genes out of the output from runnables
  Returns :   array of Bio::EnsEMBL::Gene  
  Args    :   $count: integer; 
              $genetype: string
              $analysis_obj: Bio::EnsEMBL::Analysis
              $results: reference to an array of FeaturePairs

=cut

sub make_genes {
  my ($self, $count, $genetype, $analysis_obj, $results) = @_;
  my $contig = $self->vc;
  my @tmpf   = @$results;
  my @genes;
  my $time  = time; chomp($time);
  foreach my $tmpf (@tmpf) {
    my $gene   = new Bio::EnsEMBL::Gene;
    my $tran   = new Bio::EnsEMBL::Transcript;
    my $transl = new Bio::EnsEMBL::Translation;
    
    $gene->type($genetype);
    $gene->id($self->input_id . ".$genetype.$count");
    $gene->created($time);
    $gene->modified($time);
    $gene->version(1);
    $gene->analysis($analysis_obj);

    $tran->id($self->input_id . ".$genetype.$count");
    $tran->created($time);
    $tran->modified($time);
    $tran->version(1);
    
    $transl->id($self->input_id . ".$genetype.$count");
    $transl->version(1);
    
    $count++;
    
    $gene->add_Transcript($tran);
    $tran->translation($transl);
    
    my $excount = 1;
    my @exons;
    
    foreach my $exon_pred ($tmpf->sub_SeqFeature) {
      # make an exon
      my $exon = new Bio::EnsEMBL::Exon;
      
      $exon->id($self->input_id . ".$genetype.$count.$excount");
      $exon->contig_id($contig->id);
      $exon->created($time);
      $exon->modified($time);
      $exon->version(1);
      
      $exon->start($exon_pred->start);
      $exon->end  ($exon_pred->end);
      $exon->strand($exon_pred->strand);
      
      $exon->phase($exon_pred->{_phase});
      $exon->attach_seq($self->vc->primary_seq);

      # sort out supporting evidence for this exon prediction
      foreach my $subf($exon_pred->sub_SeqFeature){
	$subf->feature1->source_tag($genetype);
	$subf->feature1->primary_tag('similarity');
	$subf->feature1->score(100);
	$subf->feature1->analysis($exon_pred->analysis);
	
	$subf->feature2->source_tag($genetype);
	$subf->feature2->primary_tag('similarity');
	$subf->feature2->score(100);
	$subf->feature2->analysis($exon_pred->analysis);
	
	$exon->add_Supporting_Feature($subf);
      }
      
      my $seq   = new Bio::Seq(-seq => $exon->seq->seq);
      
      push(@exons,$exon);
      
      $excount++;
    }
    
    if ($#exons < 0) {
      print STDERR "Odd.  No exons found\n";
    } else {
      
      push(@genes,$gene);
      
      if ($exons[0]->strand == -1) {
	@exons = sort {$b->start <=> $a->start} @exons;
      } else {
	@exons = sort {$a->start <=> $b->start} @exons;
      }
      
      foreach my $exon (@exons) {
	$tran->add_Exon($exon);
      }
      
      $transl->start_exon_id($exons[0]->id);
      $transl->end_exon_id  ($exons[$#exons]->id);
      
      if ($exons[0]->phase == 0) {
	$transl->start(1);
      } elsif ($exons[0]->phase == 1) {
	$transl->start(3);
      } elsif ($exons[0]->phase == 2) {
	$transl->start(2);
      }
      
      $transl->end  ($exons[$#exons]->end - $exons[$#exons]->start + 1);
    }
  }
  return @genes;

}

=head2 remap_genes

 Title   : remap_genes
 Usage   : $self->remap_genes($runnable, @genes)
 Function: converts the coordinates of each Bio@EnsEMBL::Gene in @genes into RawContig
           coordinates for storage.
 Example : 
 Returns : array of Bio::EnsEMBL::Gene in RawContig coordinates
 Args    : @genes: array of Bio::EnsEMBL::Gene in virtual contig coordinates


=cut

sub remap_genes {
  my ($self, @genes) = @_;
  my $contig = $self->vc;

  my @newf;
  my $trancount=1;
  foreach my $gene (@genes) {
    eval {
      my $newgene = $contig->convert_Gene_to_raw_contig($gene);
      # need to explicitly add back genetype and analysis.
      $newgene->type($gene->type);
      $newgene->analysis($gene->analysis);

      foreach my $tran ($newgene->each_Transcript) {
	foreach my $exon($tran->each_Exon) {
	  foreach my $sf($exon->each_Supporting_Feature) {
	    # this should be sorted out by the remapping to rawcontig ... strand is fine
	    if ($sf->start > $sf->end) {
	      my $tmp = $sf->start;
	      $sf->start($sf->end);
	      $sf->end($tmp);
	    }
	  }
	}
      }
      push(@newf,$newgene);

    };
    if ($@) {
      print STDERR "Couldn't reverse map gene " . $gene->id . " [$@]\n";
    }
    

  }

  return @newf;
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
    
   if(defined @genes){
     push(@{$self->{'_output'}},@genes);
   }

   return @{$self->{'_output'}};
}

1;


