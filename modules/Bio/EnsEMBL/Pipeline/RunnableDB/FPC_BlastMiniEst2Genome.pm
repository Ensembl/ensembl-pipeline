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

Bio::EnsEMBL::Pipeline::RunnableDB::FPC_BlastMiniEst2Genome

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::RunnableDB::FPC_BlastMiniEst2Genome->new(
					     -dbobj     => $db,
					     -input_id  => $id,
					     -blastdb   => $blastdb
                                             );
    $obj->fetch_input
    $obj->run

    my @genes = $obj->output;


=head1 DESCRIPTION

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::RunnableDB::FPC_BlastMiniEst2Genome;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::RootI;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::BlastMiniEst2Genome;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;

@ISA = qw( Bio::EnsEMBL::Pipeline::RunnableDB );


=head2 new

    Title   :   new
    Usage   :   $self->new(-DBOBJ       => $db
                           -INPUT_ID    => $id
                           -ANALYSIS    => $analysis);
                           
    Function:   creates a 
                Bio::EnsEMBL::Pipeline::RunnableDB::FPC_BlastMiniEst2Genome 
                object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::FPC_BlastMiniEst2Genome 
                object
    Args    :   -dbobj:      A Bio::EnsEMBL::DB::Obj (required), 
                -input_id:   Contig input id (required), 
                -seqfetcher: A Sequence Fetcher Object (required),
                -analysis:   A Bio::EnsEMBL::Pipeline::Analysis (optional) 
=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
           
    # dbobj, input_id, seqfetcher, and analysis objects are all set in
    # in superclass constructor (RunnableDB.pm)

    my( $blastdb ) = $self->_rearrange([qw(BLASTDB)], @args);
       
    $self->throw("No blast db specified") unless defined($blastdb);
    $self->blastdb($blastdb);

    if(!defined $self->seqfetcher) {
      # will look for pfetch in $PATH - change this once PipeConf up to date
      my $seqfetcher = new Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch; 
      $self->seqfetcher($seqfetcher);
    }


    return $self;
}

=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output
    Function:   Writes output data (genes) to db
    Returns :   
    Args    :   none

=cut

sub write_output {

  my($self) = @_;
  
  my @genes    = $self->output();
  my $db       = $self->dbobj();


#  my $dbout    = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
#						     -user   => 'ensadmin',
#						     -dbname => 'val_e2g_test',
#						     -host   => 'ecs1a',
#						    );
#  my $gene_obj = $dbout->gene_Obj;
  my $gene_obj = $db->gene_Obj;

  eval {
    my $genecount   = 0;
    my $transcount  = 0;
    my $translcount = 0;
    my $exoncount   = 0;
    
    my $EXON_ID_SUBSCRIPT       = "BEGE";
    my $TRANSCRIPT_ID_SUBSCRIPT = "BEGT";
    my $GENE_ID_SUBSCRIPT       = "BEGG";
    my $PROTEIN_ID_SUBSCRIPT    = "BEGP";
    
    # get counts of each type of ID we need.
    
    foreach my $gene ( @genes ) {
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
    
    my @geneids   =  $gene_obj->get_New_external_id('gene',$GENE_ID_SUBSCRIPT,$genecount);
    my @transids  =  $gene_obj->get_New_external_id('transcript',$TRANSCRIPT_ID_SUBSCRIPT,$transcount);
    my @translids =  $gene_obj->get_New_external_id('translation',$PROTEIN_ID_SUBSCRIPT,$translcount);
    my @exonsid   =  $gene_obj->get_New_external_id('exon',$EXON_ID_SUBSCRIPT,$exoncount);
    
    # database locks are over.
    
    # now assign ids. gene and transcripts are easy. Exons are harder.
    # the code currently assummes that there is one Exon object per unique
    # exon id. This might not always be the case.
    
    foreach my $gene ( @genes ) {
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
  #        $self->throw("Bailing before real write\n");
  
 GENE: foreach my $gene (@genes) {	
    # do a per gene eval...
    eval {
      
      $gene_obj->write($gene);
    }; 
    if( $@ ) {
      print STDERR "UNABLE TO WRITE GENE\n\n$@\n\nSkipping this gene\n";
    }
    
  }
  
  return 1;
}

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data for BlastMiniEst2Genome and makes runnable
    Returns :   nothing
    Args    :   none

=cut

sub fetch_input {
  my ($self) = @_;
  
  print STDERR "Fetching input \n";
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
  
  my $genseq    = $contig->get_repeatmasked_seq;
  my $blastdb   = $self->blastdb;
  print STDERR "fpc blastdb: $blastdb\n";
  
  my $runnable  = new Bio::EnsEMBL::Pipeline::Runnable::BlastMiniEst2Genome('-genomic'     => $genseq, 
									    '-queryseq'     => $blastdb,
									    '-seqfetcher'  => $self->seqfetcher);
    
  $self->runnable($runnable);
  # at present, we'll only ever have one ...
  $self->vc($contig);
}

=head2 run

    Title   :   run
    Usage   :   $self->run
    Function:   Calls run method of each runnable, & converts output into remapped genes
    Returns :   Nothing
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

=head2 _convert_output

    Title   :   _convert_output
    Usage   :   $self->_convert_output()
    Function:   Converts est2genome output into an array of genes remapped into genomic coordinates
    Returns :   Nothing, but $self->{_output} contains remapped genes
    Args    :   None
=cut

# get merged features into a form where they can be stored in the database.
sub _convert_output {
  my ($self) = @_;
  my $count  = 1;
  my $time   = time; chomp($time);
  my @genes;

  # make an array of genes for each runnable
  foreach my $runnable ($self->runnable) {
    my @g = $self->_make_genes($count, $time, $runnable);
    push(@genes, @g);
  }

  # filter out genes which are poor matches to ESTs
#  my $threshold = 95;
#  my @filtered = $self->_filter_genes($threshold, @genes);

  # map genes back to genomic coordinates
#  my @remapped = $self->_remap_genes(@filtered);	
  my @remapped = $self->_remap_genes(@genes);	
    
  if (!defined($self->{_output})) {
    $self->{_output} = [];
  }
  
  push(@{$self->{_output}},@remapped);
}


=head2 _make_genes

    Title   :   _make_genes
    Usage   :   $self->_make_gene($count, $time, $runnable)
    Function:   Converts merged exon features into a gene
    Returns :   Bio::EnsEMBL::Gene
    Args    :   

=cut

sub _make_genes {
  my ($self, $count, $time, $runnable) = @_;
  my $contig = $self->vc;
  my $genetype = 'bmeg';
  
  my @tmpf = $runnable->output; # an array of SeqFeaturesm one per gene prediction, with subseqfeatures
  print STDERR "we'll have " . scalar(@tmpf) . " genes\n";
  my @genes;
  
  foreach my $tmpf(@tmpf) {
    my $gene   = new Bio::EnsEMBL::Gene;
    my $tran   = new Bio::EnsEMBL::Transcript;
    my $transl = new Bio::EnsEMBL::Translation;
    
    $gene->type($genetype);
    $gene->id($self->input_id . ".$genetype.$count");
    $gene->created($time);
    $gene->modified($time);
    $gene->version(1);
    
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
      $exon->seqname($contig->id);
      $exon->start($exon_pred->start);
      $exon->end  ($exon_pred->end);
      $exon->strand($exon_pred->strand);
      
      #	$exon->phase($subf->feature1->{_phase});
      
      $exon->phase($exon_pred->phase);
      $exon->attach_seq($self->vc->primary_seq);
      # fix source tag and primary tag for $exon_pred - this isn;t the right place to do this.
      $exon_pred->source_tag('BME2G');
      $exon_pred->primary_tag('BME2G');

      $exon_pred->score(100); # ooooooohhhhhh
      
      print "num subf: " . scalar($exon_pred->sub_SeqFeature) . "\n";

      # sort out supporting evidence for this exon prediction
      foreach my $subf($exon_pred->sub_SeqFeature){
	$subf->feature1->source_tag($genetype);
	$subf->feature1->primary_tag('similarity');
	$subf->feature1->analysis($exon_pred->analysis);
	
	$subf->feature2->source_tag($genetype);
	$subf->feature2->primary_tag('similarity');
	$subf->feature2->analysis($exon_pred->analysis);
	
	$exon->add_Supporting_Feature($subf);
      }
      
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
      
      my $endexon = $exons[$#exons];
      
      if( $endexon->end_phase == 1 ) {
	$transl->end($endexon->length -1 );
      } elsif ( $endexon->end_phase == 2 ) {
	$transl->end($endexon->length -2 );
      } else {
	$transl->end($endexon->length);
      }
      #$transl->end  ($exons[$#exons]->end - $exons[$#exons]->start + 1);

    }
  }

  # no point in trying ttranslate these - they just won't.

  return @genes;

}

=head2 _filter_genes

    Title   :   _filter_genes
    Usage   :   $self->_filter_genes($threshold, @genes)
    Function:   Filters genes on basis of the % coverage of each exon to the EST it was predicted from.
                Only keep genes where every exon has an overall score that exceeds $threshold
    Returns :   array of Bio::EnsEMBL::Gene
    Args    :   int, array of Bio::EnsEMBL::Gene

=cut

sub _filter_genes {
  my ($self, $threshold, @genes) = @_;
  my @filtered;

  print STDERR "about to filter " . scalar(@genes) . " genes\n";

  GENE:
  foreach my $gene(@genes) {
    my @transcripts = $gene->each_Transcript;
    if(scalar(@transcripts) != 1) {
      print STDERR "eeek! got " . scalar(@transcripts) . " for " . $gene->id . " - letting it through\n";
      push (@filtered, $gene);
      next GENE;
    }

    print STDERR "processing " . $gene->id . "\n";

    foreach my $exon($transcripts[0]->each_Exon){
      my $matching_bases = 0;
      foreach my $sf($exon->each_Supporting_Feature) {
	print STDERR $sf->hseqname . " " . $sf->start . " " . $sf->end . " " . $sf->score . "\n";
	my $length = abs($sf->end - $sf->start) +1; # start-end inclusive!
	my $bases = ($sf->score * $length)/100;
	$matching_bases += $bases;
      }

      my $exon_length = abs($exon->end - $exon->start) + 1; # start-end inclusive!
      
      # sanity checks
      if ($matching_bases > $exon_length) {
	$self->warn("num matches > length! makes no sense!\n");
      }

      my $coverage = ($matching_bases / $exon_length) * 100;
      print STDERR "matching: $matching_bases / elength: $exon_length  = coverage: $coverage threshold: $threshold\n";

      next GENE unless $coverage >= $threshold;
  #    next GENE unless $coverage ge $threshold;
    }

    push(@filtered, $gene);

  }
 
   print STDERR " we are left with " . scalar (@filtered) . " genes\n";

  return (@filtered);

}

=head2 _remap_genes

    Title   :   _remap_genes
    Usage   :   $self->_remap_genes(@genes)
    Function:   Remaps predicted genes into genomic coordinates
    Returns :   array of Bio::EnsEMBL::Gene
    Args    :   Bio::EnsEMBL::Virtual::Contig, array of Bio::EnsEMBL::Gene

=cut

sub _remap_genes {
  my ($self, @genes) = @_;
  my $contig = $self->vc;
  my @remapped;

  foreach my $gene(@genes) {
    eval {
      my $newgene = $contig->convert_Gene_to_raw_contig($gene);
      $newgene->type('fpc_bme2g');
      push(@remapped,$newgene);
    };

    if ($@) {
      print STDERR "Couldn't reverse map gene " . $gene->id . " [$@]\n";
    }
  }

  return @remapped;
}

=head2 _print_FeaturePair

    Title   :   print_FeaturePair
    Usage   :   $self->_print_FeaturePair($pair)
    Function:   Prints attributes of a Bio::EnsEMBL::FeaturePair
    Returns :   Nothing
    Args    :   A Bio::EnsEMBL::FeaturePair

=cut

sub _print_FeaturePair {
  my ($self,$pair) = @_;
  
  print STDERR $pair->seqname . "\t" . $pair->start . "\t" . $pair->end . "\t" . 
               $pair->score . "\t" . $pair->strand . "\t" . $pair->hseqname . "\t" . 
	       $pair->hstart . "\t" . $pair->hend . "\t" . $pair->hstrand . "\n";
}

=head2 output

    Title   :   output
    Usage   :   $self->output
    Function:   Returns output from this RunnableDB
    Returns :   Array of Bio::EnsEMBL::Gene
    Args    :   None

=head2 vc

 Title   : vc
 Usage   : $obj->vc($newval)
 Function: 
 Returns : value of vc
 Args    : newvalue (optional)

=head2 blastdb

 Title   : blastdb
 Usage   : $obj->blastdb($newval)
 Function: 
 Returns : value of blastdb
 Args    : newvalue (optional)


=cut

sub blastdb{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'blastdb'} = $value;
    }
    return $obj->{'blastdb'};

}

1;


