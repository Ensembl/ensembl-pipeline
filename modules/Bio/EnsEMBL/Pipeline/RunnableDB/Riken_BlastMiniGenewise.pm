#
#
# Cared for by EnsEMBL  <ensembl-dev@ebi.ac.uk>
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::Riken_BlastMiniGenewise

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::RunnableDB::Riken_BlastMiniGenewise->new
    (
     '-dbobj'       => $db,
     '-input_id'    => $id,
     '-seqfetcher'  => $seqfetcher,
     '-golden_path' => $gp
     );
    $obj->fetch_input
    $obj->run

    my @newfeatures = $obj->output;


=head1 DESCRIPTION

Same basis as FPC_BlastMiniGenewise, but does not check against
previous predictions - we want to run all the Riken sequences.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::RunnableDB::Riken_BlastMiniGenewise;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::RootI;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Pipeline::SeqFetcher::getseqs;
use Data::Dumper;

# config file; parameters searched for here if not passed in as @args
require "Bio/EnsEMBL/Pipeline/GB_conf.pl";

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

=head2 new

    Title   :   new
    Usage   :   $self->new(-DBOBJ       => $db
                           -INPUT_ID    => $id
			   -SEQFETCHER  => $sf);
                           
    Function:   creates a Bio::EnsEMBL::Pipeline::RunnableDB::Riken_BlastMiniGenewise object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::Riken_BlastMiniGenewise 
                object
    Args    :   -dbobj:      A Bio::EnsEMBL::DB::Obj (required), 
                -input_id:   Contig input id (required), 
                -seqfetcher: A Bio::DB::RandomAccessI Object (required)
=cut

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);    
    
    # dbobj, input_id, seqfetcher objects are all set in
    # in superclass constructor (RunnableDB.pm)

    # check seqfetcher
    if(!defined $self->seqfetcher) {

      my $seqfetcher = $self->make_seqfetcher();
      $self->seqfetcher($seqfetcher);
    } 

print STDERR "seqfethcre is " . $self->seqfetcher . "\n";

    # check golden_path
    my ($path) = $self->_rearrange([qw(GOLDEN_PATH)], @args);
    $path = 'UCSC' unless (defined $path && $path ne '');
    $self->dbobj->static_golden_path_type($path);

    return $self;
}

=head2 dbobj

    Title   :   dbobj
    Usage   :   $self->dbobj($db)
    Function:   Get/set method for database handle
    Returns :   Bio::EnsEMBL::Pipeline::DB::ObjI
    Args    :   

=head2 input_id

    Title   :   input_id
    Usage   :   $self->input_id($input_id);
    Function:   Gets or sets the value of input_id
    Returns :   valid input id for this analysis (if set) 
    Args    :   input id for this analysis 

=head2 seqfetcher

    Title   :   seqfetcher
    Usage   :   $self->seqfetcher($seqfetcher)
    Function:   Get/set method for SeqFetcher
    Returns :   Bio::DB::RandomAccessI object
    Args    :   Bio::DB::RandomAccessI object

=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output
    Function:   Writes output data to db
    Returns :   
    Args    :   none

=cut

sub write_output {
    my($self) = @_;

    #$self->throw("exiting bfore write");

    my $gene_adaptor = $self->dbobj->get_GeneAdaptor;

  GENE: foreach my $gene ($self->output) {	
      # do a per gene eval...
      eval {
	  $gene_adaptor->store($gene);
	  print STDERR "wrote gene " . $gene->dbID . "\n";
      }; 
      if( $@ ) {
	  print STDERR "UNABLE TO WRITE GENE\n\n$@\n\nSkipping this gene\n";
      }
	    
  }
   
}

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Creates virtual contig and fetches Riken hits from the database for input to
                BlastMiniGenewise
    Returns :   nothing
    Args    :   none

=cut

sub fetch_input {
    my( $self) = @_;
    
    print STDERR "Fetching input \n";
    $self->throw("No input id") unless defined($self->input_id);

    my $chrid  = $self->input_id;
       $chrid =~ s/\.(.*)-(.*)//;

    my $chrstart = $1;
    my $chrend   = $2;

    my $stadaptor = $self->dbobj->get_StaticGoldenPathAdaptor();
    my $contig    = $stadaptor->fetch_VirtualContig_by_chr_start_end($chrid,$chrstart,$chrend);

    $contig->_chr_name($chrid);

    my $genseq    = $contig->get_repeatmasked_seq;

    print STDERR "Fetching features \n";

    my @features  = $contig->get_all_SimilarityFeatures_above_score('riken_prot',200, 0);

    
    print STDERR "Number of features = " . scalar(@features) . "\n";

    my %idhash;
    
    foreach my $f (@features) {
      if ($f->isa("Bio::EnsEMBL::FeaturePair") && 
	  defined($f->hseqname) ) {
      $idhash{$f->hseqname} = 1;
      
    }
  }
    
    my @ids = keys %idhash;

    print STDERR "Feature ids are @ids\n";

    my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise
	('-genomic'    => $genseq,
	 '-ids'        => \@ids,
	 '-seqfetcher' => $self->seqfetcher,
	 '-trim'       => 1);    
    
    $self->runnable($runnable);
    $self->{$runnable} = $contig;

}     


=head2 run

 Title   : run
 Usage   : $self->run
 Function: calls run for each stored runnable, and converts output to Bio::EnsEMBL::Gene 
           in RawContig coordinate that are stored in $slef->{'_output'}
 Example :
 Returns : 
 Args    : 


=cut

sub run {
    my ($self) = @_;

    foreach my $runnable ($self->runnable) {
	$runnable->run;
    }
    
    $self->convert_output;
}

=head2 convert_output

 Title   : convert_output
 Usage   : $self->convert_output
 Function: Converts the output of each runnable into Bio::EnsEMBL::Genes in RawContig coordinates,
           and stores these in $self->{'_output'}. Each Gene has a single Transcript, which in turn has
           Exons(with supproting features) and a Translation
 Example :
 Returns : 
 Args    :


=cut


sub convert_output {
  my ($self) =@_;

  my $genetype = 'riken_genewise';

  # get the appropriate analysis from the AnalysisAdaptor
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
       -module          => 'Riken_BlastMiniGenewise',
      );
  }

  foreach my $runnable ($self->runnable) {
    my @genes = $self->make_genes($genetype, $analysis_obj, $runnable);

    my @remapped = $self->remap_genes($runnable, @genes);

    $self->output(@remapped);

    # store the genes

    
  }
}

=head2 output

 Title   : output
 Usage   :
 Function:
 Example :
 Returns : 
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

=head2 make_genes

 Title   : make_genes
 Usage   : $self->make_genes($genetype, $analysis_obj, $runnable)
 Function: converts the output from $runnable into Bio::EnsEMBL::Genes in
           $contig(VirtualContig) coordinates. The genes have type $genetype, 
           and have $analysis_obj attached. Each Gene has a single Transcript, 
           which in turn has Exons(with supporting features) and a Translation
 Example : 
 Returns : array of Bio::EnsEMBL::Gene
 Args    : $genetype: string, $analysis_obj: Bio::EnsEMBL::Analysis, 
           $runnable: Bio::EnsEMBL::Pipeline::RunnableI


=cut

sub make_genes {

  my ($self, $genetype, $analysis_obj, $runnable) = @_;
  my $contig = $self->{$runnable};
  my @tmpf   = $runnable->output;
  my @genes;

  
  foreach my $tmpf (@tmpf) {
    my $gene       = new Bio::EnsEMBL::Gene;
    my $transcript = $self->_make_transcript($tmpf, $contig, $genetype, $analysis_obj);

    $gene->type($genetype);
    $gene->analysis($analysis_obj);
    $gene->add_Transcript($transcript);

    push (@genes, $gene);
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
 Args    : $runnable: Bio::EnsEMBL::Pipeline::RunnableI, @genes: array of Bio::EnsEMBL::Gene 
           in virtual contig coordinates


=cut

sub remap_genes {
  my ($self, $runnable, @genes) = @_;
  
  my $contig = $self->{$runnable};

  my @newf;
  my $trancount=1;
  foreach my $gene (@genes) {
    eval {
      my $newgene = $contig->convert_Gene_to_raw_contig($gene);

      # need to explicitly add back genetype and analysis.
      $newgene->type($gene->type);
      $newgene->analysis($gene->analysis);

      foreach my $tran ($newgene->each_Transcript) {
	foreach my $exon($tran->get_all_Exons) {
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

=head2 _make_transcript

 Title   : make_transcript
 Usage   : $self->make_transcript($gene, $contig, $genetype)
 Function: makes a Bio::EnsEMBL::Transcript from a SeqFeature representing a gene, 
           with sub_SeqFeatures representing exons.
 Example :
 Returns : Bio::EnsEMBL::Transcript with Bio::EnsEMBL:Exons(with supporting feature 
           data), and a Bio::EnsEMBL::translation
 Args    : $gene: Bio::EnsEMBL::SeqFeatureI, $contig: Bio::EnsEMBL::DB::ContigI,
  $genetype: string, $analysis_obj: Bio::EnsEMBL::Analysis


=cut

sub _make_transcript{
  my ($self, $gene, $contig, $genetype, $analysis_obj) = @_;
  $genetype = 'unspecified' unless defined ($genetype);

  unless ($gene->isa ("Bio::EnsEMBL::SeqFeatureI"))
    {print "$gene must be Bio::EnsEMBL::SeqFeatureI\n";}
  unless ($contig->isa ("Bio::EnsEMBL::DB::ContigI"))
    {print "$contig must be Bio::EnsEMBL::DB::ContigI\n";}

  my $transcript   = new Bio::EnsEMBL::Transcript;
  my $translation  = new Bio::EnsEMBL::Translation;    
  $transcript->translation($translation);

  my $excount = 1;
  my @exons;
    
  foreach my $exon_pred ($gene->sub_SeqFeature) {
    # make an exon
    my $exon = new Bio::EnsEMBL::Exon;
    
    $exon->contig_id($contig->id);
    $exon->start($exon_pred->start);
    $exon->end  ($exon_pred->end);
    $exon->strand($exon_pred->strand);
    
    $exon->phase($exon_pred->phase);
    $exon->attach_seq($contig);
    
    # sort out supporting evidence for this exon prediction
    foreach my $subf($exon_pred->sub_SeqFeature){
      $subf->feature1->source_tag($genetype);
      $subf->feature1->primary_tag('similarity');
      $subf->feature1->score(100);
      $subf->feature1->analysis($analysis_obj);
	
      $subf->feature2->source_tag($genetype);
      $subf->feature2->primary_tag('similarity');
      $subf->feature2->score(100);
      $subf->feature2->analysis($analysis_obj);
      
      $exon->add_Supporting_Feature($subf);
    }
    
    push(@exons,$exon);
    
    $excount++;
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
    
    $translation->start_exon($exons[0]);
    $translation->end_exon  ($exons[$#exons]);
    
    if ($exons[0]->phase == 0) {
      $translation->start(1);
    } elsif ($exons[0]->phase == 1) {
      $translation->start(3);
    } elsif ($exons[0]->phase == 2) {
      $translation->start(2);
    }
    
    $translation->end  ($exons[$#exons]->end - $exons[$#exons]->start + 1);
  }
  
  return $transcript;
}

=head2 make_seqfetcher

 Title   : make_seqfetcher
 Usage   :
 Function: get/set
 Example :
 Returns : Bio::DB::RandomAccessI
 Args    :


=cut

sub make_seqfetcher{
  my ( $self ) = @_;
  my $index = $::riken_conf{'riken_index'};
  my $seqfetcher;

  if(defined $index && $index ne ''){
    my @db = ( $index );
    $seqfetcher = new Bio::EnsEMBL::Pipeline::SeqFetcher::getseqs(
								  '-db' => \@db,
								 );
  }
  else{
    # default to Pfetch
    $seqfetcher = new Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;
  }

  return $seqfetcher;
}


1;
