#
#
# Cared for by ensembl  <ensembl-dev@ebi.ac.uk>
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::Finished_BlastMiniGenewise

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::RunnableDB::Finished_BlastMiniGenewise->new(
										  -dbobj     => $db,
										  -input_id  => $id, 
										  -golden_path => $gp,
										 );
    $obj->fetch_input;
    $obj->run;

    my @newfeatures = $obj->output;


=head1 DESCRIPTION

=head1 CONTACT

ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::RunnableDB::Finished_BlastMiniGenewise;

use vars qw(@ISA);
use strict;
require "Bio/EnsEMBL/Pipeline/GB_conf.pl";
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise;
use Bio::EnsEMBL::Pipeline::Runnable::MiniGenewise;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Getseqs;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;

use Data::Dumper;
# config file; parameters searched for here if not passed in as @args


@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB );

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);    
     
    if(!defined $self->seqfetcher) {
      my $seqfetcher =  $self->make_seqfetcher();
      $self->seqfetcher($seqfetcher);
    }
   
    my ($path, $type, $threshold) = $self->_rearrange([qw(GOLDEN_PATH TYPE THRESHOLD)], @args);
    $path = 'UCSC' unless (defined $path && $path ne '');
    $self->dbobj->assembly_type($path);
    
    

    $type = 'swall' unless (defined $type && $type ne '');
    $threshold = 200 unless (defined($threshold));

   
    $self->type($type);
    $self->threshold($threshold);


    return $self; 
}

sub type {
  my ($self,$type) = @_;

  if (defined($type)) {
    $self->{_type} = $type;
  }
  return $self->{_type};
}

sub threshold {
  my ($self,$threshold) = @_;

  if (defined($threshold)) {
    $self->{_threshold} = $threshold;
  }
  return $self->{_threshold};
}




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

    my $gene_adaptor = $self->dbobj->get_GeneAdaptor;
    #print "there are ".scalar($self->output). "genes\n";
  GENE: foreach my $gene ($self->output) {      
      # do a per gene eval..
      my @transcripts = $gene->each_Transcript;
      if(!$transcripts[0]){
	print STDERR "Gene doesn't have any transcript\n";
	next GENE;
      }else{ 
	eval {
	  $gene_adaptor->store($gene);
	  #print STDERR "wrote gene " . $gene->dbID . "\n";
	}; 
	if( $@ ) {
	  print STDERR "UNABLE TO WRITE GENE\n\n$@\n\nSkipping this gene\n";
	}
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
    
   
    $self->throw("No input id") unless defined($self->input_id);

    my $contig    = $self->dbobj->get_Contig($self->input_id);

    my $genseq    = $contig->get_repeatmasked_seq;
    
    

    my @features;
    
   
    @features  = $contig->get_all_SimilarityFeatures_above_score($self->type, $self->threshold,0);
    
    
    
    my %idhash;
    
    foreach my $f (@features) {
      if ($f->isa("Bio::EnsEMBL::FeaturePair") && 
	  defined($f->hseqname)) {
	if($f->p_value <= 0.00001){	
	  #print "seqname = ".$f->hseqname." evalue ".$f->p_value."\n";
	  if(!$idhash{$f->hseqname}){
	    $idhash{$f->hseqname} = 1;
	  }
	}
      }
    }
    my @ids = keys %idhash;
    
   
    my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise('-genomic'    => $genseq,
									   '-ids'        => \@ids,
									   '-seqfetcher' => $self->seqfetcher,
									   '-trim'       => 1);
    
    
    $self->runnable($runnable);
  

   
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

    #Now there is more than one...
    foreach my $runnable ($self->runnable) {
                if ($runnable->isa("Bio::EnsEMBL::Pipeline::Runnable::MiniGenewise")){
                        $runnable->minirun;
                }else{
                $runnable->run;
                }
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
  
  my $trancount = 1;
  my $genetype;
  foreach my $runnable ($self->runnable) {
     if ($runnable->isa("Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise")){
      $genetype = "ContigBlastMiniGenewise";
    }else{
      $self->throw("I don't know what to do with $runnable");
    }
     #$genetype = "testing";
    my $anaAdaptor = $self->dbobj->get_AnalysisAdaptor;

        #use logic name from analysis object if possible, else take $genetype;
        my $anal_logic_name = ($self->analysis->logic_name)     ?       $self->analysis->logic_name : $genetype ;       
        
    my @analyses = $anaAdaptor->fetch_by_logic_name($anal_logic_name);
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
         -db_version      => 'NULL',
         -program         => 'genewise',
         -program_version => 1,
         -gff_source      => 'genewise',
         -gff_feature     => 'gene',
         -logic_name      => $genetype,
         -module          => 'BlastMiniGenewise',
      );
    }

    my @results = $runnable->output;
    my @genes = $self->make_genes($genetype, $analysis_obj, \@results);

    $self->output(@genes);

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
  my $contig = $self->vc;
  my @tmpf   = @$results;
  my @genes;
  my $gene;
  my $transcript;
  foreach my $tmpf (@tmpf) {
    $gene = new Bio::EnsEMBL::Gene;
    $transcript = $self->_make_transcript($tmpf, $contig, $genetype, $analysis_obj);
    #my $translation = $transcript->translate;
    #if($translation->seq =~ /\*/){
    #  next;
    #}else{
    $gene->type($genetype);
    $gene->analysis($analysis_obj);
    $gene->add_Transcript($transcript);
    push (@genes, $gene);
    #}
  }
  return @genes;
}

=head2 _make_transcript

 Title   : make_transcript
 Usage   : $self->make_transcript($gene, $contig, $genetype)
 Function: makes a Bio::EnsEMBL::Transcript from a SeqFeature representing a gene, 
           with sub_SeqFeatures representing exons.
 Example :
 Returns : Bio::EnsEMBL::Transcript with Bio::EnsEMBL:Exons(with supporting feature 
           data), and a Bio::EnsEMBL::translation
 Args    : $gene: Bio::EnsEMBL::SeqFeatureI, $contig: Bio::EnsEMBL::RawContig,
  $genetype: string, $analysis_obj: Bio::EnsEMBL::Analysis


=cut

sub _make_transcript{
  my ($self, $gene, $contig, $genetype, $analysis_obj) = @_;
  $genetype = 'unspecified' unless defined ($genetype);

  unless ($gene->isa ("Bio::EnsEMBL::SeqFeatureI"))
    {print "$gene must be Bio::EnsEMBL::SeqFeatureI\n";}

  my $transcript   = new Bio::EnsEMBL::Transcript;
  my $translation  = new Bio::EnsEMBL::Translation;    
  $transcript->translation($translation);

  my $excount = 1;
  my @exons;
    
  foreach my $exon_pred ($gene->sub_SeqFeature) {
    # make an exon
    my $exon = new Bio::EnsEMBL::Exon;
    
    $exon->contig_id($contig->internal_id);
    $exon->start($exon_pred->start);
    $exon->end  ($exon_pred->end);
    $exon->strand($exon_pred->strand);
    
    $exon->phase($exon_pred->phase);
    $exon->attach_seq($contig->primary_seq);
    #print "contig mol type = ".$contig->primary_seq->moltype."\n";
    #print "exon mol type = ".$exon->entire_seq->moltype."\n";
   
    # sort out supporting evidence for this exon prediction
    foreach my $subf($exon_pred->sub_SeqFeature){
      $subf->feature1->seqname($contig->internal_id);
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
    
    #print STDERR "num exons: " . scalar(@exons) . "\n";

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

sub make_seqfetcher {
  my ($self) = @_;
  my $index = $self->analysis->db;
  my $seqfetcher;

  if(defined $index && $index ne ''){
    my @db = ( $index );
    $seqfetcher = new Bio::EnsEMBL::Pipeline::SeqFetcher::Getseqs(
                                                                  '-db' => \@db,
                                                                 );
  }
  else{
    # default to Pfetch
    $seqfetcher = new Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;
  }

  $self->seqfetcher($seqfetcher);

  return $seqfetcher;
}



=head2 _select_features
 
  Title   : _select_features
  Usage   : $self->_select_features(@features)
  Function: obtain the best scoring HSP within a certain area
  Returns : Array of FeaturePairs
  Args    : Array of selected hseqnames
 
=cut

sub _select_features {

        my ($self,@features) = @_;

        @features= sort {
                $a->strand<=> $b->strand
                       ||
                $a->start<=> $b->start
        } @features;

        my %selected_ids;

        my $best_hit = $features[0];
 
        foreach my $feat (@features){
                if ($feat->overlaps($best_hit,'strong')) {
                        if ($feat->score > $best_hit->score) {
                                $best_hit = $feat;
                        }
                        }else {
                                $selected_ids{$best_hit->hseqname} = 1;
                                $best_hit = $feat;
                        }
        }
 
        return %selected_ids;
}

1;
