

#
# Ensembl module for Bio::EnsEMBL::Pipeline::RunnableDB::DnaHmm
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::DnaHmm - DnaHmm

=head1 SYNOPSIS

  # standard runnabledb useage

=head1 DESCRIPTION

This runnabledb is really quite involved. It takes a RawContig,
predicts all ORFs and then removes those which overlap with BLAST hits
> 150bits and all repeats. The results peptides are then screened
against seg. The final datase is then searched against a library of
HMMs (this is the first expensive operation), using hmmpfam. The HMMs
which get greater than 25 bits in that round are then passed into
genewise against the original DNA sequence. The final gene fragment
predictions are written out into the gene table.

=head1 AUTHOR - Ewan Birney

This modules is part of the Ensembl project http://www.ensembl.org

Email birney@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Pipeline::RunnableDB::DnaHmm;
use vars qw(@ISA);
use strict;


use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::ORF;
use Bio::EnsEMBL::Pipeline::Runnable::Protein::Seg;
use Bio::EnsEMBL::Pipeline::Runnable::Protein::Hmmpfam;
use Bio::EnsEMBL::Pipeline::Runnable::GenewiseHmm;

# Object preamble - inherits from Bio::EnsEMBL::Root

use Bio::EnsEMBL::Root;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

# new() is written here 

=head2 new

    Title   :   new
    Usage   :   $self->new(-DBOBJ       => $db
                           -INPUT_ID    => $id
                           -ANALYSIS    => $analysis);
                           
    Function:   creates a Bio::EnsEMBL::Pipeline::RunnableDB::EPCR object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::EPCR object
    Args    :   -dbobj:     A Bio::EnsEMBL::DB::Obj, 
                input_id:   Contig input id , 
                -analysis:  A Bio::EnsEMBL::Analysis 

=cut

sub new {
    my ($class, @args) = @_;
    my $self = bless {}, $class;
    
    $self->{'_blast'}   = [];
    $self->{'_repeat'}  = [];
        
    my ( $dbobj, $input_id, $analysis, $threshold) = 
            $self->_rearrange (['DBOBJ', 'INPUT_ID', 'ANALYSIS'], @args);
    
    $self->throw('Need database handle') unless ($dbobj);
    $self->throw("[$dbobj] is not a Bio::EnsEMBL::DB::ObjI")  
                unless ($dbobj->isa ('Bio::EnsEMBL::DB::ObjI'));
    $self->dbobj($dbobj);
    
    $self->throw("No input id provided") unless ($input_id);
    $self->input_id($input_id);
    
    #$self->throw("Analysis object required") unless ($analysis);
    #$self->throw("Analysis object is not Bio::EnsEMBL::Analysis")
     #           unless ($analysis->isa("Bio::EnsEMBL::Analysis"));
    #$self->analysis($analysis);
    

    return $self;
}


=head2 fetch_input

 Title   : fetch_input
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub fetch_input{
   my ($self) = @_;

   my $rc = $self->dbobj->get_Contig($self->input_id);

   $self->vcontig($rc);

   my $seq = Bio::PrimarySeq->new( -id => $self->input_id, -seq => $rc->seq);
   $self->seq($seq);

   push(@{$self->{'_blast'}},$rc->get_all_SimilarityFeatures);
   push(@{$self->{'_repeat'}},$rc->get_all_RepeatFeatures);



   return;

}

=head2 run

 Title   : run
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub run{
   my ($self) = @_;

   my $orf = Bio::EnsEMBL::Pipeline::Runnable::ORF->new( -seq => $self->seq);
   $orf->run;
   my @orf = $orf->output;
   
   $orf = undef;

   # build up a string of seq length, setting 1 when there is 
   # repeat or Blast score > 150 bits

   my $str = '0' x $self->seq->length;

 

   foreach my $bl ( @{$self->{'_blast'}} ) {
       substr($str,$bl->start,$bl->length) = '1' x $bl->length;
   }

   foreach my $rep ( @{$self->{'_repeat'}} ) {
       substr($str,$rep->start,$rep->length) = '1' x $rep->length;
   }

   # now loop over ORF hits - if any '1s' in the string, bug out

   my @finalorf;

   foreach my $orf ( @orf ) {
       if( substr($str,$orf->start,$orf->length) =~ /1/ ) {
	   next;
       } else {
	   push(@finalorf,$orf);
       }
   }

   # loop over finalorf set, write out the files to a tempfile, seg it,
   # run hmmpfam, store the HMMs that hit in a hash

   my $count;
   my %hmmhash;

   foreach my $orf ( @finalorf ) {
       $count++;
       my $tempseq = Bio::PrimarySeq->new( -id => $self->input_id."_orf_".$count , -seq => $orf->{'_peptide'} );

       #print STDERR $orf->{'_peptide'}."\n";
       
       my $seganalysis = $self->dbobj->get_AnalysisAdaptor->fetch_by_newest_logic_name('Seg');

            
       my $seg = Bio::EnsEMBL::Pipeline::Runnable::Protein::Seg->new(-query => $tempseq, -analysis => $seganalysis);
     
       $seg->run;
    
       
       my @out = $seg->output;
       
       my $count;

       foreach my $o(@out) {
	   
	   my $low_length = $o->end - $o->start + 1;
	   
	   substr($orf->{'_peptide'},$o->start,$low_length) = 'X' x $low_length;
	   
	   #print STDERR $orf->{'_peptide'}."\n";
       }
   
       my $tempseq2 = Bio::PrimarySeq->new( -id => $self->input_id."seg".$count , -seq => $orf->{'_peptide'});
       
       
       my $hmmanalysis = $self->dbobj->get_AnalysisAdaptor->fetch_by_newest_logic_name("Pfam");
       
       my $temp = $tempseq2->seq;

       $temp =~ s/X//g;

       print STDERR "TEMPSEQ3: $temp\n";


       if (length($temp)>30 ) { 
  
	   my $pfam = Bio::EnsEMBL::Pipeline::Runnable::Protein::Hmmpfam->new(-query => $tempseq2, -analysis=> $hmmanalysis);
       
	   $pfam->run;
	   
	   foreach my $domain ( $pfam->output ) {
	       
	       if($domain->feature1->score > 25) {
	       
	       #Get the the genomic sequence corresponding to the orf + dowstream and upstream region (length to be defined)
	       #my $subseq_start = $orf->start - 100;
	       #if ($subseq_start < 1) {
	#	   $subseq_start = 1;
	       #}
	       
	       #my $subseq_end  = $orf->end + 100;
	       #if ($subseq_end > $self->seq->length) {
		#   $subseq_end = $self->seq->length;
	       #}
	       
	       #my $subseq = $self->seq->subseq($subseq_start,$subseq_end);
	       #$hmmhash{$domain->hmmname} = 1;
	       
	       #Push each genomic sequence corresponding to a given Hmm
	       #push (@{$hmmhash{$domain->feature2->seqname}},$subseq);
	       
		   my $hmmnames = $domain->feature2->seqname;

		   my $hmmtemp = "/tmp/hmmtemp".$$;
		   my $hmmfile = $hmmanalysis->db_file();
		   system("hmmfetch $hmmfile $hmmnames > $hmmtemp ");
		   
		   my $gw = Bio::EnsEMBL::Pipeline::Runnable::GenewiseHmm->new( -hmmfile => $hmmtemp, -genomic => $self->seq);
		   $gw->run;
		   
		   print STDERR "HERE1\n";
		   
		   $self->convert_output($gw);
		   print STDERR "OUTPUT".$self->output."\n";
		   
		   $self->write_output;
	       }
	      
	   }
       }
   }
   return;
}



=head2 convert_output

  Title   :   convert_output
  Usage   :   $self->convert_output
  Function:   converts output from each runnable into gene predictions
  Returns :   nothing, but $self->output contains results
  Args    :   none

=cut

sub convert_output {
  my ($self,$gw) =@_;
  
  print STDERR "HERE2\n";

  my $trancount = 1;
  my $genetype;
  my $runnable = $gw;
      print STDERR "HERE3\n";
    my $anaAdaptor = $self->dbobj->get_AnalysisAdaptor;

	#use logic name from analysis object if possible, else take $genetype;
	my $anal_logic_name = "Genewise";#($self->analysis->logic_name)	?	$self->analysis->logic_name : $genetype	;	
	
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
	 -logic_name      => 'Genewise',
	 -module          => 'BlastMiniGenewise',
      );
  }

    my @results = $runnable->output;
  print STDERR "RESULT: @results\n";

    my @genes = $self->make_genes($genetype, $analysis_obj, \@results);

  print STDERR "GENESSS: @genes\n";

    $self->output(@genes);
  print STDERR "OUTPUT0".$self->output."\n";
  
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
  my $contig = $self->vcontig;
  print STDERR "CONTIG: $contig\n";

  my @tmpf   = @$results;
  my @genes;

  foreach my $tmpf (@tmpf) {
      print STDERR "TMP: $tmpf\n";
    my $gene       = new Bio::EnsEMBL::Gene;
    my $transcript = $self->_make_transcript($tmpf, $contig, $genetype, $analysis_obj);

    $gene->type($genetype);
    $gene->analysis($analysis_obj);
    $gene->add_Transcript($transcript);

      print STDERR "GENE: ".$gene."\n";

    push (@genes, $gene);
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
    
    $exon->contig_id($contig->internal_id);
    $exon->start($exon_pred->start);


    $exon->end  ($exon_pred->end);
    $exon->strand($exon_pred->strand);
    
    $exon->phase(1);
    $exon->attach_seq($contig);
    
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
    
    print STDERR "num exons: " . scalar(@exons) . "\t".$exons[0]->phase."\n";

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

    print STDERR "OUTPUT: ".$self->output."\n";

  GENE: foreach my $gene ($self->output) {

      # do a per gene eval...
      eval {
	$gene_adaptor->store($gene);
	print STDERR "wrote gene " . $gene->dbID . "\n";
      }; 
      if( $@ ) {
	  print STDERR $@;
	  print STDERR "UNABLE TO WRITE GENE\n\n$@\n\nSkipping this gene\n";
      }
	    
  }
   
}


=head2 dbobj

 Title   : dbobj
 Usage   : $obj->dbobj($newval)
 Function: 
 Example : 
 Returns : value of dbobj
 Args    : newvalue (optional)


=cut

sub dbobj{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'dbobj'} = $value;
    }
    return $obj->{'dbobj'};

}


sub analysis {
    my ($self, $analysis) = @_;
    
    if ($analysis)
    {
        $self->throw("Not a Bio::EnsEMBL::Analysis object")
            unless ($analysis->isa("Bio::EnsEMBL::Analysis"));
        $self->{'_analysis'} = $analysis;
        $self->parameters($analysis->parameters);
    }
    return $self->{'_analysis'}
}

=head2 input_id

 Title   : input_id
 Usage   : $obj->input_id($newval)
 Function: 
 Example : 
 Returns : value of input_id
 Args    : newvalue (optional)


=cut

sub input_id{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'input_id'} = $value;
    }
    return $obj->{'input_id'};

}

=head2 seq

 Title   : seq
 Usage   : $obj->seq($newval)
 Function: 
 Example : 
 Returns : value of seq
 Args    : newvalue (optional)


=cut

sub seq{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'seq'} = $value;
    }
    return $obj->{'seq'};

}

=head2 output

 Title   : output
 Usage   : $obj->output($newval)
 Function: 
 Returns : value of output
 Args    : newvalue (optional)


=cut

sub output{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'output'} = $value;
    }
    return $obj->{'output'};

}
