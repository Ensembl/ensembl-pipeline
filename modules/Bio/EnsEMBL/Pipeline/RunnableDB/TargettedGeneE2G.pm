#
# Ensembl module for Bio::EnsEMBL::Pipeline::RunnableDB::TargettedE2G.pm
#
# Cared for by EnsEMBL  <ensembl-dev@ebi.ac.uk>
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::TargettedGeneE2G

=head1 SYNOPSIS

my $t_e2g = new Bio::EnsEMBL::Pipeline::RunnableDB::TargettedGeneE2G(
                                                                      '-db_obj'      => $dbobj,
                                                                      '-golden_path' => $gp,
                                                                      '-input_id'    => $input_id
                                                                    );

$t_e2g->fetch_input();
$t_e2g->run();
$t_e2g->output();
$t_e2g->write_output(); #writes to DB

=head1 DESCRIPTION

Run Targetted Est2Genome for building Gene objects

=head1 CONTACT

Ensembl - ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Pipeline::RunnableDB::TargettedGeneE2G;

use vars qw(@ISA);
use strict;
use Storable qw(dclone);
# Object preamble - inheriets from Bio::EnsEMBL::Root


use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Getseqs;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;
use Bio::EnsEMBL::Pipeline::Runnable::ExonerateMiniEst2Genome;
use Bio::SeqIO;
use Bio::EnsEMBL::Pipeline::GeneConf qw (
					 GB_GOLDEN_PATH
					 GB_TARGETTED_PROTEIN_INDEX
					 GB_TARGETTED_CDNA_INDEX
					 );

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  my ($path,$cdna_seqfetcher) = $self->_rearrange([qw(GOLDEN_PATH CDNA_SEQFETCHER)], @args);

  # golden path
  if(!defined $path){
    $path = $GB_GOLDEN_PATH;
  }

  $path = 'UCSC' unless (defined $path && $path ne '');
  $self->dbobj->assembly_type($path);

  # broken by test_runnableDB 
  # protein sequence fetcher
  if(!defined $self->seqfetcher) {
    my $seqfetcher = $self->make_seqfetcher($GB_TARGETTED_PROTEIN_INDEX);
    $self->seqfetcher($seqfetcher);
  }

  # cdna sequence fetcher
  if(defined $cdna_seqfetcher){
    $self->throw("[$cdna_seqfetcher] is not a Bio::DB::RandomAccessI\n") unless $cdna_seqfetcher->isa("Bio::DB::RandomAccessI");
    $self->cdna_seqfetcher($cdna_seqfetcher);
  }
  else {
    my $seqfetcher = $self->make_seqfetcher($GB_TARGETTED_CDNA_INDEX);
    $self->cdna_seqfetcher($seqfetcher);
  }

  return $self;
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
  my ( $self, $index ) = @_;
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

  return $seqfetcher;
}

=head2 cdna_seqfetcher

 Title   : cdna_seqfetcher
 Usage   :
 Function: get/set
 Example :
 Returns : 
 Args    :


=cut

sub cdna_seqfetcher {
    my( $self, $value ) = @_;    
    if ($value) {
        #need to check if passed sequence is Bio::DB::RandomAccessI object
        $value->isa("Bio::DB::RandomAccessI") || 
	    $self->throw("Input isn't a Bio::DB::RandomAccessI");
        $self->{'_cdna_seqfetcher'} = $value;
    }
    return $self->{'_cdna_seqfetcher'};
}

=head2 cdna_id

 Title   : cdna_id
 Usage   :
 Function: get/set
 Example :
 Returns : 
 Args    :


=cut

sub cdna_id {
    my( $self, $value ) = @_;    
    if ($value) {
        $self->{'_cdna_id'} = $value;
    }
    return $self->{'_cdna_id'};
}

=head2 protein_id

 Title   : protein_id
 Usage   :
 Function: get/set
 Example :
 Returns : 
 Args    :


=cut

sub protein_id {
    my( $self, $value ) = @_;    
    if ($value) {
        $self->{'_protein_id'} = $value;
    }
    return $self->{'_protein_id'};
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
  my ($self,@args) = @_;

  my $entry = $self->input_id;
  my $chrname;
  my $start;
  my $end;
  my $protein_id; 
  my $cdna_id;

  # input format: chr12:10602496,10603128:Q9UGV6:AC00012
  # or chr12:10602496,10603128:Q9UGV6 if no cDNA
  if( !(($entry =~ /(\S+):(\d+),(\d+):(\S+):(\S+)/) || ($entry =~ /(\S+):(\d+),(\d+):(\S+):/))) {
      $self->throw("Not a valid input id... $entry");
  }
  
  $chrname    = $1;
  $protein_id = $4;
  $cdna_id    = $5;
  $start   = $2;
  $end     = $3;

  if ($2 > $3) { # let blast sort it out
      $start  = $3;
      $end    = $2;
  }

  
  my $sgpa = $self->dbobj->get_StaticGoldenPathAdaptor();
  my $vc = $sgpa->fetch_VirtualContig_by_chr_start_end($chrname,$start-10000,$end+10000);
  
  $self->vcontig($vc);
  $self->cdna_id($cdna_id);
  $self->protein_id($protein_id);
  
  # genewise runnable
  # repmasking?
  my $r = Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise->new( '-genomic'    => $vc->primary_seq,
								    '-ids'        => [ $protein_id ] ,
								    '-seqfetcher' => $self->seqfetcher);
 
  $self->runnable($r);

  # est2genome runnable
  return unless defined($cdna_id);

  my $cdna;
  eval{
    $cdna = $self->cdna_seqfetcher->get_Seq_by_acc($cdna_id);
  };
  if($@) {
    $self->warn("problem fetching cdna sequence for [$cdna_id], will not be able to build UTR\n");
  }
  else{
    $self->{'_tmpfile'} = "/tmp/tge2g_" . $$ . ".fa";
    my $cdnafile = $self->{'_tmpfile'};
    
    my $seqout = new Bio::SeqIO('-file' => ">$cdnafile" , '-format' => 'Fasta');
    $seqout->write_seq($cdna);
    
    # repmasking?
    my $e2g = new Bio::EnsEMBL::Pipeline::Runnable::ExonerateMiniEst2Genome('-genomic'    => $vc->primary_seq, 
									    '-queryseq'   => $cdnafile,
									    '-seqfetcher' => $self->cdna_seqfetcher);
    
    $self->e2g_runnable($e2g);
  }
}


=head2 run

 Title   : run
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub run {
   my ($self,@args) = @_;
   
   $self->runnable->run();
   $self->convert_gw_output;

   # check to see we have genewise genes before running e2g and trying to combine results
   my @gw_genes = $self->gw_genes;
   if (scalar(@gw_genes) && defined $self->e2g_runnable) {
     $self->e2g_runnable->run;
     $self->convert_e2g_output;
     $self->combine_genes;
   }

   # clean up tmpfile
     my $tmpfile = $self->{'_tmpfile'};
     unlink $tmpfile;
   # remap to raw contig coords
   my @remapped = $self->remap_genes();
   $self->{'_output'} = \@remapped;
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
   my ($self,@args) = @_;

   return @{$self->{'_output'}};
}



=head2 e2g_runnable

 Title   : e2g_runnable
 Usage   : $obj->e2g_runnable($newval)
 Function: 
 Returns : value of e2g_runnable
 Args    : newvalue (optional)


=cut

sub e2g_runnable{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_e2g_runnable'} = $value;
    }
    return $obj->{'_e2g_runnable'};

}

=head2 runnable

 Title   : runnable
 Usage   : $obj->runnable($newval)
 Function: 
 Returns : value of runnable
 Args    : newvalue (optional)


=cut

sub runnable{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_runnable'} = $value;
    }
    return $obj->{'_runnable'};
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

    my $gene_adaptor = $self->dbobj->get_GeneAdaptor;

  GENE: foreach my $gene ($self->output) {	
      # do a per gene eval...
      eval {
	  $gene_adaptor->store($gene);
	  print STDERR "wrote gene dbID " . $gene->dbID . "\n";
      }; 
      if( $@ ) {
	  print STDERR "UNABLE TO WRITE GENE\n\n$@\n\nSkipping this gene\n";
      }
	    
  }
   
}

=head2 e2g_genes

 Title   : e2g_genes
 Usage   :
 Function: get/set for e2g gene array
 Example :
 Returns : 
 Args    :


=cut

sub e2g_genes {
  my ($self, @genes) = @_;

  if (!defined($self->{'_e2g_genes'})) {
    $self->{'_e2g_genes'} = [];
  }

  if (@genes) {
    push(@{$self->{'_e2g_genes'}},@genes);
  }
  
  return @{$self->{'_e2g_genes'}};
}

=head2 gw_genes

 Title   : gw_genes
 Usage   :
 Function: get/set for genewise gene array
 Example :
 Returns : 
 Args    :


=cut

sub gw_genes {
  my ($self, @genes) = @_;
  if (!defined($self->{'_gw_genes'})) {
    $self->{'_gw_genes'} = [];
  }

  if (scalar(@genes)) {
    push(@{$self->{'_gw_genes'}},@genes);
  }
  
  return @{$self->{'_gw_genes'}};
}

=head2 combined_genes

 Title   : combined_genes
 Usage   :
 Function: get/set for combined gene array
 Example :
 Returns : 
 Args    :


=cut

sub combined_genes {
  my ($self, @genes) = @_;

  if (!defined($self->{'_combined_genes'})) {
    $self->{'_combined_genes'} = [];
  }

  if (@genes) {
    push(@{$self->{'_combined_genes'}},@genes);
  }
  
  return @{$self->{'_combined_genes'}};
}

=head2 convert_e2g_output

 Title   : convert_e2g_output
 Usage   :
 Function: converts the output from Est2Genome into genes
 Example :
 Returns : 
 Args    :


=cut

sub convert_e2g_output {
  my ($self) = @_;
  
  my @results = $self->e2g_runnable->output;
  
  foreach my $gene(@results) {
    foreach my $ex($gene->sub_SeqFeature){
      # exonerate has no concept of phase, but remapping will fail if this is unset
#      $ex->phase(-1);
      $ex->phase(0);
      foreach my $sf($ex->sub_SeqFeature){
	# strands 
	if($sf->strand != $sf->hstrand){
	  $sf->strand(-1);
	  $sf->hstrand(1);
	  $ex->strand(-1);
	}
      }
    }
  }  

  my $count = 1;
  my $genetype = "TGE_e2g";
  
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
       -module          => 'TargettedGeneE2G',
      );
  }

  my @genes = $self->make_genes($count, $genetype, $analysis_obj, \@results);
  
  $self->e2g_genes(@genes);  
}

=head2 convert_gw_output

 Title   : convert_gw_output
 Usage   :
 Function: converts output from Genewise into genes
 Example :
 Returns : 
 Args    :


=cut

sub convert_gw_output {
  my ($self) = @_;
  my $count = 1;
  my $genetype = 'TGE_gw';
  my @results  = $self->runnable->output;

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
       -module          => 'TargettedGeneE2G',
      );
  }

  my @genes    = $self->make_genes($count, $genetype, $analysis_obj, \@results);

  my @checked  = $self->check_gw_genes(@genes);

  # check for stops?

  $self->gw_genes(@checked);
  
}

=head2 check_gw_genes

 Title   : check_gw_genes
 Usage   :
 Function: Checks coverage of parent protein, and checks all exons are sane
 Example :
 Returns : array of checked Bio::EnsEMBL::Gene
 Args    : array of Bio::EnsEMBL::Gene


=cut

sub check_gw_genes{
  my ($self, @genes) = @_;

  my @checked;

  GENE:
  foreach my $gene(@genes) {
     # single exon genes must cover 80% of the protein length; multi exon ones 
    # can get away with 25% though we might raise this later
    
    my $threshold = 80;    
    my @t = $gene->each_Transcript;
    my @exons = $t[0]->get_all_Exons;
    if(scalar(@exons) > 1){
      $threshold = 25;
    }

    my $covered = $self->_check_coverage($gene, $threshold);
    if(!$covered){
#      my $msg = "rejecting gene because < $threshold% coverage of parent protein\n";
#      $self->warn($msg);
      next GENE;
    }

    if($self->validate_gene($gene)){
      push (@checked, $gene);
    }
    else{
      my $msg = "rejecting gene\n";
      $self->warn($msg);
      next GENE;
    }
  }
  
  return @checked;
  
}

=head2 combine_genes

 Title   : combine_genes
 Usage   :
 Function: 
 Example :
 Returns : 
 Args    :


=cut

sub combine_genes{
  my ($self) = @_;

  # make array of gw_pred_genes -> merge potentially frameshifted exons together; add
  # component exons as sub_SeqFeatures so they can be retrieved later
  my @merged_gw_genes = $self->_merge_gw_genes;

  my $genetype = 'combined_gw_e2g';

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
       -module          => 'TargettedGeneE2G',
      );
  }


  my @newtrans = $self->_make_newtranscripts($genetype, $analysis_obj, @merged_gw_genes);

  $analysis_obj->gff_feature('gene');

  # make some lovely genes
  my @genes;
  my $count=0;

  foreach my $trans(@newtrans){
    $trans->sort;
    my $gene = new Bio::EnsEMBL::Gene;
    $gene->type($genetype);
    $gene->add_Transcript($trans);
    $gene->analysis($analysis_obj);

    if($self->validate_gene($gene)){
      push (@genes,$gene);
      $count++;
    }
  }

  $self->combined_genes(@genes);
  
}

=head2 _merge_gw_genes

 Title   : _merge_gw_genes
 Usage   :
 Function: merges adjacent exons if they are frameshifted; stores component exons as subSeqFeatures of the merged exon
 Example :
 Returns : 
 Args    :


=cut

sub _merge_gw_genes {
  my ($self) = @_;

  my @merged;
  my $count = 1;
  foreach my $gwg($self->gw_genes){
    my $gene = new Bio::EnsEMBL::Gene;
    $gene->type('combined');
    
    my @pred_exons;
    my $ecount = 0;
    
    # order is crucial
    my @trans = $gwg->each_Transcript;
    if(scalar(@trans) != 1) { $self->throw("expected one transcript for $gwg\n"); }
    
  EXON:      
    foreach my $exon($trans[0]->get_all_Exons){
      my $previous_exon;
      
      if ($ecount && $pred_exons[$ecount-1]){
	$previous_exon = $pred_exons[$ecount-1];
      }
      
      $ecount++;
      
      # genewise frameshift? we treat two exons separated by max 10 bases as a single exon
      if( defined($previous_exon) && abs($exon->start - $previous_exon->end) <= 10 ){
	# combine the two
	$previous_exon->end($exon->end);
	$previous_exon->add_sub_SeqFeature($exon,'');
	next EXON;
      }
      
      else{
	# make a new Exon - clone $exon
	my $cloned_exon = dclone($exon);
	$cloned_exon->attach_seq($self->vcontig->primary_seq);
	$cloned_exon->add_sub_SeqFeature($exon,'');
	push(@pred_exons, $cloned_exon);
      }
      
    }

    # transcript
    my $merged_transcript   = new Bio::EnsEMBL::Transcript;
    foreach my $pe(@pred_exons){
      $merged_transcript->add_Exon($pe);
    }
    
    $merged_transcript->sort;

    my $cloned_translation = dclone($trans[0]->translation);
    $merged_transcript->translation($cloned_translation);
    
    # and gene
    $gene->add_Transcript($merged_transcript);
    push(@merged, $gene);
    $count++;
  }

  return @merged;
}

=head2 _make_newtranscripts

 Title   : _make_newtranscripts
 Usage   :
 Function: makes new transcripts by combining the genewise and est2genome predictions. Its a monster.
 Example :
 Returns : 
 Args    :


=cut

sub _make_newtranscripts {
  my ($self, $genetype, $analysis_obj, @merged_gw_genes) = @_;
  my @combined_transcripts  = ();
  
  # we need to compare each genewise prediction with each est_genome 
  # prediction and tie them together into "combined" genes

 GENEWISE:
  foreach my $gene(@merged_gw_genes) {
    my $foundtrans = 0;  

    # should be only 1 transcript
    my @gw_tran  = $gene->each_Transcript;
    my @gw_exons = $gw_tran[0]->get_all_Exons; # ordered array of exons
    my $strand   = $gw_exons[0]->strand;

    if($gw_exons[$#gw_exons]->strand != $strand){
      $self->warn("first and last gw exons have different strands - can't make a sensible combined gene\n");
      next GENEWISE;
    }
    
  E2G:
    foreach my $eg($self->e2g_genes){
      next GENEWISE if $foundtrans;

      my @egtran = $eg->each_Transcript;
      my @e2g_exons  = $egtran[0]->get_all_Exons; # ordered array of exons
      
      # OK, let's see if we need a new gene
      # base it on the existing genewise one
      my $newtranscript = dclone($gw_tran[0]);
      my $translation   = dclone($gw_tran[0]->translation);
      $newtranscript->translation($translation);
      my $eecount = 0;
      
      $newtranscript->translation->start_exon($newtranscript->start_exon);
      $newtranscript->translation->end_exon($newtranscript->end_exon);

      # check strands are consistent
      foreach my $ee(@e2g_exons){
	if ($ee->strand != $strand){
          $self->warn("gw and e2g exons have different strands - can't combine genes\n") ;
          next GENEWISE;
        }
	
	# single exon genewise prediction?
	if(scalar(@gw_exons) == 1) {
	  $newtranscript = $self->transcript_from_single_exon_genewise( $ee, 
									$gw_exons[0], 
									$newtranscript, 
									$translation, 
									$eecount, 
									@e2g_exons);
	  
	}
	
	else {
	  $newtranscript = $self->transcript_from_multi_exon_genewise($ee, 
								      $newtranscript, 
								      $translation, 
								      $eecount,
								      $gene, 
								      $eg)
	}
	
	# increment the exon
	$eecount++;
	
      } # end foreach my $ee

      # the new transcript is made from a merged genewise gene
      # check the transcript and expand frameshifts in all but original 3' gw_exon
      if (defined($newtranscript)){

	foreach my $ex($newtranscript->get_all_Exons){

	  if(scalar($ex->sub_SeqFeature) > 1 ){
	    my @sf    = $ex->sub_SeqFeature;
	    my $first = shift(@sf);

	    $ex->end($first->end);

	    # add back the remaining component exons
	    foreach my $s(@sf){
	      $newtranscript->add_Exon($s);
	      $newtranscript->sort;
	    }
	    # flush the sub_SeqFeatures
	    $ex->flush_sub_SeqFeature;
	  }
	}
	
	# dclone messes up database handles
	foreach my $ex($newtranscript->get_all_Exons){

	  $ex->attach_seq($self->vcontig);
	  $ex->contig_id($self->vcontig->id);
	  # add new analysis object to the supporting features
	  foreach my $sf($ex->each_Supporting_Feature){
	    $sf->analysis($analysis_obj);
	    $sf->source_tag($genetype);
	  }
	}
	
      # check translation is the same as for the genewise gene we're built from
      GWG:	
	foreach my $gwg($self->gw_genes) {
	  $foundtrans = $self->compare_transcripts($gwg, $newtranscript);
 
	  if ($foundtrans == 1){
	    push (@combined_transcripts, $newtranscript); 	
	    last GWG;
	  }
	}
	# did we find any genewise gene whose translation matched our combined gene?
	if(!$foundtrans){
	  $self->warn("UTR prediction is not the same as genewise prediction - discarding it\n");
	}
	
      }
    }
  }
  return @combined_transcripts;
  
}

sub compare_transcripts{
  my ($self, $genewise_gene, $combined_transcript) = @_;
  my @genewise_transcripts = $genewise_gene->each_Transcript;
  if(scalar(@genewise_transcripts != 1)) {
    $self->warn("Panic! Got " . scalar(@genewise_transcripts) . " transcripts, expecting only 1!\n");
    return 0;
  }
	  
  my $genewise_translation;
  my $combined_translation;
  
  eval {
    $genewise_translation = $genewise_transcripts[0]->translate;
  };

  if ($@) {
    print STDERR "Couldn't translate genewise gene\n";
    return 0;
  }
  $@ = '';

#  print STDERR "translation: \n";
#  my $seqio = Bio::SeqIO->new(-fh => \*STDERR);
#  print STDERR "genewise: \n";	      
#  $seqio->write_seq($genewise_translation);


  eval{
    $combined_translation = $combined_transcript->translate;
  };


  if ($@) {
    print STDERR "Couldn't translate combined gene\n";
    return 0;
  }
	  
#  print STDERR "combined: \n";
#  $seqio->write_seq($combined_translation); 
#  print STDERR "\n ";
  
  if($genewise_translation->seq eq $combined_translation->seq) {
    return 1;
  }

  return 0;
}

sub transcript_from_single_exon_genewise {
  my ($self, $eg_exon, $gw_exon, $transcript, $translation, $exoncount, @e2g_exons) = @_;

  if ($gw_exon->start >= $eg_exon->start && $gw_exon->end <= $eg_exon->end){
#    print STDERR "single exon gene, " . $gw_exon->strand  .  " strand\n";	    

    # modify the coordinates of the first exon in $newtranscript
    my $ex = $transcript->start_exon;

    $ex->start($eg_exon->start);
    $ex->end($eg_exon->end);

    # need to explicitly set the translation start & end exons here.
    $translation->start_exon($ex);
    # end_exon may be adjusted by 3' coding exon frameshift expansion. Ouch.
    $translation->end_exon($ex);
    
    # need to add back exons, both 5' and 3'
    $self->add_5prime_exons(\$transcript, $exoncount, @e2g_exons);
    $self->add_3prime_exons(\$transcript, $exoncount, @e2g_exons);
	    
    # need to deal with translation start and end this time - varies depending on strand
    if($gw_exon->strand == 1){
      my $diff = $gw_exon->start - $ex->start;
      my $tstart = $translation->start;
      my $tend = $translation->end;
	      
      $translation->start($tstart + $diff);
      $translation->end($tend + $diff);
    }

    elsif($gw_exon->strand == -1){
      my $diff = $ex->end - $gw_exon->end;
      my $tstart = $translation->start;
      my $tend = $translation->end;
      $translation->start($tstart+$diff);
      $translation->end($tend + $diff);
    }
    
    
    # expand frameshifted exons back from one exon to multiple exons
    if(scalar($ex->sub_SeqFeature) > 1){
#      print STDERR "uh-oh frameshift\n";
      my @sf = $ex->sub_SeqFeature;
      
      # save current start and end
      my $cstart = $ex->start;
      my $cend   = $ex->end;
      
      # get first exon - this has same id as $ex
      my $first = shift(@sf);
      $ex->end($first->end);
      
      # get last exon
      my $last = pop(@sf);
      $last->end($cend);
      $transcript->add_Exon($last);
      # and adjust translation end
      $translation->end_exon($last);
      
      # get any remaining exons
      foreach my $s(@sf){
	$transcript->add_Exon($s);
	$transcript->sort;
      }
      # flush the sub_SeqFeatures
      $ex->flush_sub_SeqFeature;
    }	      
  }
  return $transcript;
}


sub transcript_from_multi_exon_genewise {
  my ($self, $current_exon, $transcript, $translation, $exoncount, $gw_gene, $eg_gene) = @_;
  
  my @gwtran  = $gw_gene->each_Transcript;
  my @gwexons = $gwtran[0]->get_all_Exons;
  
  my @egtran  = $eg_gene->each_Transcript;
  my @egexons = $egtran[0]->get_all_Exons;

  # compare to the first genewise exon
  if($gwexons[0]->strand == 1){
  FORWARD:
    if ($gwexons[0]->end == $current_exon->end && $current_exon->start <= $gwexons[0]->start){
#      print STDERR "5' exon match!\n";
      # modify the coordinates of the first exon in $newtranscript
      my $ex = $transcript->start_exon;
      $ex->start($current_exon->start);

      # add all the exons from the est2genome transcript, previous to this one
      $self->add_5prime_exons(\$transcript, $exoncount, @egexons);
      
      # fix translation start 
      # take what it was for the gw gene, and add on the extra
      my $tstart = $translation->start;
      $tstart += ($gwexons[0]->start - $ex->start);
      $translation->start($tstart);
      
    } # end 5' exon
    
    elsif ($gwexons[$#gwexons]->start == $current_exon->start && $current_exon->end >= $gwexons[$#gwexons]->end){
#      print STDERR "3' exon match\n";
      
      # modify the coordinates of the last exon in $newtranscript
      my $ex = $transcript->end_exon;
      $ex->end($current_exon->end);

      $self->expand_3prime_exon(\$ex, \$transcript);

      # need to explicitly set the translation end exon for translation to work out
      my $end_ex = $transcript->end_exon;
      $translation->end_exon($end_ex);

      $self->add_3prime_exons(\$transcript, $exoncount, @egexons);

    } # end 3' exon
    
  }
  
  elsif($gwexons[0]->strand == -1){
  REVERSE:
    if ($gwexons[0]->start == $current_exon->start && $current_exon->end >= $gwexons[0]->end){
#      print STDERR "5' exon match!\n";
      
      # modify the coordinates of the first exon in $newtranscript
      my $ex = $transcript->start_exon;
      $ex->end($current_exon->end);

      # need to explicitly set the translation start exon for translation to work out
      $translation->start_exon($ex);

      $self->add_5prime_exons(\$transcript, $exoncount, @egexons);
      
      # need to deal with translation start
      my $tstart = $translation->start;
      my $diff = $current_exon->end - $gwexons[0]->end;
      $translation->start($tstart+$diff);

    } # end 5' exon
    
    elsif ($gwexons[$#gwexons]->end == $current_exon->end && $current_exon->start <= $gwexons[$#gwexons]->start){
#      print STDERR "3' exon match\n";
      
      # modify the coordinates of the last exon in $newtranscript
      my $ex = $transcript->end_exon;
      $ex->start($current_exon->start);

      $self->expand_3prime_exon(\$ex, \$transcript);

      # need to explicitly set the translation end exon for translation to work out
      my $end_ex = $transcript->end_exon;
      $translation->end_exon($end_ex);

      $self->add_3prime_exons(\$transcript, $exoncount, @egexons);
      
    } # end 3' exon
  }  
 return $transcript; 
}

sub add_5prime_exons{
my ($self, $transcript, $exoncount, @e2g_exons);
      # add all the exons from the est2genome transcript, previous to this one
      my $c = 0;
      while($c < $exoncount){
	$$transcript->add_Exon($e2g_exons[$c]);
	$$transcript->sort;
	$c++;
      }
      
}

# $exon is the terminal exon in the genewise transcript, $transcript. We need
# to expand any frameshifts we merged in the terminal genewise exon. 
sub expand_3prime_exon{
my ($self, $exon, $transcript) = @_;
      if(scalar($$exon->sub_SeqFeature) > 1){
#	print STDERR "3' exon frameshift\n";
	my @sf = $$exon->sub_SeqFeature;
	my $last = pop(@sf);

	# sort out start, id & phase
	$$exon->start($last->start); # but don't you dare touch the end!
	$$exon->id($last->id);
	$$exon->phase($last->phase);

	# add back the remaining component exons
	foreach my $s(@sf){
	  $$transcript->add_Exon($s);
	  $$transcript->sort;
	}
	# flush the sub_SeqFeatures so we don't try to re-expand later
	$$exon->flush_sub_SeqFeature;
      }
      
}

# $exoincount tells us which position in the array 
# of e2g exons corresponds to the end of the genewise transcript so we add back 
# exons 3' to that position.
# $exon and $transcript are references to Exon and Transcript objects.
sub add_3prime_exons {
my ($self, $transcript, $exoncount, @e2g_exons) = @_;
# need to deal with frameshifts - 3' exon is a special case as its end might have changed

      # add all the exons from the est2genome transcript, subsequent to this one
      my $c = $#e2g_exons;
      while($c > $exoncount){
	$$transcript->add_Exon($e2g_exons[$c]);
	$$transcript->sort;
	$c--;
      }

}


=head2 make_genes

 Title   : make_genes
 Usage   : $self->make_genes($count, $genetype, $analysis_obj, $runnable)
 Function: converts the output from $runnable into Bio::EnsEMBL::Genes in
           $contig(VirtualContig) coordinates. The genes have type $genetype, 
           and have $analysis_obj attached. Each Gene has a single Transcript, 
           which in turn has Exons(with supporting features) and a Translation
 Example : 
 Returns : array of Bio::EnsEMBL::Gene
 Args    : $count: integer, $genetype: string, $analysis_obj: Bio::EnsEMBL::Analysis, 
           $runnable: Bio::EnsEMBL::Pipeline::RunnableI


=cut

sub make_genes {
  my ($self, $count, $genetype, $analysis_obj, $results) = @_;
  my $contig = $self->vcontig;
  my @genes;
  
#$self->throw("[$analysis_obj] is not a Bio::EnsEMBL::Pipeline::Analysis\n") unless defined($analysis_obj) && $analysis_obj->isa("Bio::EnsEMBL::Pipeline::Analysis");
$self->throw("[$analysis_obj] is not a Bio::EnsEMBL::Analysis\n") unless defined($analysis_obj) && $analysis_obj->isa("Bio::EnsEMBL::Analysis");

  foreach my $tmpf (@$results) {
    my $gene   = new Bio::EnsEMBL::Gene;
    $gene->type($genetype);

    my $transcript = $self->_make_transcript($tmpf,$self->vcontig,$genetype,$count, $analysis_obj);
	
    # add transcript to gene
    $gene->analysis($analysis_obj);
    $gene->add_Transcript($transcript);
    $count++;

    # and store it
    push(@genes,$gene);
  }
  return @genes;
}

=head2 remap_genes

 Title   : remap_genes
 Usage   :
 Function: 
 Example :
 Returns : 
 Args    :


=cut

sub remap_genes {
  my ($self) = @_;
  my @newf;  
  my $contig = $self->vcontig;

  my @genes = $self->gw_genes;
  push(@genes, $self->e2g_genes); # ??? do we want these?
  push(@genes, $self->combined_genes);

GENE:  foreach my $gene (@genes) {

    my @t = $gene->each_Transcript;
    my $tran = $t[0];

    # check that it translates - not the est2genome genes
    if($gene->type eq 'TGE_gw' || $gene->type eq 'combined_gw_e2g'){
      
      my $translates = $self->check_translation($tran);
      if(!$translates){
	my $msg = "discarding gene - tranlation has stop codons\n";
	$self->warn($msg);
	next GENE;
      }
  }
    eval {
      my $genetype = $gene->type;
      my $newgene = $contig->convert_Gene_to_raw_contig($gene);
      # need to explicitly add back genetype and analysis.
      $newgene->type($genetype);
      $newgene->analysis($gene->analysis);
      push(@newf,$newgene);

      # sort out supporting feature coordinates
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

      # is this a special case single coding exon gene with UTRS?
      if($tran->translation->start_exon() eq $tran->translation->end_exon() 
	 && $gene->type eq 'combined_gw_e2g'){
#	print STDERR "single coding exon, with UTRs\n";
	
	# problems come about when we switch from + strand on FPC contig to - strand on raw contig.
	my $fpc_strand;

# may need to compare exon objects ... with new schema
	foreach my $exon($tran->get_all_Exons) {
#	  if ($exon->id eq $tran->translation->start_exon_id()) {
	  if ($exon eq $tran->translation->start_exon()) {
	    $fpc_strand = $exon->strand;
	    last;
	  }
	}
	
	foreach my $tran ($newgene->each_Transcript) {
	  foreach my $exon($tran->get_all_Exons) {

# oh dear oh dear oh dear
#	    if ($exon->id eq $tran->translation->start_exon_id()) {
	    if ($exon eq $tran->translation->start_exon()) {
	      if($fpc_strand == 1 && $exon->strand == -1){
#		print STDERR "fpc strand 1, raw strand -1 - flipping translation start/end\n";
		$exon->end($exon->end - ($tran->translation->start -1));
		$exon->start($exon->end - ($tran->translation->end -1));
	      }
	    }
	  }
	}
      } # end special case single coding exon

    };

    # did we throw exceptions?
    if ($@) {
      print STDERR "Couldn't reverse map gene:  [$@]\n";
    }
  }

  return @newf;
}

=head2 check_translation

 Title   : check_translation
 Usage   :
 Function: checks how much of the parent protein is covered by the genewise prediction
 Example :
 Returns : 1 if transcript translates with no stops, otherwise 0
 Args    :


=cut

sub check_translation {
  my ($self, $transcript) = @_;
  my $tseq;

  eval{
    $tseq = $transcript->translate;
  };

  if((!defined $tseq) || ($@)){
    my $msg = "problem translating :\n$@\n";
    $self->warn($msg);
    return 0;
  }

  if ($tseq->seq =~ /\*/ ) {
    return 0;
  }
  else{
    return 1;
  }
}

=head2 _check_coverage

 Title   : _check_coverage
 Usage   :
 Function: checks how much of the parent protein is covered by the genewise prediction
 Example :
 Returns : 1 if > $threshold% coverage, otherwise 0
 Args    :


=cut

sub _check_coverage{
  my ($self, $gene, $threshold) = @_;
  my $pstart = 0;
  my $pend = 0;
  my $protname = $self->protein_id;
  my $plength;
  my $fetcher = new Bio::EnsEMBL::Pipeline::SeqFetcher;

  my @gw_tran = $gene->each_Transcript;
  
  my $matches = 0;

  foreach my $exon($gw_tran[0]->get_all_Exons) {
    $pstart = 0;
    $pend   = 0;
    
    foreach my $f($exon->each_Supporting_Feature){
      #      print STDERR $f->hseqname . " " . $f->hstart . " " . $f->hend . "\n";
      
      if (!defined($protname)){
	$protname = $f->hseqname;
      }
      if($protname ne $f->hseqname){
	warn("$protname ne " . $f->hseqname . "\n");
      }
      
      if((!$pstart) || $pstart > $f->hstart){
	$pstart = $f->hstart;
      }
      
      if((!$pend) || $pend < $f->hend){
	$pend= $f->hend;
      }
    }
    $matches += ($pend - $pstart + 1);
  }
  
  my $seq; 
  eval{
    $seq = $self->seqfetcher->get_Seq_by_acc($protname);
  };
  if ($@) {
    $self->throw("Error fetching sequence for [$protname]: [$@]\n");
  }
  
  $self->throw("No sequence fetched for [$protname]\n") unless defined $seq;
  
  $plength = $seq->length;

  if(!defined($plength) || $plength == 0){
    warn("no sensible length for $protname - can't get coverage\n");
    return 0;
  }

  my $coverage = $matches/$plength;
  $coverage *= 100;
  if ($coverage < $threshold){
    $self->warn ("Coverage of $protname is only $coverage - will be rejected\n");
    return 0;
  }
  
  print STDERR "Coverage of $protname is $coverage - will be accepted\n";
  return 1;
}

=head2 _make_transcript

 Title   : make_transcript
 Usage   : $self->make_transcript($gene, $contig, $genetype, $count, $analysis_obj)
 Function: makes a Bio::EnsEMBL::Transcript from a SeqFeature representing a gene, 
           with sub_SeqFeatures representing exons.
 Example :
 Returns : Bio::EnsEMBL::Transcript with Bio::EnsEMBL:Exons(with supporting feature 
           data), and a Bio::EnsEMBL::translation
 Args    : $gene: Bio::EnsEMBL::SeqFeatureI, $contig: Bio::EnsEMBL::DB::ContigI,
           $genetype: string, $count: integer
           $analysis_obj: Bio::EnsEMBL::Analysis


=cut

sub _make_transcript{
  my ($self, $gene, $contig, $genetype, $count, $analysis_obj)=@_;
  $genetype = 'unspecified' unless defined ($genetype);
  $count = 1 unless defined ($count);

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

=head2 validate_gene

 Title   : validate_gene
 Usage   : $self->validate_gene($gene)
 Function: checks start and end coordinates of each exon of each transcript are sane
 Example : 
 Returns : 1 if gene is valid, otherwise zero
 Args    : $gene: Bio::EnsEMBL::Gene


=cut

sub validate_gene{
  my ($self, $gene) = @_;

  # should be only a single transcript
  my @transcripts = $gene->each_Transcript;
  if(scalar(@transcripts) != 1) {
    my $msg = "Rejecting gene - should have one transcript, not " . scalar(@transcripts) . "\n";
    $self->warn($msg);
    return 0;
  }

  foreach my $transcript(@transcripts){
    foreach my $exon($transcript->get_all_Exons){
      if(!$self->validate_exon($exon)){
	my $msg = "Rejecting gene because of invalid exon\n";
	$self->warn($msg);
	return 0;
      }
    }
  }
  
  return 1;
}

=head2 validate_exon

 Title   : validate_exon
 Usage   : $self->validate_exon($exon)
 Function: checks start and end coordinates of exon are sane
 Example : 
 Returns : 1 if exon is valid, otherwise zero
 Args    : $exon: Bio::EnsEMBL::Exon


=cut

sub validate_exon{
  my ($self, $exon) = @_;

  if($exon->start < 0){
    my $msg = "rejecting exon, start < 0 : " . $exon->start . "\n";
    $self->warn($msg);
    return 0;
  }

  elsif($exon->start > $exon->end){
    my $msg = "rejecting exon, start > end : " . $exon->start . " > " . $exon->end . "\n";
    $self->warn($msg);
    return 0;
  }

  elsif($exon->start == $exon->end){
    my $msg = "naughty exon, start == end : " . $exon->start . " == " . $exon->end . " - letting it through\n";
    $self->warn($msg);
    return 1;
  }
  
  return 1;
}

1;
