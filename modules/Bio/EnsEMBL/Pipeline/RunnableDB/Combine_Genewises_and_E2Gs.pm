#
# Ensembl module for Bio::EnsEMBL::Pipeline::RunnableDB::Combine_Genewises_and_E2Gs
#
# Cared for by EnsEMBL  <ensembl-dev@ebi.ac.uk>
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::Combine_Genewises_and_E2Gs

=head1 SYNOPSIS

my $t_e2g = new Bio::EnsEMBL::Pipeline::RunnableDB::Combine_Genewises_and_E2Gs(
                                                                      '-db_obj'      => $dbobj,
                                                                      '-golden_path' => $gp,
                                                                      '-input_id'    => $input_id
                                                                    );

$t_e2g->fetch_input();
$t_e2g->run();
$t_e2g->output();
$t_e2g->write_output(); #writes to DB

=head1 DESCRIPTION

Combines predictions from Genewises with Est2genome predictions from cDNA alignments - ripped straight out of TargettedGeneE2G.

=head1 CONTACT

Ensembl - ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Pipeline::RunnableDB::Combine_Genewises_and_E2Gs;

use vars qw(@ISA);
use strict;
use Storable qw(dclone);
# Object preamble - inheriets from Bio::Root::RootI


use Bio::Root::RootI;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Gene;
use Bio::SeqIO;
use Bio::EnsEMBL::Pipeline::GeneConf qw (
					 GB_GOLDEN_PATH
					 GB_DBHOST
					 GB_DBNAME
					 GB_DBUSER
					 GB_DBPASS
					);
use Bio::EnsEMBL::Pipeline::ESTConf qw (
					 EST_DBHOST
					 EST_DBNAME
					 EST_REFDBUSER
					);



@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  
  my ($path) = $self->_rearrange([qw(GOLDEN_PATH)], @args);
  
  # golden path
  if(!defined $path){
    $path = $GB_GOLDEN_PATH;
  }
  
  # need 2 dbs, one for getting genewises, one for getting e2gs
  my $genedb =  new Bio::EnsEMBL::DBSQL::DBAdaptor(
						   '-host'   => $GB_DBHOST,
						   '-user'   => $GB_DBUSER,
						   '-pass'   => $GB_DBPASS,
						   '-dbname' => $GB_DBNAME,
						  );

  my $cdnadb =  new Bio::EnsEMBL::DBSQL::DBAdaptor(
						   '-host'   => $EST_DBHOST,
						   '-user'   => $EST_REFDBUSER,
						   '-dbname' => $EST_DBNAME,
						   '-dnadb'  => $genedb,
						  ); 
  
  $self->dbobj($genedb);
  $self->cdnadb($cdnadb);

  $path = 'NCBI_28' unless (defined $path && $path ne '');

  $self->dbobj->static_golden_path_type($path);
  $self->cdnadb->static_golden_path_type($path);
  return $self;
}

=head2 cdnadb

    Title   :   cdnadb
    Usage   :   $self->cdnadb($obj);
    Function:   Gets or sets the value of cdnadb. This is a handle to a database 
                containing cdna based gene predictions.
    Returns :   A Bio::EnsEMBL::DBSQL::DBAdaptor compliant object
    Args    :   A Bio::EnsEMBL::DBSQL::DBAdaptor compliant object

=cut

sub cdnadb {
  my( $self, $value ) = @_;
  
  if ($value) 
    {
      $value->isa("Bio::EnsEMBL::DBSQL::DBAdaptor")
	|| $self->throw("Input [$value] isn't a Bio::EnsEMBL::DBSQL::DBAdaptor");
      $self->{'_cdna_db'} = $value;
    }
  return $self->{'_cdna_db'};
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
  
  # input format: chrname.start-end
  if(!($entry =~ /(.*)\.(.*)\-(.*)/)) {
    $self->throw("Not a valid input id... $entry");
  }
  
  $chrname    = $1;
  $start   = $2;
  $end     = $3;
  

  my $sgpa = $self->dbobj->get_StaticGoldenPathAdaptor();
  my $vc = $sgpa->fetch_VirtualContig_by_chr_start_end($chrname,$start,$end);
  $self->vc($vc);
  print STDERR "Chromosome id : $chrname\n";
  print STDERR "Range         : $start - $end\n";
  print STDERR "Contig        : " . $vc->id . " \n";
  
  # now get vc for cdna db 
  my $tmpname = $chrname;
  $tmpname =~ s/chr//;

  $sgpa = $self->cdnadb->get_StaticGoldenPathAdaptor();
  $vc = $sgpa->fetch_VirtualContig_by_chr_start_end($tmpname,$start,$end);
  $self->cdna_vc($vc);
  
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
  
  # get genewise genes
  $self->gw_genes( $self->vc->get_Genes_by_Type('TGE_gw','evidence'));
  $self->gw_genes($self->vc->get_Genes_by_Type('similarity_genewise','evidence'));
  print STDERR "got " . scalar($self->gw_genes) . " genewise genes\n";

  # get e2g genes
  
  my @e2g = $self->cdna_vc->get_Genes_by_Type('exonerate_e2g','evidence');
  print STDERR "got " . scalar(@e2g) . " exonerate_e2g genes\n";
  my @newe2g;

  foreach my $e2g (@e2g) {
    my $found = 0;
    foreach my $tran ($e2g->each_Transcript) {
      my @exons = $tran->get_all_Exons;
      @exons = sort {$a->start <=> $b->start} @exons;
      my $i;

      for ($i = 1; $i <= $#exons; $i++) {
	my $intron = $exons[$i]->start - $exons[$i-1]->end + 1;
	if ($intron > 50000) {
	  $found = 1;
	  
	}
      }
    }
    if ($found == 0) {
      print STDERR "keeping " . $e2g->dbID . "\n";
      push(@newe2g,$e2g);
    }
  }

  $self->e2g_genes(@newe2g);

  print STDERR "got " . scalar($self->e2g_genes) . " sensible e2g genes\n";

  # find which gw matches which e2gs
  my @merged_gw_genes = $self->_merge_gw_genes;

  print STDERR "got " . scalar(@merged_gw_genes) . " merged gw geness\n";  

 GENEWISE:
  foreach my $gw(@merged_gw_genes){
    # should be only 1 transcript
    my @gw_tran  = $gw->each_Transcript;
    my @gw_exons = $gw_tran[0]->get_all_Exons; # ordered array of exons
    my $strand   = $gw_exons[0]->strand;
    
    if($gw_exons[$#gw_exons]->strand != $strand){
      $self->warn("first and last gw exons have different strands - can't make a sensible combined gene\n");
      next GENEWISE;
    }
    
    my @matching_e2gs = $self->match_gw_to_e2g($gw);
    next GENEWISE unless scalar(@matching_e2gs);

    # pick longest gene match for each gw (exon length - though this may not be the best way)
    my $chosen_e2g;
    my $longest = 0;
    foreach my $e2g(@matching_e2gs){
      my $length = 0;
      foreach my $exon ($e2g->get_all_Exons){
	$length += $exon->end - $exon->start + 1; 
      }
      if ($length > $longest){
	$chosen_e2g = $e2g;
	$longest = $length;
      }
    }
    print STDERR "combining : " . $gw->dbID . " with " . $chosen_e2g->dbID . "\n";

    # build combined genes
    $self->combine_genes($gw, $chosen_e2g);
  }
  
  # remap to raw contig coords
  my @remapped = $self->remap_genes();

  $self->output(@remapped);
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

=head2 combine_genes

 Title   : combine_genes
 Usage   :
 Function: 
 Example :
 Returns : 
 Args    :


=cut

sub combine_genes{
  my ($self, $gw, $e2g) = @_;
  
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
  
  # do the actual combining of predictions
  my @newtrans = $self->_make_newtranscript($genetype, $analysis_obj, $gw, $e2g);
  
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

=head2 match_gw_to_e2g

 Title   : match_gw_to_e2g
 Usage   :
 Function: 
 Example :
 Returns : 
 Args    :


=cut

sub match_gw_to_e2g{
  my ($self, $gw) = @_;

  my @matching_e2g;
  my @gw_tran = $gw->each_Transcript;
  my @gw_exons = $gw_tran[0]->get_all_Exons;
  
  my $strand   = $gw_exons[0]->strand;
  
  if($gw_exons[$#gw_exons]->strand != $strand){
    $self->warn("first and last gw exons have different strands - can't make a sensible combined gene\n");
    return;
  }
  
  # in order to match a starting genewise exon with an e2g exon, we need to have 
  # a. exactly coinciding exon ends
  # b. either e2g exon has start <= gw exon start, 
  # OR e2g exon start lies within $exon_slop bp of gw exon start AND the e2g transcript will add extra UTR exons. 
  # previously we had required e2g start to be strictly <= gw start, but this will lose us some valid UTRs
  # substitute "end" for "start" for 3' ends of transcripts
  # BUT don't allow through any e2gs that will result just in a shortened prediction without additional UTR exons.
  my $exon_slop = 20;
  
 E2G:
  foreach my $e2g($self->e2g_genes){
    my @egtran  = $e2g->each_Transcript;
    my @eg_exons = $egtran[0]->get_all_Exons;

    $strand   = $eg_exons[0]->strand;
    
    if($eg_exons[$#eg_exons]->strand != $strand){
#      $self->warn("first and last e2g exons have different strands - skip it\n");
      next E2G;
    }

    my $fiveprime_match = 0;
    my $threeprime_match = 0;
    

    # Lets deal with single exon genes first

    if ($#gw_exons == 0) {
      foreach my $current_exon (@eg_exons) {
	
	if($current_exon->strand != $gw_exons[0]->strand){
	  next E2G;
	}

# don't yet deal with genewise leakage for single exon genes
	if ($gw_exons[0]->end   <= $current_exon->end && 
	    $gw_exons[0]->start >= $current_exon->start){
	  $fiveprime_match = 1;
	  $threeprime_match = 1;
	}
	
      }
      if($fiveprime_match && $threeprime_match){
	push(@matching_e2g, $e2g);
      }
      
      # Now the multi exon genewises
    } else {
      foreach my $current_exon (@eg_exons) {
	if($current_exon->strand != $gw_exons[0]->strand){
	  next E2G;
	}
	
	if($gw_exons[0]->strand == 1){

	FORWARD:
	  if ($gw_exons[0]->end == $current_exon->end && 
	      # either e2g exon starts before genewise exon
	      ($current_exon->start <= $gw_exons[0]->start || 
	       # or e2g exon is a bit shorter but there are spliced UTR exons as well
	       (abs($current_exon->start - $gw_exons[0]->start) <= $exon_slop && 
		$current_exon != $eg_exons[0]))){
#	    print STDERR "fiveprime match\n";
	    $fiveprime_match = 1;
	  }
	  
	  elsif($gw_exons[$#gw_exons]->start == $current_exon->start &&
		# either e2g exon ends after genewise exon
		($current_exon->end >= $gw_exons[$#gw_exons]->end ||
		 # or there are UTR exons to be added
		 (abs($current_exon->end - $gw_exons[$#gw_exons]->end) <= $exon_slop && 
		 $current_exon != $eg_exons[$#eg_exons]))){
#	    print STDERR "threeprime match\n";
	    $threeprime_match = 1;
	  }
	}
	
	elsif($gw_exons[0]->strand == -1){
	REVERSE:
	  if ($gw_exons[0]->start == $current_exon->start &&
	      # either e2g exon ends after gw exon
	      ($current_exon->end >= $gw_exons[0]->end ||
	       # or there are UTR exons to be added
	       (abs($current_exon->end - $gw_exons[0]->end) <= $exon_slop &&
	       $current_exon != $eg_exons[0]))){

	    $fiveprime_match = 1;
	  }
	  elsif ($gw_exons[$#gw_exons]->end == $current_exon->end &&
		 # either e2g exon starts before gw exon
		 ($current_exon->start <= $gw_exons[$#gw_exons]->start ||
		  # or there are UTR exons to be added
		  (abs($current_exon->start - $gw_exons[$#gw_exons]->start) <= $exon_slop &&
		  $current_exon != $eg_exons[$#eg_exons]))){

	    $threeprime_match = 1;
	  }
	}
      }
      if($fiveprime_match && $threeprime_match){
	push(@matching_e2g, $e2g);
      }
    }
  }
  return @matching_e2g;
  
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
    $gene->dbID($gwg->dbID);
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
	foreach my $suppfeat($exon->each_Supporting_Feature){
	  $previous_exon->add_Supporting_Feature($suppfeat);
	}
	next EXON;
      }
      
      else{
	# make a new Exon - clone $exon
	my $cloned_exon = new Bio::EnsEMBL::Exon;
	$cloned_exon->start($exon->start);
	$cloned_exon->end($exon->end);
	$cloned_exon->strand($exon->strand);
	$cloned_exon->phase($exon->phase);
	$cloned_exon->end_phase($exon->end_phase);
	$cloned_exon->contig_id($exon->contig_id);
	
	$cloned_exon->attach_seq($self->vc->primary_seq);
	$cloned_exon->add_sub_SeqFeature($exon,'');

	foreach my $suppfeat($exon->each_Supporting_Feature){
	  $cloned_exon->add_Supporting_Feature($suppfeat);
	}
	push(@pred_exons, $cloned_exon);
      }
      
    }

    # transcript
    my $merged_transcript   = new Bio::EnsEMBL::Transcript;
    foreach my $pe(@pred_exons){
      $merged_transcript->add_Exon($pe);
    }
    
    $merged_transcript->sort;

    my $cloned_translation = new Bio::EnsEMBL::Translation;
    $cloned_translation->start($trans[0]->translation->start);
    $cloned_translation->start_exon($trans[0]->translation->start_exon);
    $cloned_translation->end($trans[0]->translation->end);
    $cloned_translation->end_exon($trans[0]->translation->end_exon);

    $merged_transcript->translation($cloned_translation);
    
    # and gene
    $gene->add_Transcript($merged_transcript);
    push(@merged, $gene);
    $count++;
  }

  return @merged;
}

=head2 _make_newtranscript

 Title   : _make_newtranscript
 Usage   :
 Function: makes new transcript by combining the genewise and est2genome predictions. Its a monster.
 Example :
 Returns : 
 Args    :


=cut

sub _make_newtranscript {
  my ($self, $genetype, $analysis_obj, $gw, $e2g) = @_;
  my @combined_transcripts  = ();
  
  # should be only 1 transcript
  my @gw_tran  = $gw->each_Transcript;
#  $gw_tran[0]->sort;
  my @gw_exons = $gw_tran[0]->get_all_Exons; # ordered array of exons
  my @egtran = $e2g->each_Transcript;
#  $egtran[0]->sort;
  my @e2g_exons  = $egtran[0]->get_all_Exons; # ordered array of exons
  
  # OK, let's see if we need a new gene
  # base it on the existing genewise one
  my $newtranscript = new Bio::EnsEMBL::Transcript;
  foreach my $exon(@gw_exons){
    $newtranscript->add_Exon($exon);
  }

  my $translation   = new Bio::EnsEMBL::Translation;
  $translation->start($gw_tran[0]->translation->start);
  $translation->end($gw_tran[0]->translation->end);
  $translation->start_exon($gw_tran[0]->translation->start_exon);
  $translation->end_exon($gw_tran[0]->translation->end_exon);

 $newtranscript->translation($translation);
  my $eecount = 0;
  
  $newtranscript->translation->start_exon($newtranscript->start_exon);
 $newtranscript->translation->end_exon($newtranscript->end_exon);

  # check strands are consistent
  foreach my $ee(@e2g_exons){
    if ($ee->strand != $gw_exons[0]->strand){
      $self->warn("gw and e2g exons have different strands - can't combine genes\n") ;
      next GENEWISE;
    }
    
    # single exon genewise prediction?
    if(scalar(@gw_exons) == 1) {
      print STDERR "Making new transcript from single exon genewise\n";

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
								  $gw, 
								  $e2g)
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
      $ex->attach_seq($self->vc);
      $ex->contig_id($self->vc->id);
      # add new analysis object to the supporting features
      foreach my $sf($ex->each_Supporting_Feature){
	$sf->analysis($analysis_obj);
	$sf->source_tag($genetype);
      }
    }

    # double check sane translation start & end
    if($newtranscript->translation->start < 1){
      print STDERR "dodgy translation start - defaulting to 1\n";
      $newtranscript->translation->start(1);
    }

    if($newtranscript->translation->end < 1 || 
       $newtranscript->translation->end > $newtranscript->translation->end_exon->length ){
      print STDERR "dodgy translation end " . $newtranscript->translation->end . " - defaulting to end_exon length " . $newtranscript->translation->end_exon->length. "\n";
      $newtranscript->translation->end($newtranscript->translation->end_exon->length);
    }
    
    # check translation is the same as for the genewise gene we're built from
    # not the merged gene you idiot!
    my $foundtrans = 0;
    GWG:
	foreach my $gwg($self->gw_genes) {
	  if($gw->dbID == $gwg->dbID){
	    print STDERR "comparing transcripts\n";
	    my @tran = $gwg->each_Transcript;
	    $foundtrans = $self->compare_transcripts($gwg, $newtranscript);
	    
	    if ($foundtrans == 1){
	      push (@combined_transcripts, $newtranscript); 	
	      last GWG;
	    }
	  }
	}

 
    if(!$foundtrans){
      $self->warn("UTR prediction is not the same as genewise prediction - discarding it\n");
    }
    
  }
  
  
  return @combined_transcripts;
  
}

sub compare_transcripts{
  my ($self, $genewise_gene, $combined_transcript) = @_;
  my @genewise_transcripts = $genewise_gene->each_Transcript;

  my $seqout = new Bio::SeqIO->new(-fh => \*STDERR);

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
  }
  else{
#    print STDERR "genewise: \n";             
#    $seqout->write_seq($genewise_translation);
  }

  $@ = '';

  eval{
    $combined_translation = $combined_transcript->translate;
  };


  if ($@) {
    print STDERR "Couldn't translate combined gene:[$@]\n";
    return 0;
  }
  else{
#    print STDERR "combined: \n";             
#    $seqout->write_seq($combined_translation);
  }	 

  my $gwseq  = $genewise_translation->seq;
  my $comseq = $combined_translation->seq;
 
  if($gwseq eq $comseq) {
    print STDERR "combined translation is identical to genewise translation\n";
    return 1;
  }
  elsif($gwseq =~ /$comseq/){
    print STDERR "combined translation is a truncated version of genewise translation\n";
    return 1;
  }

  return 0;
}

sub transcript_from_single_exon_genewise {
  my ($self, $eg_exon, $gw_exon, $transcript, $translation, $exoncount, @e2g_exons) = @_;

  # save out current translation end - we will need this if we have to unmerge frameshifted exons later
  my $orig_tend = $translation->end;

  # stay with being strict about gw vs e2g coords - may change this later ...
  if ($gw_exon->start >= $eg_exon->start && $gw_exon->end <= $eg_exon->end){

    my $egstart = $eg_exon->start;
    my $egend   = $eg_exon->end;
    my $gwstart = $gw_exon->start;
    my $gwend   = $gw_exon->end;

#    print STDERR "single exon gene, " . $gw_exon->strand  .  " strand\n";	    

    # modify the coordinates of the first exon in $newtranscript
    
    my $ex = $transcript->start_exon;

    $ex->start($eg_exon->start);
    $ex->end($eg_exon->end);

    # need to explicitly set the translation start & end exons here.
    $translation->start_exon($ex);

    # end_exon may be adjusted by 3' coding exon frameshift expansion. Ouch.
    $translation->end_exon($ex);
    
    # need to deal with translation start and end this time - varies depending on strand
    if($gw_exon->strand == 1){
      my $diff = $gwstart - $egstart;
      my $tstart = $translation->start;
      my $tend = $translation->end;
	      
      $translation->start($tstart + $diff);
      $translation->end($tend + $diff);

      $self->throw("setting very dodgy translation start: " . $translation->start.  "\n") unless $translation->start > 0;
      $self->throw("setting dodgy translation end: " . $translation->end . " exon_length: " . $translation->end_exon->length . "\n") unless $translation->end <= $translation->end_exon->length;
    }

    elsif($gw_exon->strand == -1){
      my $diff = $egend - $gwend;
      my $tstart = $translation->start;
      my $tend = $translation->end;
      $translation->start($tstart+$diff);
      $translation->end($tend + $diff);


      $self->throw("setting very dodgy translation start: " . $translation->start.  "\n") unless $translation->start > 0;
      $self->throw("setting dodgy translation end: " . $translation->end . " exon_length: " . $translation->end_exon->length . "\n") unless $translation->end <= $translation->end_exon->length;
    }
    
    
    # expand frameshifted single exon genewises back from one exon to multiple exons
    if(scalar($ex->sub_SeqFeature) > 1){
#      print STDERR "uh-oh frameshift\n";
      my @sf = $ex->sub_SeqFeature;
      
      # save current start and end of modified exon
      my $cstart = $ex->start;
      my $cend   = $ex->end;
      my $exlength = $ex->length;
      
      # get first exon - this has same id as $ex
      my $first = shift(@sf);
      $ex->end($first->end); # NB end has changed!!!
      
      # get last exon
      my $last = pop(@sf);
      $last->end($cend);
      $transcript->add_Exon($last);
      # and adjust translation end - the end is still relative to the merged gw exon
      $translation->end_exon($last);
      $translation->end($orig_tend);

      # get any remaining exons
      foreach my $s(@sf){
	$transcript->add_Exon($s);
	$transcript->sort;
      }
      # flush the sub_SeqFeatures
      $ex->flush_sub_SeqFeature;
    }

    # need to add back exons, both 5' and 3'
    $self->add_5prime_exons(\$transcript, $exoncount, @e2g_exons);
    $self->add_3prime_exons(\$transcript, $exoncount, @e2g_exons);
	      
  }
  return $transcript;
}


sub transcript_from_multi_exon_genewise {
  my ($self, $current_exon, $transcript, $translation, $exoncount, $gw_gene, $eg_gene) = @_;
  # save out current translation->end - we'll need it if we have to expand 3prime exon later
  my $orig_tend = $translation->end;

  my @gwtran  = $gw_gene->each_Transcript;
  my @gwexons = $gwtran[0]->get_all_Exons;
  
  my @egtran  = $eg_gene->each_Transcript;
  my @egexons = $egtran[0]->get_all_Exons;

  # in order to match a starting genewise exon with an e2g exon, we need to have 
  # a. exactly coinciding exon ends
  # b. exon starts lying within $exon_slop bp of each other. 
  # previously we had required e2g start to be strictly <= gw start, but this will lose us some valid UTRs
  # substitute "end" for "start" for 3' ends of transcripts
  my $exon_slop = 20;

  # compare to the first genewise exon
  if($gwexons[0]->strand == 1){
  FORWARD:
    if ($gwexons[0]->end == $current_exon->end && 
	# either e2g exon starts before genewise exon
	($current_exon->start <= $gwexons[0]->start || 
	 # or e2g exon is a bit shorter but there are spliced UTR exons as well
	 (abs($current_exon->start - $gwexons[0]->start) <= $exon_slop && 
	  $current_exon != $egexons[0]))){

      my $current_start = $current_exon->start;
      my $gwstart = $gwexons[0]->start;
      
      # modify the coordinates of the first exon in $newtranscript
      my $ex = $transcript->start_exon;
      $ex->start($current_exon->start);

      # add all the exons from the est2genome transcript, previous to this one
      $self->add_5prime_exons(\$transcript, $exoncount, @egexons);
      
      # fix translation start 

      if($gwstart >= $current_start){
	# take what it was for the gw gene, and add on the extra
	my $tstart = $translation->start;
	$tstart += ($gwstart - $current_start);
	$translation->start($tstart);

      }
      else{
	# genewise has leaked over the start. Tougher call - we need to take into account the 
	# frame here as well
	my $diff = $current_start - $gwstart;
	my $tstart = $translation->start;

	if    ($diff % 3 == 0) { $translation->start(1); }
	elsif ($diff % 3 == 1) { $translation->start(3); }
	elsif ($diff % 3 == 2) { $translation->start(2); }
	else {
	  $translation->start(1);
	  $self->warn("very odd - $diff mod 3 = " . $diff % 3 . "\n");}
      }
      $self->throw("setting very dodgy translation start: " . $translation->start.  "\n") unless $translation->start > 0;

    } # end 5' exon
    
    elsif ($gwexons[$#gwexons]->start == $current_exon->start && 
	# either e2g exon ends after genewise exon
	   ($current_exon->end >= $gwexons[$#gwexons]->end ||
	    # or there are UTR exons to be added
	    (abs($current_exon->end - $gwexons[$#gwexons]->end) <= $exon_slop && 
	     $current_exon != $egexons[$#egexons]))){   
      
      #      print STDERR "3' exon match\n";
      
      # modify the coordinates of the last exon in $newtranscript
      my $ex = $transcript->end_exon;
      $ex->end($current_exon->end);

      my $expanded = $self->expand_3prime_exon(\$ex, \$transcript);

      # need to explicitly set the translation end exon for translation to work out
      my $end_ex = $transcript->end_exon;
      $translation->end_exon($end_ex);

      if($expanded){
	# set translation end to what it originally was in the unmerged genewise gene
#	print STDERR "setting translation end to $orig_tend\n";
	$translation->end($orig_tend);
      }

      # fix translation end iff genewise has leaked over - will need truncating
      if($current_exon->end < $gwexons[$#gwexons]->end){
#	print STDERR "FORWARD exon length: " . $current_exon->length . "\n";
	$translation->end($current_exon->length);
      }

      
      # finally add any 3 prime e2g exons
      $self->add_3prime_exons(\$transcript, $exoncount, @egexons);

    } # end 3' exon
    
  }
  
  elsif($gwexons[0]->strand == -1){
  REVERSE:
    if ($gwexons[0]->start == $current_exon->start && 
	# either e2g exon ends after gw exon
	($current_exon->end >= $gwexons[0]->end ||
	 # or there are UTR exons to be added
	 (abs($current_exon->end - $gwexons[0]->end) <= $exon_slop &&
	  $current_exon != $egexons[0]))){
      
      # sort out translation start
      if($current_exon->end >= $gwexons[0]->end){
	# take what it was for the gw gene, and add on the extra
	my $tstart = $translation->start;
	$tstart += $current_exon->end - $gwexons[0]->end;
	$translation->start($tstart);

      }
      else{
	# genewise has leaked over the start. Tougher call - we need to take into account the 
	# frame here as well
	my $diff = $gwexons[0]->end - $current_exon->end;
	my $gwstart = $gwexons[0]->end;
	my $current_start = $current_exon->end;
	my $tstart = $translation->start;

	if    ($diff % 3 == 0) { $translation->start(1); }
	elsif ($diff % 3 == 1) { $translation->start(3); }
	elsif ($diff % 3 == 2) { $translation->start(2); }
	else {
	  $translation->start(1);
	  $self->warn("very odd - $diff mod 3 = " . $diff % 3 . "\n");}
      }
      $self->throw("setting very dodgy translation start: " . $translation->start.  "\n") unless $translation->start > 0;

      # modify the coordinates of the first exon in $newtranscript
      my $ex = $transcript->start_exon;
      $ex->end($current_exon->end);
     
      # need to explicitly set the translation start exon for translation to work out
      $translation->start_exon($ex);

      $self->add_5prime_exons(\$transcript, $exoncount, @egexons);
      

    } # end 5' exon
    
    elsif ($gwexons[$#gwexons]->end == $current_exon->end && 
	   # either e2g exon starts before gw exon
	   ($current_exon->start <= $gwexons[$#gwexons]->start ||
	    # or there are UTR exons to be added
	    (abs($current_exon->start - $gwexons[$#gwexons]->start) <= $exon_slop &&
	     $current_exon != $egexons[$#egexons]))){
      
      #      print STDERR "3' exon match\n";
  
      # modify the coordinates of the last exon in $newtranscript
      my $ex = $transcript->end_exon;

      $ex->start($current_exon->start);

      my $expanded = $self->expand_3prime_exon(\$ex, \$transcript);

      # need to explicitly set the translation end exon for translation to work out
      my $end_ex = $transcript->end_exon;
      $translation->end_exon($end_ex);

      if($expanded){
	# set translation end to what it originally was in the unmerged genewise gene
#	print STDERR "setting translation end to $orig_tend\n";
	$translation->end($orig_tend);
      }

      # adjust translation end iff genewise has leaked
      if($current_exon->start > $gwexons[$#gwexons]->start){
#	print STDERR "REVERSE exon length: " . $current_exon->length . "\n";
	$translation->end($current_exon->length);
      }

      $self->add_3prime_exons(\$transcript, $exoncount, @egexons);
      
    } # end 3' exon
  }  

 return $transcript; 
}

sub add_5prime_exons{
my ($self, $transcript, $exoncount, @e2g_exons) = @_;

      # add all the exons from the est2genome transcript, previous to this one
      # db handle will be screwed up, need to make new exons from these
      my $c = 0;
      while($c < $exoncount){
	my $newexon = new Bio::EnsEMBL::Exon;
	my $oldexon = $e2g_exons[$c];
	$newexon->start($oldexon->start);
	$newexon->end($oldexon->end);
	$newexon->strand($oldexon->strand);
	$newexon->phase($oldexon->phase);
	$newexon->contig_id($oldexon->contig_id);
	$newexon->attach_seq($self->vc);
	foreach my $sf($oldexon->each_Supporting_Feature){
	  $newexon->add_Supporting_Feature($sf);
	}
#	print STDERR "Adding 5prime exon " . $newexon->start . " " . $newexon->end . "\n";
	$$transcript->add_Exon($newexon);
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
	return 1;
      }
      
  # no expansion
  return 0;
}

# $exoncount tells us which position in the array 
# of e2g exons corresponds to the end of the genewise transcript so we add back 
# exons 3' to that position.
# $exon and $transcript are references to Exon and Transcript objects.
sub add_3prime_exons {
my ($self, $transcript, $exoncount, @e2g_exons) = @_;
# need to deal with frameshifts - 3' exon is a special case as its end might have changed

      # add all the exons from the est2genome transcript, subsequent to this one
      my $c = $#e2g_exons;
      while($c > $exoncount){
	my $newexon = new Bio::EnsEMBL::Exon;
	my $oldexon = $e2g_exons[$c];
	$newexon->start($oldexon->start);
	$newexon->end($oldexon->end);
	$newexon->strand($oldexon->strand);
	$newexon->phase($oldexon->phase);
	$newexon->contig_id($oldexon->contig_id);
	$newexon->attach_seq($self->vc);
	foreach my $sf($oldexon->each_Supporting_Feature){
	  $newexon->add_Supporting_Feature($sf);
	}
#	print STDERR "Adding 3prime exon " . $newexon->start . " " . $newexon->end . "\n";
	$$transcript->add_Exon($newexon);
	$$transcript->sort;
	$c--;
      }

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
  my $contig = $self->vc;

  my @genes = $self->combined_genes;

GENE:  foreach my $gene (@genes) {
    my @t = $gene->each_Transcript;
    my $tran = $t[0];

    # check that it translates - not the est2genome genes
    if($gene->type eq 'TGE_gw' || $gene->type eq 'combined_gw_e2g'){
      
      my $translates = $self->check_translation($tran);
      next GENE unless $translates;
  }

    eval {
      my $genetype = $gene->type;
      my $newgene = $contig->convert_Gene_to_raw_contig($gene);
      # need to explicitly add back genetype and analysis.
      $newgene->type($genetype);
      $newgene->analysis($gene->analysis);

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
	
	foreach my $exon($tran->get_all_Exons) {
	  if ($exon eq $tran->translation->start_exon()) {
	    $fpc_strand = $exon->strand;
	    last;
	  }
	}
	
	foreach my $tran ($newgene->each_Transcript) {
	  foreach my $exon($tran->get_all_Exons) {
	    
	    # oh dear oh dear oh dear
	    # this is still giving some problems
	    if ($exon eq $tran->translation->start_exon()) {
	      if($fpc_strand == 1 && $exon->strand == -1){
		print STDERR "fpc strand 1, raw strand -1 - flipping translation start/end\n";
		$exon->end($exon->end - ($tran->translation->start -1));
		
		$exon->start($exon->end - ($tran->translation->end -1));
	      }
	    }
	  }
	}
	
      } # end special case single coding exon
      
      # final exon coord sanity check
      foreach my $exon($newgene->get_all_Exons){
	# make sure we deal with stickies!
	if($exon->isa("Bio::EnsEMBL::StickyExon")){
	  foreach my $ce($exon->each_component_Exon){
	    # exon start and end must both be within the raw contig!!!
	    if($ce->start < 1){
	      $self->throw("can't set exon->start < 1 (" . $ce->start . ") - discarding gene\n");
	    }
	    
	    if($ce->end > $ce->contig->primary_seq->length){
	      $self->throw("exon extends beyond end of contig - discarding gene\n");
	    }
	  }
	}
	else{
	  # regular exon
	  # exon start and end must both be within the raw contig!!!
	  if($exon->start < 1){
	    $self->throw("can't set exon->start < 1 (" . $exon->start . ") - discarding gene\n");
	  }
	  
	  if($exon->end > $exon->contig->primary_seq->length){
	    $self->throw("exon extends beyond end of contig - discarding gene\n");
	  }
	}
      }
      # if we get to here, the gene is fine, so push it onto the array to be returned
      push(@newf,$newgene);

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
 Function: checks that transcript translates
 Example :
 Returns : 1 if transcript translates with no stops; otherwise 0
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
    $self->warn("discarding gene because of stops\n");
    return 0;
  }

  # if we get to here, all is well
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

=head2 cdna_vc

 Title   : cdna_vc
 Usage   : $obj->cdna_vc($newval)
 Function: 
 Returns : value of cdna_vc
 Args    : newvalue (optional)


=cut

sub cdna_vc {
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_cdna_vc'} = $value;
    }
    return $obj->{'_cdna_vc'};

}

1;
