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
use Bio::EnsEMBL::Pipeline::GeneConf qw (EXON_ID_SUBSCRIPT
					 TRANSCRIPT_ID_SUBSCRIPT
					 GENE_ID_SUBSCRIPT
					 PROTEIN_ID_SUBSCRIPT
					 );
use Storable qw(dclone);
# Object preamble - inheriets from Bio::Root::RootI


use Bio::Root::RootI;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;
use Bio::EnsEMBL::Pipeline::Runnable::ExonerateMiniEst2Genome;
use Bio::SeqIO;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  if(!defined $self->seqfetcher) {
    # will look for pfetch in $PATH - change this once PipeConf up to date
    my $seqfetcher = new Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch; 
    $self->seqfetcher($seqfetcher);
  }

  my ($path) = $self->_rearrange([qw(GOLDEN_PATH)], @args);
  $path = 'UCSC' unless (defined $path && $path ne '');
  $self->dbobj->static_golden_path_type($path);

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
  my ($self,@args) = @_;

  my $entry = $self->input_id;
  my $fpc;
  my $start;
  my $end;
  my $protein_id; 
  my $cdna_id;

  # input format: ctg1234:10602496,10603128:Q9UGV6:AC00012
  # or ctg1234:10602496,10603128:Q9UGV6 if no cDNA
  if( !(($entry =~ /(\S+):(\d+),(\d+):(\S+):(\S+)/) || ($entry =~ /(\S+):(\d+),(\d+):(\S+):/))) {
      $self->throw("Not a valid input id... $entry");
  }
  
  print STDERR "input: ".$entry . "\n";
  
  $fpc = $1;
  $protein_id = $4;
  $cdna_id = $5;

  $start   = $2;
  $end     = $3;

  if ($2 > $3) { # let blast sort it out
      $start  = $3;
      $end    = $2;
  }
  
  my $sgpa = $self->dbobj->get_StaticGoldenPathAdaptor();

  print STDERR "$fpc $start $end\n";
  my ($chrname,$chrstart,$chrend) = $sgpa->convert_fpc_to_chromosome($fpc,$start-10000,$end+10000);
  print STDERR "$chrname $chrstart $chrend\n";
  my $vc = $sgpa->fetch_VirtualContig_by_chr_start_end($chrname,$chrstart,$chrend);

#  my $vc = $sgpa->fetch_VirtualContig_by_chr_start_end($chrname,$start-10000,$end+10000);
  
  $self->vc($vc);
  
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
    $cdna = $self->seqfetcher->get_Seq_by_acc($cdna_id);
  };
  if($@) {
    $self->throw("problem fetching [$cdna_id]: [$@]\n");
  }
  
  $self->{'_tmpfile'} = "/tmp/tge2g_" . $$ . ".fa";
  my $cdnafile = $self->{'_tmpfile'};

  my $seqout = new Bio::SeqIO('-file' => ">$cdnafile" , '-format' => 'Fasta');
  $seqout->write_seq($cdna);

  # repmasking?
  my $e2g = new Bio::EnsEMBL::Pipeline::Runnable::ExonerateMiniEst2Genome('-genomic'    => $vc->primary_seq, 
									  '-queryseq'   => $cdnafile,
									  '-seqfetcher' => $self->seqfetcher);

  $self->e2g_runnable($e2g);
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

   if (defined $self->e2g_runnable) {
     $self->e2g_runnable->run;
     $self->convert_e2g_output;
     $self->combine_genes;

     # clean up tmpfile
     my $tmpfile = $self->{'_tmpfile'};
     unlink $tmpfile;
   }

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

#    $self->throw("exiting bfore write");

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
	  
#	  print STDERR "genetype: " . $gene->type . "\n";
	  
	  
	  foreach my $trans ( $gene->each_Transcript ) {
	    $transcount++;
	    $translcount++;
	  }

	  foreach my $exon ( $gene->each_unique_Exon() ) {
	    $exoncount++;
	    }
	}
	
#	$self->throw("exiting bfore write");
	
	# get that number of ids. This locks the database
	
	my @geneids  =  $gene_obj->get_New_external_id('gene',$GENE_ID_SUBSCRIPT,$genecount);
	my @transids =  $gene_obj->get_New_external_id('transcript',$TRANSCRIPT_ID_SUBSCRIPT,$transcount);
	my @translids=  $gene_obj->get_New_external_id('translation',$PROTEIN_ID_SUBSCRIPT,$translcount);
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
	if( scalar(@geneids)  != 0 || scalar(@exonsid)   != 0 || 
	    scalar(@transids) != 0 || scalar(@translids) != 0 ) {
	    $self->throw("In id assignment, left with unassigned ids ".
			 scalar(@geneids)  . " " .
			 scalar(@transids) . " " .
			 scalar(@translids)." " .
			 scalar(@exonsid));
	}

    };
    if( $@ ) {
	$self->throw("Exception in getting new ids. Exiting befor write\n\n$@" );
    }


    # this now assummes that we are building on a single VC.

#    $self->throw("Bailing before real write\n");
    
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
  my $time  = time; chomp($time);
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
  
  print STDERR "e2g genes: " . scalar(@genes) . "\n";

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

  # check for stops?
  print STDERR "gw genes: " . scalar(@genes) . "\n";
  $self->gw_genes(@genes);
  
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

  # merge the genewise and est2genome predictions. Only one new transcript per genewise gene ...
  my @newtrans = $self->_make_newtranscripts($genetype, $analysis_obj, @merged_gw_genes);

  $analysis_obj->gff_feature('gene');

  # make some lovely genes
  my @genes;
  my $count=0;
  my $time = time;
  chomp($time);

  foreach my $trans(@newtrans){
    $trans->sort;
    my $gene = new Bio::EnsEMBL::Gene;
    $gene->type($genetype);
    $gene->id($self->input_id . "$genetype.$count"); 
    $gene->version(1);
    $gene->add_Transcript($trans);
    $gene->analysis($analysis_obj);
    $gene->created($time);
    $gene->modified($time);
    push (@genes,$gene);
    $count++;
  }

  $self->combined_genes(@genes);
  
}

=head2 _merge_gw_genes

 Title   : _merge_gw_genes
 Usage   :
 Function: merges adjacent exons if they are frameshifted; stores component exons
 Example :
 Returns : 
 Args    :


=cut

sub _merge_gw_genes {
  my ($self) = @_;

  my @merged;
  my $count = 1;
  my $contig = $self->vc;
  foreach my $gwg($self->gw_genes){
    my $gene = new Bio::EnsEMBL::Gene;
    $gene->type('combined');
    $gene->id($gwg->id);
    $gene->version(1);
    
    my @pred_exons;
    my $ecount = 0;
    
    # order is crucial
    my @trans = $gwg->each_Transcript;
    if(scalar(@trans) != 1) { $self->throw("expected one transcript for $gwg\n"); }
    
  EXON:      
    foreach my $exon($trans[0]->each_Exon){
      my $prev_exon;
      
      if ($ecount && $pred_exons[$ecount-1]){
	$prev_exon = $pred_exons[$ecount-1];
      }
      
      $ecount++;
      
      # frameshift? we treat two exons separated by max 10 bases as a single exon
      if( defined($prev_exon) && abs($exon->start - $prev_exon->end) <= 10 ){
	# combine the two
	$prev_exon->end($exon->end);
	$prev_exon->add_sub_SeqFeature($exon,'');
	next EXON;
      }
      
      else{
	# make a new Exon - clone $exon
	my $ne = dclone($exon);
	$ne->attach_seq($self->vc->primary_seq);
	$ne->add_sub_SeqFeature($exon,'');
	push(@pred_exons, $ne);
      }
      
    }

    # transcript
    my $transcript   = new Bio::EnsEMBL::Transcript;
    $transcript->id($contig->id . ".combined.$count");
    $transcript->version(1);
    foreach my $pe(@pred_exons){
      $transcript->add_Exon($pe);
    }
    
    my $gw_translation = dclone($trans[0]->translation);
    $transcript->translation($gw_translation);
    
    # and gene
    $gene->add_Transcript($transcript);
    push(@merged, $gene);
    $count++;

    # phase check??
    print STDERR "Phase check\n";
    $transcript->sort;
    $trans[0]->sort;
    print STDERR "new transcript\n";
    foreach my $e($transcript->each_Exon) {
      foreach my $sf($e->sub_SeqFeature){
	print STDERR $sf->start . " - " . $sf->end . " " .$sf->phase . "\t" . $sf->end_phase . "\n";
      }
    }
    print STDERR "original transcript\n";
    foreach my $e($trans[0]->each_Exon) {
      print STDERR $e->start . " - " . $e->end . " " . $e->phase . "\t" . $e->end_phase . "\n";
    }

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
  my @gw_genes  = $self->gw_genes;
  my @e2g_genes = $self->e2g_genes;
  my @newtrans  = ();



 GENE:
  foreach my $gene(@merged_gw_genes) {
    print "\nGENE\n";
    my @gw_tran = $gene->each_Transcript;
    my @gw_ex = $gw_tran[0]->each_Exon; # need ordered array
    my $strand = $gw_ex[0]->strand;
    if($gw_ex[$#gw_ex]->strand != $strand){
      $self->warn("first and last gw exons have different strands - can't make a sensible combined gene\n");
      next GENE;
    }
    
    my $foundtrans = 0;  
    if(scalar(@gw_ex) == 1){
      my $covered = $self->_check_coverage($gene);
      print STDERR "Single exon\n";
      next GENE unless $covered;
    }

 E2G:
    foreach my $eg(@e2g_genes){
      # yuk yuk yuk only want 1 prediction per cDNA - SP:refseq duplicates
      next GENE if $foundtrans;
      my @egtran = $eg->each_Transcript;
      my @eg_ex = $egtran[0]->each_Exon; # need ordered array again
      
      print STDERR "comparing " . $gene->id . " with " . $eg->id . "\n";
      # OK, let's see if we need a new gene
      # base it on the existing genewise one
      my $newtranscript = dclone($gw_tran[0]);
      my $tmpid = $newtranscript->id;
      $tmpid .= "_combined";
      $newtranscript->id($tmpid);
      my $translation  = dclone($gw_tran[0]->translation);
      $tmpid = $translation->id;
      $tmpid .= "_combined";
      $translation->id($tmpid);
      $newtranscript->translation($translation);
      my $eecount = 0;

      print "e2g exons: " . scalar(@eg_ex) . "\n";

      foreach my $ee(@eg_ex){
	if ($ee->strand != $strand){
          $self->warn("gw and e2g exons have different strands - can't combine genes\n") ;
          next GENE;
        }

	# single exon genewise prediction?
	if(scalar(@gw_ex) == 1) {# eeeep

	  $newtranscript = $self->transcript_from_single_exon_genewise( $ee, $gw_ex[0], $newtranscript, $translation, $eecount, @eg_ex);
	  
	}
	
	# multiple exon genewise prediction
	else {#multi
	  $newtranscript = $self->transcript_from_multi_exon_genewise($ee, $newtranscript, 
								      $translation, $eecount, $gene, $eg)
	} # end multi exon predictions
  
	# increment the exon
	$eecount++;
	
      } # end foreach my $ee
      
      # check the transcript and expand frameshifts in all but original 3' gw_exon
      # check translation does not change from original genewise prediction
      if (defined($newtranscript)){
	foreach my $ex($newtranscript->each_Exon){
	  if(scalar($ex->sub_SeqFeature) > 1 ){
	    print STDERR "frameshift: " . $ex->id . "\n";
	    my @sf = $ex->sub_SeqFeature;
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
	foreach my $ex($newtranscript->each_Exon){
	  $ex->attach_seq($self->vc);
	  $ex->contig_id($self->vc->id);
	  # add new analysis object to the supporting features
	  foreach my $sf($ex->each_Supporting_Feature){
	    $sf->analysis($analysis_obj);
	    $sf->source_tag($genetype);
	  }
	}
	
	# compare UTR-gene translation with genewise translation NB NOT the merged gw gene 
	# we have been using for building!!!
      GWG:	foreach my $gwg($self->gw_genes) {
	  if ($gwg->id eq $gene->id){
	    my @gwgtran = $gwg->each_Transcript;
	    if(scalar(@gwgtran != 1)) {
	      $self->warn("Panic! Got " . scalar(@gwgtran) . " transcripts from " . $gwg->id . "\n");
	      next GWG;
	    }
	    
	    my $genewise;
	    my $combined;
	    
	    eval {
	      $genewise = $gwgtran[0]->translate;
	      $combined = $newtranscript->translate;
	    };
	    
	    if ($@) {
	      print STDERR "Couldn't translate: " . $gene->id . " plus " . $eg->id  . "[$@]\n";
	    }
	    
	    print STDERR "translation: \n";
	    my $seqio = Bio::SeqIO->new(-fh => \*STDERR);
	    print STDERR "genewise: \n";	      
	    $seqio->write_seq($genewise);	      
	    print STDERR "combined: \n";
	    $seqio->write_seq($combined); 
	    print STDERR "\n ";
	    
	    if($genewise->seq eq $combined->seq) {
	      $foundtrans = 1;
	      push (@newtrans, $newtranscript); 	
	    }
	    else {
	      $self->warn("UTR prediction is not the same as genewise prediction - discarding it\n");
	    }
	  }
	}
      }
    }
  }
  return @newtrans;
  
}

sub transcript_from_single_exon_genewise {
  my ($self, $eg_exon, $gw_exon, $transcript, $translation, $exoncount, @e2g_exons) = @_;
  if ($gw_exon->start >= $eg_exon->start && $gw_exon->end <= $eg_exon->end){
    print STDERR "single exon gene, " . $gw_exon->strand  .  " strand\n";	    
    # modify the coordinates of the first exon in $newtranscript
    my $ex = $transcript->start_exon;
    
    $ex->start($eg_exon->start);
    $ex->end($eg_exon->end);
    
    #	    print STDERR "eecount: $eecount\n";
    
    # need to add back exons, both 5' and 3'
    my $c = 0;
    while($c < $exoncount){
      print STDERR "adding 5' exon\n";
      $transcript->add_Exon($e2g_exons[$c]);
      $transcript->sort;
      $c++;
    }
	    
    # add all the exons from the est2genome transcript, subsequent to this one
    $c = $#e2g_exons;
    while($c > $exoncount){
      print STDERR "adding 3' exon\n";
      $transcript->add_Exon($e2g_exons[$c]);
      $transcript->sort;
      $c--;
    }
	    
    # need to deal with translation start and end this time - varies depending on strand
    if($gw_exon->strand == 1){
      my $diff = $gw_exon->start - $ex->start;
      my $tstart = $translation->start;
      my $tend = $translation->end;
	      
      #	      print STDERR "***gw  " . $gw_ex[0]->start . " " . $gw_ex[0]->end . "\n";
      $translation->start($tstart + $diff);
      $translation->end($tend + $diff);
    }

    elsif($gw_exon->strand == -1){
      #	      print STDERR "***reverse\n";
      #	    my $diff = $e2g_exon->end - $gw_ex[0]->end;
#      my $diff = $gw_exon->start - $eg_exon->start;
      my $diff = $ex->end - $gw_exon->end;
      print STDERR "diff $diff\n";
      my $tstart = $translation->start;
      my $tend = $translation->end;
      print STDERR "***gw  " . $gw_exon->start . " " . $gw_exon->end . "\n";
      print STDERR "***e2g  " . $ex->start . " " . $ex->end . "\n";
      print STDERR "*** translation - before  " . $translation->start . " " . $translation->end . "\n";
      
      $translation->start($tstart+$diff);
      $translation->end($tend + $diff);
 #    $translation->start(18);
 #     $translation->end(1008);
      print STDERR "*** translation - after  " . $translation->start . " " . $translation->end . "\n";
    }
    
    
    # frameshifts - if > 1 frameshift we may just be buggered. My brain hurts.
    if(scalar($ex->sub_SeqFeature) > 1){
      print STDERR "uh-oh frameshift\n";
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
      
      # get any remaining exons
      foreach my $s(@sf){
	$transcript->add_Exon($s);
	$transcript->sort;
      }
      # flush the sub_SeqFeatures
      $ex->flush_sub_SeqFeature;
    }	      
  }
  return $transcript
}


sub transcript_from_multi_exon_genewise {
  my ($self, $current_exon, $transcript, $translation, $exoncount, $gw_gene, $eg_gene) = @_;
  
  my @gwtran  = $gw_gene->each_Transcript;
  my @gwexons = $gwtran[0]->each_Exon;
  
  my @egtran  = $eg_gene->each_Transcript;
  my @egexons = $egtran[0]->each_Exon;

  
  # compare to the first genewise exon
  if($gwexons[0]->strand == 1){
    if ($gwexons[0]->end == $current_exon->end && $current_exon->start <= $gwexons[0]->start){
      print STDERR "5' exon match!\n";
      # modify the coordinates of the first exon in $newtranscript
      my $ex = $transcript->start_exon;
      $ex->start($current_exon->start);
	      
      # add all the exons from the est2genome transcript, previous to this one
      my $c = 0;
      while($c < $exoncount){
	$transcript->add_Exon($egexons[$c]);
	$transcript->sort;
	$c++;
      }
      
      # fix translation start 
      # take what it was for the gw gene, and add on the extra
      my $tstart = $translation->start;
      $tstart += ($gwexons[0]->start - $ex->start);
      $translation->start($tstart);
      
    } # end 5' exon
    
    elsif ($gwexons[$#gwexons]->start == $current_exon->start && $current_exon->end >= $gwexons[$#gwexons]->end){
      print STDERR "3' exon match\n";
      
      # modify the coordinates of the last exon in $newtranscript
      my $ex = $transcript->end_exon;
      $ex->end($current_exon->end);
      
      # is it a frameshifted one?	3' exon is a special case as its end might have changed
      if(scalar($ex->sub_SeqFeature) > 1){
	print STDERR "3' exon frameshift + strand\n";
	my @sf = $ex->sub_SeqFeature;
	my $last = pop(@sf);
	
	# sort out start, id & phase
	$ex->start($last->start); # but don't you dare touch the end!
	$ex->id($last->id);
	$ex->phase($last->phase);

	# add back the remaining component exons
	foreach my $s(@sf){
	  $transcript->add_Exon($s);
	  print STDERR "added exon: " . $s->id . " : " . $s->start . "-" . $s->end . " end phase " . $s->end_phase . "\n";
	  $transcript->sort;
	}
	# flush the sub_SeqFeatures so we don't try to add this one again later
	$ex->flush_sub_SeqFeature;
      }

      # add all the exons from the est2genome transcript, subsequent to this one ($current_exon)
      my $c = $#egexons;
      while($c > $exoncount){
	$transcript->add_Exon($egexons[$c]);
	$transcript->sort;
	$c--;
      }
      
    } # end 3' exon
    
  }
  
  elsif($gwexons[0]->strand == -1){
    #	    print STDERR "***reverse strand\n";
    # first is last and last is first
    if ($gwexons[0]->start == $current_exon->start && $current_exon->end >= $gwexons[0]->end){
      print STDERR "5' exon match!\n";
      
      # modify the coordinates of the first exon in $newtranscript
      my $ex = $transcript->start_exon;
      $ex->end($current_exon->end);
      #	      print STDERR " here\n";
      # need to add back exons
      my $c = 0;
      while($c < $exoncount){
	$transcript->add_Exon($egexons[$c]);
	$transcript->sort;
	$c++;
      }
      
      # need to deal with translation start
      my $tstart = $translation->start;
      my $diff = $current_exon->end - $gwexons[0]->end;
      $translation->start($tstart+$diff);
      
      
      
    } # end 5' exon
    
    elsif ($gwexons[$#gwexons]->end == $current_exon->end && $current_exon->start <= $gwexons[$#gwexons]->start){
      print STDERR "3' exon match\n";
      
      # modify the coordinates of the last exon in $newtranscript
      my $ex = $transcript->end_exon;
      $ex->start($current_exon->start);
      
      # need to deal with frameshifts - 3' exon is a special case as its end might have changed
      if(scalar($ex->sub_SeqFeature) > 1){
	print STDERR "3' exon frameshift - strand\n";
	my @sf = $ex->sub_SeqFeature;
	my $last = pop(@sf);

	# sort out start, id & phase
	$ex->start($last->start); # but don't you dare touch the end!
	$ex->id($last->id);
	$ex->phase($last->phase);

	# add back the remaining component exons
	foreach my $s(@sf){
	  $transcript->add_Exon($s);
	  $transcript->sort;
	}
	# flush the sub_SeqFeatures?
	$ex->flush_sub_SeqFeature;
      }
      
      # add all the exons from the est2genome transcript, subsequent to this one
      my $c = $#egexons;
      while($c > $exoncount){
	$transcript->add_Exon($egexons[$c]);
	$transcript->sort;
	$c--;
      }
      
    } # end 3' exon
  }  
 return $transcript; 
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
  my $contig = $self->vc;
  my @genes;
  
$self->throw("[$analysis_obj] is not a Bio::EnsEMBL::Analysis\n") unless defined($analysis_obj) && $analysis_obj->isa("Bio::EnsEMBL::Analysis");

print STDERR "***analysis dbID: " . $analysis_obj->dbID . "\n";

  my $time = time; chomp($time);

  foreach my $tmpf (@$results) {
    my $gene   = new Bio::EnsEMBL::Gene;
    $gene->type($genetype);
    $gene->id($self->input_id . ".$genetype.$count");
    $gene->version(1);
    $gene->created($time);
    $gene->modified($time);

    my $transcript = $self->_make_transcript($tmpf,$self->vc,$genetype,$count);
	
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
  my $contig = $self->vc;

  my @genes = $self->gw_genes;
  push(@genes, $self->e2g_genes); # ??? do we want these?
  push(@genes, $self->combined_genes);

  foreach my $gene (@genes) {
    print STDERR "about to remap " . $gene->id . "\n";
    my @t = $gene->each_Transcript;
    my $tran = $t[0];

    eval {
      my $genetype = $gene->type;
      my $newgene = $contig->convert_Gene_to_raw_contig($gene);
      # need to explicitly add back genetype and analysis.
      $newgene->type($genetype);
      $newgene->analysis($gene->analysis);
      push(@newf,$newgene);

      # sort out supporting feature coordinates
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

      # is this a special case single coding exon gene with UTRS?
      if($tran->translation->start_exon_id() eq $tran->translation->end_exon_id() 
	 && $gene->type eq 'combined_gw_e2g'){
	print STDERR "single coding exon, with UTRs\n";
	
	# problems come about when we switch from + strand on FPC contig to - strand on raw contig.
	my $fpc_strand;

	foreach my $exon($tran->each_Exon) {
	  if ($exon->id eq $tran->translation->start_exon_id()) {
	    $fpc_strand = $exon->strand;
	    last;
	  }
	}
	
	foreach my $tran ($newgene->each_Transcript) {
	  foreach my $exon($tran->each_Exon) {
	    if ($exon->id eq $tran->translation->start_exon_id()) {
	      if($fpc_strand == 1 && $exon->strand == -1){
		print STDERR "fpc strand 1, raw strand -1 - flipping translation start/end\n";
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
      print STDERR "Couldn't reverse map gene " . $gene->id . " [$@]\n";
    }
  }

  return @newf;
}



=head2 _check_coverage

 Title   : _check_coverage
 Usage   :
 Function: checks how much of the parent protein is covered by the genewise prediction
 Example :
 Returns : 1 if > 80% coverage, otherwise 0
 Args    :


=cut

sub _check_coverage{
  my ($self, $gene) = @_;
  my $pstart = 0;
  my $pend = 0;
  my $protname;
  my $plength;
  my $fetcher = new Bio::EnsEMBL::Pipeline::SeqFetcher;

  my @gw_tran = $gene->each_Transcript;
  my @gw_ex = $gw_tran[0]->each_Exon;
  return 0 unless scalar(@gw_ex) == 1;

  foreach my $f($gw_ex[0]->each_Supporting_Feature){
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

  
  my $coverage = $pend - $pstart;
  $coverage /= $plength;
  $coverage *= 100;
  if ($coverage < 80){
    warn "Coverage of $protname by " . $gene->id . " is only $coverage\n";
    return 0;
  }
  
  print STDERR "***Coverage of $protname by " . $gene->id . " is $coverage\n";
  return 1;
}

=head2 _make_transcript

 Title   : make_transcript
 Usage   : $self->make_transcript($gene, $contig, $genetype, $count)
 Function: makes a Bio::EnsEMBL::Transcript from a SeqFeature representing a gene, 
           with sub_SeqFeatures representing exons.
 Example :
 Returns : Bio::EnsEMBL::Transcript with Bio::EnsEMBL:Exons(with supporting feature 
           data), and a Bio::EnsEMBL::translation
 Args    : $gene: Bio::EnsEMBL::SeqFeatureI, $contig: Bio::EnsEMBL::DB::ContigI,
           $genetype: string, $count: integer


=cut

sub _make_transcript{
  my ($self, $gene, $contig, $genetype, $count)=@_;
  $genetype = 'unspecified' unless defined ($genetype);
  $count = 1 unless defined ($count);

  unless ($gene->isa ("Bio::EnsEMBL::SeqFeatureI"))
    {print "$gene must be Bio::EnsEMBL::SeqFeatureI\n";}
  unless ($contig->isa ("Bio::EnsEMBL::DB::ContigI"))
    {print "$contig must be Bio::EnsEMBL::DB::ContigI\n";}

  my $time  = time; 
  chomp($time);

  my $transcript   = new Bio::EnsEMBL::Transcript;
  $transcript->id($contig->id . ".$genetype.$count");
  $transcript->version(1);

  my $translation  = new Bio::EnsEMBL::Translation;    
  $translation->id($contig->id . ".$genetype.$count");
  $translation->version(1);

  $transcript->translation($translation);

  my $excount = 1;
  my @exons;
    
  foreach my $exon_pred ($gene->sub_SeqFeature) {
    # make an exon
    my $exon = new Bio::EnsEMBL::Exon;
    
    $exon->id($contig->id . ".$genetype.$count.$excount");
    $exon->contig_id($contig->id);
    $exon->created($time);
    $exon->modified($time);
    $exon->version(1);
      
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
      $subf->feature1->analysis($exon_pred->analysis);
	
      $subf->feature2->source_tag($genetype);
      $subf->feature2->primary_tag('similarity');
      $subf->feature2->score(100);
      $subf->feature2->analysis($exon_pred->analysis);
      
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
  }
  
  return $transcript;
}

1;
