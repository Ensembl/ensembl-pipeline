#
# Ensembl module for Bio::EnsEMBL::Pipeline::RunnableDB::TargettedGeneWise.pm
#
# Cared for by Val Curwen <vac@sanger.ac.uk>
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::TargettedGeneDNA

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

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
use Bio::EnsEMBL::Pipeline::SeqFetcher;
use Bio::EnsEMBL::Pipeline::Runnable::Est2Genome;

@ISA = qw(Bio::Root::RootI Bio::EnsEMBL::Pipeline::RunnableDB);

sub new {
  my ($class,@args) = @_;
  my $self = bless {}, $class;

           
  my( $dbobj, $input_id ) = $self->_rearrange(['DBOBJ',
					       'INPUT_ID'], @args);
  
  $self->throw("No database handle input")           unless defined($dbobj);
  
#  $self->throw("[$dbobj] is not a Bio::EnsEMBL::DBSQL::Obj") unless $dbobj->isa("Bio::EnsEMBL::DBSQL::Obj");
  $self->throw("[$dbobj] is not a Bio::EnsEMBL::DBSQL::DBAdaptor") unless $dbobj->isa("Bio::EnsEMBL::DBSQL::DBAdaptor");
  $self->dbobj($dbobj);
  $dbobj->static_golden_path_type('UCSC');
  
  $self->throw("No input id input") unless defined($input_id);
  $self->input_id($input_id);

  return $self; # success - we hope!
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
  my @fps;
  my $fpc;
  my $pid; 
  
  my $start;
  my $end;
  my $seqfetcher = new Bio::EnsEMBL::Pipeline::SeqFetcher;


  
  # input format: chr22:10602496,10603128;Q9UGV6:AC00012
  if( !($entry =~ /(\S+):(\d+),(\d+):(\S+):(\S+)/) ) {
      $self->throw("Not a valid input id... $entry");
  }
  
  print STDERR "input: ".$entry . "\n";
  
  if ($fpc) { $self->throw("mixed fpc contigs") unless $fpc = $1; }
  if ($pid) { $self->throw("mixed protein hits") unless $pid = $4; }
 
  $fpc = $1;
  $pid = $4;
  my $dnaid = $5;

  my $fpcstart   = $2;
  my $fpcend     = $3;
  my $fpcstrand  = 1;
  
  if ($2 > $3) { # let blast sort it out
      $fpcstart  = $3;
      $fpcend    = $2;
      $fpcstrand = -1;
  }
  
  $start = $fpcstart unless (defined $start && $start < $fpcstart);
  $end = $fpcend unless (defined $end && $end > $fpcend);

    
  
  my $sgpa = $self->dbobj->get_StaticGoldenPathAdaptor();

  print STDERR "$fpc $start $end\n";

  my ($chrname,$chrstart,$chrend) = $sgpa->convert_fpc_to_chromosome($fpc,$start-10000,$end+10000);
  print STDERR "$chrname $chrstart $chrend\n";
  my $vc = $sgpa->fetch_VirtualContig_by_chr_start_end($chrname,$chrstart,$chrend);
  
  $self->vc($vc);
  
  my $r = Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise->new( -genomic => $vc->primary_seq,
								    -ids => [ $pid ] );
 
# ARGH"!!!! 
#  my $cdna = $seqfetcher->run_efetch($dnaid);
  my $cdna = $seqfetcher->run_pfetch($dnaid);

  my $e2g = Bio::EnsEMBL::Pipeline::Runnable::Est2Genome->new( -genomic => $vc->primary_seq,
							       -est => $cdna
							     );
    
  $self->runnable($r);
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
   $self->e2g_runnable->run();

   $self->convert_gw_output;
   $self->convert_e2g_output;

   $self->combine_genes;

   # remap to raw contig coords
   my @remapped = $self->remap_genes();
   
   # check translations
   foreach my $gene(@remapped){
     next unless $gene->type eq 'combined_gw_e2g';
     foreach my $trans ( $gene->each_Transcript ) {
       my $tseq;
       eval{
	 $tseq = $trans->translate();
       };
       
       if ($@) {
	 print STDERR "Couldn't translate: " . $gene->id . "[$@]\n";
       } 
       
       # check for stops
       if ( $tseq->seq =~ /\*/ ) {
	 $self->throw("UTR gene translation has stop codons, something is wrong\n");
       }

       # if we get here, all is well, so print out translation
       print STDERR "translation: \n";
       my $seqio = Bio::SeqIO->new(-fh => \*STDERR);
       print STDERR "remapped: ";
       $seqio->write_seq($trans->translate); 
       print STDERR "\n ";
     }
   }

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
      $obj->{'e2g_runnable'} = $value;
    }
    return $obj->{'e2g_runnable'};

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

    #$self->throw("exiting bfore write");

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
	  
	  print STDERR "genetype: " . $gene->type . "\n";
	  
	  
	  foreach my $trans ( $gene->each_Transcript ) {
	    $transcount++;
	    $translcount++;
	    if($gene->type eq 'combined_gw_e2g') {
	      eval {
		print STDERR "translation: \n";
		my $seqio = Bio::SeqIO->new(-fh => \*STDERR);
		print STDERR "checktrans: ";
		$seqio->write_seq($trans->translate); 
		print STDERR "\n ";
	      };
	      
	      if ($@) {
		print STDERR "Couldn't translate: " . $gene->id . "[$@]\n";
	      }
	    }
	  }
	  foreach my $exon ( $gene->each_unique_Exon() ) {
	    $exoncount++;
	    foreach my $sf($exon->each_Supporting_Feature) {
	      print STDERR "***sub_align: " . 
		           $sf->seqname  . "\t" .
		           $sf->start    . "\t" .
		           $sf->end      . "\t" .
		           $sf->strand   . "\t" .
			   $sf->score   . "\t" .
			   $sf->hseqname . "\t" .
			   $sf->hstart   . "\t" .
			   $sf->hend     . "\n";
	    }
	    
	    }
	}
	
	$self->throw("exiting bfore write");
	
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
	  
	  $gene_obj->write($gene);
      }; 
      if( $@ ) {
	  print STDERR "UNABLE TO WRITE GENE\n\n$@\n\nSkipping this gene\n";
      }
	    
  }
   
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
  my @genes = $self->make_genes($count, $genetype, \@results);
  
  if (!defined($self->{_e2g_genes})) {
    $self->{_e2g_genes} = [];
  }

  print STDERR "e2g genes: " . scalar(@genes) . "\n";

  push(@{$self->{_e2g_genes}},@genes);  
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
  my @genes    = $self->make_genes($count, $genetype, \@results);

  if (!defined($self->{_gw_genes})) {
    $self->{_gw_genes} = [];
  }

  print STDERR "gw genes: " . scalar(@genes) . "\n";

  push(@{$self->{_gw_genes}},@genes);
  
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
  
  # merge the genewise and est2genome predictions. Only one new transcript per genewise gene ...
  my @newtrans = $self->_make_newtranscripts(@merged_gw_genes);

  # make some lovely genes
  my @genes;
  my $count=0;
  my $genetype = 'combined_gw_e2g';
  foreach my $trans(@newtrans){
    my $gene = new Bio::EnsEMBL::Gene;
    $gene->type($genetype);
    $gene->id($self->input_id . "$genetype.$count"); 
    $gene->version(1);
    $gene->add_Transcript($trans);
    push (@genes,$gene);
    $count++;
  }

  if (!defined($self->{_combined_genes})) {
    $self->{_combined_genes} = [];
  }
  push(@{$self->{_combined_genes}},@genes);
  
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
  foreach my $gwg(@{$self->{_gw_genes}}){
    my $gene = new Bio::EnsEMBL::Gene;
    $gene->type('combined');
    #  $gene->id($self->input_id . ".combined.$count");
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
      
      print STDERR "exon id: " . $exon->id . "\n";

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
  }
  return @merged;
}

=head2 _merge_gw_genes

 Title   : _merge_gw_genes
 Usage   :
 Function: makes new transcripts by combining the genewise and est2genome predictions. Its a monster.
 Example :
 Returns : 
 Args    :


=cut

sub _make_newtranscripts {
  my ($self, @merged_gw_genes) = @_;
  my @gw_genes  = @{$self->{_gw_genes}};
  my @e2g_genes = @{$self->{_e2g_genes}};
  my @newtrans  = ();

 GENE:
  foreach my $gene(@merged_gw_genes) {
    print "\nGENE\n";
    my @gw_tran = $gene->each_Transcript;
    my @gw_ex = $gw_tran[0]->each_Exon; # need ordered array
    my $strand = $gw_ex[0]->strand;
    $self->warn("first and last gw exons have different strands - odd things will happen\n") if ($gw_ex[$#gw_ex]->strand != $strand);
    
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
      my $newtranscript = dclone($gw_tran[0]);;
      my $translation  = dclone($gw_tran[0]->translation);
      $newtranscript->translation($translation);
      my $eecount = 0;

      print "e2g exons: " . scalar(@eg_ex) . "\n";

      foreach my $ee(@eg_ex){
	$self->warn("gw and e2g exons have different strands - odd things will happen\n") 
	  if ($ee->strand != $strand);

	# single exon genewise prediction?
	if(scalar(@gw_ex) == 1) {# eeeep
	  if ($gw_ex[0]->start >= $ee->start && $gw_ex[0]->end <= $ee->end){
	    print STDERR "single exon gene\n";	    
	    # modify the coordinates of the first exon in $newtranscript
	    my $ex = $newtranscript->start_exon;
	    
	    $ex->start($ee->start);
	    $ex->end($ee->end);

	    print STDERR "eecount: $eecount\n";
    
	    # need to add back exons, both 5' and 3'
	    my $c = 0;
	    while($c < $eecount){
	      print STDERR "adding 5' exon\n";
	      $newtranscript->add_Exon($eg_ex[$c]);
	      $newtranscript->sort;
	      $c++;
	    }
	    
	    # add all the exons from the est2genome transcript, subsequent to this one
	    my $c = $#eg_ex;
	    while($c > $eecount){
	      print STDERR "adding 3' exon\n";
	      $newtranscript->add_Exon($eg_ex[$c]);
	      $newtranscript->sort;
	      $c--;
	    }
	    
	    # need to deal with translation start and end this time - varies depending on strand
	    if($strand == 1){
	      my $diff = $gw_ex[0]->start - $ex->start;
	      my $tstart = $translation->start;
	      my $tend = $translation->end;
	      
	      print STDERR "***gw  " . $gw_ex[0]->start . " " . $gw_ex[0]->end . "\n";
	      $translation->start($tstart + $diff);
	      $translation->end($tend + $diff);
	    }
	    elsif($strand == -1){
	      print STDERR "***reverse\n";
	      #	    my $diff = $ee->end - $gw_ex[0]->end;
	      my $diff = $gw_ex[0]->start - $ee->start;
	      my $tstart = $translation->start;
	      my $tend = $translation->end;
	      print STDERR "***gw  " . $gw_ex[0]->start . " " . $gw_ex[0]->end . "\n";
	      
	      $translation->start($tstart+$diff);
	      $translation->end($tend + $diff);
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
	      $newtranscript->add_Exon($last);
	      
	      # get any remaining exons
	      foreach my $s(@sf){
		$newtranscript->add_Exon($s);
		$newtranscript->sort;
	      }
	      # flush the sub_SeqFeatures
	      $ex->flush_sub_SeqFeature;
	    }	      
	  }
	}
	
	# multiple exon genewise prediction
	else {

	  # compare to the first genewise exon
	  if($strand == 1){
	    if ($gw_ex[0]->end == $ee->end && $ee->start <= $gw_ex[0]->start){
	      print STDERR "5' exon match!\n";
	      # modify the coordinates of the first exon in $newtranscript
	      my $ex = $newtranscript->start_exon;
	      $ex->start($ee->start);
	      
	      # add all the exons from the est2genome transcript, previous to this one
	      my $c = 0;
	      while($c < $eecount){
		$newtranscript->add_Exon($eg_ex[$c]);
		$newtranscript->sort;
		$c++;
	      }
	      
	      # fix translation start 
	      # take what it was for the gw gene, and add on the extra
	      my $tstart = $translation->start;
	      $tstart += ($gw_ex[0]->start - $ex->start);
	      $translation->start($tstart);
	      
	    } # end 5' exon
	    
	    elsif ($gw_ex[$#gw_ex]->start == $ee->start && $ee->end >= $gw_ex[$#gw_ex]->end){
	      print STDERR "3' exon match\n";
	      
	      # modify the coordinates of the last exon in $newtranscript
	      my $ex = $newtranscript->end_exon;
	      $ex->end($ee->end);
	      
	      # is it a frameshifted one?	3' exon is a special case as its end might have changed
	      if(scalar($ex->sub_SeqFeature) > 1){
		print STDERR "3' exon frameshift\n";
		print STDERR $ex->id . "\n";
		print STDERR " before: " . $ex->id . " : " .$ex->start . "-" . $ex->end . "\n";
		my @sf = $ex->sub_SeqFeature;
		my $last = pop(@sf);
		
		# ids will be messed up
		
		$ex->start($last->start); # but don't you dare touch the end!
		$ex->id($last->id);
		print STDERR " after: " . $ex->id . " : " .$ex->start . "-" . $ex->end . "\n";
		# add back the remaining component exons
		foreach my $s(@sf){
		  $newtranscript->add_Exon($s);
		  print STDERR "added exon: " . $s->id . " : " . $s->start . "-" . $s->end . "\n";
		  $newtranscript->sort;
		}
		# flush the sub_SeqFeatures so we don't try to add this one again later
		$ex->flush_sub_SeqFeature;
	      }
	      
	      # add all the exons from the est2genome transcript, subsequent to this one ($ee)
	      my $c = $#eg_ex;
	      while($c > $eecount){
		$newtranscript->add_Exon($eg_ex[$c]);
		$c--;
	      }
	      
	    } # end 3' exon
	    
	  }

	  elsif($strand == -1){
	    print STDERR "***reverse strand\n";
	    # first is last and last is first
	    if ($gw_ex[0]->start == $ee->start && $ee->end >= $gw_ex[0]->end){
	      print STDERR "5' exon match!\n";
	      
	      # modify the coordinates of the first exon in $newtranscript
	      my $ex = $newtranscript->start_exon;
	      $ex->end($ee->end);
	      print STDERR " here\n";
	      # need to add back exons
	      my $c = 0;
	      while($c < $eecount){
		$newtranscript->add_Exon($eg_ex[$c]);
		$newtranscript->sort;
		$c++;
	      }
	      
	      # need to deal with translation start
	      my $tstart = $translation->start;
	      my $diff = $ee->end - $gw_ex[0]->end;
	      $translation->start($tstart+$diff);
	      
	      
	      
	    } # end 5' exon
	    
	    elsif ($gw_ex[$#gw_ex]->end == $ee->end && $ee->start <= $gw_ex[$#gw_ex]->start){
	      print STDERR "3' exon match\n";
	      
	      # modify the coordinates of the last exon in $newtranscript
	      my $ex = $newtranscript->end_exon;
	      $ex->start($ee->start);
	      
	      # need to deal with frameshifts - 3' exon is a special case as its end might have changed
	      if(scalar($ex->sub_SeqFeature) > 1){
		print STDERR "3' exon frameshift\n";
		my @sf = $ex->sub_SeqFeature;
		my $last = pop(@sf);
		$ex->start($last->start); # but don't you dare touch the end!
		# add back the remaining component exons
		foreach my $s(@sf){
		  $newtranscript->add_Exon($s);
		  $newtranscript->sort;
		}
		# flush the sub_SeqFeatures?
		$ex->flush_sub_SeqFeature;
	      }
	      
	      # add all the exons from the est2genome transcript, subsequent to this one
	      my $c = $#eg_ex;
	      while($c > $eecount){
		$newtranscript->add_Exon($eg_ex[$c]);
		$c--;
	      }
	      
	    } # end 3' exon
	  }
	}
  
	# increment the exon
	$eecount++;
	
      } # end foreach my $ee
      
      # check the transcript and expand frameshifts in all but original 3' gw_exon
      print STDERR "NEWTRANSCRIPT $newtranscript\n";
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
	}
	
	eval {
	  print STDERR "translation: \n";
	  my $seqio = Bio::SeqIO->new(-fh => \*STDERR);
	  print STDERR "grargh: ";
	  $seqio->write_seq($newtranscript->translate); 
	  print STDERR "\n ";
	};
	
	 if ($@) {
	   print STDERR "Couldn't translate: " . $gene->id . " plus " . $eg->id  . "[$@]\n";
	 }
	
      }
      $foundtrans = 1;
      push (@newtrans, $newtranscript); 
    }
  }
  return @newtrans;

}


=head2 make_genes

 Title   : make_genes
 Usage   :
 Function: 
 Example :
 Returns : 
 Args    :


=cut

sub make_genes {
  my ($self,$count,$genetype,$results) = @_;
  my $contig = $self->vc;
  my @genes;
  
  foreach my $tmpf (@$results) {
    my $gene   = new Bio::EnsEMBL::Gene;
    $gene->type($genetype);
    $gene->id($self->input_id . ".$genetype.$count");
    $gene->version(1);

    my $transcript = $self->_make_transcript($tmpf,$self->vc,$genetype,$count);

    # add transcript to gene
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

  my @genes = @{$self->{_gw_genes}};
  push(@genes, @{$self->{_e2g_genes}});
  push(@genes, @{$self->{_combined_genes}});

  foreach my $gene (@genes) {
    print STDERR "about to remap\n";
    my @t = $gene->each_Transcript;
    my $tran = $t[0];
    eval {
      my $genetype = $gene->type;
      my $newgene = $contig->convert_Gene_to_raw_contig($gene);
      $newgene->type($genetype);

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
      print STDERR $f->hseqname . " " . $f->hstart . " " . $f->hend . "\n";

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
  
  my $seq = $fetcher->run_pfetch($protname);
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
 Usage   :
 Function: 
 Example :
 Returns : 
 Args    :


=cut

sub _make_transcript{
  my ($self,$gene,$contig,$genetype,$count)=@_;
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
    
    $exon->phase($exon_pred->{_phase});
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

=head2 input_id

 Title   : input_id
 Usage   : $obj->input_id($newval)
 Function: 
 Returns : value of input_id
 Args    : newvalue (optional)


=cut

sub input_id{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'input_id'} = $value;
    }
    return $obj->{'input_id'};

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
      $obj->{'runnable'} = $value;
    }
    return $obj->{'runnable'};

}

=head2 vc

 Title   : vc
 Usage   : $obj->vc($newval)
 Function: 
 Returns : value of vc
 Args    : newvalue (optional)


=cut

sub vc{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'vc'} = $value;
    }
    return $obj->{'vc'};

}

1;
