#!/usr/local/bin/perl

#
#
# Cared for by Michele Clamp  <michele@sanger.ac.uk>
#
# Copyright Michele Clamp
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::Riken_BlastMiniGenewise

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::RunnableDB::Riken_BlastMiniGenewise->new(
					     -dbobj     => $db,
					     -input_id  => $id
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

# Object preamble - inherits from Bio::Root::Object;
use Bio::EnsEMBL::Pipeline::RunnableDBI;
use Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise;
use Bio::EnsEMBL::Pipeline::GeneConf qw (EXON_ID_SUBSCRIPT
					 TRANSCRIPT_ID_SUBSCRIPT
					 GENE_ID_SUBSCRIPT
					 PROTEIN_ID_SUBSCRIPT
					 );

use Data::Dumper;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDBI Bio::Root::Object );

sub _initialize {
    my ($self,@args) = @_;
    my $make = $self->SUPER::_initialize(@_);    
           
    my( $dbobj,$input_id ) = $self->_rearrange(['DBOBJ',
						'INPUT_ID'], @args);
       
    $self->throw("No database handle input")                 unless defined($dbobj);
    $self->dbobj($dbobj);

    $self->throw("No input id input") unless defined($input_id);
    $self->input_id($input_id);
    
    return $self; # success - we hope!
}
sub input_id {
	my ($self,$arg) = @_;

   if (defined($arg)) {
      $self->{_input_id} = $arg;
   }

   return $self->{_input_id};
}

=head2 dbobj

    Title   :   dbobj
    Usage   :   $self->dbobj($db)
    Function:   Get/set method for database handle
    Returns :   Bio::EnsEMBL::Pipeline::DB::ObjI
    Args    :   

=cut

sub dbobj {
    my( $self, $value ) = @_;    
    if ($value) {

        $value->isa("Bio::EnsEMBL::DB::ObjI") || $self->throw("Input [$value] isn't a Bio::EnsEMBL::DB::ObjI");
        $self->{'_dbobj'} = $value;
    }
    return $self->{'_dbobj'};
}

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

#$self->throw("Bailing before real write\n");

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
    Function:   Fetches input data for est2genome from the database
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

    my @features  = $contig->get_all_SimilarityFeatures_above_score('riken_prot',200);

    
    print STDERR "Number of features = " . scalar(@features) . "\n";

    my %idhash;
    
    foreach my $f (@features) {
#        print "Feature " . $f . " " . $f->seqname . " " . $f->source_tag . "\n";
      if ($f->isa("Bio::EnsEMBL::FeaturePair") && 
	  defined($f->hseqname) ) {
      $idhash{$f->hseqname} = 1;
      
    }
  }
    
    my @ids = keys %idhash;

    print STDERR "Feature ids are @ids\n";

    my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise('-genomic'  => $genseq,
									   '-ids'      => \@ids,
									   '-trim'     => 1);
    
    
    $self->add_Runnable($runnable);
    $self->{$runnable} = $contig;

}
     
sub add_Runnable {
    my ($self,$arg) = @_;

    if (!defined($self->{_runnables})) {
	$self->{_runnables} = [];
    }

    if (defined($arg)) {
	if ($arg->isa("Bio::EnsEMBL::Pipeline::RunnableI")) {
	    push(@{$self->{_runnables}},$arg);
	} else {
	    $self->throw("[$arg] is not a Bio::EnsEMBL::Pipeline::RunnableI");
	}
    }
}
sub get_Runnables {
    my ($self) = @_;

    if (!defined($self->{_runnables})) {
	$self->{_runnables} = [];
    }
    
    return @{$self->{_runnables}};
}

sub run {
    my ($self) = @_;

    foreach my $runnable ($self->get_Runnables) {
	$runnable->run;
    }
    
    $self->convert_output;

}

sub convert_output {
  my ($self) =@_;

#$self->throw("exiting before convert_output\n"); 
 
  my $count = 1;
  my $time  = time; chomp($time);
  
  my $trancount = 1;
  
  foreach my $runnable ($self->get_Runnables) {
    my @genes = $self->make_genes($count,$time,$runnable);

    my @remapped = $self->remap_genes($runnable,@genes);

    # change of plan - store the genes
    #    # now get out the feature data we want and store THAT in the output array
    #    $self->generate_outfeat(@remapped);
    if (!defined($self->{_output})) {
      $self->{_output} = [];
    }
    
    push(@{$self->{_output}},@remapped);
    
  }
}

sub check_splice {
    my ($self,$f1,$f2) = @_;
    
    my $splice1 = substr($self->{_genseq}->seq,$f1->end,2);
    my $splice2 = substr($self->{_genseq}->seq,$f2->start-3,2);
    
    if (abs($f2->start - $f1->end) > 50) {
	print ("Splices are " . $f1->hseqname . " [" . 
	                        $splice1      . "][" . 
	                        $splice2      . "] " . 
	       ($f2->start - $f1->end)        . "\n");
    }
}


sub output {
    my ($self) = @_;
   
    if (!defined($self->{_output})) {
      $self->{_output} = [];
    } 
    return @{$self->{_output}};
}

sub make_genes {

  my ($self,$count,$time,$runnable) = @_;
  my $contig = $self->{$runnable};
  my @tmpf   = $runnable->output;
  
  my @genes;
    
  foreach my $tmpf (@tmpf) {
    my $gene   = new Bio::EnsEMBL::Gene;
    my $tran   = new Bio::EnsEMBL::Transcript;
    my $transl = new Bio::EnsEMBL::Translation;
    
    $gene->type('riken_genewise');
    $gene->id($self->input_id . ".riken_genewise.$count");
    $gene->created($time);
    $gene->modified($time);
    $gene->version(1);
    
    $tran->id($self->input_id . ".riken_genewise.$count");
    $tran->created($time);
    $tran->modified($time);
    $tran->version(1);
    
    $transl->id($self->input_id . ".riken_genewise.$count");
    $transl->version(1);
    
    $count++;
    
    $gene->add_Transcript($tran);
    $tran->translation($transl);
    
    my $excount = 1;
    my @exons;
    
    foreach my $subf ($tmpf->sub_SeqFeature) {
      $subf->feature1->source_tag('riken_genewise');
      $subf->feature1->primary_tag('similarity');
      $subf->feature1->score(100);
      #	$subf->feature1->analysis($analysis);
      
      $subf->feature2->source_tag('riken_genewise');
      $subf->feature2->primary_tag('similarity');
      $subf->feature2->score(100);
      #	$subf->feature2->analysis($analysis);
      
      my $exon = new Bio::EnsEMBL::Exon;
      
      $exon->id($self->input_id . ".riken_genewise.$count.$excount");
      $exon->contig_id($contig->id);
      $exon->created($time);
      $exon->modified($time);
      $exon->version(1);
      
      $exon->start($subf->start);
      $exon->end  ($subf->end);
      $exon->strand($subf->strand);
      
      print STDERR "\tFeaturePair " . $subf->gffstring . "\n";
      
      $exon->phase($subf->feature1->{_phase});
      $exon->attach_seq($self->{$runnable}->primary_seq);
      $exon->add_Supporting_Feature($subf);
      
      my $seq   = new Bio::Seq(-seq => $exon->seq->seq);
      
      my $tran0 =  $seq->translate('*','X',0)->seq;
      my $tran1 =  $seq->translate('*','X',2)->seq;
      my $tran2 =  $seq->translate('*','X',1)->seq;
      
      print STDERR "\n\t exon phase 0 : " . $tran0 . " " . $exon->phase . "\n";
      print STDERR "\t exon phase 1 : " . $tran1 . "\n";
      print STDERR "\t exon phase 2 : " . $tran2 . "\n";
      
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

sub remap_genes {
  my ($self,$runnable,@genes) = @_;
  
  my $contig = $self->{$runnable};
  print STDERR "genes: " . scalar(@genes) . "\n";
  print STDERR "contig: $contig\n";

  my @newf;
  my $trancount=1;
  foreach my $gene (@genes) {
    foreach my $tran ($gene->each_Transcript) {
      print STDERR " Translation is " . $tran->translate->seq . "\n";
      foreach my $exon ($tran->each_Exon) {
	my $strand = "+";
	if ($exon->strand == -1) {
	  $strand = "-";
	}
	print STDERR $exon->contig_id . "\tgenewise\texon\t" . $exon->start . "\t" . $exon->end . "\t100\t" . $strand .  "\t" . $exon->phase . "\t" . $tran->id . ".$trancount\n";
      }
      $trancount++;
    }
    
    eval {
      my $newgene = $contig->convert_Gene_to_raw_contig($gene);
      $newgene->type('riken_genewise');
      push(@newf,$newgene);

    };
    if ($@) {
      print STDERR "Couldn't reverse map gene " . $gene->id . " [$@]\n";
    }
    

  }

  return @newf;

}

sub generate_outfeat {
  my ($self, @genes) = @_;

  my @newf;

  foreach my $gene(@genes) {
    foreach my $tran($gene->each_Transcript){
      foreach my $ex($tran->each_Exon){
	foreach my $sf($ex->each_Supporting_Feature){
	  # already is a feature pair
	  # need to sort out prot coords ...
	  print STDERR "SF: " . 
	    $sf->seqname . "\t" .
	    $sf->start . "\t" .
	    $sf->end . "\t" .
	    $sf->strand . "\t" .
	    $sf->hseqname . "\t" .
	    $sf->hstart . "\t" .
	    $sf->hend . "\n" ;
	  push (@newf, $sf);
	}
      }
    }
  }

  return @newf;

}

1;


