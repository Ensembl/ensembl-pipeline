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

Bio::EnsEMBL::Pipeline::RunnableDB::AlignFeature

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::RunnableDB::MiniGenewise->new(
					     -dbobj     => $db,
					     -input_id  => $id
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

# Object preamble - inherits from Bio::Root::Object;
use Bio::EnsEMBL::Pipeline::RunnableDBI;
use Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise;

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

    my $db = $self->dbobj();
    my $gene_obj = $db->gene_Obj;


    my $EXON_ID_SUBSCRIPT       = "BGWE";
    my $TRANSCRIPT_ID_SUBSCRIPT = "BGWT";
    my $GENE_ID_SUBSCRIPT       = "BGWG";
    my $PROTEIN_ID_SUBSCRIPT    = "BGWP";

    my $sth = $db->prepare("lock tables gene write, genetype write, exon write, transcript write, exon_transcript write, translation write,dna read,contig read,clone read,feature read,analysis read");
    $sth->execute;


    eval {
	(my $gcount = $gene_obj->get_new_GeneID($GENE_ID_SUBSCRIPT))
	    =~ s/$GENE_ID_SUBSCRIPT//;
	(my $tcount = $gene_obj->get_new_TranscriptID($TRANSCRIPT_ID_SUBSCRIPT))
	    =~ s/$TRANSCRIPT_ID_SUBSCRIPT//;
	(my $pcount = $gene_obj->get_new_TranslationID($PROTEIN_ID_SUBSCRIPT))
	    =~ s/$PROTEIN_ID_SUBSCRIPT//;
	print STDERR "Weiiiiird $EXON_ID_SUBSCRIPT\n";
	(my $ecount = $gene_obj->get_new_ExonID($EXON_ID_SUBSCRIPT))
	    =~ s/$EXON_ID_SUBSCRIPT//;
	

	foreach my $gene (@features) {
	    print STDERR "Feature is $gene\n";

	    print STDERR "Exon stub is $EXON_ID_SUBSCRIPT\n";
	    
	    
	    $gene->id($GENE_ID_SUBSCRIPT . $gcount);
	    $gcount++;
	    print (STDERR "Writing gene " . $gene->id . "\n");
	    
            # Convert all exon ids and save in a hash
            my %namehash;
	    
            foreach my $ex ($gene->each_unique_Exon) {
		print STDERR "Exon id " . $ex . " " . $ex->id . " " . $EXON_ID_SUBSCRIPT . "\n";
		$namehash{$ex->id} = $EXON_ID_SUBSCRIPT.$ecount;
		$ex->id($EXON_ID_SUBSCRIPT.$ecount);
		$ecount++;
            }
	    
	    print (STDERR "Transcripts are\n");
	    foreach my $tran ($gene->each_Transcript) {
		$tran->id             ($TRANSCRIPT_ID_SUBSCRIPT . $tcount);
		$tran->translation->id($PROTEIN_ID_SUBSCRIPT . $pcount);
		
		my $translation = $tran->translation;
		
		$tcount++;
		$pcount++;
		
		print (STDERR "Transcript  " . $tran->id . "\n");
		print (STDERR "Translation " . $tran->translation->id . "\n");
		
		foreach my $ex ($tran->each_Exon) {
                    my @sf = $ex->each_Supporting_Feature;
                    print STDERR "Supporting features are " . scalar(@sf) . "\n";
		    
                    if ($namehash{$translation->start_exon_id} ne "") {
			$translation->start_exon_id($namehash{$translation->start_exon_id});
                    }
                    if ($namehash{$translation->end_exon_id} ne "") {
			$translation->end_exon_id  ($namehash{$translation->end_exon_id});
                    }
		    print(STDERR "Exon         " . $ex->id . "\n");
		}
		
	    }
	    
	    $gene_obj->write($gene);
	}
	};
	if ($@) {
	    $sth = $db->prepare("unlock tables");
	    $sth->execute;
	    
	    $self->throw("Error writing gene for " . $self->input_id . " [$@]\n");
	} else {
	    $sth = $db->prepare("unlock tables");
	    $sth->execute;
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

    my $contigid  = $self->input_id;
       $contigid =~ s/\.(.*)//;
    my $contignum = $1;

    print STDERR "Contig id = $contigid , number $contignum\n";

    $self->dbobj->static_golden_path_type('SANGER');

    my $stadaptor = $self->dbobj->get_StaticGoldenPathAdaptor();

    my $contig    = $stadaptor->fetch_VirtualContig_by_chr_start_end('chr20',1000000,2000000);
    #my $contig    = $contig[$contignum];
    $contig->_chr_name('chr20');
    print STDERR "Analysing contig " . $contig->id . " " . $contignum . "\n";

    foreach my $rc ($contig->_vmap->each_MapContig) {
	my $strand = "+";
	if ($rc->orientation == -1) {
	    $strand = "-";
	}
	
	print STDERR $rc->contig->id . "\tsequence\t" . $rc->contig->id . "\t" . $rc->start . "\t" . $rc->end . "\t100\t" . $strand . "\t0\n";
    }
    my $genseq    = $contig->primary_seq;

    print STDERR "Length is " . $genseq->length . "\n";
    print STDERR "Fetching features " . $contig . "\n";

    my @features  = $contig->get_all_SimilarityFeatures_above_score('swir',200);

    print STDERR "Number of features = " . scalar(@features) . "\n";

    my %idhash;
    
    foreach my $f (@features) {
        print "Feature " . $f . "\n";
	if ($f->isa("Bio::EnsEMBL::FeaturePair") && defined($f->hseqname)) {
	    $idhash{$f->hseqname} = 1;
	}
    }
    
    my @ids = keys %idhash;
    
    print STDERR "Feature ids are @ids\n";
    
    my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise('-genomic'  => $genseq,
									   '-ids'      => \@ids);
    
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

    my @tmpf;
    my $count = 1;
    my $time  = time; chomp($time);

    # This BAD! Shouldn't be using internal ids.
    # <sigh> no time to change it now
    my $analysis = $self->dbobj->get_OldAnalysis(7);
    my $trancount = 1;
    foreach my $runnable ($self->get_Runnables) {
	my $contig = $self->{$runnable};
	my @tmpf   = $runnable->output;

	my @genes;
	
	foreach my $tmpf (@tmpf) {

	    my $gene   = new Bio::EnsEMBL::Gene;
	    my $tran   = new Bio::EnsEMBL::Transcript;
	    my $transl = new Bio::EnsEMBL::Translation;
	    
	    $gene->type('genewise');
	    $gene->id($self->input_id . ".genewise.$count");
	    $gene->created($time);
	    $gene->modified($time);
	    $gene->version(1);
	    
	    $tran->id($self->input_id . ".genewise.$count");
	    $tran->created($time);
	    $tran->modified($time);
	    $tran->version(1);
	    
	    $transl->id($self->input_id . ".genewise.$count");
	    $transl->version(1);
	    
	    $count++;
	    
	    $gene->add_Transcript($tran);
	    $tran->translation($transl);
	    
	    
	    my $excount = 1;
	    my @exons;
	    
	    foreach my $subf ($tmpf->sub_SeqFeature) {
		$subf->feature1->source_tag('genewise');
		$subf->feature1->primary_tag('similarity');
		$subf->feature1->score(100);
		$subf->feature1->analysis($analysis);

		$subf->feature2->source_tag('genewise');
		$subf->feature2->primary_tag('similarity');
		$subf->feature2->score(100);
		$subf->feature2->analysis($analysis);

		my $exon = new Bio::EnsEMBL::Exon;

		$exon->id($self->input_id . ".genewise.$count.$excount");
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
		
		if ($exons[0]->strand == 1) {
		    $transl->start($exons[0]->start);
		    $transl->end  ($exons[$#exons]->end);
		} else {
		    $transl->start($exons[0]->end);
		    $transl->end  ($exons[$#exons]->start);
		} 
		
		my @newf;

		foreach my $gene (@genes) {
		    foreach my $tran ($gene->each_Transcript) {
			print STDERR " Translation is " . $tran->translate->seq . "\n";
			foreach my $exon ($tran->each_Exon) {
			    my $strand = "+";
			    if ($exon->strand == -1) {
				$strand = "-";
			    }
			    print STDERR $exon->contig_id . "\tgenewise\tsexon\t" . $exon->start . "\t" . $exon->end . "\t100\t" . $strand .  "\t" . $exon->phase . "\t" . $tran->id . ".$trancount\n";
			}
			$trancount++;
		    }
		    
		    eval {
			my $newgene = $contig->convert_Gene_to_raw_contig($gene);
			$newgene->type('genewise');
			push(@newf,$newgene);
		    };
		    if ($@) {
			print STDERR "Couldn't reverse map gene " . $gene->id . " [$@]\n";
		    }
		}
		
		
		if (!defined($self->{_output})) {
		    $self->{_output} = [];
		}
		
		push(@{$self->{_output}},@newf);
	    }
	}
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


1;


