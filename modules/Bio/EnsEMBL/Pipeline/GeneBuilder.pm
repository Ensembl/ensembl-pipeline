#
# Object for submitting jobs to and querying the LSF queue
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

Bio::EnsEMBL::Pipeline::GeneBuilder

=head1 SYNOPSIS

# This is the main analysis database

my $db = new Bio::EnsEMBL::DBSQL::Obj(-host   => 'obi-wan',
				      -user   => 'ensro',
				      -dbname => 'ens500',
				      );

# Fetch a clone and its contigs from the database
my $clone       = $db   ->get_Clone($clone);
my @contigs     = $clone->get_all_Contigs;

# The genebuilder object will fetch all the features from the contigs and use them
# to first construct exons, then join those exons into exon pairs.  These exon apris are
# then made into transcripts and finally all overlapping transcripts are put together into one gene.


my $genebuilder = new Bio::EnsEMBL::Pipeline::GeneBuilder(-contigs => \@contigs);

my @genes       = $genebuilder->build_Genes;

# After the genes are built they can be used to order the contigs they are on.

my @contigs     = $genebuilder->order_Contigs;


=head1 DESCRIPTION

Takes in contigs and returns genes.  The procedure is currently reimplementing the TimDB method of
building genes where genscan exons are confirmed by similarity features which are then joined together
into exon pairs.  An exon pair is constructed as follows :

  ---------          --------    genscan exons
    -------          ----->      blast hit which spans an intron
    1     10        11    22        

For an exon pair to make it into a gene there must be at least 2 blast hits (features) that span across 
an intron.  This is called the coverage of the exon pair.

After all exon pairs have been generated for all the genscan exons there is a recursive routine 
(_recurseTranscripts) that looks for all exons that are the start of an exon pair with no 
preceding exons.  The exon pairs are followed recursively (including alternative splices) to build up 
full set of transcripts.

To generate the genes the transcripts are grouped together into sets with overlapping exons.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::GeneBuilder;

use Bio::EnsEMBL::Pipeline::ExonPair;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::MappedExon;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Pipeline::GeneConf qw (EXON_ID_SUBSCRIPT
					 TRANSCRIPT_ID_SUBSCRIPT
					 GENE_ID_SUBSCRIPT
					 PROTEIN_ID_SUBSCRIPT
					 );
use FreezeThaw;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object;

@ISA = qw(Bio::Root::Object);

sub _initialize {
    my ($self,@args) = @_;

    my $make = $self->SUPER::_initialize(@args);

    my ($contig) = $self->_rearrange([qw(CONTIG
					 )],@args);

    $self->throw("Must input a contig to GeneBuilder") unless defined($contig);
    
    $self->contig($contig);
    $self->get_all_Features;

    return $make; # success - we hope!
}

=head2 build_Genes

 Title   : 
 Usage   : my @genes = $self->build_Genes
 Function: 
 Example : 
 Returns : 
 Args    : 
           


=cut

sub build_Genes {
    my ($self) = @_;

    print STDERR "Building genes\n";

    $self->get_Features;
    $self->make_Exons;
    
    $self->make_ExonPairs;
    $self->print_ExonPairs;
    
    $self->link_ExonPairs;

    $self->filter_Transcripts;
    $self->make_Genes;
    $self->print_Genes;

    print STDERR "Finished printing genes...\n";

    $self->print_gff;

    print STDERR "Out of build Genes...\n";

}

=head2 get_Features

 Title   : get_Features
 Usage   : my ($features,$genscan) = $self->get_Features;
 Function: Gets all the features from all the contigs and sorts them into
           similarity features and genscan exons.
 Example : 
 Returns : Array ref of Bio::EnsEMBL::SeqFeature,array ref of Bio::EnsEMBL::SeqFeature
 Args    : none

=cut

sub get_Features {
    my ($self) = @_;

    if (defined($self->{_got_features})) {
	return;
    }
    
    my @features;
    my @genscan;
    my @genewise;

    my $contig = $self->contig;
    my @tmp    = $self->get_all_Features;

    @tmp = sort {$a->start <=> $b->start} @tmp;

    foreach my $f (@tmp) {

	if (defined($f->sub_SeqFeature)) {

	    print STDERR "\nGot feature set " . $f->seqname    . "\t" . 
		                                $f->start      . "\t" . 
						$f->end        . "\t" . 
						$f->strand     . "\t" . 
						$f->score      . "\t" . 
						$f->source_tag . "\t" . 
						$f->primary_tag ."\n\n";

	    my @fset = $self->set_phases($f,$contig);

	    if (defined(@fset) && $f->source_tag eq "genewise") {
		push(@genewise,@fset);

	    } elsif (defined(@fset) && $f->source_tag eq "genscan") {
		push(@genscan,@fset);

	    } elsif (defined(@fset))  {
		print STDERR "Unrecognised feature set : " . $f->source_tag . "\n";

	    } else {
		print STDERR "Problem with validating feature set : skipping\n";
	    }

	} elsif ($f->primary_tag eq "similarity") {
	    
	    $f->seqname  ($contig->id);

	    if ($f->source_tag eq "est2genome") {
#		push(@extras,$f);
	    } elsif ($f->source_tag eq "hmmpfam" && $f->score > 25) {
		push(@features,$f);
            } elsif ($f->score >= 100) {
		push(@features,$f);
            }
	}
    }

    $self->genscan (@genscan);
    $self->genewise(@genewise);
    $self->feature (@features);
    
    print STDERR "\nNumber of similarity features " . scalar(@features) . "\n";
    print STDERR "Number of genscans              " . scalar(@genscan)  . "\n";
    print STDERR "Number of genewise features     " . scalar(@genewise) . "\n\n";

}
    
=head2 make_Exons

 Title   : make_Exons
 Usage   : my @exons = $self->make_Exons(\@features,\@genscan);
 Function: Turns features into exons with the help of genscan predictions
 Example : 
 Returns : Array of Bio::EnsEMBL::Exon
 Args    : array ref of Bio::EnsEMBL::SeqFeature
           array ref of Bio::EnsEMBL::SeqFeature


=cut

sub make_Exons {
    my ($self) = @_;


    $self->make_genscanExons;
    $self->make_genewiseExons;

    print STDERR "\nFinal list of non redundant exons is \n";

    for my $ex ($self->genscan_exons) {
	print STDERR "Genscan exon is " . $ex->id . "\t" . $ex->start . "\t" . $ex->end . "\t" . $ex->strand  . "\n";
    }

    print STDERR "\n";

    for my $ex ($self->genewise_exons) {
	print STDERR "Genewise exon is " . $ex->id . "\t" . $ex->start . "\t" . $ex->end . "\t" . $ex->strand  . "\n";
    }
    
}

sub make_genscanExons {
    my ($self) = @_;

    my @exons;
    my @features = $self->feature;
    my $gscount  = 1;
    
    foreach my $gs ($self->genscan) {
	
	my $excount    = 1;

      EXON: foreach my $f ($gs->sub_SeqFeature) {
	  print STDERR "Looking at " . $f->id . "\t" . $f->start . "\t" . $f->end . "\n";  
	  # Don't include any genscans that are inside a genewise
	  foreach my $gw ($self->genewise) {
	      if (!(($gw->end < $f->start) || $gw->start > $f->end)) {
		  print STDERR "Ignoring genscan exon\n";
		  next EXON;
	      }
	  }
	  
	  my $newexon = $self->_make_Exon($f,$excount,"genscan." . $gscount . "." . $excount );
	  
	                $newexon->find_supporting_evidence(\@features);
	  my @support = $newexon->each_Supporting_Feature;
	  
	  print ("Supporting evidence for " . $newexon->id .  " " . @support . "\n");
	  if ($#support >= 0) {
	      push(@exons,$newexon);
	      $excount++;
	  }
      }
	
	$gscount++;
    }

    $self->genscan_exons(@exons);
}

sub make_genewiseExons {
    my ($self) = @_;

    my @newexons;


    my $gwcount = 1;

    my @fset = $self->genewise;

    @fset = sort {$b->sub_SeqFeature <=> $a->sub_SeqFeature} @fset;

    my @newfset;

     # Remove overlapping genewise sets - this confuses the exon pair building process
     # 
     #FSET:  foreach my $f (@fset) {
     # foreach my $f2 (@newfset) {
#	  if (!($f->start > $f2->end || $f->end < $f2->start)) {
#	      next FSET;
#	  } else {
#	      push(@newfset,$f);
#	  }
#      }
#  }

      foreach my $f (@fset) {
	my @gwexons;

	my $excount = 1;
      SUBF: foreach my $subf ($f->sub_SeqFeature) {

	  my $found   = 0;
	  my $foundex;

	  # Make sure non redundant with all the new exons
	  foreach my $ex (@newexons) {
	      if ($subf->start == $ex->start && $subf->end == $ex->end) {
		  $found = 1;
		  $foundex = $ex;
	      }
	  }

	  # Make sure non redundant with current sub seqfeatures
	  foreach my $ex (@gwexons) {
	      if ($subf->start == $ex->start && $subf->end == $ex->end) {
		  $found = 1;
		  $foundex = $ex;
	      }
	  }
	  
	  if ($found == 1) {
	      print STDERR "\nFound duplicate exon - adding evidence\n";
	      $foundex->add_Supporting_Feature($subf);
	  } else {
	      print STDERR "Making new exon from " .  $subf->gffstring . "\n";
	      
	      my $newexon = $self->_make_Exon($subf,"genewise." . $gwcount . "." . $excount);
	      
	      push(@gwexons,$newexon);
	      
	      $excount++;
	  }
      }

	$gwcount++;


	if ($#gwexons >= 0) {
	    if ($gwexons[0]->strand == 1) {
		@gwexons = sort {$a->start <=> $b->start} @gwexons;
	    } else {
		@gwexons = sort {$b->start <=> $a->start} @gwexons;
	    }
	    push(@newexons,@gwexons);
	}


    }

    $self->genewise_exons(@newexons);
}

sub make_genewise_ExonPairs {
    my ($self) = @_;

    my $count = 0;
    my @exons = $self->genewise_exons;

    foreach my $ex (@exons) {
	eval {
	if ($count > 0) {
	    my $makepair = 0;

	    my $f1 = $exons[$count-1];
	    my $f2 = $exons[$count];

	    my $spliceseq = $f1->{_3splice} . $f2->{_5splice};

	    if ($spliceseq eq "GTAG") {
		$makepair = 1;
	    } elsif ($f1->strand == $f2->strand) {
		if ($f1->strand == 1) {
		    if (abs($f1->end - $f2->start) <= 20) {
			$makepair = 1;
		    }
		} elsif ($f1->strand == -1) {
		    if (abs($f1->start - $f2->end) <= 20) {
			$makepair = 1;
		    }
		} 
	    } 
	    if ($makepair == 1) {
		my $tmppair = new Bio::EnsEMBL::Pipeline::ExonPair(-exon1 => $f1,
								   -exon2 => $f2,
								   -type  => 'ABUTTING',
								   );
		my $found = 0;

		foreach my $p ($self->get_all_ExonPairs) {
		    if ($p->compare($tmppair) == 1) {
			$p->add_coverage;
			$found = 1;
		    }
		} 
		if ($found == 0) {
		    $self->add_ExonPair($tmppair);
		    $tmppair->add_Evidence($f1->each_Supporting_Feature);
		    $tmppair->add_Evidence($f2->each_Supporting_Feature);
		    $tmppair->splice_seq(new Bio::Seq(-id => "splice", 
						      -seq => $f1->{_3splice} . $f2->{_5splice}));
		}
	    }
	}
    };
	
	if ($@) {
	    print STDERR "ERROR making exon pair $@\n";
	}
	$count++;
    }
  
}
	

=head2 make_ExonPairs

 Title   : make_ExonPairs
 Usage   : my @pairs = $self->make_ExonPairs(@exons);
 Function: Links exons with supporting evidence into ExonPairs
 Example : 
 Returns : Array of Bio::EnsEMBL::Exon
 Args    : Array of Bio::EnsEMBL::Pipeline::ExonPair

=cut

sub  make_ExonPairs {
    my ($self) = @_;

    $self->make_genewise_ExonPairs;

    my $gap = 5;

    my %pairhash;

    my @exons = $self->genscan_exons;

    @exons = sort { $a->start <=> $b->start } @exons;

    #  We can cut down the number of searches we do if we
    #  only search forward exons within one contig - I
    #  have yet to implement this.

    print STDERR "Making exon pairs \n";

    for (my $i = 0; $i < scalar(@exons)-1; $i++) {

	print ("Looking at exon $i\n");
	my %idhash;
	my $exon1 = $exons[$i];
	
	my $jstart = $i - 10;  if ($jstart < 0) {$jstart = 0;}
	my $jend   = $i + 10;  if ($jend > scalar(@exons)) {$jend    = scalar(@exons);}

	J: for (my $j = $jstart ; $j < $jend; $j++) {
	    next J if ($i == $j);
            next J if ($exons[$i]->strand != $exons[$j]->strand);
	    next J if ($exons[$i]->id     eq $exons[$j]->id);

	    my $exon2 = $exons[$j];

	    my %doneidhash;

	    print ("EXONS : " . $exon1->id . "\t" . $exon1->start . "\t" . $exon1->end . "\t" . $exon1->strand . "\n");
	    print ("EXONS : " . $exon2->id . "\t" . $exon2->start . "\t" . $exon2->end . "\t" . $exon2->strand . "\n");

	    # For the two exons we compare all of their supporting features.
	    # If any of the supporting features of the two exons
            # span across an intron a pair is made.
	  F1: foreach my $f1 ($exon1->each_Supporting_Feature) {
	      next F1 if (!$f1->isa("Bio::EnsEMBL::FeaturePair"));
	#      print ("FEATURE1 " . $f1->gffstring . "\n");
	      my @f = $exon2->each_Supporting_Feature;
	    F2: foreach my $f2 (@f) {
		next F2 if (!$f2->isa("Bio::EnsEMBL::FeaturePair"));
#		print ("FEATURE2 " . $f1->gffstring . "\n");
		my @pairs = $self->get_all_ExonPairs;
		
		next F1 if (!($f1->isa("Bio::EnsEMBL::FeaturePairI")));
		next F2 if (!($f2->isa("Bio::EnsEMBL::FeaturePairI")));
		
		# Do we have hits from the same sequence
		# n.b. We only allow each database hit to span once
		# across the intron (%idhash) and once the pair coverage between
		# the two exons reaches $minimum_coverage we 
		# stop finding evidence. (%pairhash)
		
		if ($f1->hseqname eq $f2->hseqname &&
		    $f1->strand   == $f2->strand   &&
		    !(defined($idhash{$f1->hseqname}))) {# &&
#		    !(defined($pairhash{$exon1}{$exon2}))) {
							     
		    my $ispair = 0;
		    my $thresh = $self->threshold;

		    if ($f1->strand == 1) {
			if (abs($f2->hstart - $f1->hend) < $gap) {

			    if (!(defined($doneidhash{$f1->hseqname}))) {
				$ispair = 1;
			    }
			}
		    } elsif ($f1->strand == -1) {
			if (abs($f1->hend - $f2->hstart) < $gap) {
			    if (!(defined($doneidhash{$f1->hseqname}))) {
				$ispair = 1;
			    }
			}
		    }
		    
		    # This checks if the coordinates are consistent if the 
		    # exons are on the same contig
		    if ($ispair == 1) {
			if ($exon1->contig_id eq $exon2->contig_id) {
			    if ($f1->strand == 1) {
				if ($f1->end >  $f2->start) {
				    $ispair = 0;
				}
			    } else {
				if ($f2->end >  $f1->start) {
				    $ispair = 0;
				}
			    }
			}
		    }

		    # We finally get to make a pair
		    if ($ispair == 1) {
			eval {
			    my $check = $self->check_link($exon1,$exon2,$f1,$f2);
			    print STDERR "\nPossible pair - checking link - $check ". 
				$exon1->start . "\t" . $exon1->end . "\n";
			    
			    next J unless $check;
			    
			    print STDERR "Making new pair " . $exon1->start . " " . 
				                              $exon1->end   . " " . 
							      $exon2->start . " " . 
							      $exon2->end . "\n";

			    my $pair = $self->makePair($exon1,$exon2,"ABUTTING");
			    
			    if (defined $pair) {
				$idhash    {$f1->hseqname} = 1;
				$doneidhash{$f1->hseqname} = 1;
				
				$pair->add_Evidence($f1);
				$pair->add_Evidence($f2);
				
				if ($pair->is_Covered == 1) {
				    $pairhash{$exon1}{$exon2}  = 1;
				}
			    }
			    
			};
			if ($@) {
			    warn("Error making ExonPair from [" . $exon1->id . "][" .$exon2->id ."] $@");
			}
		    }
		}
	    }
	  }
	}
    }
    return $self->get_all_ExonPairs;
}

sub check_link {
    my ($self,$exon1,$exon2,$f1,$f2) = @_;

    my @pairs = $self->get_all_ExonPairs;
    print STDERR "Checking link for " . $f1->hseqname . " " . $f1->hstart . " " . $f1->hend . " " . $f2->hstart . " " . $f2->hend . "\n";
    foreach my $pair (@pairs) {

	if ($exon1->strand == 1) {
	    if ($exon1 == $pair->exon1) {
		my @linked_features = $pair->get_all_Evidence;
		
		foreach my $f (@linked_features) {
		    
		    if ($f->hseqname eq $f2->hseqname && $f->hstrand == $f2->hstrand) {
			return 0;
		    }
		}
	    }
	} else {
	     if ($exon2 == $pair->exon2) {
		my @linked_features = $pair->get_all_Evidence;
		
		foreach my $f (@linked_features) {
		    
		    if ($f->hseqname eq $f2->hseqname && $f->hstrand == $f2->hstrand) {
			return 0;
		    }
		}
	    }
	 }
    }
    return 1;
}





=head2 link_ExonPairs

 Title   : link_ExonPairs
 Usage   : my @transcript = $self->make_ExonPairs(@exons);
 Function: Links ExonPairs into Transcripts
 Example : 
 Returns : Array of Bio::EnsEMBL::Pipeline::ExonPair
 Args    : Array of Bio::EnsEMBL::Transcript

=cut

sub link_ExonPairs {
    my ($self) = @_;

    my @genewise = $self->genewise_exons;
    my @genscan  = $self->genscan_exons;

    my @exons;

    push(@exons,@genewise);
    push(@exons,@genscan);

  EXON: foreach my $exon (@exons) {
	$self->throw("[$exon] is not a Bio::EnsEMBL::Exon") unless $exon->isa("Bio::EnsEMBL::Exon");

	if ($self->isHead($exon) == 1) {
	    
	    # We have a higher score threshold for single exons
	    # and we need a protein hit

	    if ($self->isTail($exon)) {
		my $found = 0;
		foreach my $f ($exon->each_Supporting_Feature) {
		    if ($f->analysis->db eq "swir" && $f->score > 200) {
			$found = 1;
		    } 
		}
		next EXON unless ($found == 1);

	    }

	    print STDERR "Found new transcript start [" . $exon->id . "]\n";

	    my $transcript = new Bio::EnsEMBL::Transcript;

	    $self      ->add_Transcript($transcript);
	    $transcript->add_Exon       ($exon);

	    $self->_recurseTranscript($exon,$transcript);
	}
    }
    my $count = 1;



    foreach my $tran ($self->get_all_Transcripts) {
	$tran->id ($TRANSCRIPT_ID_SUBSCRIPT . "." . $self->contig->id . "." .$count);
	$tran->version(1);
	$self->make_Translation($tran,$count);
	$count++;
    }
    return $self->get_all_Transcripts;

}


=head2 _recurseTranscript

 Title   : _recurseTranscript
 Usage   : $self->_recurseTranscript($exon,$transcript)
 Function: Follows ExonPairs to form a new transcript
 Example : 
 Returns : nothing
 Args    : Bio::EnsEMBL::Exon Bio::EnsEMBL::Transcript

=cut


sub _recurseTranscript {
    my ($self,$exon,$tran) = @_;
    
    print STDERR "\nRecursing transcript \n";

    if (defined($exon) && defined($tran)) {
	$self->throw("[$exon] is not a Bio::EnsEMBL::Exon")       unless $exon->isa("Bio::EnsEMBL::Exon");
	$self->throw("[$tran] is not a Bio::EnsEMBL::Transcript") unless $tran->isa("Bio::EnsEMBL::Transcript");
    } else {
	$self->throw("Wrong number of arguments [$exon][$tran] to _recurseTranscript");
    }

    # Checks for circular genes here.
    my %exonhash;

    foreach my $exon ($tran->each_Exon) {
	$exonhash{$exon->id}++;
    }

    foreach my $exon (keys %exonhash) {
	if ($exonhash{$exon} > 1) {
	    $self->warn("Eeeek! Found exon " . $exon . " more than once in the same gene. Bailing out");

	    $tran = undef;

	    return;
	}
    }

    # First copy all the exons into a new gene
    my $tmptran = new Bio::EnsEMBL::Transcript;

    foreach my $ex ($tran->each_Exon) {
	$tmptran->add_Exon($ex);
    }

    my $count = 0;

    my @pairs = $self->_getPairs($exon);

    foreach my $pair (@pairs) {
	$self->print_ExonPair($pair);
    }

    my @exons = $tran->each_Exon;;
    
    if ($exons[0]->strand == 1) {
	@exons = sort {$a->start <=> $b->start} @exons;
	
    } else {
	@exons = sort {$b->start <=> $a->start} @exons;
    }

    
    PAIR: foreach my $pair (@pairs) {
#	print STDERR "Comparing " . $exons[$#exons]->id . "\t" . $exons[$#exons]->end_phase . "\t" . 
#	    $pair->exon2->id . "\t" . $pair->exon2->phase . "\n";
#	next PAIR if ($exons[$#exons]->end_phase != $pair->exon2->phase);

	# Only use multiple pairs once
	#if ($#pairs > 0) {
	#    next PAIR if ($self->{_usedPairs}{$pair} == 1);
	#    print ("Already used pair = skipping\n");
	#}

	$self->{_usedPairs}{$pair} = 1;

	if ($count > 0) {
	    my $newtran = new Bio::EnsEMBL::Transcript;
	    $self->add_Transcript($newtran);

	    foreach my $tmpex ($tmptran->each_Exon) {
		$newtran->add_Exon($tmpex);
	    }

	    $newtran->add_Exon($pair->exon2);
	    $self->_recurseTranscript($pair->exon2,$newtran);
	} else {
	    $tran->add_Exon($pair->exon2);
	    $self->_recurseTranscript($pair->exon2,$tran);
	}
	$count++;
    }
}


=head2 add_Transcript

 Title   : add_Transcript
 Usage   : $self->add_Transcript
 Function:  
 Example : 
 Returns : nothing
 Args    : Bio::EnsEMBL::Transcript

=cut

sub add_Transcript {
    my ($self,$transcript) = @_;

    $self->throw("No transcript input") unless defined($transcript);
    $self->throw("Input must be Bio::EnsEMBL::Transcript") unless $transcript->isa("Bio::EnsEMBL::Transcript");

    if (!defined($self->{_transcripts})) {
	$self->{_transcripts} = [];
    }

    push(@{$self->{_transcripts}},$transcript);
}


=head2 get_all_Transcripts

 Title   : get_all_Transcripts
 Usage   : my @tran = $self->get_all_Transcripts
 Function:  
 Example : 
 Returns : @Bio::EnsEMBL::Transcript
 Args    : none

=cut

sub get_all_Transcripts {
    my ($self) = @_;

    if (!defined($self->{_transcripts})) {
	$self->{_transcripts} = [];
    }

    return (@{$self->{_transcripts}});
}

=head2 _getPairs

 Title   : _getPairs
 Usage   : my @pairs = $self->_getPairs($exon)
 Function: Returns an array of all the ExonPairs 
           in which this exon is exon1
 Example : 
 Returns : @Bio::EnsEMBL::Pipeline::ExonPair
 Args    : Bio::EnsEMBL::Exon

=cut

sub _getPairs {
    my ($self,$exon) = @_;

    my $minimum_coverage = 2;
    my @pairs;

    $self->throw("No exon input") unless defined($exon);
    $self->throw("Input must be Bio::EnsEMBL::Exon") unless $exon->isa("Bio::EnsEMBL::Exon");

    foreach my $pair ($self->get_all_ExonPairs) {

	if (($pair->exon1->id eq $exon->id) && ($pair->is_Covered == 1)) {
	    push(@pairs,$pair);
	}
    }

    @pairs = sort { $a->exon2->start <=> $b->exon2->start} @pairs;
    return @pairs;
}
	
	
=head2 isHead

 Title   : isHead
 Usage   : my $foundhead = $self->isHead($exon)
 Function: checks through all ExonPairs to see whether this
           exon is connected to a preceding exon.
 Example : 
 Returns : 0,1
 Args    : Bio::EnsEMBL::Exon

=cut

sub isHead {
    my ($self,$exon) = @_;

    my $minimum_coverage = 2;

    foreach my $pair ($self->get_all_ExonPairs) {

	my $exon2 = $pair->exon2;
	if (($exon->id  eq $exon2->id) && ($pair->is_Covered == 1)) {
	    return 0;
	}
    }

    return 1;
}


=head2 isTail

 Title   : isTail
 Usage   : my $foundtail = $self->isTail($exon)
 Function: checks through all ExonPairs to see whether this
           exon is connected to a following exon.
 Example : 
 Returns : 0,1
 Args    : Bio::EnsEMBL::Exon

=cut

sub isTail {
    my ($self,$exon) = @_;

    my $minimum_coverage = 2;

    foreach my $pair ($self->get_all_ExonPairs) {
	my $exon1 = $pair->exon1;

	if ($exon == $exon1 && $pair->is_Covered == 1) {
	    return 0;
	}
    }
    
    return 1;
}


=head2 make_Genes

 Title   : make_Genes
 Usage   : my @genes = $self->make_Genes(@transcript)
 Function: Turns a set of transcripts into an array of genes
           Transcripts with shared exons go into the same gene
 Example : 
 Returns : Array of Bio::EnsEMBL::Gene
 Args    : Array of Bio::EnsEMBL::Transcript

=cut

sub make_Genes {
    my ($self) = @_;
    
    my @genes;
    
    my $trancount = 1;
    my $genecount = 1;
    my $contigid = $self->contig->id;

    my @transcripts = $self->get_all_Transcripts;

    foreach my $tran (@transcripts) {
	$trancount++;
	
	my $found = undef;

      GENE: foreach my $gene (@genes) {

	  EXON: foreach my $gene_exon ($gene->each_unique_Exon) {

	      foreach my $exon ($tran->each_Exon) {
		  next EXON if ($exon->contig_id ne $gene_exon->contig_id);
		
		  if ($exon->overlaps($gene_exon) && $exon->strand == $gene_exon->strand) {
#		      $self->print_Exon($exon);
#		      $self->print_Exon($gene_exon);
		      $found = $gene;
		      last GENE;
		  }
	      }
	  }
	}

	my $time = time; chomp($time);

	if (defined($found)) {
	    $found->add_Transcript($tran);
	} else {
	    my $gene = new Bio::EnsEMBL::Gene;
	    my $geneid = "TMPG_$contigid.$genecount";
	    $gene->id($geneid);
	    $gene->created($time);
	    $gene->modified($time);
	    $gene->version(1);
	    $genecount++;
	    $gene->add_Transcript($tran);
	    push(@genes,$gene);
	}
    }
    foreach my $gene (@genes) {
        $self->prune_Exons($gene);
	$self->add_Gene($gene);
    }
}

sub add_Gene {
    my ($self,$gene) = @_;

    if (!defined($self->{_genes})) {
	$self->{_genes} = [];
    }
    push(@{$self->{_genes}},$gene);
}

sub each_Gene {
    my ($self) = @_;

    if (!defined($self->{_genes})) {
	$self->{_genes} = [];
    }

    return (@{$self->{_genes}});
}
sub prune_Exons {
    my ($self,$gene) = @_;

    my @unique_Exons; 

    foreach my $tran ($gene->each_Transcript) {
       my @newexons;

       foreach my $exon ($tran->each_Exon) {
           my $found;
           foreach my $uni (@unique_Exons) {
              if ($uni->start  == $exon->start &&
                  $uni->end    == $exon->end   &&
                  $uni->strand == $exon->strand ) {
                  $found = $uni;
                  last $uni;
              }
           }
           if (defined($found)) {
              push(@newexons,$found);
           } else {
              push(@newexons,$exon);
           }
                 
         }          
      $tran->flush_Exon;
      foreach my $exon (@newexons) {
         $tran->add_Exon($exon);
      }
   }
}
sub make_id_hash {
    my ($self,@features) = @_;

    my %id;

    foreach my $f (@features) {
	print STDERR "Feautre $f\n";
	if (!defined($id{$f->hseqname})) {
	    $id{$f->hseqname} = [];
	}
	push(@{$id{$f->hseqname}},$f);
    }

    return \%id;
}




sub translate_Exon {
    my ($self,$exon) = @_;

    my $dna = $exon->seq;
    
    print ("Tran 0 " . $dna->translate('*','X',0)->seq . "\n");
    print ("Tran 1 " . $dna->translate('*','X',1)->seq . "\n");
    print ("Tran 2 " . $dna->translate('*','X',2)->seq . "\n");
}

sub make_Translation{
    my ($self,$transcript,$count) = @_;

    my $translation = new Bio::EnsEMBL::Translation;

    my @exons = $transcript->each_Exon;
    my $exon  = $exons[0];

    $translation->id("TMPP_" . $exon->contig_id . "." . $count);
    $translation->version(1);
#    print("Starting exon is " . $exon->id . "\t" . $exon->start . "\t" . $exon->end . "\t" . $exon->strand . "\n");

    if ($exon->phase != 0) {
	my $tmpphase = $exon->phase;
#	print ("Old coords are " . $exon-> start . "\t" . $exon->end . "\t" . $exon->phase . "\t" . $exon->end_phase . "\n");

	print("Starting phase is not 0 " . $tmpphase . "\t" . $exon->strand ."\n");
	if ($exon->strand == 1) {
	    my $tmpstart = $exon->start;
	    $exon->start($tmpstart + 3 - $tmpphase);
	    $exon->phase(0);
	} else {
	    my $tmpend= $exon->end;
	    $exon->end($tmpend - 3 + $tmpphase);
	    print ("New start end " . $exon->start . "\t" . $exon->end . "\n");
	    $exon->phase(0);
	}
	print ("New coords are " . $exon-> start . "\t" . $exon->end . "\t" . $exon->phase . "\t" . $exon->end_phase . "\n");
    }   


    print ("Transcript strand is " . $exons[0]->strand . "\n");

    if ($exons[0]->strand == 1) {
	@exons = sort {$a->start <=> $b->start} @exons;
	$translation->start        ($exons[0]->start);
	$translation->end          ($exons[$#exons]->end);
	
    } else {
	@exons = sort {$b->start <=> $a->start} @exons;
	print("Start exon is " . $exons[0]->id . "\n");
	$translation->start        ($exons[0]->end);
	$translation->end          ($exons[$#exons]->start);
	
    }
    
    $translation->start_exon_id($exons[0]->id);
    $translation->end_exon_id  ($exons[$#exons]->id);
    
    my $endphase = $exons[0]->end_phase
;
    print ("Making translation " . $endphase . "\n");

    if ($exons[0]->phase != 0) {
	print STDERR "ERROR: start exon phase is not 0";
    }
    shift @exons;

#    foreach my $exon (@exons) {
#	$exon->phase         ($endphase);
#	$endphase = $exon->end_phase;
	
#    }

    $transcript->translation($translation);
}   

sub threshold {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_threshold} = $arg;
    }

    return $self->{_threshold} || 100;
}


sub set_ExonEnds {
    my ($self,$exon) = @_;

    my @genscan = @{$self->{_pred}};
    my $contig  = $self->contig;

    # find genscan ends if poss
    my $leftend;
    my $rightend;

    my $gap = 10;
    my $fover;

    foreach my $gs (@genscan) {

	foreach my $genscan ($gs->sub_SeqFeature) {

	    if ($genscan->overlaps($exon)) {
		print (STDERR "Found overlap : " . $genscan->id . "\t" . $exon->id . "\n");

		$fover = $genscan;

		if (abs($genscan->start - $exon->start) < $gap) {
		    $leftend = $genscan->start;
		}
		
		if (abs($genscan->end - $exon->end) < $gap) {
		    $rightend = $genscan->end;
		}
		
	    }
	}
    }
    return unless defined($fover);
    # else find average ends
    if (!defined($leftend)) {
	
	my %lefthash;

	foreach my $f ($exon->each_Supporting_Feature) {
	    $lefthash{$f->start}++;
	}

	my $maxcount = 0;
	my $maxpos   = 0;

	foreach my $start (keys %lefthash) {
	    if ($lefthash{$start} >= $maxcount)  {
		$maxcount = $lefthash{$start};
		$maxpos   = $start;
	    }
	}
    } 
       

    if (!defined($rightend)) {
	my %righthash;

	foreach my $f ($exon->each_Supporting_Feature) {
	    $righthash{$f->end}++;
	}

	my $maxcount = 0;
	my $maxpos   = 0;

	foreach my $end (keys %righthash) {
	    if ($righthash{$end} >= $maxcount) {
		$maxcount = $righthash{$end};
		$maxpos   = $end;
	    }
	}
    }

    $leftend  = $fover->start;
    $rightend = $fover->end;

    $exon->phase($fover->{phase});

    print(STDERR "New exon start end = $leftend\t$rightend\t" . $exon->translate->seq . "\n");

    $exon->start($leftend);
    $exon->end  ($rightend);

    return;

    my $seq = $exon->seq;

    my $prim = new Bio::PrimarySeq(-id => $exon->id ,
				   -seq => $seq->seq);

    my %stophash;

    my $tran1 = $prim->translate('*','X',0);
    if ($tran1->seq !~ /\*/) {
	$self->add_ExonFrame($exon,1);
	print ("found 1\n");
    }

    my $tran2 = $prim->translate('*','X',2);
    if ($tran2->seq !~ /\*/) {
	$self->add_ExonFrame($exon,2);
	print ("found 2\n");
    }

    my $tran3 = $prim->translate('*','X',1);
    if ($tran3->seq !~ /\*/) {
	$self->add_ExonFrame($exon,3);
	print ("found 3\n");
    }

    my $trancount = 0;
    my $framehash = $self->each_ExonFrame($exon);
    $trancount = scalar(keys(%{$framehash}));
    print ("trancount $trancount for " . $exon->id . "\n");
    print(STDERR "translation 1 " . $tran1->seq . "\n");
    print(STDERR "translation 2 " . $tran2->seq . "\n");
    print(STDERR "translation 3 " . $tran3->seq . "\n");


    if ($trancount == 0) {
	print STDERR "ERROR: No translatable frame\n";
    } 
	

}

sub find_frame_from_evidence {
    my ($self,$exon) = @_;
}

sub check_ExonPair {
    my ($self,$pair) = @_;

    my $exon1  = $pair->exon1;
    my $exon2  = $pair->exon2;

    my $frame1 = $self->each_ExonFrame($exon1);
    my $frame2 = $self->each_ExonFrame($exon2);

    my $trans1 = $pair->exon1->translate();
    my $trans2 = $pair->exon2->translate();

    my $splice1;
    my $splice2;

    my $spliceseq;

    if ($pair->exon1->strand == 1) {
	$splice1 = $exon1->{_gsf_seq}->subseq($exon1->end+1,$exon1->end+2);
	$splice2 = $exon2->{_gsf_seq}->subseq($exon2->start-2,$exon2->start-1);
	$spliceseq = new Bio::Seq(-id => "splice",-seq => "$splice1$splice2");
    } else {
	$splice1 = $exon1->{_gsf_seq}->subseq($exon1->start-2,$exon1->start-1);
	$splice2 = $exon2->{_gsf_seq}->subseq($exon2->end+1,$exon2->end+2);
	$spliceseq = new Bio::Seq(-id => "splice",-seq => "$splice2$splice1");
	$spliceseq = $spliceseq->revcom;
    }

    $pair->splice_seq($spliceseq);

    print (STDERR "Splice " . $spliceseq->seq ."\n");

    return 1;
    return 0 if ($spliceseq ->seq ne "GTAG");
    return 1 if ($pair->exon1->end_phase == $pair->exon2->phase);
    return 0;

    my $match  = 0;
    my $oldphase1 = $exon1->phase;
    my $oldphase2 = $exon2->phase;
    
    foreach my $frame (keys %$frame1) {
	
	$exon1->phase($frame-1);
	my $endphase = $exon1->end_phase;
	print (STDERR "Looking for exon2 phase of " . $endphase+1 . "\n");

	if ($frame2->{$endphase+1} == 1) {
	    print STDERR "Hooray! Found matching phases\n";
	    $match = 1;
	    $exon2->phase($endphase);

	    my $trans1 = $exon1->seq->translate('*','X',(3-$exon1->phase)%3)->seq;
	    my $trans2 = $exon2->seq->translate('*','X',(3-$exon2->phase)%3)->seq;

	    print(STDERR "exon 1 " . $exon1->id . " translation $frame : " . $trans1. "\n");
	    print(STDERR "exon 2 " . $exon2->id . " translation " . ($endphase+1) . " : " . $trans2. "\n");

	    if ($self->add_ExonPhase($exon1) && $self->add_ExonPhase($exon2)) {
		return $match;
	    } else {
		$exon1->phase($oldphase1);
		$exon2->phase($oldphase2);
	    }
	} else {
	    $exon1->phase($oldphase1);
	}
    }




    return $match;
}

sub find_PairSpliceSites {
    my ($self,$pair) = @_;
}

sub find_ExonSpliceSites {
    my ($self,$exon) = @_;

}

sub find_FrameSplices {
    my ($self,$exon,$frame);
}


sub add_ExonFrame {
    my ($self,$exon,$frame) = @_;

    $self->{_framehash}{$exon}{$frame} = 1;
}

sub each_ExonFrame {
    my ($self,$exon) = @_;

    return $self->{_framehash}{$exon};
}


sub add_ExonPhase {
    my ($self,$exon) = @_;

    if (defined($self->{_exonphase}{$exon})) {
	print STDERR "Already defined phase : old phase " . $self->{_exonphase}{$exon} . " new " . $exon->phase . "\n";
	if ($self->{_exonphase}{$exon} != $exon->phase) {
	    return 0;
	}
    } else {
	$self->{_exonphase}{$exon} = $exon->phase;
	return 1;
    }


}

sub set_phases {
    my ($self,$fset) = @_;

    my $strand;
    my $mrna     = "";
    
    $fset->attach_seq($self->contig->primary_seq);

    my @subf = $fset->sub_SeqFeature;
    my @nrf;

    if ($#subf >= 0) {
	$strand = $subf[0]->strand;
	
	if ($strand == 1) {
	    @subf = sort { $a->start <=> $b->start } @subf;
	} else {
	    @subf = sort { $b->start <=> $a->start } @subf;
	}
    }

    foreach my $ex (@subf) {

	my $found = 0;
	foreach my $nex (@nrf) {
	    if ($ex->start == $nex->start && $ex->end == $nex->end) {
		$found = 1;
	    }
	}
	
	if ($found == 1) {
	#    print STDERR "Duplicate sub feature  - skipping\n";
	} else {
	    
	    if ($ex->strand == 1) {
		my $splice3 = $self->contig->primary_seq->subseq($ex->end+1,$ex->end+2);
		my $splice5 = $self->contig->primary_seq->subseq($ex->start-2,$ex->start-1);

		$ex->{_3splice} = $splice3;
		$ex->{_5splice} = $splice5;

	    } else {
		my $splice3 = $self->contig->primary_seq->subseq($ex->start-2,$ex->start-1);
		my $splice5 = $self->contig->primary_seq->subseq($ex->end+1  ,$ex->end+2);

		$splice3 = new Bio::Seq(-id => "splice",-seq => "$splice3");
		$splice5 = new Bio::Seq(-id => "splice",-seq => "$splice5");

		$ex->{_3splice} = $splice3->revcom->seq;
		$ex->{_5splice} = $splice5->revcom->seq;

	    }

	    push(@nrf,$ex);
	}
    }

    my @fset;

    my $count = 0;
    my $sf    = new Bio::EnsEMBL::SeqFeature;

    foreach my $f (@nrf) {
	if ($count > 0) {
	    my $spliceseq = $nrf[$count-1]{_3splice} . $nrf[$count]{_5splice};

	    if ($spliceseq ne "GTAG") {
		push(@fset,$sf);
		$sf = new Bio::EnsEMBL::SeqFeature;
		$sf->add_sub_SeqFeature($f,'EXPAND');
	    } else {
		$sf->add_sub_SeqFeature($f,'EXPAND');
	    }
	} else {
	    $sf->add_sub_SeqFeature($f,'EXPAND');
	}
	$count++;
    }

    push(@fset,$sf);
    
    my @newfset;

    foreach my $fset (@fset) {
	my $phase = $self->find_phase($fset->sub_SeqFeature);
	
	if ($phase == -1) {
	    $self->warn("Couldn't find correct phase for feature set - ignoring feature set\n");
	} else {
	    push(@newfset,$fset);
	}


	foreach my $f ($fset->sub_SeqFeature) {

	    $f->{phase} = $phase;

	    my $len   = $f->end() - $f->start() + 1;
	    my $left_overhang = (3 - $phase)%3;
	    
	    $phase = ($len - $left_overhang)%3;
	    
	    $f->{'end_phase'} = $phase;
	}
	
	foreach my $ex ($fset->sub_SeqFeature) {
	 
	    print STDERR "\tSub feature " . $ex->id          . "\t" .
		$ex->start       . "\t" . 
		$ex->end         . "\t" . 
		$ex->strand      . "\t" . 
		$ex->{phase}     . "\t" . 
		$ex->{end_phase} . "\t" . 
		$ex->{_5splice}   . "\t" . 
		$ex->{_3splice}   . "\n";	
	}
    }
	
    return @newfset;
}

# This is a somewhat kludgy attempt to deal with gene repeats
# If, out of a set of transcripts, one transcript encompasses both
# the start and end of other transcripts it is suspected of leaping 
# between gene repeats and is deleted.
#
# I suspect this will introduce extra fragmentation but it should deal with 
# some of the pathological cases

sub filter_Transcripts {
    my ($self,@transcripts) = @_;

    my @new;

    TRAN: foreach my $tran1 (@transcripts) {
	my $foundstart = 0;
	my $foundend   = 0;

	TRAN2: foreach my $tran2 (@transcripts) {
	    
	    next TRAN2 if ($tran1 == $tran2);
	    
	    if (($tran2->first_exon->start > $tran1->first_exon->start) && 
		($tran2->first_exon->start < $tran1->last_exon->end)) {
		$foundstart = 1;
	    }
	    if (($tran2->last_exon->end > $tran1->first_exon->start) && 
		($tran2->last_exon->end < $tran1->last_exon->end)) {
		$foundend = 1;
	    }
	}
	if ($foundstart == 0 || $foundend == 0) {
	    push(@new,$tran1);
	}
    }
    return @new;
}


sub  genscan {
    my ($self,@genscan) = @_;

    if (!defined($self->{_genscan})) {
        $self->{_genscan} = [];
    }

    if (defined(@genscan)) {
	push(@{$self->{_genscan}},@genscan);
    }

    return @{$self->{_genscan}};
}

sub  feature {
    my ($self,@features) = @_;

    if (!defined($self->{_feature})) {
        $self->{_feature} = [];
    }

    if (defined(@features)) {
	push(@{$self->{_feature}},@features);
    }

    return @{$self->{_feature}};
}

sub  genewise {
    my ($self,@genewise) = @_;

    if (!defined($self->{_genewise})) {
        $self->{_genewise} = [];
    }

    if (defined(@genewise)) {
	push(@{$self->{_genewise}},@genewise);
    }

    return @{$self->{_genewise}};
}

sub genewise_exons {
    my ($self,@exons) = @_;

    if (!defined($self->{_genewise_exons})) {
	$self->{_genewise_exons} = [];
    }
    if (defined(@exons)) {
	push(@{$self->{_genewise_exons}},@exons);
    }
    return @{$self->{_genewise_exons}};
}

sub genscan_exons {
    my ($self,@exons) = @_;

    if (!defined($self->{_genscan_exons})) {
	$self->{_genscan_exons} = [];
    }
    if (defined(@exons)) {
	push(@{$self->{_genscan_exons}},@exons);
    }
    return @{$self->{_genscan_exons}};
}

sub get_all_Features {
    my ($self) = @_;


    if (!defined($self->{_all_Features})) {
	my @tmp = $self->contig->get_all_SimilarityFeatures;
	print STDERR "Got features " . @tmp . "\n";
	$self->{_all_Features} = [];
	push(@{$self->{_all_Features}},@tmp);
    }

    return @{$self->{_all_Features}};
}

=head2 contig

 Title   : 
 Usage   : $self->contig($contig);
 Function: 
 Example : 
 Returns : nothing
 Args    : Bio::EnsEMBL::DB::ContigI
           


=cut

sub contig {
    my ($self,$contig) = @_;
    
    if (defined($contig)) {
	if ($contig->isa("Bio::EnsEMBL::DB::ContigI")) {
	    $self->{_contig} = $contig;
	} else {
	    $self->throw("[$contig] is not a Bio::EnsEMBL::DB::ContigI");

	}
    }
    return $self->{_contig};
}

sub _make_Exon { 
    my ($self,$subf,$stub) = @_;

    my $contigid = $self->contig->id;
    my $exon     = new Bio::EnsEMBL::MappedExon;
    my $time     = time; chomp($time);		

    $exon->id       ("TMPE_" . $contigid . "." . $subf->id . "." . $stub);
    $exon->seqname  ($exon->id);
    $exon->contig_id($contigid);
    $exon->start    ($subf->start);
    $exon->end      ($subf->end  );
    $exon->strand   ($subf->strand);
    $exon->phase    ($subf->{phase});
    $exon->created  ($time);
    $exon->modified ($time);
    $exon->version  (1);
    $exon->attach_seq($self->contig->primary_seq);
    $exon->add_Supporting_Feature($subf);
    
    $exon->{_5splice} = $subf->{_5splice};
    $exon->{_3splice} = $subf->{_3splice};

    return $exon;
}

sub find_phase {
    my ($self,@exons) = @_;

    my $seq = "";
    
    foreach my $exon (@exons) {
	$seq .= $exon->seq->seq;
    }

    my $dna = new Bio::PrimarySeq(-id => "cdna" ,-seq => $seq);
	 
    my $tran0 =  $dna->translate('*','X',0)->seq;
    my $tran1 =  $dna->translate('*','X',2)->seq;
    my $tran2 =  $dna->translate('*','X',1)->seq;

    my $phase;

    if ($tran0 !~ /\*/) {
	$phase = 0;
    } elsif ($tran1 !~ /\*/) {
	$phase = 1;
    } elsif ($tran2 !~ /\*/) {
	$phase = 2;
    } else {
	$phase = -1;
    }
    
    if ($phase == -1) {
	print STDERR "Couldn't find correct phase for transcript. Translations were :\n";

	print STDERR "Phase 0 : $tran0\n";
	print STDERR "Phase 1 : $tran0\n";
	print STDERR "Phase 2 : $tran0\n";
    } else {
	print(STDERR "\n\tTranslation ($phase) =  " . $dna->translate('*','X',(3-$phase)%3)->seq . " " . $phase . "\n\n");
    }

    return $phase;
}

=head2 makePair

 Title   : makePair
 Usage   : my $pair = $self->makePair($exon1,$exon2)
 Function:  
 Example : 
 Returns : Bio::EnsEMBL::Pipeline::ExonPair
 Args    : Bio::EnsEMBL::Exon,Bio::EnsEMBL::Exon

=cut

sub makePair {
    my ($self,$exon1,$exon2,$type) = @_;

    if (!defined($exon1) || !defined($exon2)) {
	$self->throw("Wrong number of arguments [$exon1][$exon2] to makePair");
    }

    $self->throw("[$exon1] is not a Bio::EnsEMBL::Exon") unless $exon1->isa("Bio::EnsEMBL::Exon");
    $self->throw("[$exon2] is not a Bio::EnsEMBL::Exon") unless $exon2->isa("Bio::EnsEMBL::Exon");

    my $tmppair = new Bio::EnsEMBL::Pipeline::ExonPair(-exon1 => $exon1,
						       -exon2 => $exon2,
						       -type  => $type,
						       );
    my $found = 0;

    foreach my $p ($self->get_all_ExonPairs) {
	if ($p->compare($tmppair) == 1) {
	    $p->add_coverage;
	    $tmppair = $p;
	    $found = 1;
	}
    }

    if ($found == 0 && $self->check_ExonPair($tmppair)) {
	$self->add_ExonPair($tmppair);
	return $tmppair;
    }

    return;

}

=head2 get_all_ExonPairs

 Title   : 
 Usage   : 
 Function: 
 Example : 
 Returns : 
 Args    : 

=cut


sub get_all_ExonPairs {
    my ($self) = @_;

    if (!defined($self->{_exon_pairs})) {
	$self->{_exon_pairs} = [];
    }
    return @{$self->{_exon_pairs}};
}


=head2 add_ExonPair

 Title   : add_ExonPair
 Usage   : 
 Function: 
 Example : 
 Returns : 
 Args    : 

=cut

sub add_ExonPair {
    my ($self,$arg) = @_;


    if (!defined($self->{_exon_pairs})) {
	$self->{_exon_pairs} = [];
    }

    if (defined($arg) && $arg->isa("Bio::EnsEMBL::Pipeline::ExonPair")) {
	push(@{$self->{_exon_pairs}},$arg);
    } else {
	$self->throw("[$arg] is not a Bio::EnsEMBL::Pipeline::ExonPair");
    }
}

#############################################################################
# Printing routines
#############################################################################

sub print_Exon {
    my ($self,$exon) = @_;

    print STDERR $exon->seqname . "\t" . $exon->id . "\t" . $exon->start . "\t" . $exon->end . "\t" . $exon->strand . "\n"; #"\t" . $exon->phase . "\t" . $exon->end_phase ."\n";
}

sub print_ExonPairs {
    my ($self) = @_;

    foreach my $pair ($self->get_all_ExonPairs) {
	$self->print_ExonPair($pair);
    }
}

sub print_ExonPair {
    my ($self,$pair) = @_;

    print(STDERR "\nExon Pair (splice - " . $pair->splice_seq->seq . ")\n");
    $self->print_Exon($pair->exon1);
    $self->print_Exon($pair->exon2);
    foreach my $ev ($pair->get_all_Evidence) {
	print(STDERR "   -  " . $ev->hseqname . "\t" . $ev->hstart . "\t" . $ev->hend . "\t" . $ev->strand . "\n");
    }
}
sub print_Genes {
    my ($self) = @_;

    my @genes = $self->each_Gene;

    foreach my $gene (@genes) {
	print(STDERR "\nNew gene - " . $gene->id . "\n");
	my $contig;

	foreach my $tran ($gene->each_Transcript) {
	    print STDERR "\nTranscript - " . $tran->id . "\n";
	    my $cdna;
	    foreach my $exon ($tran->each_Exon) {
		$contig = $exon->contig_id;
#		$self->print_Exon($exon);
#		$self->translate_Exon($exon);
		$cdna .= $exon->seq->seq;
	    }
	    print STDERR "\nTranslation is " . $contig . " " . $tran->translate->seq . "\n";

#	    my $dna = new Bio::PrimarySeq(-id => "cdna" ,-seq => $cdna);
	    
#	    print ("Tran 0 " . $dna->translate('*','X',0)->seq . "\n");
#	    print ("Tran 1 " . $dna->translate('*','X',1)->seq . "\n");
#	    print ("Tran 2 " . $dna->translate('*','X',2)->seq . "\n");
	}
    }
}

sub print_gff {
    my ($self) = @_;
    
    open (POG,">pog.gff");
    
    foreach my $gene ($self->each_Gene) {
	foreach my $tran ($gene->each_Transcript) {
	    foreach my $exon ($tran->each_Exon) {
		print POG $exon->id . "\tSPAN\texon\t" . 
		    $exon->start . "\t" . $exon->end . "\t100\t" ;
		if ($exon->strand == 1) {
		    print POG "+\t.\t";
		} else {
		    print POG ("-\t.\t");
		}
		print POG $tran->id . "\n";
	    }
	}
    }

    foreach my $f ($self->feature) {
	print POG $f->seqname . "\t" . $f->source_tag . "\tsimilarity\t" .
	    $f->start . "\t" . $f->end . "\t" . $f->score . "\t";
	if ($f->strand == 1) {
	    print POG "+\t.\t";
	} else {
	    print POG ("-\t.\t");
	}
	if (ref($f) =~ "FeaturePair") {
	    print (POG $f->hseqname . "\t" . $f->hstart . "\t" . $f->hend . "\n");
	}

    }
    foreach my $gen ($self->genscan) {
	foreach my $ex ($gen->sub_SeqFeature) {
	    print POG  $ex->id . "\t" . $ex->source_tag . "\texon\t" . 
		$ex->start . "\t" . $ex->end . "\t100\t";
	    if ($ex->strand == 1) {
		print POG "+\t.\t";
	    } else {
		print POG ("-\t.\t");
	    }
	    print POG $gen->id . "\n";
	}
    }
    foreach my $gen ($self->genewise) {
	foreach my $ex ($gen->sub_SeqFeature) {
	    print POG  $ex->id . "\t" . $ex->source_tag . "\texon\t" . 
		$ex->start . "\t" . $ex->end . "\t100\t";
	    if ($ex->strand == 1) {
		print POG "+\t.\t";
	    } else {
		print POG ("-\t.\t");
	    }
	    print POG $gen->id . "\n";
	}
    }

    close(POG);

    open (POG,">pog.fa");
    print POG ">". $self->contig->id . "\n";

    my $seq = $self->contig->primary_seq->seq;
    $seq =~ s/(.{72})/$1\n/g;
    
    print POG $seq . "\n";
    close(POG);
}

1;


























