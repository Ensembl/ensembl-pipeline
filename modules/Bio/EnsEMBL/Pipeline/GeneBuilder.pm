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

sub get_all_Features {
    my ($self) = @_;


    if (!defined($self->{_all_Features})) {
	my @tmp = $self->contig->get_all_SimilarityFeatures;

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

    print STDERR "Buildling genes\n";

    my ($features,$genscan,$repeats) = $self->get_Features;
    
    my @exons               = $self->make_Exons    ($features,$genscan);
    my @pairs               = $self->make_ExonPairs(@exons);
    
    $self->print_ExonPairs;
    
    my @transcripts         = $self->link_ExonPairs(@exons);
    my @contigoverlaps      = $self->order_Contigs (@transcripts);
    my @genes               = $self->make_Genes    (@transcripts);

    push(@{$self->{_genes}},@genes);

    $self->print_Genes(@genes);
    $self->print_gff  (\@genes,$features,$genscan,$repeats);

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

    if (defined($self->{_features}) && defined($self->{_genscan})) {
           return ($self->{_features},$self->{_genscan});
    }
    
    my @features;
    my @genscan;


    my $contig = $self->contig;

    # This is just working with all the similarity features here
    # When EST2genome features are in we need to be more selective
    
    my @tmp     = $self->get_all_Features;

    my @extras;

    foreach my $f (@tmp) {

	if ($f->source_tag eq "genscan") {
	    my $newgs = $self->set_genscan_phases($f,$contig);

	    if (defined $newgs) {
		push(@genscan,$newgs);
	    } else {
		print STDERR "Problem with genscan prediction : skipping\n";
	    }

	} elsif ($f->primary_tag eq "similarity") {
	    $f->seqname($contig->id);
	    if ($f->source_tag eq "vert_est2genome") {
		push(@extras,$f);
	    } elsif ($f->source_tag eq "hmmpfam" && $f->score > 25) {
		push(@features,$f);
            } elsif ($f->score > 100) {
		push(@features,$f);
            }
	}
    }

#    foreach my $x (@extras) {
#	print ($x->gffstring . "\n");
#    }

    my @merged = $self->merge_features(@extras);


#    foreach my $x (@merged) {
#	print ($x->gffstring . "\n");
#    }

    $self->predictions(@genscan);

    $self->{_features} = [];
    $self->{_genscan}  = [];
    $self->{_extras}   = [];

    push(@{$self->{_features}},@features);
    push(@{$self->{_features}},@merged);
    push(@{$self->{_merged}},@merged);

    push(@{$self->{_genscan}},@genscan);
    
    print STDERR "Number of similarity features " . scalar(@features) . "\n";
    print STDERR "Number of genscans "            . scalar(@genscan)  . "\n";
#    print STDERR "Number of repeats "            . scalar(@repeats)  . "\n";
    return (\@features,\@genscan);
}

sub predictions {
    my ($self,@pred) = @_;

    if (defined(@pred)) {
	$self->{_pred} = [];
	push(@{$self->{_pred}},@pred);
    }
    if (!defined($self->{_pred})) {
        $self->{_pred} = [];
    }
    return @{$self->{_pred}};
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
    my ($self,$features,$genscan) = @_;

    my @exons;
    my @supported_exons;


    # ********** TEST ONLY **************
    # To try and replicate what TimDB is doing I verify genscan exons here
    
    my $gscount = 1;

    foreach my $gs (@$genscan) {

	my $excount    = 1;
	my $contigid   = $gs->seqname;
	my $contig     = $self->contig;

	foreach my $f ($gs->sub_SeqFeature) {

	    my $exon  = new Bio::EnsEMBL::MappedExon;
	    $exon->id       ("TMPE_" . $contigid . ".$gscount.$excount");
	    $exon->seqname  ($contigid . ".$gscount.$excount");
	    $exon->contig_id($contigid);
	    $exon->start    ($f->start);
	    $exon->end      ($f->end  );
	    $exon->strand   ($f->strand);
	    $exon->attach_seq($contig->primary_seq);
	    push(@exons,$exon);
	    $excount++;
	    
	}
	$gscount++;
    }
    
#    my @extraexons = $self->find_extra_exons(@exons,@{$self->{_merged}});
#    my $excount = 1;
    
#    foreach my $ex (@extraexons) {
#	my $contig = $self->contig;
#	my $contigid = $contig->id;
	
#	print (STDERR "Extra exon " . $ex->id . " "  . $ex->start . " " . $ex->end . "\n");

#	my $exon  = new Bio::EnsEMBL::MappedExon;

#	$exon->id       ("TMPE_" . $contigid . ".extra.$excount");
#	$exon->seqname  ($contigid . ".extra.$excount");
#	$exon->contig_id($contigid);
#	$exon->start    ($ex->start);
#	$exon->end      ($ex->end  );
#	$exon->strand   ($ex->strand);
#	$exon->attach_seq($self->contig->primary_seq);
#	push(@exons,$exon);
#	$excount++;
#    }


    foreach my $ex (@exons) {
	$ex->contig_id($self->contig->id);
	$ex->find_supporting_evidence($features);
	my @support = $ex->each_Supporting_Feature;
	
	foreach my $s (@support) {
	    $self->add_toExonHash($s,$ex);
	}

	my $time = time; chomp($time);
	if ($#support >= 0) {
	    push(@supported_exons,$ex);
	    $ex->version(1);
	    $ex->created($time);
	    $ex->modified($time);
	    $self->set_ExonEnds($ex);
	}
    }

    return @supported_exons;

}

sub find_extra_exons {
    my ($self,@exons,@features) = @_;

    my @newexons;

    foreach my $f (@features) {
	my $found = 1;
	print (STDERR "Feature = " . $f->gffstring . "\n");

	foreach my $ex (@exons) {
	    print (STDERR "Comparing to " . $ex->start . " " . $ex->end . " " . $ex->strand . "\n");
	    if ($ex->overlaps($f)) {
		$found = 0;
	    }
	}
	if ($found == 1) {
	    push(@newexons,$f);
	}
    }
    return @newexons;
}


sub add_toExonHash {
    my ($self,$f,$exon) = @_;

    if (!defined($self->{_exonhash}{$f->hseqname})) {
	$self->{_exonhash}{$f->hseqname} = [];
    }

    push(@{$self->{_exonhash}{$f->hseqname}},$exon);
}

sub get_ExonHash {
    my ($self,$hid) = @_;

    return $self->{_exonhash}{$hid};
}

=head2 cluster_Features

 Title   : cluster_Features
 Usage   : my @features = $self->cluster_Features(@features)
 Function: clusters overlapping features into one feature
 Example : 
 Returns : Array of Bio::EnsEMBL::SeqFeature
 Args    : Array of Bio::EnsEMBL::SeqFeature

=cut

sub  cluster_Features {
    my ($self,@features) = @_;

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
    my ($self,@exons) = @_;

    my $gap = 5;

    my %pairhash;

    @exons = sort { $a->start <=> $b->start } @exons;

    #  We can cut down the number of searches we do if we
    #  only search forward exons within one contig - I
    #  have yet to implement this.

    print STDERR "Making exon pairs\n";

    for (my $i = 0; $i < scalar(@exons)-1; $i++) {

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

#	    print ("EXONS : " . $exon1->id . "\t" . $exon1->start . "\t" . $exon1->end . "\t" . $exon1->strand . "\n");
#	    print ("EXONS : " . $exon2->id . "\t" . $exon2->start . "\t" . $exon2->end . "\t" . $exon2->strand . "\n");

	    # For the two exons we compare all of their supporting features.
	    # If any of the supporting features of the two exons
            # span across an intron a pair is made.
	  F1: foreach my $f1 ($exon1->each_Supporting_Feature) {
#	      print ("FEATURE1 " . $f1->gffstring . "\n");
	    F2: foreach my $f2 ($exon2->each_Supporting_Feature) {

#		print ("FEATURE2 " . $f1->gffstring . "\n");
		my @pairs = $self->get_all_ExonPairs;

		next F1 if (!($f1->isa("Bio::EnsEMBL::FeaturePairI")));
		next F2 if (!($f2->isa("Bio::EnsEMBL::FeaturePairI")));

#		print ("pog\n");
		# Do we have hits from the same sequence
		# n.b. We only allow each database hit to span once
		# across the intron (%idhash) and once the pair coverage between
		# the two exons reaches $minimum_coverage we 
		# stop finding evidence. (%pairhash)

		if ($f1->hseqname eq $f2->hseqname &&
		    $f1->strand   == $f2->strand   &&

		    
		    !(defined($idhash{$f1->hseqname})) &&
		    !(defined($pairhash{$exon1}{$exon2}))) {
		    
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
#			    print STDERR "Making new pair\n";
			    my $pair = $self->makePair($exon1,$exon2,"ABUTTING");
			    
			    if (defined $pair) {
				$idhash    {$f1->hseqname} = 1;
				$doneidhash{$f1->hseqname} = 1;
				
				$pair->add_Evidence($f1);
				$pair->add_Evidence($f2);
				
				if ($pair->coverage > 1) {
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

sub find_ExonSet {
    my ($self,$exon) = $_;

    my %set;
    
    foreach my $f ($exon->each_SupportingFeature) {
	my @ex = $self->get_ExonHash($f->hseqname);
	
	foreach my $exon (@ex) {
	    $set{$exon->id} = $exon;
	}
    }

    my @out = values %set;

    return @out;
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
	print STDERR "Adding new pair\n";
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



=head2 link_ExonPairs

 Title   : link_ExonPairs
 Usage   : my @transcript = $self->make_ExonPairs(@exons);
 Function: Links ExonPairs into Transcripts
 Example : 
 Returns : Array of Bio::EnsEMBL::Pipeline::ExonPair
 Args    : Array of Bio::EnsEMBL::Transcript

=cut

sub link_ExonPairs {
    my ($self,@exons) = @_;

    if (!defined(@exons)) {
	$self->warn("No array of exons input to link_ExonPairs");
	return $self->get_all_Transcripts;
    }
    my $contigid = $self->contig->id;

    my @newexons;
    foreach my $exon (@exons) {
	if (defined($exon->phase)) {
	    push(@newexons,$exon);
	}
    }

    
  EXON: foreach my $exon (@newexons) {
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
	$tran->id ($TRANSCRIPT_ID_SUBSCRIPT . $contigid.$count);
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

    my @exons = $tran->each_Exon;;
    
    if ($exons[0]->strand == 1) {
	@exons = sort {$a->start <=> $b->start} @exons;
	
    } else {
	@exons = sort {$b->start <=> $a->start} @exons;
    }


    PAIR: foreach my $pair (@pairs) {
	print STDERR "Comparing " . $exons[$#exons]->id . "\t" . $exons[$#exons]->end_phase . "\t" . 
	    $pair->exon2->id . "\t" . $pair->exon2->phase . "\n";
	next PAIR if ($exons[$#exons]->end_phase != $pair->exon2->phase);

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

	if (($pair->exon1->id eq $exon->id) && ($pair->coverage >= $minimum_coverage)) {
	    push(@pairs,$pair);
	}
    }

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

	if (($exon->id  eq $exon2->id) && ($pair->coverage >= $minimum_coverage)) {
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

	if ($exon == $exon1 && $pair->coverage >= $minimum_coverage) {
	    return 0;
	}
    }
    
    return 1;
}

=head2 order_Contigs

 Title   : order_Contigs
 Usage   : my @contigs = $self->order_Contigs(@transcript);
 Function: Orders contigs based on exon order in transcripts
 Example : 
 Returns : Array of Bio::EnsEMBL::Transcript
 Args    : Array of Bio::EnsEMBL::DB::ContigI

=cut

sub order_Contigs {
    my ($self,@contigs) = @_;

    return @contigs;
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
    my ($self,@transcripts) = @_;
    
    my @genes;
    
    my $trancount = 1;
    my $genecount = 1;
    my $contigid = $self->contig->id;

    foreach my $tran (@transcripts) {
	$trancount++;
	
	my $found = undef;

      GENE: foreach my $gene (@genes) {

	  EXON: foreach my $gene_exon ($gene->each_unique_Exon) {

	      foreach my $exon ($tran->each_Exon) {
		  next EXON if ($exon->contig_id ne $gene_exon->contig_id);
		
		  if ($exon->overlaps($gene_exon) && $exon->strand == $gene_exon->strand) {
		      $self->print_Exon($exon);
		      $self->print_Exon($gene_exon);
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

    return @genes;
}

sub make_id_hash {
    my ($self,@features) = @_;

    my %id;

    foreach my $f (@features) {
	if (!defined($id{$f->hseqname})) {
	    $id{$f->hseqname} = [];
	}
	push(@{$id{$f->hseqname}},$f);
    }

    return \%id;
}

sub merge_features {
    my ($self,@features) = @_;

    my $idhash = $self->make_id_hash(@features);

    my @newfeatures;

    foreach my $hid (keys %$idhash) {
	my @f = @{$idhash->{$hid}};
	my @mergef = $self->merge_feature_array(@f);

	push(@newfeatures,@mergef);
    }

    return @newfeatures;
}

sub merge_feature_array {
    my ($self,@features) = @_;
    
    my $gap = 5;
    my @newfeatures;
    
    @features = sort { $a->start <=> $b->start } @features;

    FEAT: foreach my $f (@features) {
	my $found = 0;

	foreach my $newf (@newfeatures) {

	    if ($newf->strand == $f->strand) {

		if ($f->strand == 1) {
		    if (abs($f->start  - $newf->end) <= $gap &&
			abs($f->hstart - $newf->hend) <= $gap) {

			$newf->end($f->end);
			$newf->hend($f->hend);
			$found = 1;
			next FEAT;
		    }
		} elsif ($f->strand == -1) {
		    if (abs($f->start - $newf->end) <= $gap && 
			abs($f->hend - $newf->hstart) <= $gap) {
			$newf->end($f->end);
			$newf->hstart($f->hstart);
			$found = 1;
			next FEAT;
		    }
		}
	    }
	}
	if ($found == 0) {
	    my $newf = new Bio::EnsEMBL::FeaturePair(-feature1 => $f->feature1,
						     -feature2 => $f->feature2);
	    push(@newfeatures,$newf);
	}
    }
    return @newfeatures;
}



sub print_Exon {
    my ($self,$exon) = @_;

    print STDERR $exon->seqname . "\t" . $exon->id . "\t" . $exon->start . "\t" . $exon->end . "\t" . $exon->strand . "\t" . $exon->phase . "\t" . $exon->end_phase ."\n";
}

sub print_ExonPairs {
    my ($self) = @_;

    foreach my $pair ($self->get_all_ExonPairs) {
	print(STDERR "\nExon Pair (splice - " . $pair->splice_seq->seq . ")\n");
	$self->print_Exon($pair->exon1);
	$self->print_Exon($pair->exon2);
	foreach my $ev ($pair->get_all_Evidence) {
	    print(STDERR "   -  " . $ev->hseqname . "\t" . $ev->hstart . "\t" . $ev->hend . "\t" . $ev->strand . "\n");
	}
	
    }
}

sub print_Genes {
    my ($self,@genes) = @_;

    foreach my $gene (@genes) {
	print(STDERR "\nNew gene - " . $gene->id . "\n");
	my $contig;

	foreach my $tran ($gene->each_Transcript) {
	    print STDERR "\nTranscript - " . $tran->id . "\n";
	    my $cdna;
	    foreach my $exon ($tran->each_Exon) {
		$contig = $exon->contig_id;
		$self->print_Exon($exon);
		$self->translate_Exon($exon);
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
    my ($self,$genes,$features,$genscan,$repeats) = @_;
    
    open (POG,">pog.gff");
    
    foreach my $gene (@$genes) {
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

    foreach my $f (@$features) {
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
    foreach my $f (@$repeats) {
	print POG $f->seqname . "\t" . $f->source_tag . "\tRepeatMasker\t" .
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

    foreach my $gen (@$genscan) {
	foreach my $ex ($gen->sub_SeqFeature) {
	    print POG  $ex->id . "\tgenscan\texon\t" . 
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
    print("Starting exon is " . $exon->id . "\t" . $exon->start . "\t" . $exon->end . "\t" . $exon->strand . "\n");

    if ($exon->phase != 0) {
	my $tmpphase = $exon->phase;
	print ("Old coords are " . $exon-> start . "\t" . $exon->end . "\t" . $exon->phase . "\t" . $exon->end_phase . "\n");

#	print ("Old translation is " . $exon->seq->translate('*','X',0)->seq . "\n");
#	print ("Old translation is " . $exon->seq->translate('*','X',1)->seq . "\n");
#	print ("Old translation is " . $exon->seq->translate('*','X',2)->seq . "\n");
	
	print("Starting phase is not 0 " . $tmpphase . "\t" . $exon->strand ."\n");
	if ($exon->strand == 1) {
	    my $tmpstart = $exon->start;
	    $exon->start($tmpstart + 3 - $tmpphase);
	    $exon->phase(0);
#	    print ("New translation is " . $exon->seq->translate('*','X',0)->seq . "\n");
#	    print ("New translation is " . $exon->seq->translate('*','X',1)->seq . "\n");
#	    print ("New translation is " . $exon->seq->translate('*','X',2)->seq . "\n");
	} else {
	    my $tmpend= $exon->end;
	    $exon->end($tmpend - 3 + $tmpphase);
	    print ("New start end " . $exon->start . "\t" . $exon->end . "\n");
	    $exon->phase(0);
#	    print ("New translation is " . $exon->seq->translate('*','X',0)->seq . "\n");
#	    print ("New translation is " . $exon->seq->translate('*','X',1)->seq . "\n");
#	    print ("New translation is " . $exon->seq->translate('*','X',2)->seq . "\n");
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

    foreach my $exon (@exons) {
	$exon->phase         ($endphase);
	$endphase = $exon->end_phase;
	
    }

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

sub set_genscan_phases {
    my ($self,$gs,$contig) = @_;
    my $strand;
    
    $gs->attach_seq($contig->primary_seq);
    print STDERR "Got contig primary_seq\n";
    my $mrna = "";

    foreach my $f ($gs->sub_SeqFeature) {
       print STDERR "Fetching seq for " . $f->id . " "  . $f->start  . " " . $f->end . "\n";
	$mrna .= $f->seq->seq;
	$strand = $f->strand;
    }
    
    my $dnaseq = new Bio::PrimarySeq(-id  => "mrna",
				     -seq => $mrna);
    
    my $phase;
    if ($dnaseq->translate('*','X',0)->seq !~ /\*/) {
	$phase = 0;
    } elsif ($dnaseq->translate('*','X',2)->seq !~ /\*/) {
	$phase = 1;
    } elsif ($dnaseq->translate('*','X',1)->seq !~ /\*/) {
	$phase = 2;
    } else {
	$self->warn("No translateable frame for genscan prediction".$strand.":".$gs->seqname);
	return;
    }


    print(STDERR "genscan translation = " . $gs->id . " " . $dnaseq->translate('*','X',(3-$phase)%3)->seq . "\n");
    foreach my $f ($gs->sub_SeqFeature) {
	print(STDERR "Genscan exon ids " . $f->id . ":". $f->seq->seq . "\n");
	$f->{phase} = $phase;
	my $len   = $f->end() - $f->start() + 1;
	my $left_overhang = (3 - $phase)%3;

	$phase = ($len - $left_overhang)%3;

	$f->{'end_phase'} = $phase;

	my $splice1;
	my $splice2;
	
	my $spliceseq;

	if ($f->strand == 1) {
	$splice1 = $contig->primary_seq->subseq($f->end+1,$f->end+2);
	$splice2 = $contig->primary_seq->subseq($f->start-2,$f->start-1);
	$spliceseq = new Bio::Seq(-id => "splice",-seq => "$splice1$splice2");
    } else {
	$splice1 = $contig->primary_seq->subseq($f->start-2,$f->start-1);
	$splice2 = $contig->primary_seq->subseq($f->end+1,$f->end+2);
	$spliceseq = new Bio::Seq(-id => "splice",-seq => "$splice2$splice1");
	$spliceseq = $spliceseq->revcom;
    }

    print (STDERR "Genscan splices : " . $spliceseq->seq ."\n");

	
    }

    return $gs;
}

1;


























