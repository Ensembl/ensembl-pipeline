
# Spangle version of the ensembl Transcript object

# POD documentation - main docs before the code

=head1 NAME

MiniSeq - an artificially constructed cDNA on a genomic sequence

=head1 SYNOPSIS

This module is used when we only want to run an analysis over
a part or multiple parts of a sequence. 

=head1 DESCRIPTION

Contains details of coordinates of all exons that make
up a gene transcript.

Creation:
   
     my $mini = new Bio::EnsEMBL::Pipeline::MiniSeq();
     my $mini = new Bio::EnsEMBL::Pipeline::MiniSeq($id,\@exons,$seq);

Manipulation:

    $mini->add_Exon  ($exon);
    $mini->contig_dna($seq);

    my @exons = $mini->get_all_Exons      # Returns an array of Exon objects
    my $dna   = $mini->dna_seq            # Returns the dna sequence
     
    $tran->sort()                        # Sorts exons into order (forward for + strand, reverse fo - strand)

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::Pipeline::MiniSeq;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object

use Bio::Root::Object;
use Bio::EnsEMBL::Exon;


@ISA = qw(Bio::Root::Object);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;

  $self->{'_exons'} = [];

  my $id         = $args[0];
  my $exons      = $args[1];
  my $contig_dna = $args[2];

  $self->throw("Wrong number of args [@args] for MiniSeq") unless $#args = 2;
  
  $self->id($id);

  if (defined($exons) && ref($exons) eq 'ARRAY') {
      foreach my $exon (@$exons) {
	  $self->add_Exon($exon);
      }
  }

  if (defined($contig_dna)) {
      $self->contig_dna($contig_dna);
  }

  return $self; # success - we hope!
}

=head2 id

 Title   : id
 Usage   : $obj->id($newval)
 Function: 
 Returns : value of id
 Args    : newvalue (optional)


=cut

sub id {
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'id'} = $value;
    }
    return $self->{'id'};

}

=head2 add_Exon

 Title   : add_Exon
 Usage   : $trans->add_Exon($exon)
 Returns : Nothing
 Args    :


=cut

sub add_Exon{
   my ($self,$exon) = @_;

   if( ! $exon->isa("Bio::EnsEMBL::Exon") ) {
       $self->throw("$exon is not a Bio::EnsEMBL::Exon!");
   }
   push(@{$self->{'_exons'}},$exon);
   
}

=head2 get_all_Exons

 Title   : get_all_Exons
 Usage   : foreach $exon ( $trans->get_all_Exons)
 Function: Returns an array of exons in the transcript
 Example : my @exons = $tr->get_all_Exons
 Returns : An array of exon objects
 Args    : none


=cut

sub get_all_Exons {
   my ($self) = @_;
   my @sub;
   my @ret;

   return @{$self->{'_exons'}};
}

=head2 dna_seq

  Title   : dna_seq
  Usage   : $dna = $feat->dna_seq
   Function: Returns the dna sequence of the gene
  Returns : Bio::Seq
  Args    : none

=cut

sub dna_seq {
  my ($self) = @_;

  my $mrna   = "";
  my @exons  = $self->get_all_Exons;
  my $strand = $exons[0]->strand;

  
  foreach my $exon ($self->get_all_Exons) {
      my $tmp = $self->contig_dna->subseq($exon->start,$exon->end);

    if ($strand == -1) {
      $tmp =~ tr/ATGCatgc/TACGtacg/;
      $tmp = reverse($tmp);
    }

    $mrna  .= $tmp;
  }


  my $seq = new Bio::PrimarySeq(-id => $self->id,
				-seq => $mrna);

  return $seq;
}

=head2 contig_dna

  Title   : contig_dna
  Usage   : $tran->contig_dna($dna);
   Function: Sets the dna sequence of the contig
  Returns : Bio::Seq
  Args    : Bio::Seq

=cut

sub contig_dna {
  my ($self,$dna) = @_;

  if (defined($dna)) {
      if ($dna->isa("Bio::PrimarySeq")) {
	  $self->{_contig_dna} = $dna;
      }
  }

  return $self->{_contig_dna};
}

=head2 sort

 Title   : sort
 Usage   : $feat->sort()
 Function: Sorts the exon features by start coordinate
           Sorts forward for forward strand and reverse for reverse strand
 Returns : none
 Args    : none

=cut

sub sort {
  my $self = shift;

  # Fetch all the features
  my @exons = $self->get_all_Exons;

  # Empty the feature table
  $self->{_exons} = [];

  # Now sort the exons and put back in the feature table
  my $strand = $exons[0]->strand;

  if ($strand == 1) {
    @exons = sort { $a->start <=> $b->start } @exons;
  } elsif ($strand == -1) {
    @exons = sort { $b->start <=> $a->start } @exons;
  }

  foreach my $e (@exons) {
    $self->add_Exon($e);
  }
}


sub convert_feature {
    my ($self,$f) = @_;

    if ($f->start > $f->end) {
	my $tmp = $f->end;
	$f->end  ($f->start);
	$f->start($tmp);
    }
    
    my $start = $f->start;
    my $end   = $f->end;
    
    my @out;
    my @exons = $self->get_all_Exons;
    
    if ($exons[0]->strand == 1) {
	@exons = sort {$a->start <=> $b->start} @exons;
	
	my $foundstart = 0;
	my $foundend   = 0;
	
	my $count = 0;
	my $dna   = $exons[0]->start - 1;
	
	foreach my $ex (@exons) {
	    my $tmpstart;
	    my $tmpend;
	    
	    if ($foundstart == 0 ) {
		if ($start >= ($ex->start - $dna) && $start <= ( $ex->end - $dna)) {
		    $foundstart = 1;
		    $tmpstart = $dna + $start;
		}
	    }
	    
	    if (!$foundend && $foundstart) {
		$tmpstart = $ex->start unless $tmpstart;
		$tmpend   = $ex->end;
		
		# Before making a new exon check we haven't also
		# Got the end point in the same exon
		
		if ($end <= $ex->end - $dna ) {
		    $tmpend = $end + $dna;
		    $foundend = 1;
		}
		
	    }
	    
	    if (defined($tmpstart) && defined($tmpend)) {
		my $newex = Bio::EnsEMBL::Exon->new;
		   $newex->start ($tmpstart);
		   $newex->end   ($tmpend);
		push(@out,$newex);
		
		if ($foundend == 1) {
		    return @out;
		}
		
	    }
	    
	    $dna = $dna + $exons[$count+1]->start - $ex->end - 1;
	    $count++;
	}
	return @out; 
    } else {
	
	@exons = sort {$b->start <=> $a->start} @exons;
	
	my $foundstart = 0;
	my $foundend   = 0;
	
	my $count = 0;
	my $dna   = $exons[0]->end;
	

	foreach my $ex (@exons) {
	    my $tmpstart;
	    my $tmpend;
	    
	    if ($foundstart == 0 ) {
		if ($start >= ($dna - $ex->end) && $start <= ($dna-$ex->start)) {
		    $foundstart = 1;
		    $tmpstart = $dna - $start + 1;
		}
	    }
	    
	    if (!$foundend && $foundstart) {
		$tmpstart = $ex->end unless $tmpstart;   ### is this right?
		$tmpend   = $ex->start;
		
		if ($end <= ($dna - $ex->start)) {
		    $tmpend = $dna - $end + 1;
		    $foundend = 1;
		}
		
	    }
	    
	    if (defined($tmpstart) && defined($tmpend)) {
		my $newex = Bio::EnsEMBL::Exon->new;
		$newex->start ($tmpend);
		$newex->end   ($tmpstart);
		push(@out,$newex);
		
		if ($foundend == 1) {
		    return @out;
		}
		
	    }
	    
	    if ($count < $#exons) {
		$dna = $dna - ($ex->start - $exons[$count+1]->end) + 1;
	    }

	    $count++;
	}
	return @out; 
    }
}


sub exon_dna {
  my ($self,$exon) = @_;

  my $tmpseq = $self->contig_dna->str($exon->start,$exon->end);
  
  if ($exon->strand == -1) {
    $tmpseq =~ tr/ATGCatgc/TACGtacg/;
    $tmpseq = reverse($tmpseq);
  }
  return new Bio::Seq(-seq => $tmpseq);
}

sub pep_coords {
  my $self = shift;

  # for mapping the peptide coords back onto the dna sequence
  # it would be handy to have a list of the peptide start end coords
  # for each exon
  
  my @starts;
  my @ends;
  
  my $fullpep = $self->translate()->seq;

  foreach my $ex ($self->get_all_Exons) {
    
    my @tmp = $self->translate_exon($ex);
    my $pep = $tmp[$ex->phase]->seq;
    
    my $start = index($fullpep,$pep) + 1;
    my $end = $start + length($pep) - 1;
    
    push(@starts,$start);
    push(@ends,$end);
    
  }

  return \@starts,\@ends;
}


sub translate_region {
  my ($self,$start,$end) = @_;
  
  my $mrna = "";
  my $count = 0;
  
  # Loop through the exons until we find the start
  # We adjust the start coordinate so we translate
  # in the right frame

  my $foundstart;    # We have found the start of the region to translate
  my $foundend;      # We have found the end of the regino to translate

  my @exons  = $self->get_all_Exons;
  my $strand = $exons[0]->strand;

  # Separate the forward and reverse strands
  if ($strand eq 1) {

    foreach my $exon (@exons) {
      
      my $tmpstart;
      my $tmpend;
      
      # Find the coords in the current exon that we wish to translate
      if (!$foundstart) {
	if ($start >= $exon->start && $start <= $exon->end) {
	  $foundstart = 1;
	  $tmpstart = $start;
	  
	  # Adjust the start coordinate so we are in the
	  # right frame
	  $tmpstart = $tmpstart + (3 - ($tmpstart - $exon->phase - $exon->start)%3) % 3;
	  
	}
      }
      
      if ($foundstart && !$foundend) {
	$tmpend   = $exon->end;
	$tmpstart = $exon->start unless $tmpstart;
	
	# Check to see if we have the end coord as well
	if ($end <= $exon->end) {
	  $foundend = 1;
	  $tmpend = $end;
	} elsif  ($count < $#exons && $end <= $exons[$count+1]->start) {
	  # Or does the end occur in the intron (shouldn't do!!!)
	  $tmpend = $exon->end;
	  $foundend = 1;
	} 
      }

      
      # Only tack on sequence to the mrna if we are in the middle of the translated region
      
      if (defined($tmpstart) && defined($tmpend)) {
	
	my $newstart = $tmpstart  - $exon->start + 1;
	my $newend   = $tmpend    - $exon->start + 1;
	
	my $seq    = $self->exon_dna($exon)->seq();
	my $tmpseq = substr($seq,$newstart-1,($newend-$newstart+1));
	
	$mrna = $mrna . $tmpseq;
      }
      $count++;
      
    }

  } else {
    foreach my $exon (@exons) {
      
      my $tmpstart;
      my $tmpend;
      
      if (!$foundstart) {
	if ($end <= $exon->end && $end >= $exon->start) {
	  $foundstart = 1;
	  $tmpstart = $end;

	  # Adjust the start coordinate so we are in the
	  # right frame

	  $tmpstart = $tmpstart  - (3 - ($exon->end - $end - $exon->phase)%3)%3;
#	  print("Adjusting by " . (3-($exon->end - $end - $exon->phase)%3)%3  . "\n");

	}
      }

      if (!$foundend && $foundstart) {
	$tmpend   = $exon->start;
	$tmpstart = $exon->end  unless $tmpstart;

	if ($start >= $exon->start) {
	  $foundend = 1;
	  $tmpend = $start;
	} elsif  ($count < $#exons && $start >= $exons[$count+1]->end) {
	  # Or does the end occur in the intron (shouldn't do!!!)
	  $tmpend = $exon->start;
	  $foundend = 1;
	} 
	
	if (defined($tmpstart) && defined($tmpend)) {
#	  print("Exon $count\t" . $exon->start . "\t" . $exon->end . "\t" .  $tmpstart . "\t" . $tmpend . "\n");
	  
	  my $newstart = $exon->end - $tmpstart + 1;
	  my $newend   = $exon->end - $tmpend   + 1;
	  
	  my $seq    = $self->exon_dna($exon)->seq();
	  my $tmpseq = substr($seq,$newstart-1,($newend-$newstart+1));
	  
	  $mrna .= $tmpseq;
	}
      }
      $count++;
    }
  }

  
  my $seq = new Bio::Seq(-seq => $mrna);
  my $i = 0;
  my @out;
  
  for ($i=0; $i < 3; $i++) {
    my $subseq = new Bio::Seq(-seq => substr($mrna,$i));
    my $trans = $subseq->translate();
    push(@out,$trans);
  }

  return \@out;
}


sub find_coord {
  my ($self,$coord,$type) = @_;
 
  my $count = 0;
  my @exons = $self->get_all_Exons;
  my $end   = $#exons;
  my $dna;

  my ($starts,$ends) = $self->pep_coords;
  my $strand = $exons[0]->strand;

  # $starts and $ends are array refs containing the _peptide_ coordinates
  # of each exon. We may have 1 missing residue that spans an intron.
  # We ignore these.

  if ($strand == 1) {
    foreach my $ex ($self->get_all_Exons) {
      
      if ($coord >= $starts->[$count] && $coord <= $ends->[$count]) {
	my $dna   = $ex->start + $ex->phase;
	my $nopep = $coord - $starts->[$count];
	
	$dna += 3 * $nopep;

	if ($type eq "end") {
	  $dna += 2;
	}
	
	return $dna;
	
      } elsif ($count < $end) {
	my $endpep = $ends->[$count]+1;
	if ($endpep == $coord) {

	  my $dna;

	  if ($type eq "end") {
	    my $end_phase = $ex->end_phase;
	    $dna = $ex->end - 3 + $end_phase;
	  } else {
	    $dna = $exons[$count+1]->start + $exons[$count+1]->phase;
	  }
	  return $dna;
	}
      }
      $count++;
    }
  } else {

    foreach my $ex ($self->get_all_Exons) {
      
      if ($coord >= $starts->[$count] && $coord <= $ends->[$count]) {
	
	my $dna   = $ex->end - $ex->phase;
	my $nopep = $coord - $starts->[$count];

	$dna -= 3*$nopep;

	if ($type eq "end") {
	  $dna -= 2;
	}
	
	return $dna;
	
      } elsif ($count < $end) {
	my $endpep = $ends->[$count]+1;

	if ($endpep == $coord) {
	  my $dna;

	  if ($type eq "end") {
	    my $end_phase = $ex->end_phase;
	    $dna = $ex->start + 3 - $end_phase;
	  } else {
	    $dna = $exons[$count+1]->end - $exons[$count+1]->phase;
	  }
	  return $dna;
	}
      }
      $count++;
    } 
  }
}



1;







