#
#
# BioPerl module for Bio::EnsEMBL::Pipeline::ExonPair
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Analysis::ExonPair - stores 2 exons surrounding a confirmed intron

=head1 SYNOPSIS

    
    $pair = new Bio::EnsEMBL::Pipeline::ExonPair(-exon1 => $exon1,
						 -exon2 => $exon2,
						 -type  => $type);
    
    $exon1 and $exon2 are Bio::EnsEMBL::Exon;


=head1 DESCRIPTION

This object stores 2 exon objects which surround a confirmed intron.  A set of these pair
objects are generated from similarity matches and gene structures are then generated.

    my $exon1 = $pair->exon1;
    my $exon2 = $pair->exon2; 

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Pipeline::ExonPair;

use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;

use Bio::EnsEMBL::Exon;

@ISA = qw(Bio::Root::Object);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;

  my ($exon1,$exon2,$type) =     $self->_rearrange([qw(EXON1
						       EXON2
						       TYPE
						 )],@args);

  
  $exon1 || $self->throw("First exon not defined in ExonPair constructor");
  $exon2 || $self->throw("Second exon not defined in ExonPair constructor");
  $type  || $self->throw("No type entered for exon pair");

  $self->throw("Exon1 not Bio::EnsEMBL::Exon") unless $exon1->isa("Bio::EnsEMBL::Exon");
  $self->throw("Exon2 not Bio::EnsEMBL::Exon") unless $exon2->isa("Bio::EnsEMBL::Exon");

  $self->exon1($exon1);
  $self->exon2($exon2);
  $self->type ($type);
  
  $self->verifyExons();
  $self->add_coverage();

  return $make; # success - we hope!

}


=head2 exon1

 Title   : exon1
 Usage   : $pair->exon1($exon)
 Function: get/set method for the first exon
 Example : 
 Returns : Bio::EnsEMBL::Exon
 Args    : 

=cut

sub exon1 {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->throw("Argument must be Bio::EnsEMBL::Exon") unless $arg->isa("Bio::EnsEMBL::Exon");
	$self->{_exon1} = $arg;
    }

    return $self->{_exon1};
}


=head2 exon2

 Title   : exon2
 Usage   : $pair->exon2($exon)
 Function: get/set method for the second exon
 Example : 
 Returns : Bio::EnsEMBL::Exon
 Args    : 

=cut

sub exon2 {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->throw("Argument must be Bio::EnsEMBL::Exon") unless $arg->isa("Bio::EnsEMBL::Exon");
	$self->{_exon2} = $arg;
    }
    return $self->{_exon2};
}


=head2 type

 Title   : type
 Usage   : $pair->type($type)
 Function: get/set method for the type of exon pair
 Example : 
 Returns : string
 Args    : string

=cut

sub type {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_type} = $arg;
    }
    return $self->{_type};
}


=head2 verifyExons

 Title   : verifyExons
 Usage   : $pair->verifyExons
 Function: Consistency checks for the 2 exons in the pair
 Example : 
 Returns : nothing
 Args    : 

=cut

sub verifyExons {
    my ($self) = @_;

    $self->throw("No exon1 defined") unless $self->exon1;
    $self->throw("No exon2 defined") unless $self->exon2;

    if ($self->exon1->contig_id eq $self->exon2->contig_id) {
	# Should we check contig version here?
	if ($self->exon1->strand != $self->exon2->strand) {
	    $self->throw("Exons are on opposite strands");
	}

	if ($self->exon1->strand == 1) {
	    if ($self->exon1->end >  $self->exon2->start) {
		$self->throw("Inconsistent coordinates for exons [" . $self->exon1->end . "][" . $self->exon2->start ."]");
	    }
	} else {
	    if ($self->exon2->end >  $self->exon1->start) {
		$self->throw("Inconsistent coordinates for exons (-1)[" . $self->exon2->end . "][" . $self->exon1->start ."]");
	    }
	}
    }
}


=head2 add_coverage

 Title   : add_coverage
 Usage   : $pair->add_coverage
 Function: Increases the coverage level by 1
 Example : 
 Returns : nothing
 Args    : none

=cut

sub add_coverage {
    my ($self) = @_;

    $self->{_coverage}++;
}


=head2 coverage

 Title   : coverage
 Usage   : $pair->coverage
 Function: Returns the coverage level for the pair
 Example : 
 Returns : int
 Args    : none

=cut

sub coverage {
    my ($self) = @_;

    return $self->{_coverage};
}

=head2 is_Covered

 Title   : is_Covered
 Usage   : $pair->is_Covered
 Function: Returns 1 if this pair is validated by enough evidence
 Example : 
 Returns : int
 Args    : none

=cut

sub is_Covered {
    my ($self) = @_;

    my $est = 0;

    for my $f ($self->get_all_Evidence) {
	print STDERR "Evidence " . $f->source_tag . "\n";
	
	if ($f->source_tag eq "est2genome" || $f->source_tag eq "genewise") {
	    $est = 1;
	}
    }
    print STDERR "EST " . $est . " " . $self->coverage . "\n";
    if ($est == 1 || $self->coverage >= 2) {
	return 1;
    } else {
	return 0;
    }
}

=head2 compare

 Title   : compare
 Usage   : $pair->compare($pair2);
 Function: Compares self to another pair to 
           see if they\'re the same
 Example : 
 Returns : 0,1
 Args    : Bio::EnsEMBL::Pipeline::ExonPair

=cut

sub compare {
    my ($self,$pair) = @_;

    if ($self->exon1 == $pair->exon1 && 
	$self->exon2 == $pair->exon2) {
	return 1;
    } else {
	return 0;
    }
}



=head2 add_Evidence

 Title   : add_Evidence
 Usage   : $pair->add_Evidence($f);
 Function: Adds a sequence feature that was
           used to create this pair
 Example : 
 Returns : nothing
 Args    : Bio::EnsEMBL::FeaturePair

=cut

sub add_Evidence {
    my ($self,$feature) = @_;


    if (!defined($self->{_evidence})) {
	$self->{_evidence} = [];
    }

    $self->throw("Argument to add_Evidence must be Bio::EnsEMBL::FeaturePair") unless defined($feature);
    $self->throw("Argument to add_Evidence must be Bio::EnsEMBL::FeaturePair") unless $feature->isa("Bio::EnsEMBL::FeaturePair");

    push(@{$self->{_evidence}},$feature);

}



=head2 get_all_Evidence

 Title   : get_all_Evidence
 Usage   : my @features = $pair->get_all_Evidence;
 Function: Returns all the feature pairs that were
           used to create this pair
 Example : 
 Returns : Array of Bio::EnsEMBL::FeaturePair
 Args    : none

=cut

sub get_all_Evidence {
    my ($self) = @_;


    if (!defined($self->{_evidence})) {
	$self->{_evidence} = [];
    }

    return @{$self->{_evidence}};
}


sub splice_seq {
    my ($self,$seq) = @_;

    if (defined($seq)) {
	$self->{_spliceseq} = $seq;
    }

    return $self->{_spliceseq};
}

1;

