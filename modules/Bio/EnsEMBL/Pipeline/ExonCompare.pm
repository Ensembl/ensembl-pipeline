
#
# BioPerl module for Bio::EnsEMBL::Pipeline::ExonCompare
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::EoxnCompare - Comparison of old and new Gene/Transcript/Exons

=head1 SYNOPSIS

    # $dbobj is the analysis database
    # $timdb is the tim database implementing get_old_Exons on the clone
    # (in the future $dbobj and $timdb could be the same object)
    # @newexons is a set of exons with temporary ids assigned
    # @mappedexons is the set of exons with olds ids, version numbers and new ids etc...

    ($mapped,$new,$untransfered) = 
             Bio::EnsEMBL::Pipeline::GeneComp::map_temp_Exons_to_real_Exons($dbobj,
									    $timdb,
									    @newexons);


    # $mapped - reference to array of tempexons mapped to their new ids, with versions
    # and modified stamps correctly placed

    # $new - reference to array of new exons, with new ids, versions set to 1 and
    # modified/created time stamps.



=head1 DESCRIPTION

This is a methods bag, not a real object. It deals with mapping exons,
transcript and genes from old versions through to new versions. This
is where calls to get_new_ExonID etc are actually made, and where the
version logic happens. To do the mapping we need to get the old exons
out in the new coordinates (remapping). This currently is hidden
behind the method call get_old_Exons() on a contig object. This call
returns old exon objects in the new coordinates, with the method
->identical_dna set to true or not.


=head1 CONTACT

Ensembl - ensembl-dev@ebi.ac.uk

=head1 APPENDIX

=cut



package Bio::EnsEMBL::Pipeline::ExonCompare;

use strict;
use vars qw(@ISA);
use Bio::Root::RootI;
use Bio::EnsEMBL::Exon;

@ISA = qw(Bio::Root::RootI);



=head2 new

 Title   : new
 Usage   : GeneCompare->new()
 Function: 
 Example : 
 Returns : Reference to an object
 Args    : 

=cut

sub new {
    my ($class, @args) = @_;
    my $self = bless {}, $class;

    @{$self->{'_predictorExons'}} = $self->_rearrange([qw(PREDICTOREXONS)],@args);

    return $self;
}



=head2 setStandardExon

 Title   : setStandardExon
 Usage   : $obj->setStandardExon($newval)
 Function: 
 Example : 
 Returns : 
 Args    : newvalue


=cut

sub setStandardExon {
    my ($self, $exon) = @_;
    
    $self->throw("$exon is not an Exon") unless ($exon->isa('Bio::EnsEMBL::Exon'));
    $self->{'_standardExon'} = $exon;
}


=head2 hasOverlap

 Title   : hasOverlap
 Usage   : $overlap = $obj->hasOverlap()
 Function: Checks if any of the predictor exons have an overlap with the standard exon.
 Example : 
 Returns : 0 if none of the of the predictor exons overlap the standard, otherwise 1
 Args    : 


=cut

sub hasOverlap {
    my ($self) = @_;

    my $exon1 = $self->{'_standardExon'};
    
    foreach my $exon2 (@{$self->{'_predictorExons'}}) {               
         if (($exon2->end   > $exon1->start && $exon2->start < $exon1->end) ||
             ($exon2->start < $exon1->end   && $exon2->end   > $exon1->start)) {                
                 return 1;
        }
    }
    
    return 0;
}



=head2 hasExactOverlap

 Title   : hasExactOverlap
 Usage   : $obj->hasExactOverlap()
 Function: 
 Example : 
 Returns : 1 if one of the predictor exons exactly matches the standard exon, otherwise 0.
 Args    : 


=cut

sub hasExactOverlap {
    my ($self) = @_;
    
    my $exon1 = $self->{'_standardExon'};
   
    foreach my $exon2 (@{$self->{'_predictorExons'}}) { 
        if (($exon2->start == $exon1->start) && ($exon2->end == $exon1->end)) {
            return 1;
        }
    }
    
    return 0;
}



=head2 getOverlaps

 Title   : getOverlaps
 Usage   : $obj->getOverlaps()
 Function: 
 Example : 
 Returns : The number of predictor exons with which the standard exon overlaps.
 Args    : 


=cut

sub getOverlaps {
    my ($self) = @_;
    
    my $overlaps = 0;
    my $exon1 = $self->{'_standardExon'};
    
    foreach my $exon2 (@{$self->{'_predictorExons'}}) { 
         if (($exon2->end   > $exon1->start && $exon2->start < $exon1->end) ||
             ($exon2->start < $exon1->end   && $exon2->end   > $exon1->start)) {
                 $overlaps++;
        }
    }
    
    return $overlaps;
}
   


=head2 getBaseOverlaps

 Title   : getBaseOverlaps
 Usage   : $obj->getBaseOverlaps()
 Function: This just considers the case of the first predictor exon found which has an overlap.
 Example : 
 Returns : The number of bases on the standard exon that have true positive, true negative
            and false positive overlaps with a predictor exon. 
 Args    : 


=cut

sub getBaseOverlaps {
    my ($self) = @_;
        
    my $exon1 = $self->{'_standardExon'};
    
    foreach my $exon2 (@{$self->{'_predictorExons'}}) { 
         if (($exon2->end   > $exon1->start && $exon2->start < $exon1->end) ||
             ($exon2->start < $exon1->end   && $exon2->end   > $exon1->start)) {
     
            # Note that the end should always be greater than the start because
            # if it's not the values are swapped when a new Exon is created.            
            my $truePositive = $exon1->end - $exon1->start;
            my $falsePositive = 0;
            my $falseNegative = 0;
                         
            my $leftEnd = $exon1->start - $exon2->start;
            if ($leftEnd > 0) {
                $falsePositive += $leftEnd;
            } else {
                $falseNegative += -$leftEnd;
                $truePositive += $leftEnd;
            }
            
            my $rightEnd = $exon2->end - $exon1->end;
            if ($rightEnd > 0) {
                $falsePositive += $rightEnd;
            } else {
                $falseNegative += -$rightEnd;
                $truePositive += $rightEnd;
            }   
            
            return ($truePositive, $falsePositive, $falseNegative);         
        }
    }
    
    # If no overlaping exon is found there are no true or false positives and
    # the whole length of the exon is false negative
    return (0, 0, $exon1->length());
}


1;
