
#
# BioPerl module for Bio::EnsEMBL::Pipeline::GeneCompare
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::GeneCompare - Comparison of old and new Gene/Transcript/Exons

=head1 SYNOPSIS

    # $dbobj is the analysis database
    # $timdb is the tim database implementing get_old_Exons on the clone
    # (in the future $dbobj and $timdb could be the same object)
    # @newexons is a set of exons with temporary ids assigned
    # @mappedexons is the set of exons with olds ids, version numbers and new ids etc...

    # this module *does not* deal with writing the new, mapped, exons into the database.


    ($mapped,$new,$untransfered) = 
             Bio::EnsEMBL::Pipeline::GeneComp::map_temp_Exons_to_real_Exons($dbobj,
									    $timdb,
									    @newexons);


    # $mapped - reference to array of tempexons mapped to their new ids, with versions
    # and modified stamps correctly placed

    # $new - reference to array of new exons, with new ids, versions set to 1 and
    # modified/created time stamps.

    # $untransfered - reference to array of old exons, with old ids which although were 
    # remapped did not have exons in the new database to map to.


=head1 DESCRIPTION

This is a methods bag, not a real object. It deals with mapping exons,
transcript and genes from old versions through to new versions. This
is where calls to get_new_ExonID etc are actually made, and where the
version logic happens. To do the mapping we need to get the old exons
out in the new coordinates (remapping). This currently is hidden
behind the method call get_old_Exons() on a contig object. This call
returns old exon objects in the new coordinates, with the method
->identical_dna set to true or not.

This module is complex. I dont think there is anyway around
this. There are two basic pieces of logic - rules for exon migration
and rules for gene/transcript migration.

For exon migration, if the start/end/strand in new coordinates of the
exons are the same then it gets the old exon id. If the dna sequence
has changed, this increments the version number.  If not it stays the
same.

=head1 CONTACT

Ensembl - ensembl-dev@ebi.ac.uk

=head1 APPENDIX

=cut


package Bio::EnsEMBL::Pipeline::GeneCompare;

use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::Pipeline::ExonCompare;
use Bio::Root::RootI;

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

    my @predictorGenes = $self->_rearrange([qw(PREDICTORGENES)],@args);
                                          
    # Note that if different genes have the same exon, @predictorGenes will
    # include multiple copies of this exon but this is precisely what we want!
    foreach my $gene (@predictorGenes) {
        foreach my $exon ($gene->each_unique_Exon() ) {
            push(@{$self->{'_predictorExons'}}, $exon);
        }
    }   

    return $self;
}



=head2 setStandardGene

 Title   : setStandardGene
 Usage   : $obj->setStandardGene($standardGene)
 Function: 
 Example : 
 Returns : 
 Args    : $standardGene - reference to a Gene object


=cut

sub setStandardGene {
    my ($self, $standardGene) = @_;
    
    $self->{'_standardExons'} = [$standardGene->each_unique_Exon()];
}



=head2 isMissed

 Title   : isMissed
 Usage   : $missed = $obj->isMissed()
 Function: A gene is considered missed if none of its exons are overlapped by a predicted gene.
 Example : 
 Returns : 1 or 0
 Args    : None


=cut

sub isMissed {
    my ($self) = @_;

    my $comparer = new Bio::EnsEMBL::Pipeline::ExonCompare(-predictorExons => @{$self->{'_predictorExons'}});

    foreach my $standardExon (@{$self->{'_standardExons'}}) {
        $comparer->setStandardExon($standardExon);            
        if ($comparer->hasOverlap()) {
            return 0;
        }
    }
    
    return 1;
}



=head2 isExactlyMatched

 Title   : isExactlyMatched
 Usage   : $obj->isExactlyMatched()
 Function: A gene is considered exactly matched if all the exons 
            on the standard gene are exactly identified
 Example : 
 Returns : 1 or 0 
 Args    : None


=cut

sub isExactlyMatched {
    my ($self) = @_;
       
    my $comparer = new Bio::EnsEMBL::Pipeline::ExonCompare(-predictorExons => @{$self->{'_predictorExons'}});

    foreach my $standardExon (@{$self->{'_standardExons'}}) {
        $comparer->setStandardExon($standardExon);
        unless ($comparer->hasExactOverlap() == 1) {
            return 0;
        } 
    }
    
    return 1;
}



=head2 getOverlaps

 Title   : getOverlaps
 Usage   : $obj->getOverlaps()
 Function: Calculates the number of standard exons that are overlapped by predictor exons. 
 Example : 
 Returns : int
 Args    : None


=cut

sub getOverlaps {
    my ($self) = @_;
    
    my $comparer = new Bio::EnsEMBL::Pipeline::ExonCompare(-predictorExons => @{$self->{'_predictorExons'}});
    
    my $overlaps = 0;
    foreach my $standardExon (@{$self->{'_standardExons'}}) {
        $comparer->setStandardExon($standardExon);            
        $overlaps += $comparer->getOverlaps();
    }
    
    return $overlaps;
}



=head2 getExactOverlapRatio

 Title   : getExactOverlaps
 Usage   : $obj->getExactOverlapRatio()
 Function: Calculates the number of standard exons that are exactly overlapped and the number that are not.
 Example : 
 Returns : (int, int) - $overlaps, $nonOverlaps
 Args    : None


=cut

sub getExactOverlapRatio {
    my ($self) = @_;

    my $comparer = new Bio::EnsEMBL::Pipeline::ExonCompare(-predictorExons => @{$self->{'_predictorExons'}});
    
    my $overlaps = 0;
    my $nonOverlaps = 0;
    foreach my $standardExon (@{$self->{'_standardExons'}}) {
    
        $comparer->setStandardExon($standardExon);            
        if ($comparer->hasExactOverlap()) {
            $overlaps++;
        } else {
            $nonOverlaps++;
        }
    }
    
    return ($overlaps, $nonOverlaps);
}



=head2 getMissed

 Title   : getMissed
 Usage   : $missed = $obj->getMissed()
 Function: 
 Example : 
 Returns : The number of standard exons which are not overlapped by predictor exons.
 Args    : None


=cut

sub getMissed {
    my ($self) = @_;

    my $comparer = new Bio::EnsEMBL::Pipeline::ExonCompare(-predictorExons => @{$self->{'_predictorExons'}});

    my $missed = 0;
    foreach my $standardExon (@{$self->{'_standardExons'}}) {
        $comparer->setStandardExon($standardExon);            
        unless ($comparer->hasOverlap()) {
            $missed++;
        }
    }
    
    return $missed;
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
        
    my $comparer = new Bio::EnsEMBL::Pipeline::ExonCompare(-predictorExons => @{$self->{'_predictorExons'}});
    
    my $truePositive = 0;
    my $falsePositive = 0;
    my $falseNegative = 0;
    
    foreach my $standardExon (@{$self->{'_standardExons'}}) {
        $comparer->setStandardExon($standardExon);            
        my ($tP, $fP, $fN) = $comparer->getBaseOverlaps();
  
        $truePositive += $tP;
        $falsePositive += $fP;
        $falseNegative += $fN;
    }  
        
    return ($truePositive, $falsePositive, $falseNegative);         
}

1;
