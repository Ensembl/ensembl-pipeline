
#
# BioPerl module for Bio::EnsEMBL::Pipeline::GeneComparisonStats
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::GeneComparisonStats - Comparison of old and new Gene/Transcript/Exons

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



package Bio::EnsEMBL::Pipeline::GeneComparisonStats;

use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::DB::CloneI;
use Bio::Root::RootI;
use Bio::EnsEMBL::Pipeline::GeneCompare;

@ISA = qw(Bio::Root::RootI);

=head2 new

 Title   : new
 Usage   : GeneComparisonStats->new()
 Function: 
 Example : 
 Returns : Reference to an object
 Args    : 

=cut

sub new {
    my ($class, @args) = @_;
    my $self = bless {}, $class;

    my ($standard, $predictor) = $self->_rearrange([qw(STANDARD PREDICTOR)],@args);

    ($standard && $standard->isa('Bio::EnsEMBL::DB::CloneI')) 
        || $self->throw("GeneComparisonStats requires a standard object which implements CloneI");
    ($predictor && $predictor->isa('Bio::EnsEMBL::DB::CloneI')) 
        || $self->throw("GeneComparisonStats requires a predictor object which implements CloneI");

    @{$self->{'_standardGenes'}} = $standard->get_all_Genes();
    @{$self->{'_predictorGenes'}} = $predictor->get_all_Genes();

    return $self;
}



=head2 getGeneSpecificity

 Title   : getGeneSpecificity
 Usage   : $obj->getGeneSpecificity()
 Function: 
 Example : 
 Returns : The specificity at which the predictor is able to correctly identify
            and asemble all of a gene's exons. 
 Args    : 

=cut

sub getGeneSpecificity {
    my ($obj) = @_;
    
    unless ($obj->{'_geneSpecificity'}) {
        $obj->_genePredictions();       
    }
    
    return $obj->{'_geneSpecificity'};
}



=head2 getGeneSensitivity

 Title   : getGeneSensitivity
 Usage   : $obj->getGeneSensitivity()
 Function: 
 Example : 
 Returns : The sensitivity at which the predictor is able to correctly identify
            and asemble all of a gene's exons. 
 Args    : 

=cut

sub getGeneSensitivity {
    my ($obj) = @_;
    
    unless ($obj->{'_geneSensitivity'}) {
        $obj->_genePredictions();
    }
    
    return $obj->{'_geneSensitivity'};
}



=head2 _genePredictions

 Title   : _genePredictions
 Usage   : $obj->_genePredictions()
 Function: Calculates the specificity and sensitivity at which the predictor is able to correctly identify
            and asemble all of a gene's exons. They are both calculated at the same time because almost
            certainly they will both be required and half the data for each calculation is the same.
            If all the exons of a standard gene are identified and every intron-exon boundary is 
            correct, i.e. each exon has an exact overlap, the gene is a true positive; otherwise
            a false negative. The false positives are then found separately by counting any predictor genes which 
            are missing from the standard set.
 Example : 
 Returns : Nothing
 Args    : None

=cut

sub _genePredictions {
    my ($self) = @_;
    
    my $truePositive = 0;
    my $falsePositive = 0;
    my $falseNegative = 0;    
    my $comparer = new Bio::EnsEMBL::Pipeline::GeneCompare(-predictorGenes => @{$self->{'_predictorGenes'}});
    
    foreach my $standardGene (@{$self->{'_standardGenes'}}) {
        $comparer->setStandardGene($standardGene);
        if ($comparer->isExactlyMatched()) {
            $truePositive++;
        } else {
            $falseNegative++;
        }
    }
    
    $comparer = new Bio::EnsEMBL::Pipeline::GeneCompare(-predictorGenes => @{$self->{'_standardGenes'}});
    
    foreach my $predictorGene (@{$self->{'_predictorGenes'}}) {
        $comparer->setStandardGene($predictorGene);
        $falsePositive += $comparer->isMissed();
    }    
     
    if ($truePositive == 0) {
        $self->{'_geneSpecificity'} = 0;  
        $self->{'_geneSensitivity'} = 0;
    } else {           
        $self->{'_geneSpecificity'} = $truePositive / ($truePositive + $falsePositive);  
        $self->{'_geneSensitivity'} = $truePositive / ($truePositive + $falseNegative); 
    }      
}



=head2 getMissedGeneScore

 Title   : getMissedGeneScore
 Usage   : $obj->getMissedGeneScore()
 Function: A gene is considered missed if none of its exons are overlapped by a predicted gene.
            If necessary this is calculated by a call to _getMissed with each standard gene being 
            compared against all the predictor genes.
 Example : 
 Returns : The frequency at which the predictor completely fails to 
            identify an gene. 
 Args    : 

=cut

sub getMissedGeneScore {
    my ($self) = @_;
    
    unless ($self->{'_missedGeneScore'}) { 
        my @array1 = @{$self->{'_standardGenes'}};
        my @array2 = @{$self->{'_predictorGenes'}};            
        $self->{'_missedGeneScore'} = $self->_getMissed(@array1, @array2, "gene"); 
    }
    
    return $self->{'_missedGeneScore'};
}



=head2 getWrongGeneScore

 Title   : getWrongGeneScore
 Usage   : $obj->getWrongGeneScore()
 Function: A prediction is considered wrong if none of its exons are overlapped by a gene from the standard set.
            If necessary this is calculated by a call to _getMissed with each predictor gene being 
            compared against all the standard genes. 
 Example : 
 Returns : The frequency at which the predictor incorrectly identifies a gene.  
 Args    : 

=cut

sub getWrongGeneScore {
    my ($self) = @_;
    
    unless ($self->{'_wrongGeneScore'}) { 
        my @array1 = @{$self->{'_predictorGenes'}};
        my @array2 = @{$self->{'_standardGenes'}};           
        $self->{'_wrongGeneScore'} = $self->_getMissed(@array1, @array2, "gene");
    }
    
    return $self->{'_wrongGeneScore'};
}



=head2 _getMissed

 Title   : _getMissed
 Usage   : $obj->_getMissed()
 Function: 
 Example : 
 Returns : The frequency at which the genes in array2 completely fails to 
            identify a gene in array1. A gene is considered missed if none of
            its exons are overlapped by a predicted gene.
 Args    : Two arrays of genes

=cut

sub _getMissed {
    my ($self, $array1, $array2, $type) = @_;
    
    my $comparer = new Bio::EnsEMBL::Pipeline::GeneCompare(-predictorGenes => @$array2);
    my $missed = 0;
    my $count = 0;
        
    foreach my $gene (@$array1) {
        $comparer->setStandardGene($gene);
        if ($type eq "gene") {
            $missed += $comparer->isMissed();
        } 
        else {
            $missed += $comparer->getMissed(); 
        }   
        $count++;
    }
            
    return $missed / $count;
}



=head2 getSplitGeneScore

 Title   : getSplitGeneScore
 Usage   : $obj->getSplitGeneScore()
 Function: 
 Example : 
 Returns : The frequency at which the predictor incorrectly splits a gene's 
            exons into multiple genes. A gene from the standard set is 
            considered split if it overlaps more than one predicted gene.
            The score is defined as the sum of the number of predicted genes
            that overlap each standard gene divided by the number of 
            standard genes that were split.
 Args    : None

=cut

sub getSplitGeneScore {
    my ($self) = @_;
    
    unless ($self->{'_splitGeneScore'}) {
        my @array1 = @{$self->{'_standardGenes'}};
        my @array2 = @{$self->{'_predictorGenes'}};        
        $self->{'_splitGeneScore'} = $self->_getMulitipleOverlaps(@array1, @array2);
    }
    
    return $self->{'_splitGeneScore'};
}



=head2 getJoinedGeneScore

 Title   : getJoinedGeneScore
 Usage   : $obj->getJoinedGeneScore()
 Function: 
 Example : 
 Returns : The frequency at which the predictor incorrectly assembles multiple
            genes' exons into a single gene. A predicted gene is considered 
            joined if it overlaps more than one gene in the standard set.
            The score is defined as the sum of the number of standard genes that
            overlap each predicted genes divided by the number of predicted genes
            that were joined.
 Args    : 

=cut

sub getJoinedGeneScore {
    my ($self) = @_;
    
    unless ($self->{'_joinedGeneScore'}) {
        my @array1 = @{$self->{'_predictorGenes'}};
        my @array2 = @{$self->{'_standardGenes'}};     
        $self->{'_joinedGeneScore'} = $self->_getMulitipleOverlaps(@array1, @array2);
    }
    
    return $self->{'_joinedGeneScore'};
}



=head2 _getMulitipleOverlaps

 Title   : _getMulitipleOverlaps
 Usage   : $obj->_getMulitipleOverlaps()
 Function: 
 Example : 
 Returns : The ratio of the number of overlaps between the genes in array1 and array2
            to the number of genes in array2 that overlap more than one gene in array1.
 Args    : 

=cut

sub _getMulitipleOverlaps {
    my ($self, $array1, $array2) = @_;
    
    my $comparer = new Bio::EnsEMBL::Pipeline::GeneCompare(-predictorGenes => @$array2);
    my $multiples = 0;
    my $overlaps = 0;
        
    foreach my $gene (@$array1) {
       $comparer->setStandardGene($gene);
       my $count = $comparer->getOverlaps();
       if ($overlaps > 1) {
           $multiples++;
       }
       $overlaps += $count;
    }
    
    return $overlaps / $multiples;
}


=head2 getExonSpecificity

 Title   : getExonSpecificity
 Usage   : $obj->getExonSpecificity()
 Function: 
 Example : 
 Returns : The specificity at which the predictor identifies exons 
            and correctly recognises their boundaries.
 Args    : 

=cut

sub getExonSpecificity {
    my ($obj) = @_;
    
    unless ($obj->{'_exonSpecificity'}) {
        $obj->_exonPredictions();
    }
    
    return $obj->{'_exonSpecificity'};
}


=head2 getExonSensitivity

 Title   : getExonSensitivity
 Usage   : $obj->getExonSensitivity()
 Function: 
 Example : 
 Returns : The sensitivity at which the predictor identifies exons 
            and correctly recognises their boundaries.
 Args    : 


=cut

sub getExonSensitivity {
    my ($obj) = @_;
    
    unless ($obj->{'_exonSensitivity'}) {
        $obj->_exonPredictions();
    }
    
    return $obj->{'_exonSensitivity'};
}



=head2 _exonPredictions

 Title   : _exonPredictions
 Usage   : $obj->_exonPredictions()
 Function: Calculates the specificity and sensitivity at which the predictor is able to correctly identify
            exons and correctly recognises their boundaries. They are both calculated at the same time because almost
            certainly they will both be required and half the data for each calculation is the same.
            If a standard exon is exactly overlapped by a predictor exon it is a true positive; otherwise
            a false negative. The false positives are then found separately by counting any predictor exons which 
            are missing from the standard set.
 Example : 
 Returns : Nothing
 Args    : None

=cut

sub _exonPredictions {
    my ($self) = @_;
    
    my $truePositive = 0;
    my $falsePositive = 0;
    my $falseNegative = 0;    
    my $comparer = new Bio::EnsEMBL::Pipeline::GeneCompare(-predictorGenes => @{$self->{'_predictorGenes'}});
    
    foreach my $standardGene (@{$self->{'_standardGenes'}}) {
        $comparer->setStandardGene($standardGene);
        my ($tP, $fN) = $comparer->getExactOverlapRatio(); 
        $truePositive += $tP;
        $falseNegative += $fN;        
    }
    
    $comparer = new Bio::EnsEMBL::Pipeline::GeneCompare(-predictorGenes => @{$self->{'_standardGenes'}});
    
    foreach my $predictorGene (@{$self->{'_predictorGenes'}}) {
        $comparer->setStandardGene($predictorGene);
        $falsePositive += $comparer->getMissed(); 
    }    
            
    $self->{'_exonSpecificity'} = $truePositive / ($truePositive + $falsePositive);  
    $self->{'_exonSensitivity'} = $truePositive / ($truePositive + $falseNegative);       
}



=head2 getMissedExonScore

 Title   : getMissedExonScore
 Usage   : $obj->getMissedExonScore()
 Function: 
 Example : 
 Returns : The frequency at which the predictor completely fails to 
            identify an exon (no prediction or overlap). 
 Args    : 


=cut

sub getMissedExonScore {
    my ($self) = @_;
    
    unless ($self->{'_missedExonScore'}) {
        my @array1 = @{$self->{'_standardGenes'}};
        my @array2 = @{$self->{'_predictorGenes'}};            
        $self->{'_missedGeneScore'} = $self->_getMissed(@array1, @array2, "exon"); 
    }
    
    return $self->{'_missedExonScore'};
}



=head2 getWrongExonScore

 Title   : getWrongExonScore
 Usage   : $obj->getWrongExonScore()
 Function: 
 Example : 
 Returns : The frequency at which the predictor identifies an exon
            that has no overlap with any exon from the standard. 
 Args    : 


=cut

sub getWrongExonScore {
    my ($self) = @_;
    
    unless ($self->{'_wrongExonScore'}) {
        my @array1 = @{$self->{'_predictorGenes'}};
        my @array2 = @{$self->{'_standardGenes'}};           
        $self->{'_wrongGeneScore'} = $self->_getMissed(@array1, @array2, "exon");
    }
    
    return $self->{'_wrongExonScore'};
}



=head2 getBaseSpecificity

 Title   : getBaseSpecificity
 Usage   : $obj->getBaseSpecificity()
 Function: 
 Example : 
 Returns : The specificity at which the predictor correctly labels
            a base in the genomic sequence as being part of each gene.
            This rewards predictors that get most of each gene correct
            and penalises those that miss large parts.
 Args    : 


=cut

sub getBaseSpecificity {
    my ($obj) = @_;
    unless ($obj->{'_baseSpecificity'}) {
        $obj->_basePredictions();
    }
    return $obj->{'_baseSpecificity'};
}



=head2 getBaseSensitivity

 Title   : getBaseSensitivity
 Usage   : $obj->getBaseSensitivity()
 Function: 
 Example : 
 Returns : The sensitivity at which the predictor correctly labels
            a base in the genomic sequence as being part of each gene.
            This rewards predictors that get most of each gene correct
            and penalises those that miss large parts.
 Args    : 


=cut

sub getBaseSensitivity {
    my ($obj) = @_;
    unless ($obj->{'_baseSensitivity'}) {
        $obj->_basePredictions();
    }
    return $obj->{'_baseSensitivity'};
}



=head2 _basePredictions

 Title   : _basePredictions
 Usage   : $obj->_basePredictions()
 Function: Calculates the specificity and sensitivity at which the predictor is able to correctly 
            identify exons and correctly recognises their boundaries. They are both calculated at 
            the same time because almost certainly they will both be required and half the data for 
            each calculation is the same. If a standard exon is exactly overlapped by a predictor 
            exon it is a true positive; otherwise a false negative. The false positives are then 
            found separately by counting any predictor exons which are missing from the standard set.
 Example : 
 Returns : Nothing
 Args    : None

=cut

sub _basePredictions {
    my ($self) = @_;
    
    my $truePositive = 0;
    my $falsePositive = 0;
    my $falseNegative = 0; 
             
    my $comparer = new Bio::EnsEMBL::Pipeline::GeneCompare(-predictorGenes => @{$self->{'_predictorGenes'}});
    
    foreach my $standardGene (@{$self->{'_standardGenes'}}) {
        $comparer->setStandardGene($standardGene);
        my ($tP, $fP, $fN) = $comparer->getBaseOverlaps(); 
        $truePositive += $tP;
        $falsePositive += $fP;
        $falseNegative += $fN;        
    } 
            
    $self->{'_baseSpecificity'} = $truePositive / ($truePositive + $falsePositive);  
    $self->{'_baseSensitivity'} = $truePositive / ($truePositive + $falseNegative);       
}

1;
