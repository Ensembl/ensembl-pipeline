
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

For gene/transcript migration, the following things happen. Old and
New genes are clustered into 4 sets on the basis of shared exons (this
occurs after exon mapping, done outside of this module)

There is the possibility of >1 old gene with >1 new gene. Depending on
the order of discovery, this will be classified as a split or
merge. This is a known bug/feature.

For each cluster, old transcripts are sorted by length and then fitted
to new transcripts, with the best fit taking a win (fit on the number
of co-linear identical id'd exons). Perfect matches (all exons the
same id) trigger a direct assignment.

Versioning for transcripts is that any addition/removal of an exon, or
any update in sequence of an exon rolls up the transcript version. The
gene version clicks up on any transcript version or any transcript
addition/deletion.

This is thick code. It is less thick than tims code. There are comments,
but probably not enough. I cant see many other areas for subroutines...


=head1 CONTACT

Ensembl - ensembl-dev@ebi.ac.uk

=head1 APPENDIX

=cut


package Bio::EnsEMBL::Pipeline::GeneCompare;

use strict;
use Bio::EnsEMBL::Pipeline::GeneCompare;

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

    my (@predictorGenes) = $self->_rearrange([qw(
					  PREDICTORGENES
					  )],@args);
                                          
    # Note that if different genes have the same exon @predictorGenes will
    # include multiple copies of this exon but this is precisely what we want!
    my @predictorExons = ();
    foreach $gene (@predictorGenes) {
        @predictorExons = @predictorExons, [$gene->each_unique_Exon()];
    }
    
    $self->{'_predictorExons'} = @predictorExons;


    return $self;
}


=head2 setStandardGene

 Title   : setStandardGene
 Usage   : $obj->setStandardGene($newval)
 Function: 
 Example : 
 Returns : 
 Args    : newvalue


=cut

sub setStandardGene {
    my ($self, $standardGene) = @_;
    
    $self->{'_standardExons'} = [$standardGene->each_unique_Exon()];
}


=head2 isMissed

 Title   : isMissed
 Usage   : $missed = $obj->isMissed()
 Function: A gene is considered missed if none of
            its exons are overlapped by a predicted gene.
 Example : 
 Returns : 1 or 0
 Args    : None


=cut

sub isMissed {
    my ($self) = @_;

    my $comparer = new Bio::EnsEMBL::Pipeline::ExonCompare(-predictorExons => $self->{'_predictorExons'});

    my $missed = 1;
    foreach $standardExon (@{$self->{'_standardExons'}}) {
        $comparer->setStandardExon($standardExon);            
        if ($comparer->isOverlapped()) {
            $missed = 0;
            break;
        }
    }
    
    return $missed;
}


=head2 getOverlaps

 Title   : getOverlaps
 Usage   : $obj->getOverlaps()
 Function: 
 Example : 
 Returns : The number of predictor genes that overlap the standard gene. 
 Args    : None


=cut


sub getOverlaps {
    my ($self) = @_;
    
    my $comparer = new Bio::EnsEMBL::Pipeline::ExonCompare(-predictorExons => $self->{'_predictorExons'});
    
    my $overlaps = 0;
    foreach $standardExon (@{$self->{'_standardExons'}}) {
        $comparer->setStandardExon($standardExon);            
        $overlaps += $comparer->getOverlaps();
    }
    
    return $overlaps;
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
     # write implementation here
    }
    return $obj->{'_geneSensitivity'};
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
     # write implementation here
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
     # write implementation here
    }
    return $obj->{'_exonSensitivity'};
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
    my ($obj) = @_;
    unless ($obj->{'_missedExonScore'}) {
     # write implementation here
    }
    return $obj->{'_missedExonScore'};
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
    my ($obj) = @_;
    unless ($obj->{'_wrongExonScore'}) {
     # write implementation here
    }
    return $obj->{'_wrongExonScore'};
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
     # write implementation here
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
     # write implementation here
    }
    return $obj->{'_baseSensitivity'};
}

1;
