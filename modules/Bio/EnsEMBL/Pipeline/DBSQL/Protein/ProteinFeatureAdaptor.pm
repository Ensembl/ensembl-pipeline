# Author: Marc Sohrmann (ms2@sanger.ac.uk)
# Copyright (c) Marc Sohrmann, 2001
# You may distribute this code under the same terms as perl itself
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code


=head1 NAME

  Bio::EnsEMBL::Pipeline::DBSQL::ProteinFeatureAdaptor - Object representing an adaptor
                                                         for protein features to an EnsEMBL DB

=head1 SYNOPSIS

  use Bio::EnsEMBL::DBSQL::Protein::DBAdaptor;
  use Bio::EnsEMBL::Pipeline::DBSQL::Protein::ProteinFeatureAdaptor;

  my $db = Bio::EnsEMBL::Pipeline::DBSQL::Protein::DBAdaptor->new ( -host   => 'ics1e',
                                                                    -dbname => 'test',
                                                                    -user   => 'marc',
                                                                    -pass   => 'xyz',
                                                         );               

  my $proteinFeatureAdaptor = $db->get_ProteinFeatureAdaptor;
  my $proteinFeature = $proteinFeatureAdaptor->fetch_Feature_by_dbid ($id);


=head1 DESCRIPTION

  This Object inherits from BaseAdaptor, following the new adaptor rules.
  The main method is fetch_Feature_by_dbid, which returns a complete feature object.


=head1 CONTACT

  Marc Sohrmann (ms2@sanger.ac.uk)

=head1 APPENDIX

  The rest of the documentation details each of the object methods.
  Internal methods are usually preceded with a _.

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::DBSQL::Protein::ProteinFeatureAdaptor;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::SeqFeature;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


=head2 write_Protein_feature

 Title    : write_Protein_feature
 Usage    : $self->write_Protein_feature (@features)
 Function :
 Example  : 
 Returns  : 
 Args     : array of Bio::EnsEMBL::SeqFeatureI
 Throws   :

=cut

sub write_Protein_feature {
    my ($self, @features) = @_;

    # Get AnalysisAdaptor
    my $analysisAdaptor = $self->db->get_AnalysisAdaptor;

    FEATURE:foreach my $feature (@features) { 
        # check that we have got features
        if( ! $feature->isa('Bio::EnsEMBL::SeqFeatureI') ) {
            $self->throw("Feature $feature is not a SeqFeatureI");
        }
        # check that we have Analysis
        my $analysisid;

        if (!defined ($feature->analysis)) {
            $self->throw ("Feature ".$feature->seqname." ". 
                                     $feature->source_tag." ".
                                     $feature->logic_name. 
                          " doesn't have analysis. Can't write to database");
        }
        else {
	    if ($analysisid = $analysisAdaptor->exists ($feature->analysis)) {
	    }
            else {
                $analysisid = $analysisAdaptor->store ($feature->analysis);
	    }
        }
        # check what kind of feature we're dealing with
        if ($feature->isa ('Bio::EnsEMBL::FeaturePair')) {
            $self->_store_FeaturePair ($feature, $analysisid);
        }
        else {
            # it is a Bio::EnsEMBL::SeqFeatureI
            $self->_store_Feature ($feature, $analysisid);
        }
    }
}


=head2 _store_Feature

 Title    : _store_Feature
 Usage    : $self-> _store_Feature ($feature, $analysisid);
 Function : writes a SeqFeature to the database
 Example  : 
 Returns  : 
 Args     : a Bio::EnsEMBL::SeqFeatureI, an AnalysisId
 Throws   :

=cut

sub _store_Feature {
    my ($self, $feature, $analysisid) = @_;

    $feature->validate_prot_feature;

    my $sth = $self->prepare ( q{ INSERT INTO protein_feature
                                              (protein_featureId,
                                               analysis,
                                               proteinId, start, end,
                                               score, evalue, percent_id, strand)
                                       VALUES ('NULL', ?, ?, ?, ?, ?, ?, ?, ?)
				} );

    $sth->execute ($analysisid,
                   $feature->seqname, $feature->start, $feature->end,
                   $feature->score,
                   $feature->p_value,
                   $feature->percent_id,
                   $feature->strand);
    $sth->finish;
}


=head2 

 Title    : write_Protein_Blast_feature
 Usage    : $self-> write_Protein_Blast_feature ($proteinFeature)
 Function : writes a protein FeaturePair (Blast) to the database,
            including the gapped_alignment coordinates
 Example  : 
 Returns  : 
 Args     : a Bio::EnsEMBL::FeaturePair
 Throws   :

=cut _store_FeaturePair

sub _store_FeaturePair{
    my ($self, $featurepair, $analysisid) = @_;
        
    $featurepair->feature1->validate_prot_feature;
    $featurepair->feature2->validate_prot_feature;

    my $cigar_string;
    if ($featurepair->feature1->has_tag ('cigar')) {
        my @cigar_tags = $featurepair->feature1->each_tag_value ('cigar');
        $cigar_string = $cigar_tags[0];

    }


    ################
    # a little temporary hack to modify gadfly id's (Fly database)

    my $target_seqname;

    if ($featurepair->hseqname =~ /^\S+\|FB\S+\|(CT\d+)\|FB\S+/) {
        $target_seqname = "$1";
    }
    else {
        $target_seqname = $featurepair->hseqname;
    }
    ################

    my $sth = $self->prepare ( q{ INSERT INTO protein_feature
                                              (protein_featureId,
                                               analysis,
                                               proteinId, start, end,
                                               hId, hstart, hend,
                                               score, evalue, percent_id, strand, cigar)
                                       VALUES ('NULL', ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
				 } );

    $sth->execute ($analysisid,
                   $featurepair->seqname, $featurepair->start, $featurepair->end,
                   $target_seqname, $featurepair->hstart, $featurepair->hend,
                   $featurepair->score,
                   $featurepair->p_value,
                   $featurepair->percent_id,
                   $featurepair->strand,
                   $cigar_string);

    $sth->finish;
}

1;
