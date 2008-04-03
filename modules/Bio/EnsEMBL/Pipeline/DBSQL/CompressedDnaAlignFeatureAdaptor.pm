#
# EnsEMBL module for Bio::EnsEMBL::Pipeline::DBSQL::CompressedDnaAlignFeatureAdaptor
#
# Copyright (c) 2003 EnsEMBL
#
# You may distribute this module under the same terms as perl itself

=head1 NAME

Bio::EnsEMBL::Pipeline::DBSQL::CompressedDnaAlignFeatureAdaptor - Adaptor for compressed DnaAlignFeatures

=head1 SYNOPSIS


    $cdafa = Bio::EnsEMBL::Pipeline::DBSQL::CompressedDnaAlignFeatureAdaptor->new
      (
       $db_adaptor
      );

    @features = @{$cdafa->fetch_by_Slice($slice)};

    $cdafa->store(@features);

=head1 DESCRIPTION

This is an adaptor responsible for the retrieval and storage of 
CompressedDnaDnaAlignFeatures from the database. This adaptor inherits most of its 
functionality from the DnaAlignFeatureAdaptor and BaseFeatureAdaptor 
superclasses. It just overrides the store method to allow a the features to be
compressed

=head1 CONTACT

Post questions to the EnsEMBL development list <ensembl-dev@ebi.ac.uk>

=cut


package Bio::EnsEMBL::Pipeline::DBSQL::CompressedDnaAlignFeatureAdaptor;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::DBSQL::BaseAlignFeatureAdaptor;
use Bio::EnsEMBL::DBSQL::DnaAlignFeatureAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

@ISA = qw(Bio::EnsEMBL::DBSQL::DnaAlignFeatureAdaptor);

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);
  return $self;
}


sub store {
  my ($self, @feats) = @_;

  throw("Must call store with features") if( scalar(@feats) == 0 );

  my @tabs = $self->_tables;
  my ($tablename) = @{$tabs[0]};

  my $db = $self->db();
  my $analysis_adaptor = $db->get_AnalysisAdaptor();

  my $sth = $self->prepare
    (
     "UPDATE IGNORE $tablename set evalue = evalue+1 
      WHERE seq_region_id = ? 
      AND seq_region_start = ?
      AND seq_region_end  = ?
      AND seq_region_strand = ?
      AND analysis_id = ?");
  
 FEATURE: foreach my $feat ( @feats ) {
    if( !ref $feat || !$feat->isa("Bio::EnsEMBL::DnaDnaAlignFeature") ) {
      throw("feature must be a Bio::EnsEMBL::DnaDnaAlignFeature,"
            . " not a [".ref($feat)."].");
    }

    if($feat->is_stored($db)) {
      warning("DnaDnaAlignFeature [".$feat->dbID."] is already stored" .
              " in this database.");
      next FEATURE;
    }

    my $hstart = $feat->hstart();
    my $hend   = $feat->hend();
    my $hstrand = $feat->hstrand();
    $self->_check_start_end_strand($hstart,$hend, $hstrand);

    my $cigar_string = $feat->cigar_string();
    if(!$cigar_string) {
      $cigar_string = $feat->length() . 'M';
      warning("DnaDnaAlignFeature does not define a cigar_string.\n" .
              "Assuming ungapped block with cigar_line=$cigar_string .");
    }

    my $hseqname = $feat->hseqname();
    if(!$hseqname) {
      throw("DnaDnaAlignFeature must define an hseqname.");
    }

    if(!defined($feat->analysis)) {
      throw("An analysis must be attached to the features to be stored.");
    }

    #store the analysis if it has not been stored yet
    if(!$feat->analysis->is_stored($db)) {
      $analysis_adaptor->store($feat->analysis());
    }

    my $original = $feat;
    my $seq_region_id;
    ($feat, $seq_region_id) = $self->_pre_store($feat);
    $sth->bind_param(1,$seq_region_id,SQL_INTEGER);
    $sth->bind_param(2,$feat->start,SQL_INTEGER);
    $sth->bind_param(3,$feat->end,SQL_INTEGER);
    $sth->bind_param(4,$feat->strand,SQL_TINYINT);
    $sth->bind_param(5,$feat->analysis->dbID,SQL_INTEGER);
    my $result = $sth->execute();
 #   print " $result \n";
    if ($result eq "0E0"){
      $self->SUPER::store($original);
      return;
    }
    $original->dbID($sth->{'mysql_insertid'});
    $original->adaptor($self);
  }

  $sth->finish();
}

1;
