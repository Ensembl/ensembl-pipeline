=head1 LICENSE

# Copyright [1999-2013] Genome Research Ltd. and the EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

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

=cut


package Bio::EnsEMBL::Pipeline::DBSQL::CompressedDnaAlignFeatureAdaptor;
use warnings ;
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
