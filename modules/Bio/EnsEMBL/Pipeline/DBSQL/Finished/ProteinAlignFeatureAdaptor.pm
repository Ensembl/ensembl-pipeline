### Bio::EnsEMBL::Pipeline::DBSQL::Finished::ProteinAlignFeatureAdaptor

package Bio::EnsEMBL::Pipeline::DBSQL::Finished::ProteinAlignFeatureAdaptor;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DnaPepAlignFeature;
use Bio::EnsEMBL::DBSQL::ProteinAlignFeatureAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

@ISA = qw(Bio::EnsEMBL::DBSQL::ProteinAlignFeatureAdaptor);

sub store{
 my ($self, @feats) = @_;

  throw("Must call store with features") if( scalar(@feats) == 0 );

  my $seq_region_id;
  my $analysis_id;
  my $first_dbID;
  my $last_dbID;

  my @tabs = $self->_tables;
  my ($tablename) = @{$tabs[0]};
  my $history_table = $tablename.'_history';

  my $db = $self->db();
  my $analysis_adaptor = $db->get_AnalysisAdaptor();

  my $sth = $self->prepare(
     "INSERT INTO $tablename (seq_region_id, seq_region_start, seq_region_end,
                             seq_region_strand, hit_start, hit_end,
                             hit_name, cigar_line,
                             analysis_id, score, evalue, perc_ident,external_db_id, hcoverage)
     VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)");

  my $sth_history = $self->prepare(
  	 "INSERT INTO $history_table (seq_region_id, analysis_id, align_feature_id_start,
	 							align_feature_id_end, db_version , date)
	 VALUES (?,?,?,?,?,NOW())");

 FEATURE: foreach my $feat ( @feats ) {
   if( !ref $feat || !$feat->isa("Bio::EnsEMBL::DnaPepAlignFeature") ) {
     throw("feature must be a Bio::EnsEMBL::DnaPepAlignFeature,"
           . " not a [".ref($feat)."].");
   }

   if($feat->is_stored($db)) {
     warning("PepDnaAlignFeature [".$feat->dbID."] is already stored" .
             " in this database.");
     next FEATURE;
   }

   #sanity check the hstart and hend
   my $hstart  = $feat->hstart();
   my $hend    = $feat->hend();
   $self->_check_start_end_strand($hstart,$hend,1);

   my $cigar_string = $feat->cigar_string();
   if(!$cigar_string) {
     $cigar_string = $feat->length() . 'M';
     warning("DnaPepAlignFeature does not define a cigar_string.\n" .
             "Assuming ungapped block with cigar_string=$cigar_string\n");
   }

   my $hseqname = $feat->hseqname();
   if(!$hseqname) {
     throw("DnaPepAlignFeature must define an hseqname.");
   }

   if(!defined($feat->analysis)) {
     throw("An analysis must be attached to the features to be stored.");
   }

   #store the analysis if it has not been stored yet
   if(!$feat->analysis->is_stored($db)) {
     $analysis_adaptor->store($feat->analysis());
   }

   my $original = $feat;
    $analysis_id = $feat->analysis->dbID;
   ($feat, $seq_region_id) = $self->_pre_store($feat);

   $sth->bind_param(1,$seq_region_id,{ TYPE => 'SQL_INTEGER' });
   $sth->bind_param(2,$feat->start,{ TYPE => 'SQL_INTEGER' });
   $sth->bind_param(3,$feat->end,{ TYPE => 'SQL_INTEGER' });
   $sth->bind_param(4,$feat->strand,{ TYPE => 'SQL_TINYINT' });
   $sth->bind_param(5,$feat->hstart,{ TYPE => 'SQL_INTEGER' });
   $sth->bind_param(6,$feat->hend,{ TYPE => 'SQL_INTEGER' });
   $sth->bind_param(7,$feat->hseqname,{ TYPE => 'SQL_VARCHAR' });
   $sth->bind_param(8,$feat->cigar_string,{ TYPE => 'SQL_LONGVARCHAR' });
   $sth->bind_param(9,$analysis_id,{ TYPE => 'SQL_INTEGER' });
   $sth->bind_param(10,$feat->score,{ TYPE => 'SQL_DOUBLE' });
   $sth->bind_param(11,$feat->p_value,{ TYPE => 'SQL_DOUBLE' });
   $sth->bind_param(12,$feat->percent_id,{ TYPE => 'SQL_FLOAT' });
   $sth->bind_param(13,$feat->external_db_id,{ TYPE => 'SQL_INTEGER' });
   $sth->bind_param(14,$feat->hcoverage,{ TYPE => 'SQL_DOUBLE' });

   $sth->execute();
   $original->dbID($sth->{'mysql_insertid'});
   $original->adaptor($self);

    $first_dbID ||= $original->dbID;
    $last_dbID = $original->dbID;
  }
  # save dbIDs, time and db version into history table
  $sth_history->execute($seq_region_id,
  						$analysis_id,
  						$first_dbID,
  						$last_dbID,
  						$self->db_version);
  $sth_history->finish();
 $sth->finish();
}

sub db_version {
  my $self = shift;
  $self->{'db_version'} = shift if(@_);
  return $self->{'db_version'};
}

1;
