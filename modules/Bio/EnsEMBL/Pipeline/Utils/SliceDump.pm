package Bio::EnsEMBL::Pipeline::Utils::SliceDump;

use strict;
use warnings;

use Exporter;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Root);


use Bio::EnsEMBL::Utils::Exception qw(throw warning verbose);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );


sub new{
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  &verbose('WARNING');
  my ($db, $slice, $output_dir ) = rearrange ([ 'DB', 
                                                'SLICE', 'OUTPUT_DIR'], @args);
  if(!$db){
    throw("Can't run SliceDump without a database adaptor");
  }
  $self->db($db);
  $self->slice($slice);
  $self->output_dir($output_dir);
  return $self;
}


###################
#container methods#
###################


=head2 db

  Arg [1]   : Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
  Function  : stores the DBadaptor for the object
  Returntype: Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
  Exceptions: throws if argument passed in isn't a 
              Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor'
  Example   : my $rules_adaptor = $self->db->get_RulesAdaptor;

=cut

sub db{
  my ($self, $db) = @_;

  if($db){
    if(!$db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor')){
      throw("Can't run the RuleManager with $db you need a ".
            "Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor");
    }
    $self->{'dbadaptor'} = $db;
  }

  return $self->{'dbadaptor'};
}


=head2 slice

  Arg [1]    : (optional) Bio::EnsEMBL::Slice $slice
  Example    : $seqname = $feature->slice()->name();
  Description: Getter/Setter for the Slice that data is desired from
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : thrown if an invalid argument is passed
  Caller     : general

=cut

sub slice {
  my $self = shift;

  if(@_) {
    my $sl = shift;
    if(defined($sl) && (!ref($sl) || !$sl->isa('Bio::EnsEMBL::Slice'))) {
      throw('slice argument must be a Bio::EnsEMBL::Slice');
    }

    $self->{'slice'} = $sl;
  }

  return $self->{'slice'};
}



sub output_dir{
  my ($self, $output_dir) = @_;
  if($output_dir){
    throw("Output dir ".$output_dir." must exist or we can't use it") 
      unless (-d $output_dir);
    $self->{'output_dir'} = $output_dir;
  }
  return $self->{'output_dir'};
}


#################
#Utility methods#
#################



=head2 dump_table

  Arg [1]   : 
  Function  : method to specified table
  into the desired directory
  Returntype: filename
  Exceptions: 
  Example   : my $filename = $slice_dump->dump_table()


=cut

#sub dump_table{
#  my ($self, $table_name, $filename, $where_clause) = @_;
#  #print STDERR "Dumping ".$table_name."\n";
#  #print STDERR "with where clause ".$where_clause."\n" if($where_clause);
#  if(!$table_name){
#    throw("Can't dump table without tablename");
#  }
#  if(!$filename){
#    $filename = $self->output_dir."/".$table_name;
#    #warning("No filename provided so using ".$filename);
#  }
#  if(-e $filename){
#    throw($filename." exists mysql can't dump into an existing file ");
#  }
#  my $sql = "select * from ".$table_name;
#  $sql .= " where ".$where_clause if($where_clause);
#  $sql .= " into outfile '".$filename."'";
#  my $sth = $self->db->prepare($sql);
#  $sth->execute();
#  return $filename;
#}

sub dump_table{
  my ($self, $table_name, $filename, $where, $select, $from, $out) = @_;
  if(!$table_name){
    throw("Can't dump table without tablename");
  }
  if(!$filename){
    $filename = $self->output_dir."/".$table_name;
  }
  if(-e $filename){
    throw($filename." exists mysql can't dump into an existing file ");
  }
  $select = "select * " unless($select);
  $from = "from ".$table_name." " unless($from);
  $out = " into outfile '".$filename."'" unless($out);
  my $sql = $select;
  $sql .= $from;
  $sql .= $where if($where);
  $sql .= $out;
  print $sql."\n";
  my $sth = $self->db->prepare($sql);
  $sth->execute;
  return $out;
}


sub generate_where_clause{
  my ($self, $id, $start, $end, $overlaps_boundaries) = @_;
  
  if(!$id){
    throw("Can't generate a where clause without seq_region_id");
  }
  my $where_clause = "where seq_region_id = ".$id;
  if($overlaps_boundaries){
    $where_clause .= " and seq_region_start <= ".$end if($end);
    $where_clause .= " and seq_region_end >= ".$start if($start);
  }else{
    $where_clause .= " and seq_region_start >= ".$start if($start);
    $where_clause .= " and seq_region_end <= ".$end if($end);
  }
  return $where_clause;
}


sub can_dump{
  my ($self, $table_name, $id, $column_name) = @_;

  if(!$table_name || !$id){
    throw("Can't check if can dump for a table with a specific seq_region_id ".
          " if table_name $table_name or id $id isn't defined");
  }
  if(!$column_name){
    $column_name = 'seq_region_id';
  }
  my $sql = "select count(*) from ".$table_name." where ".$column_name.
    " = ".$id;
  my $sth = $self->db->prepare($sql);
  $sth->execute;
  my ($count) = $sth->fetchrow;
  return $count;
}



sub get_filename{
  my ($self, $table_name, $slice) = @_;
  return $self->output_dir."/".$table_name.".".$slice->seq_region_name.".".
    $slice->start."-".$slice->end;
}


#
#Dump methods for tables which require partial dumps
#All these dumps will require a slice object and will take either the
#slice held by the SliceDump object or a slice which is passed in
#

sub dump_seq_region_table{
  my ($self, $slice) = @_;
  if(!$slice){
    $slice = $self->slice;
  }
  if(!$slice || !$slice->isa("Bio::EnsEMBL::Slice")){
    throw("Can't dump partial table without slice information");
  }
  my $id = $slice->get_seq_region_id;
  if($self->can_dump('seq_region', $id)){
    my $filename = $self->get_filename('seq_region', $slice);
    my $where_clause = $self->generate_where_clause($id);
    $self->dump_table('seq_region', $filename, $where_clause);
    return 1;
  }
  return 0;
}



sub dump_dna_table{
  my ($self, $slice) = @_;
  if(!$slice){
    $slice = $self->slice;
  }
  if(!$slice || !$slice->isa("Bio::EnsEMBL::Slice")){
    throw("Can't dump partial table without slice information");
  }
  if(!$slice->coord_system->is_sequence_level){
    return 1;
  }
  my $id = $slice->get_seq_region_id;
  if($self->can_dump('dna', $id)){
    my $filename = $self->get_filename('dna', $slice);
    my $where_clause = $self->generate_where_clause($id);
    $self->dump_table('dna', $filename, $where_clause);
    return 1;
  }
  return 0;
}

sub dump_seq_region_attrib_table{
  my ($self, $slice) = @_;
  if(!$slice){
    $slice = $self->slice;
  }
  if(!$slice || !$slice->isa("Bio::EnsEMBL::Slice")){
    throw("Can't dump partial table without slice information");
  }
  my $id = $slice->get_seq_region_id;
  if($self->can_dump('seq_region_attrib', $id)){
    my $filename = $self->get_filename('seq_region_attrib', $slice);
    my $where_clause = $self->generate_where_clause($id);
    $self->dump_table('seq_region_attrib', $filename, $where_clause);
    return 1;
  }
  return 0;
}


sub dump_assembly_exception_table{
  my ($self, $slice) = @_;
  if(!$slice){
    $slice = $self->slice;
  }
  if(!$slice || !$slice->isa("Bio::EnsEMBL::Slice")){
    throw("Can't dump partial table without slice information");
  }
  my $id = $slice->get_seq_region_id;
  if($self->can_dump('assembly_exception', $id)){
    my $filename = $self->get_filename('assembly_exception', $slice);
    my $where_clause = $self->generate_where_clause($id);
    $self->dump_table('assembly_exception', $filename, $where_clause);
    return 1;
  }
  return 0;
}


sub dump_repeat_feature_table{
  my ($self, $slice) = @_;
  if(!$slice){
    $slice = $self->slice;
  }
  if(!$slice || !$slice->isa("Bio::EnsEMBL::Slice")){
    throw("Can't dump partial table without slice information");
  }
  my $id = $slice->get_seq_region_id;
  if($self->can_dump('repeat_feature', $id)){
    my $filename = $self->get_filename('repeat_feature', $slice);
    my $where_clause = $self->generate_where_clause($id, $slice->start, 
                                                    $slice->end);
    $self->dump_table('repeat_feature', $filename, $where_clause);
    return 1;
  }
  return 0;
}

sub dump_simple_feature_table{
  my ($self, $slice) = @_;
  if(!$slice){
    $slice = $self->slice;
  }
  if(!$slice || !$slice->isa("Bio::EnsEMBL::Slice")){
    throw("Can't dump partial table without slice information");
  }
  my $id = $slice->get_seq_region_id;
  if($self->can_dump('simple_feature', $id)){
    my $filename = $self->get_filename('simple_feature', $slice);
    my $where_clause = $self->generate_where_clause($id, $slice->start, 
                                                    $slice->end);
    $self->dump_table('simple_feature', $filename, $where_clause);
    return 1;
  }
  return 0;
}

sub dump_protein_align_feature_table{
  my ($self, $slice) = @_;
  if(!$slice){
    $slice = $self->slice;
  }
  if(!$slice || !$slice->isa("Bio::EnsEMBL::Slice")){
    throw("Can't dump partial table without slice information");
  }
  my $id = $slice->get_seq_region_id;
  if($self->can_dump('protein_align_feature', $id)){
    my $filename = $self->get_filename('protein_align_feature', $slice);
    my $where_clause = $self->generate_where_clause($id, $slice->start, 
                                                    $slice->end);
    $self->dump_table('protein_align_feature', $filename, $where_clause);
    return 1;
  }
  return 0;
}

sub dump_dna_align_feature_table{
  my ($self, $slice) = @_;
  if(!$slice){
    $slice = $self->slice;
  }
  if(!$slice || !$slice->isa("Bio::EnsEMBL::Slice")){
    throw("Can't dump partial table without slice information");
  }
  my $id = $slice->get_seq_region_id;
  if($self->can_dump('dna_align_feature', $id)){
    my $filename = $self->get_filename('dna_align_feature', $slice);
    my $where_clause = $self->generate_where_clause($id, $slice->start, 
                                                    $slice->end);
    $self->dump_table('dna_align_feature', $filename, $where_clause);
    return 1;
  }
  return 0;
}

sub dump_prediction_exon_table{
  my ($self, $slice) = @_;
  if(!$slice){
    $slice = $self->slice;
  }
  if(!$slice || !$slice->isa("Bio::EnsEMBL::Slice")){
    throw("Can't dump partial table without slice information");
  }
  my $id = $slice->get_seq_region_id;
  if($self->can_dump('prediction_exon', $id)){
    my $filename = $self->get_filename('prediction_exon', $slice);
    my $where_clause = $self->generate_where_clause($id, $slice->start, 
                                                    $slice->end);
    $self->dump_table('prediction_exon', $filename, $where_clause);
    return 1;
  }
  return 0;
}

sub dump_prediction_transcript_table{
  my ($self, $slice) = @_;
  if(!$slice){
    $slice = $self->slice;
  }
  if(!$slice || !$slice->isa("Bio::EnsEMBL::Slice")){
    throw("Can't dump partial table without slice information");
  }
  my $id = $slice->get_seq_region_id;
  if($self->can_dump('prediction_transcript', $id)){
    my $filename = $self->get_filename('prediction_transcript', $slice);
    my $where_clause = $self->generate_where_clause($id, $slice->start, 
                                                    $slice->end);
    $self->dump_table('prediction_transcript', $filename, $where_clause);
    return 1;
  }
  return 0;
}

sub dump_gene_table{
  my ($self, $slice) = @_;
  if(!$slice){
    $slice = $self->slice;
  }
  if(!$slice || !$slice->isa("Bio::EnsEMBL::Slice")){
    throw("Can't dump partial table without slice information");
  }
  my $id = $slice->get_seq_region_id;
  if($self->can_dump('gene', $id)){
    my $filename = $self->get_filename('gene', $slice);
    my $where_clause = $self->generate_where_clause($id, $slice->start, 
                                                    $slice->end);
    $self->dump_table('gene', $filename, $where_clause);
    return 1;
  }
  return 0;

}

sub dump_exon_table{
  my ($self, $slice) = @_;
  if(!$slice){
    $slice = $self->slice;
  }
  if(!$slice || !$slice->isa("Bio::EnsEMBL::Slice")){
    throw("Can't dump partial table without slice information");
  }
  my $id = $slice->get_seq_region_id;
  if($self->can_dump('exon', $id)){
    my $filename = $self->get_filename('exon', $slice);
    my $where_clause = $self->generate_where_clause($id, $slice->start, 
                                                    $slice->end);
    $self->dump_table('exon', $filename, $where_clause);
    return 1;
  }
  return 0;
}

sub dump_transcript_table{
  my ($self, $slice) = @_;
  if(!$slice){
    $slice = $self->slice;
  }
  if(!$slice || !$slice->isa("Bio::EnsEMBL::Slice")){
    throw("Can't dump partial table without slice information");
  }
  my $id = $slice->get_seq_region_id;
  if($self->can_dump('transcript', $id)){
    my $filename = $self->get_filename('transcript', $slice);
    my $where_clause = $self->generate_where_clause($id, $slice->start, 
                                                    $slice->end);
    $self->dump_table('transcript', $filename, $where_clause);
    return 1;
  }
  return 0;
}
##
##methods for tables which don't contain straight ref to seq_region_id
##


sub dump_assembly_table{
  my ($self, $slice) = @_;
  if(!$slice){
    $slice = $self->slice;
  }
  if(!$slice || !$slice->isa("Bio::EnsEMBL::Slice")){
    throw("Can't dump partial table without slice information");
  }
  my $id = $slice->get_seq_region_id;
  if($self->can_dump('assembly', $id, 'asm_seq_region_id')){
    my $filename = $self->get_filename('assembly', $slice);
    my $where_clause = "where asm_seq_region_id = ".$id.
      " and asm_start <= ".$slice->end.
        " and asm_end >= ".$slice->start;
    $self->dump_table('assembly', $filename, $where_clause);
    return 1;
  }
  return 0;
}

sub dump_repeat_consensus_table{
  my ($self, $slice) = @_;
  if(!$slice){
    $slice = $self->slice;
  }
  if(!$slice || !$slice->isa("Bio::EnsEMBL::Slice")){
    throw("Can't dump partial table without slice information");
  }
  my $id = $slice->get_seq_region_id;
  if($self->can_dump('repeat_feature', $id)){
    my $filename = $self->get_filename('repeat_consensus', $slice);
    my $select = "select rc.* ";
    my $from = "from repeat_consensus rc, repeat_feature rf ";
    my $where = "where rf.repeat_consensus_id = rc.repeat_consensus_id ".
      "and rf.seq_region_id = ".$id." ".
        "and rf.seq_region_start = ".$slice->start." ".
          "and rf.seq_region_end = ".$slice->end." ";
    my $out = "into outfile '".$filename."'";
    $self->dump_table('repeat_consensus', $filename, $where, $select, 
                      $from, $out);
    return 1;
  }
  return 0;
}


sub dump_translation_table{
  my ($self, $slice) = @_;
  if(!$slice){
    $slice = $self->slice;
  }
  if(!$slice || !$slice->isa("Bio::EnsEMBL::Slice")){
    throw("Can't dump partial table without slice information");
  }
  my $id = $slice->get_seq_region_id;
  if($self->can_dump('transcript', $id)){
    my $filename = $self->get_filename('translation', $slice);
    my $select = "select translation.* ";
    my $from = "from translation, transcript ";
    my $where = "where ".
      "translation.transcript_id = transcript.transcript_id and ".
        "seq_region_id = ".$id." and ".
          "seq_region_start >= ".$slice->start." and ".
            "seq_region_end <= ".$slice->end." ";
    my $out = "into outfile '".$filename."'";
    $self->dump_table('translation', $filename, $where, $select, 
                      $from, $out);
    return 1;
  }
  return 0;
}

sub dump_exon_transcript_table{
  my ($self, $slice) = @_;
  if(!$slice){
    $slice = $self->slice;
  }
  if(!$slice || !$slice->isa("Bio::EnsEMBL::Slice")){
    throw("Can't dump partial table without slice information");
  }
  my $id = $slice->get_seq_region_id;
  if($self->can_dump('transcript', $id)){
    my $filename = $self->get_filename('exon_transcript', $slice);
    my $select = "select exon_transcript.* ";
    my $from = "from exon_transcript, transcript ";
    my $where = "where ".
      "exon_transcript.transcript_id = transcript.transcript_id and ".
        "seq_region_id = ".$id." and ".
          "seq_region_start >= ".$slice->start." and ".
            "seq_region_end <= ".$slice->end." ";
    my $out = "into outfile '".$filename."'";
    $self->dump_table('exon_transcript', $filename, $where, $select, 
                      $from, $out);
    return 1;
  }
  return 0;
}

sub dump_supporting_feature_table{
  my ($self, $slice) = @_;
  if(!$slice){
    $slice = $self->slice;
  }
  if(!$slice || !$slice->isa("Bio::EnsEMBL::Slice")){
    throw("Can't dump partial table without slice information");
  }
  my $id = $slice->get_seq_region_id;
  if($self->can_dump('exon', $id)){
    my $filename = $self->get_filename('supporting_feature', $slice);
    my $select = "select supporting_feature.* ";
    my $from = "from supporting_feature, exon ";
    my $where = "where ".
      "supporting_feature.exon_id = exon.exon_id and ".
        "seq_region_id = ".$id." and ".
          "seq_region_start >= ".$slice->start." and ".
            "seq_region_end <= ".$slice->end." ";
    my $out = "into outfile '".$filename."'";
    $self->dump_table('supporting_feature', $filename, $where, $select, 
                      $from, $out);
    return 1;
  }
  return 0;
}



sub dump_transcript_stable_id_table{
  my ($self, $slice) = @_;
  if(!$slice){
    $slice = $self->slice;
  }
  if(!$slice || !$slice->isa("Bio::EnsEMBL::Slice")){
    throw("Can't dump partial table without slice information");
  }
  my $id = $slice->get_seq_region_id;
  if($self->can_dump('transcript', $id)){
    my $filename = $self->get_filename('transcript_stable_id', $slice);
    my $select = "select transcript_stable_id.* ";
    my $from = "from transcript_stable_id, transcript ";
    my $where = "where ".
      "transcript_stable_id.transcript_id = transcript.transcript_id ".
        " and seq_region_id = ".$id." and ".
          "seq_region_start >= ".$slice->start." and ".
            "seq_region_end <= ".$slice->end." ";
    my $out = "into outfile '".$filename."'";
    $self->dump_table('transcript_stable_id', $filename, $where, $select,
                      $from, $out);
    return 1;
  }
  return 0;
}

sub dump_exon_stable_id_table{
  my ($self, $slice) = @_;
  if(!$slice){
    $slice = $self->slice;
  }
  if(!$slice || !$slice->isa("Bio::EnsEMBL::Slice")){
    throw("Can't dump partial table without slice information");
  }
  my $id = $slice->get_seq_region_id;
  if($self->can_dump('exon', $id)){
    my $filename = $self->get_filename('exon_stable_id', $slice);
    my $select = "select exon_stable_id.* ";
    my $from = "from exon_stable_id, exon ";
    my $where = "where ".
      "exon_stable_id.exon_id = exon.exon_id and ".
        "seq_region_id = ".$id." and ".
          "seq_region_start >= ".$slice->start." and ".
            "seq_region_end <= ".$slice->end." ";
    my $out = "into outfile '".$filename."'";
    $self->dump_table('exon_stable_id', $filename, $where, $select, 
                      $from, $out);
    return 1;
  }
  return 0;
}

sub dump_gene_stable_id_table{
  my ($self, $slice) = @_;
  if(!$slice){
    $slice = $self->slice;
  }
  if(!$slice || !$slice->isa("Bio::EnsEMBL::Slice")){
    throw("Can't dump partial table without slice information");
  }
  my $id = $slice->get_seq_region_id;
  if($self->can_dump('gene', $id)){
    my $filename = $self->get_filename('gene_stable_id', $slice);
    my $select = "select gene_stable_id.* ";
    my $from = "from gene_stable_id, gene ";
    my $where = "where ".
      "gene_stable_id.gene_id = gene.gene_id and ".
        "seq_region_id = ".$id." and ".
          "seq_region_start >= ".$slice->start." and ".
            "seq_region_end <= ".$slice->end." ";
    my $out = "into outfile '".$filename."'";
    $self->dump_table('gene_stable_id', $filename, $where, $select, 
                      $from, $out);
    return 1;
  }
  return 0;
}





sub dump_protein_feature_table{
  my ($self, $slice) = @_;
  if(!$slice){
    $slice = $self->slice;
  }
  if(!$slice || !$slice->isa("Bio::EnsEMBL::Slice")){
    throw("Can't dump partial table without slice information");
  }
  my $id = $slice->get_seq_region_id;
  if($self->can_dump('transcript', $id)){
    my $filename = $self->get_filename('protein_feature', $slice);
    if(-e $filename){
      throw($filename." exists mysql can't dump into an existing file ");
    }
    my $sql = "select protein_feature.* from protein_feature, translation, ".
      "transcript where translation.transcript_id = transcript.transcript_id ".
        "and seq_region_id = ".$id." and seq_region_start >= ".
          $slice->start." and seq_region_end <= ".$slice->end." ".
            "and protein_feature.translation_id = translation.translation_id ".
              "into outfile '".$filename."'";
    my $sth = $self->db->prepare($sql);
    $sth->execute;
    return 1;
  }
  return 0;
}


1;
