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



=head2 output_dir

  Arg [1]   : string, path to directory
  Function  : container for path to directoy
  Returntype: string,
  Exceptions: throws if directory doesn't exist'
  Example   : my $filename = $self->output_dir'/filename'

=cut


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
#these are methods which provide generic functionality for the
#object


=head2 dump_table

  Arg [1]   : string, name of table
  Arg [2]   : string, name of file to dump to
  Arg [3]   : string, where clause for sql statement if wanted
  Arg [4]   : string, select clause for sql if wanted
  Arg [5]   : string, from clause of sql if wanted
  Arg [6]   : string, into outfile clause of wanted
  Function  : construction the sql to dump a particular tables
  contents for a particular slice if desire. The last four arguments
  are options to allow for various different statements to be generated
  as not all tables fit the select * from table where seq_region_id = $id
  model
  Returntype: filename
  Exceptions: throws if not passed a table name or if the filename
  already exists as mysql won't dump to an existing file'
  Example   : my $filename = $slice_dump->dump_table()


=cut


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
  return $filename;
}



=head2 genereate_where_clause

  Arg [1]   : int, seq_region_id
  Arg [2]   : int, seq_region_start
  Arg [3]   : int, seq_region_end
  Arg [4]   : int, boolean toggle 
  Function  : This method will construct a standard where clause
  working on the assumption that the column names are seq_region_id,
  seq_region_start and seq_region_end the last toggle is to whether to
  get features which overlap the specified boundaries
  Returntype: string, the where clause
  Exceptions: throws if not passed a seq_region_id
  Example   : 

=cut



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



=head2 can_dump

  Arg [1]   : string table_name
  Arg [2]   : int seq_region_id
  Arg [3]   : string column name
  Function  : This checks if the id passed in has any entries in the
  table specified for the column name specified, the default column name
  is seq_region_id
  Returntype: int, count
  Exceptions: throws if not passed either a table name or a seq_region_id
  Example   : if($self->can_dump('seq_region', 1))

=cut



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



=head2 get_filename

  Arg [1]   : string table name
  Arg [2]   : Bio::EnsEMBL::Slice
  Function  : generate a filename useable by mysql import on the basis of
  the tablename and info from the slice, can't use Slice::name as the 
  format upsets mysqlimport
  Returntype: string, filename
  Exceptions: none
  Example   : my $filename = $self->get_filename(seq_region', $slice)

=cut


sub get_filename{
  my ($self, $table_name, $slice) = @_;
  return $self->output_dir."/".$table_name.".".$slice->seq_region_name.
    ".".$slice->start."-".$slice->end;
}


#
#Dump methods for tables which require partial dumps
#All these dumps will require a slice object and will take either the
#slice held by the SliceDump object or a slice which is passed in
#

#the following tables have structures which fit a standard model
#they all have columns seq_region_id, seq_region_start and seq_region_end
#all these methods have the name which follows the structure
#dump_table_name_table to allow for on the fly calling of methods
#for appropriate tables

=head2 dump_name_table

  Arg [1]   : Bio::EnsEMBL::Slice
  Function  : check if any data can be dumped from the defined table
  and if it can call the dump table method with the appropriate
  arguments
  Returntype: returns 1 if a dump was made 0 if not 
  Exceptions: throws if no slice is defined either in the method args or
  SliceDump::Slice or if that slice isn't a slice object'
  Example   : $SliceDump->dump_seq_region_table

=cut



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



#dumps entries from the assembly table
#can use standard model as must take entries whose column labelled
#asm_seq_region_id match the given id and the start and end are asm_start
#and asm_end respectively

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

#dumps entries from the repeat_consensus table
#this must generate sql which joins to the repeat feature table
# in order to dump its entries

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


#dumps entries from the translation table
#must join to the transcript table to do this

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


#dumps entries from the exon_transcript table
#must join to the transcript table to do this

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

#dumps entries from the supporting_feature table
#must join to the exon table to do this

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

#dumps entries from the transcript_stable_id table
#must join to the transcript table to do this

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

#dumps entries from the exon_stable_id table
#must join to the exon table to do this

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

#dumps entries from the gene_stable_id table
#must join to the gene table to do this 

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



#dumps entries from the protein_feature table
#must join to the translation and transcript tables to do this

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
