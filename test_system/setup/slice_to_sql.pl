#!/usr/local/ensembl/bin/perl -w

use strict;
use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::ProjectionSegment;
use Bio::EnsEMBL::Pipeline::Utils::SliceDump;

$| = 1;

my $host   = '';
my $port   = '';
my $dbname = '';
my $dbuser = '';
my $dbpass = '';
my $cs_name;
my $cs_version;
my $seq_region_name;
my $start;
my $end;
my $strand = 1;
my $output_dir;
my $verbose;
my @whole_tables;
my @partial_tables;
my @whole_standard_tables = ('meta', 'meta_coord', 
                             'coord_system', 'analysis', 'attrib_type');
my @whole_pipeline_tables = ('rule_goal', 'rule_conditions', 
                             'input_id_type_analysis') ;
my @partial_standard_tables = ('seq_region', 'assembly', 'dna', 
                               'seq_region_attrib', 'assembly_exception');
my @raw_compute_tables = ('repeat_consensus', 'repeat_feature',
                          'prediction_exon', 'prediction_transcript',
                          'dna_align_feature', 'protein_align_feature',
                          'simple_feature');
my @gene_tables = ('gene', 'exon', 'transcript', 'translation',
                   'exon_transcript', 'supporting_feature',
                   'protein_align_feature', 'dna_align_feature',
                   'gene_stable_id', 'exon_stable_id',
                   'translation_stable_id', 'transcript_stable_id');
my @protein_annotation_tables = ('protein_feature');
my $raw_computes;
my $genes;
my $protein_annotation;
my @whole_commandline_tables;
my @partial_commandline_tables;
my $pipeline_tables = 1;
my $sequence_tables = 1;
my $no_defaults;
&GetOptions(
            'dbhost:s'   => \$host,
            'dbport:n'   => \$port,
            'dbname:s'   => \$dbname,
            'dbuser:s'   => \$dbuser,
            'dbpass:s'   => \$dbpass,
            'coord_system_name:s' => \$cs_name,
            'coord_system_version:s' => \$cs_version,
            'seq_region_name:s' => \$seq_region_name,
            'seq_region_start:s' => \$start,
            'seq_region_end:s' => \$end,
            'seq_region_strand:s' => \$strand,
            'verbose!' => \$verbose,
            'output_dir:s' => \$output_dir,
            'whole_table:s@' => \@whole_commandline_tables,
            'partial_table:s@' => \@partial_commandline_tables,
            'pipeline!' => \$pipeline_tables,
            'raw_computes!' => \$raw_computes,
            'genes!' => \$genes,
            'protein_annotation!' => \$protein_annotation,
           ) or throw("Can't get options");

unless ($host && $dbname && $dbuser) {
  throw("Can't run without database argument");
}
my %whole_tables;
my %partial_tables;
if(!$no_defaults){
  %whole_tables = %{setup_tablelist(\%whole_tables, 
                                     \@whole_commandline_tables)};
  %partial_tables = %{setup_tablelist(\%partial_tables, 
                                       \@partial_commandline_tables)};
  %partial_tables = %{setup_tablelist(\%partial_tables, 
                                       \@partial_standard_tables)};
  %whole_tables = %{setup_tablelist(\%whole_tables, 
                                     \@whole_standard_tables)};
  if($pipeline_tables){
    %whole_tables = %{setup_tablelist(\%whole_tables, 
                                       \@whole_pipeline_tables)};
  }
  if($raw_computes){
    %partial_tables = %{setup_tablelist(\%partial_tables, 
                                         \@raw_compute_tables)};
  }
  if($genes){
    %partial_tables = %{setup_tablelist(\%partial_tables, 
                                        \@gene_tables)};
  }
  if($protein_annotation){
    %partial_tables = %{setup_tablelist(\%partial_tables, 
                                         \@protein_annotation_tables)};
  }
}else{
  %whole_tables = %{setup_tablelist(\%whole_tables, 
                                     \@whole_commandline_tables)};
  %partial_tables = %{setup_tablelist(\%partial_tables, 
                                       \@partial_commandline_tables)};
}
@whole_tables = keys(%whole_tables);
@partial_tables = keys(%partial_tables);
my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -dbname => $dbname,
    -host   => $host,
    -user   => $dbuser,
    -port   => $port,
    -pass   => $dbpass
);

my $dump_dir = $output_dir."/".$dbname;
if(! -e $dump_dir){
  my $mkdir_cmd = "mkdir $dump_dir";
  print $mkdir_cmd."\n";
  my $chmod_cmd = "chmod 777 $dump_dir";
  print $chmod_cmd."\n";
  system($mkdir_cmd);
  throw("Couldn't make ".$dump_dir." $! ") unless(-d $dump_dir);
  system($chmod_cmd); 
}
my %dumped_whole;
my $slicedump = Bio::EnsEMBL::Pipeline::Utils::SliceDump->new
  (
   -DB => $db,
   -OUTPUT_DIR => $dump_dir,
  );

my $sa = $db->get_SliceAdaptor();
my $ca = $db->get_CoordSystemAdaptor();

my $slice = $sa->fetch_by_region($cs_name, $seq_region_name, $start,
                                 $end, 1, $cs_version);
if(!$slice){
  throw("Have been unable to fetch slice for ".$seq_region_name." ".
        $cs_name);
}
foreach my $table(@whole_tables){
  if(!$dumped_whole{$table} || -e ($dump_dir."/".$table)){
    $slicedump->dump_table($table);
    $dumped_whole{$table} = 1;
  }else{
    print "Have already dumped table ".$table."\n" if($verbose);
  }
}



my $coord_system = $ca->fetch_by_name($cs_name);
my $seq_coord_system = $ca->fetch_by_name('seqlevel');


my @paths = @{$ca->get_mapping_path($coord_system, $seq_coord_system)};

if(!@paths || @paths == 0){
  throw("Can't produce a database dump of ".$cs_name." ".$seq_region_name.
        ":".$start.":".$end.":".$strand." if there is no mapping path ".
        "to the seqlevel coord system ");
}

my @pieces;
push(@pieces, $slice);
my %coord_system;
$coord_system{$cs_name} = [];
my $projection = bless([$start, $end, $slice],
                       "Bio::EnsEMBL::ProjectionSegment");
push(@{$coord_system{$cs_name}}, $projection);

PATH:foreach my $path(@paths){
  if($path->name eq $cs_name){
    next PATH;
  }
  my @projections = @{$slice->project($path->name, $path->version)};
  if(!$coord_system{$path->name}){
    $coord_system{$path->name} = [];
  }
  push(@{$coord_system{$path->name}}, @projections);
  foreach $projection(@projections){
    push(@pieces, $projection->to_Slice);
    print "Slice name ".$projection->to_Slice->name."\n" 
      if($verbose);
  }
}

#print STDERR "Have ".@pieces." slices\n";
#print STDERR "Have hash with ".keys(%coord_system)." keys\n";
my %partial_dumped;
SLICE:foreach my $consituent_slice(@pieces){
  print "Replacing ".$slicedump->slice->name." with ".
    $consituent_slice->name."\n" if($verbose);
  $slicedump->slice($consituent_slice);
  if(!$partial_dumped{$consituent_slice->name}){
    $partial_dumped{$consituent_slice->name} = {};
  }
  TABLE:foreach my $table(@partial_tables){
      if( $partial_dumped{$consituent_slice->name}->{$table}){
        next TABLE;
      }else{
        $partial_dumped{$consituent_slice->name}->{$table} = 1;
      }
      #print "Dumping table ".$table."\n";
      my $method = "dump_".$table."_table";
      print "Trying to use ".$method." method\n";
      if($slicedump->can($method)){
        $slicedump->$method;
      }else{
        #if(! -e ($dump_dir."/".$table)){
        if(! -e ($dump_dir."/".$table) || !$dumped_whole{$table}){
          #if(!$dumped_whole{$table}){
          $slicedump->dump_table($table);
          $dumped_whole{$table} = 1;
        }else{
          print "Have already dumped table ".$table."\n" if($verbose);
        }
      }
    }
}


sub setup_tablelist{
  my ($hash, $table_list) = @_;
  foreach my $table(@$table_list){
    if(!$hash->{table}){
      $hash->{$table} = 1;
    }
  }
  return $hash;
}
