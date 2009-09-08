=pod

=head1 NAME slice_to_sql.pl

 mail ensembl-dev@ebi.ac.uk

=head1 SYNOPSIS

This is a script which given information about a slice and a
list of tables will dump information about that slice in the form
of text files which can then be reloaded into a DB using mysqlimport.

=head1 OPTIONS

DB Connection Details (compulsory)


   -dbhost     The host where the pipeline database is.
   -dbport     The port.
   -dbuser     The user to connect as.
   -dbpass     The password to use.
   -dbname     The database name.


Slice details

  -coord_system_name    (compulsory) the name of the coordinate system 
                         you want the slice from
  -coord_system_version the version of the coordinate system you want the
                        slice from.
                        This has to be specified in some cases where two 
                        parallel versions are present in the same DB, e.g.
                        in homo_sapiens_core_51_36m DB which contains
                        both NCBI35 and NCBI36 seq_region data, or else
                        the system may not fetch the correct slice and
                        it causes problems when getting the mapping path
                        between e.g. contigs vs a chromosome slice.
  -seq_region_name      (compulsory)
  -seq_region_start
  -seq_region_end
  -seq_region_strand    (default set to '1', i.e. forward/sense strand)

Other details

  -output_dir the directory the files will be dumped into
  -verbose prints information about what the script is doing

Tables to dump

The script requires a list of tables to dump.
There is a list of standard tables which involve sequence storage.
There are other groupings of tables like pipeline, raw compute and genes
which mean you don't need to specify every table required on the
commandline but just a single option. Those tables which are to be dumped
partially need to have a method implemented in SliceDump in the form
dump_TABLENAME_table otherwise the whole table will just be dumped.

  -whole_table           name of a table to dump the entire contents of. If specifying
                         multiple tables, separate table names with a comma with no extra space, e.g.
                         "-whole table table_A,table_B,table_C".

  -partial_table         name of the table to dump contents based on the
                         specified seq_region

  -pipeline              dump the pipeline tables (incl. 'rule_goal', 'rule_conditions', 
                        'input_id_type_analysis')

  -raw_computes          dump the raw compute tables (incl.'repeat_consensus', 'repeat_feature', 
                        'prediction_exon', 'prediction_transcript', 'dna_align_feature',
                         'protein_align_feature', 'simple_feature')

  -genes                 dump the gene tables (incl. 'gene', 'exon', 'transcript', 'translation', 
                        'exon_transcript', 'supporting_feature', 'protein_align_feature',       
                         'dna_align_feature', 'gene_stable_id', 'exon_stable_id', 'translation_stable_id',
                         'transcript_stable_id')

  -protein_annotation    dump the protein_annotation tables ('protein_feature' table)

  -no_defaults           this option means only the tables specified on the commandline with the 
                         options -whole_table or -partial_table will be dumped or those from the 
                         table_group flags ('pipeline', 'raw_computes', etc) but none of the 
                         standard tables.

  -all                   dump all the tables in the database as defined by "show tables;"


Each table which is dumped "whole" will be written into a single tab-delimited text file.

For a table which is dumped "partial" (e.g. where the requested slice is a specified region 
on a single chromosome), data for the requested slice will be dumped into 1 file, and data for
the constituent/component slices projected from the requested slice (e.g. the contigs which
have been projected from the chromosomal region) will be dumped into individual text files
(1 text file per component slice). Filename convention follows: table_name.slice_name.

For example, the 'dna' table dumped partially will generate files like these:

dna.AL645608.30.1.186759.955-186759
dna.AL645703.18.1.47568.2001-47568
dna.AL645728.31.1.87105.1-85105 

Concatenate these dna.* files will give the entire partially-dumped 'dna' table, ready for
importing into a new DB using mysqlimport.

The complusory options are the database arguments, slice name and slice coord system name.
All other arguments are optional.

=head1 EXAMPLES

To dump out data for positions 1-10,000,000 on human chromosome 1,  Partial tables 
for sequence/assembly-related data:

perl path/to/script/slice_to_sql2.pl -dbhost xxxx  -dbuser xxxx -dbport xxxx \
-dbpass xxxx -dbname homo_sapiens_core_51_36m -coord_system_name chromosome \
-coord_system_version NCBI36 -seq_region_name 1 -seq_region_start 1 \
-seq_region_end 10000000 -output_dir output/directory/location \
-whole_table meta,meta_coord,coord_system,attrib_type \
-partial_table seq_region,assembly,dna,seq_region_attrib \
-no_defaults -verbose


=head1 SEE ALSO

  Bio::EnsEMBL::Pipeline::Utils::SliceDump
 
=cut


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
my @whole_standard_tables =
  ( 'meta', 'meta_coord', 'coord_system', 'analysis', 'attrib_type' );

my @whole_pipeline_tables =
  ( 'rule_goal', 'rule_conditions', 'input_id_type_analysis' );

my @pmatch_tables = ( 'pmatch_feature', 'protein' );
my @partial_standard_tables = ( 'seq_region', 'assembly',
                                'dna',        'seq_region_attrib',
                                'assembly_exception' );

my @raw_compute_tables = ( 'repeat_consensus',  'repeat_feature',
                           'prediction_exon',   'prediction_transcript',
                           'dna_align_feature', 'protein_align_feature',
                           'simple_feature' );

my @gene_tables = ( 'gene',                  'exon',
                    'transcript',            'translation',
                    'exon_transcript',       'supporting_feature',
                    'protein_align_feature', 'dna_align_feature',
                    'gene_stable_id',        'exon_stable_id',
                    'translation_stable_id', 'transcript_stable_id' );

my @protein_annotation_tables = ('protein_feature');
my $raw_computes;
my $genes;
my $protein_annotation;
my @whole_commandline_tables;
my @partial_commandline_tables;
my $pipeline_tables;
my $pmatch;
my $no_defaults;
my $help;
my $all;
&GetOptions( 'dbhost:s'               => \$host,
             'dbport:n'               => \$port,
             'dbname:s'               => \$dbname,
             'dbuser:s'               => \$dbuser,
             'dbpass:s'               => \$dbpass,
             'coord_system_name:s'    => \$cs_name,
             'coord_system_version:s' => \$cs_version,
             'seq_region_name:s'      => \$seq_region_name,
             'seq_region_start:s'     => \$start,
             'seq_region_end:s'       => \$end,
             'seq_region_strand:s'    => \$strand,
             'verbose!'               => \$verbose,
             'output_dir:s'           => \$output_dir,
             'whole_table:s'          => \@whole_commandline_tables,
             'partial_table:s'        => \@partial_commandline_tables,
             'pipeline!'              => \$pipeline_tables,
             'raw_computes!'          => \$raw_computes,
             'genes!'                 => \$genes,
             'protein_annotation!'    => \$protein_annotation,
             'pmatch!'                => \$pmatch,
             'no_defaults!'           => \$no_defaults,
             'all!'                   => \$all,
             'help'                   => \$help
) or throw("Can't get options");

if ($help) {
    exec('perldoc', $0);
}

unless ($host && $dbname && $dbuser) {
  throw("Can't run without database argument, use -help option ".
        "to see commandline args");
}

unless($output_dir){
  throw("You must specify an output directory on the command line with the ".
        "option -output_dir, use -help to get docs");
}

if(!$seq_region_name || !$cs_name){
  throw("You must specify a seq_region_name and a coord_system_name on the ".
        "commandline use -help for more information about commandline options");
}


if($all){
  print "All tables will be dumped regardless of your commandline ".
    "options\n";
}

if($all && $no_defaults){
  throw("You cannot specify both -all and -no_defaults as they are mutually ". 
        "exclusive options. See -help to read docs");
}



@whole_commandline_tables = split(/,/,join(',', @whole_commandline_tables));

@partial_commandline_tables = split (/,/, join(',', @partial_commandline_tables));

# check if the table names have been split up properly:

print "The whole_commandline_tables are:\n";

foreach (@whole_commandline_tables) {
        print "$_\n" if($verbose);
}

print "\nThe partial_commandline_tables are:\n";

foreach (@partial_commandline_tables) {
        print "$_\n\n" if($verbose);
}


my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -dbname => $dbname,
    -host   => $host,
    -user   => $dbuser,
    -port   => $port,
    -pass   => $dbpass,
);

# An output directory named after the dumped DB will be created if it doesn't exist. 
# A "mkdir" statement will be printed out as a reminder.

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


# Creating slicedump object:
# Here the exact slice required for partial dump has not be specified.
# The slicedump object created will be mainly for dumping whole tables.
# Exact "slice" info will be specified for the slicedump object later, just before dumping partial tables.

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
  throw("Have been unable to fetch slice for ". $seq_region_name." ".$cs_name);
       
}

# Check which whole/partial tables are to be dumped, whether it's specified by
# the user, dumping all tables, or dumping default tables, etc.

my %whole_tables;  
my %partial_tables;

# The %whole_tables and %partial_tables hashes are set up by &setup_tablelist. This subroutine
# takes a list ref of table names (e.g. \@whole_commandline_tables or \@raw_compute_tables defined at 
# the beginning of this script), then uses each of the table names as individual keys in an 
# anonymous hash and assigns the value of "1" for each key.  The value is a dummy, what matters is the
# list of keys which will be eventually used in sqldump.

if($no_defaults){  #If all required tables (whole and/or partial) are explicitly specified...
  %whole_tables = %{setup_tablelist(\%whole_tables, 
                                    \@whole_commandline_tables)};
  %partial_tables = %{setup_tablelist(\%partial_tables, 
                                      \@partial_commandline_tables)};
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
  if($pmatch){
     %partial_tables = %{setup_tablelist(\%partial_tables, 
                                         \@pmatch_tables)};
  }
 
  @whole_tables = keys(%whole_tables);
  @partial_tables = keys(%partial_tables);
}

elsif($all){  # All tables as defined by "show tables" will be specified, some dumped whole, some partially. (see below)
  my @tables = @{get_tables($db)}; # &get_tables retrieves a list of table names from a DB like "show tables;" in mysql
  foreach my $table(@tables){
    my $method = "dump_".$table."_table";  #slicedump object can call methods like "dump_gene_table"
    if($slicedump->can($method)){
      push(@partial_tables, $table);  #dump partial tables if slice-specific data are available (e.g. seq_region table)
    }
    else{
      push(@whole_tables, $table);  #dump out whole tables if there's no slice-specific data (e.g. meta table)
    }
  }
}

# The final "else" loop below is executed when the user hasn't sepcified either "all" nor "no_defaults"
# on the command line. If there are command line arguments for "whole_table" option, they will be
# taken, or else the script uses the set of "standard" whole tables as hardcoded at the beginning
# of the script. The same logic applies to the partial tables.

else {
  %whole_tables = %{setup_tablelist(\%whole_tables, 
                                    \@whole_commandline_tables)};

  %whole_tables = %{setup_tablelist(\%whole_tables, 
                                    \@whole_standard_tables)};
 
  %partial_tables = %{setup_tablelist(\%partial_tables, 
                                      \@partial_commandline_tables)};

  %partial_tables = %{setup_tablelist(\%partial_tables, 
                                      \@partial_standard_tables)};
  
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
  if($pmatch){
     %partial_tables = %{setup_tablelist(\%partial_tables, 
                                         \@pmatch_tables)};
  }

  @whole_tables = keys(%whole_tables);
  @partial_tables = keys(%partial_tables);
}

# After confirming the names of the tables to be dumped, dump out whole tables as 
# specified in @whole_tables unless the table has already been dumped previously:

my %dumped_whole;

foreach my $table(@whole_tables){
  print STDERR "Looking at table ".$table."\n";  
  if(!$dumped_whole{$table} || -e ($dump_dir."/".$table)){  
    $slicedump->dump_table($table);
    $dumped_whole{$table} = 1;
  }else{
    print "Have already dumped table ".$table."\n" if($verbose);
  }
}


# Before dumping out partial tables, one has to first get the mapping path between toplevel 
# (chromosome) seq_regions and e.g. contigs:

my $coord_system = $ca->fetch_by_name($cs_name,$cs_version); 

my $seq_coord_system = $ca->fetch_by_name('seqlevel');

my @paths = @{$ca->get_mapping_path($coord_system, $seq_coord_system)};# @paths is a CoordSystem object.

# Next, check if the paths are defined properly.

# This check is done because in CoordSystemAdaptor.pm, an undef element was 
# spliced into the 2nd element of @coord_systems. As $coord_system is passed 
# to the method get_mapping_path,it results in the 2nd element of @paths always 
# being undef.  This will have a knock-on effect later on when we call methods 
# foreach (@paths) as the script aborts.  

print "The number of mapping path objects is:";
print scalar(@paths)."\n\n";

foreach my $p(@paths) {
   if(defined $p) { 
     print "$p ".ref($p)." ".$p->name."\n" if($verbose); 
# Should print out sth like --- Bio::EnsEMBL::CoordSystem=HASH(0xfb9f00) Bio::EnsEMBL::CoordSystem chromosome 
   } 
   else {  
     print "path is not defined ...\n" if($verbose) ; 
   } 
}

if(!@paths || scalar(@paths) < 1){   # throw error if there are no mapping paths
         throw("Can't produce a database dump of ".$cs_name." ".$seq_region_name.
        ":".$start.":".$end.":".$strand." if there is no mapping path ".
        "to the seqlevel coord system ");
} 
elsif (scalar(@paths) >2) {
  print STDERR "Have more than 2 elements in the paths list ref.\n" if($verbose);
} 
else {
  print STDERR "Expecting 2 elements in the list ref, but have ".scalar(@paths)."  instead.\n"
  if($verbose);
}


my @pieces;
push(@pieces, $slice);  #  $slice was created at the beginning of the code 
my %coord_system;
$coord_system{$cs_name} = [];  
my $projection = bless([$start, $end, $slice],
                       "Bio::EnsEMBL::ProjectionSegment");
push(@{$coord_system{$cs_name}}, $projection);  

PATH:foreach my $path(@paths){ 

  unless ( $path ) {
    next PATH ;
  }

  if($path->name eq $cs_name){   #$path->name will be sth like "chromosome", "contig".
  next PATH; # Nothing to project if path name = coord_sys name!
  }
 
  # project the "$slice" (what we provided at the start of the script, e.g. "chromosome") 
  # onto the coord_sys of the $path (e.g. "contig"):
  my @projections = @{$slice->project($path->name, $path->version)};

  if(!$coord_system{$path->name}){
    $coord_system{$path->name} = []; 
  }

  push(@{$coord_system{$path->name}}, @projections);  

  foreach $projection(@projections){
    push(@pieces, $projection->to_Slice);  # @pieces collects projected slice objects
    print "Slice name ".$projection->to_Slice->name."\n"   #print out e.g."Slice name contig...."
      if($verbose);
  }
}

print "Have ".@pieces." slices\n" if($verbose);  
# When projecting a chromosome slice onto contigs, @pieces should contain 
# all commandline-specified slice (1 piece) +  constituent slices. 

# Next, for each "partial" table, work out what coord_system's data to be dumped out. 
# (E.g. chr-level data only when contig-level data don't exist?)
# The table (key) and "relevant" coord_sys (value) pairs will be stored in %table_coord_systems.

my %table_coord_systems;
my $meta_container = $db->get_MetaCoordContainer;

TABLE:foreach my $table(@partial_tables){
  my @coord_systems = @{$meta_container-> fetch_all_CoordSystems_by_feature_type($table)};
  if(@coord_systems == 0){
    next TABLE;
  }
  foreach my $cs(@coord_systems){
    if(!$table_coord_systems{$table}){
      $table_coord_systems{$table} = {};
    }
    $table_coord_systems{$table}->{$cs->name} = 1;
  }
}

my %partial_dumped;

SLICE:foreach my $consituent_slice(@pieces){

  $slicedump->slice($consituent_slice); #$slicedump previously defined only by DB and dump_dir, now it knows its slice too
  print "Replacing ".$slicedump->slice->name." with " if($verbose); 
  print $consituent_slice->name."\n" if($verbose); 

  if(!$partial_dumped{$consituent_slice->name}){
    $partial_dumped{$consituent_slice->name} = {};
  }

 #first check if the table name - coord_sys relationship has been stored
 #then check if table can be dumped out for const.slice's coord_sys
 
   TABLE:foreach my $table(@partial_tables){
      if($table_coord_systems{$table}){  
        if(!$table_coord_systems{$table}->{$consituent_slice->coord_system->name}){
          print "Can't dump ".$table." for ".
            $consituent_slice->coord_system->name."\n" if($verbose);
          next TABLE;
        }
      }
      
      # The rest of the code follows the same logic as dumping whole tables 
 
      if( $partial_dumped{$consituent_slice->name}->{$table}){
        print STDERR "Have already dumped ".$table." for ".
          $consituent_slice->name."\n" if($verbose);
        next TABLE;
      }
      else{
        $partial_dumped{$consituent_slice->name}->{$table} = 1;
      }
      
      my $method = "dump_".$table."_table";
      
      if($slicedump->can($method)){
        $slicedump->$method;
      }
      
      else{
          if(! -e ($dump_dir."/".$table) || !$dumped_whole{$table}){ 
          $slicedump->dump_table($table);
          $dumped_whole{$table} = 1;
        }
          else{
          print "Have already dumped table ".$table."\n" if($verbose);
           }
      }
    }
}

# A subroutine used to set up the list of required table names. What matters are the keys in $hash.
sub setup_tablelist{
  my ($hash, $table_list) = @_;
  foreach my $table(@$table_list){
    if(!$hash->{table}){
      $hash->{$table} = 1;
    }
  }
  return $hash;
}

# A subroutine used to retrieve a list of all table names from a given DB
sub get_tables {
  my ($db) = @_;

 	my $query = "show tables";
  
	my $sth = $db->prepare($query);
	my $res = $sth->execute;
  
	my @tables;
  
	while (my $ref = $sth->fetchrow_arrayref) {
    push(@tables,$ref->[0]);
	}

  return \@tables;
}

