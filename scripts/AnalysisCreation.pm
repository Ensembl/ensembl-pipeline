# AnalysisCreation, class for creating analysis tables from config and 
# writing config from analysis tables and blantanly plager from code
# written by Glenn Procter while designig a pipeline alternative
#
# Cared for by Laura Clarke
#
# Copyright Laura Clarke
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

AnalysisCreation

=head1 SYNOPSIS

This class will parse config files to produce analysis objects and store them
and will take analysis tables and write a config file

[RepeatMask]
db=repbase
db_version=020713
db_file=repbase
program=RepeatMasker
program_version=1
program_file=RepeatMasker
parameters=-low, -lib, /acari/work5a/lec/briggsae.lib 
module=RepeatMasker
module_version=1
gff_source=RepeatMasker
gff_feature=Repeat
type=CONTIG


=head1 DESCRIPTION



=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

the class itself obviously doesn't' need to be instantiated but the either
the script which uses it should be in the same directory as it or the 
directory which contains it should be in you PERL5LIB

the analysis_setup script which should be found in the directory can
perform both functions for you if you have the appropriate database
and config files

=cut
package AnalysisCreation;

use strict;
use warnings;
use Bio::EnsEMBL::Pipeline::Analysis;
require Exporter;


our @ISA = qw(Exporter);
our @EXPORT = qw(parse_files write_into_db read_db write_file);




=head2 parse_files

  Arg [1]   : array of filenames
  Function  : parse a analysis config file and produce analysis objects
  Returntype: an array ref for an array of analyses objects
  Exceptions: if file doesn't exist
              if config format is in correct
              if key already exists for a particular header'
  Caller    : 
  Example   : my @analyses = @{&parse_files($file)};

=cut



sub parse_files {

  my @files = shift;

  my %headers;     # will store names of headers and number of keys for each
  my $config = {};
  # read each file

  foreach my $file (@files) {

    if (! -e $file) {
      die("analysis file $file not found\n");
    }
    my $header = "";

    open (FILE, $file) or die "Couldn't open file $file";
    while (<FILE>) {
      chomp();

      # Comment or blank line
      next if (/^\s$/ || /^\#/);

      # [HEADER]
      if (/^\[(.*)\]$/) {         # $1 will be the header name, without the []
	$header = $1;
	$headers{$header} = 0;
	#print "Reading stanza $header\n";
      }

      # key=value
      if (/^([^=\s]+)\s*=\s*(.+)/) {   # $1 = key, $2 = value

	my $key = lc($1);           # keys stored as all lowercase, values have case preserved
	my $value = $2;

	if (length($header) == 0) {
	  die("Found key/value pair $key/$value outside stanza");
	}

	#print "Key: $key Value: $value\n";

	# Check if this header/key is already defined
	if (exists($config->{$header}->{$key})) {
	  die("$key is already defined for [$header]; cannot be redefined");
	} else {
	  # store them in the config hash
	  $config->{$header}->{$key} = $value;
	  #print "$header:$key=$value\n";
	  $headers{$header}++;  # will be used to check for headers with no keys
	}

      }

    } # while <FILE>

    close FILE;
  }

  my @analyses;
  # add a blank key/value for any headers that have no keys
  foreach my $h (keys (%headers)) {
    if ($headers{$h} == 0) {
      print STDERR "you seem to have no values for $h will create a".
	"analysis object with just a logic_name";
    }
    my $analysis = Bio::EnsEMBL::Pipeline::Analysis->new
      (
       -db              => $config->{$h}->{db},
       -db_file         => $config->{$h}->{db_file},
       -db_version      => $config->{$h}->{db_version},
       -program         => $config->{$h}->{program},
       -program_version => $config->{$h}->{program_version},
       -program_file    => $config->{$h}->{program_file},
       -gff_source      => $config->{$h}->{gff_source},
       -gff_feature     => $config->{$h}->{gff_feature},
       -module          => $config->{$h}->{module},
       -module_version  => $config->{$h}->{module_version},
       -parameters      => $config->{$h}->{parameters},
       -logic_name      => $h,
       -input_id_type   => $config->{$h}->{type},
    );
    push(@analyses, $analysis);
  }

  return \@analyses;

}



=head2 write_into_db

  Arg [1]   : Bio::EnsEMBL::DBSQL::DBAdaptor
      [2]   : Ref to an array of analysis objects
  Function  : Write the analysis objects into the database
  Returntype: N/A
  Exceptions: if dbadaptor is the wrong type of object
  Caller    : 
  Example   : &write_into_db($db, \@analyses);

=cut



sub write_into_db{
  my $db = shift;
  my $analyses = shift;
  
  #print "have analysis adaptor ".$analysis_adaptor."\n";
  if(!($db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor'))){
    die("need a Pipeline::DBAdaptor not ".$db);
  }

  my $analysis_adaptor = $db->get_AnalysisAdaptor;
  my $sql = "select analysis_id from analysis where logic_name = ?";
  my $sth = $db->prepare($sql);
  ANALYSIS:foreach my $a(@$analyses){ 
    $sth->execute($a->logic_name);
    my ($analysis_id)= $sth->fetchrow;
    if($analysis_id){
      if($a->input_id_type){
	if($db->isa('Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor')){
	  my $sql = "select input_id_type from input_id_type_analysis ".
	    "where analysis_id = ?";
	  my $sth = $db->prepare($sql);
	  $sth->execute($analysis_id);
	  my ($type) = $sth->fetchrow;
	  if($type){
	    die("need ".$type." to be the same as ".$a->type) 
	      unless($type eq $a->type);
	  }else{
	    my $stored_sql = "insert into input_id_type_analysis ".
	      "(analysis_id, input_id_type) values(?, ?)";
	    my $stored_sth = $db->prepare($stored_sql);
	    $stored_sth->execute($analysis_id, $a->type);
	  }
	}
      }
      next ANALYSIS;
    }else{
      $analysis_adaptor->store($a);
    }
  }
}



=head2 read_db

  Arg [1]   : Bio::EnsEMBL::DBAdaptor
  Function  : Read the analysis objects from the database
  Returntype: array ref of analysis objects
  Exceptions: if db isn't the correct type'
  Caller    : 
  Example   : my $analyses = &read_db($db);

=cut



sub read_db{
  my $db = shift;

  if(!($db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor'))){
    die("need a DBAdaptor not ".$db);
  }
 
  my $analysis_adaptor = $db->get_AnalysisAdaptor;

  return $analysis_adaptor->fetch_all;
}


=head2 write_file

  Arg [1]   : filename
      [2]   : arrayref of analysis objects
  Function  : write a config file for the objects given
  Returntype: N/A
  Exceptions: if file doesnt exist
  Caller    : 
  Example   : &write_file($file, $analyses);

=cut



sub write_file{
  my $file = shift;
  my $analyses = shift;
  print STDERR "opening ".$file."\n";
  open (FH, '>'.$file) or die ("couldn't open $file to write to");
  foreach my $a(@$analyses){
    print FH "[".$a->logic_name."]\n";
    print FH "db=".$a->db."\n" if($a->db);
    print FH "db_version=".$a->db_version."\n" if($a->db_version);
    print FH "db_file=".$a->db_file."\n" if($a->db_file);
    print FH "program=".$a->program."\n" if($a->program);
    print FH "program_version=".$a->program_version."\n" if($a->program_version);
    print FH "program_file=".$a->program_file."\n" if($a->program_file);
    print FH "parameters=".$a->parameters."\n" if($a->parameters);
    print FH "module=".$a->module."\n" if($a->module);
    print FH "gff_source=".$a->gff_source."\n" if($a->gff_source);
    print FH "gff_feature=".$a->gff_feature."\n" if($a->gff_feature);
    if($a->can("input_id_type")){
      print FH "type=".$a->input_id_type."\n" if($a->input_id_type);
    }
    print FH "\n\n";
  }
  
}
