#!/usr/bin/env perl

#Written in June 2010, based on make_input_ids, to make input IDs for working with the 
#patches introduced into the human genome in patch release 1 of GRCh37

=head1 NAME

make_input_ids  - handles insertion, deletion and listing of rules in a database

=head1 SYNOPSIS

make_input_ids -dbhost ecs1a -dbuser ensadmin -dbpass **** -dbname pipeline_db -contig

=head1 DESCRIPTION

this script allows input_ids to be generated and written to the input_id_analysis table for

=head1 OPTIONS

    -dbhost    host name for database (gets put as host= in locator)

    -dbport    For RDBs, what port to connect to (port= in locator)

    -dbname    For RDBs, what name to connect to (dbname= in locator)

    -dbuser    For RDBs, what username to connect as (user= in locator)

    -dbpass    For RDBs, what password to use (pass= in locator)

    -help      Displays script documentation with PERLDOC


    -logic_name the logic_name of the analysis object which needs to be 
                associated with these entries 
    
    -coord_system the coordinate system you want slices in
    
    -coord_system_version the version of the coord system you want

    -input_id_type if you want to specific an input_id_type not already
     used by the analysis object
    -insert_analysis if you want to insert an analysis object if it doesn't
     already exist in the database'
    -top_level this will fetch all the non_redundant pieces in
     the database this may produce ids which are a mixture of different
     coordinate systems, if -coord_system_version is specified it will
     be ignored

=head1 EXAMPLES


perl  make_input_ids -dbhost <HOST> -dbname <DBNAME> -dbport <PORT> 
         -dbuser <USER> -dbpass ***** 
           -slice -slice_size 150000 
             -coord_system chromosome 
              -logic_name SubmitSlice150k 
                -input_id_type Slice150k

This will create input-id's in the format  'chromosome:BDGP4.1:2L:1:149386:1' of size 150k
Remember to create an analysis 'SubmitSlice150k'(module 'Dummy') and an entry in 
input_id_type_analysis as well !



./make_input_ids -dbhost host -dbuser user -dbpass *** -dbport 3306 
  -dbname my_database -contig

this will use all the contig names are input_ids


=cut




use warnings ;
use Bio::EnsEMBL::Pipeline::Analysis;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Utils::InputIDFactory;
use Bio::EnsEMBL::Pipeline::DBSQL::StateInfoContainer;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use strict;
use Getopt::Long qw(:config no_ignore_case);

my $host;
my $user;
my $pass;
my $port = 3306;
my $dbname;
my $logic_name;
my $coord_system;
my $coord_system_version;
my $input_id_type;
my $help = 0;
my $insert_analysis = 0;
my $date;
my $patch_novel;
my $patch_fix;
my $all_patch_types;
my $on_or_after;
my $all_dates;
my $slice_size = undef;
GetOptions(
            'host|dbhost|h:s'     => \$host,
            'port|dbport|P:n'     => \$port,
            'user|dbuser|u:s'     => \$user,
            'pass|dbpass|p:s'     => \$pass,
            'dbname|db|D:s'     => \$dbname,
            'coord_system|cs_name:s'       => \$coord_system,
            'coord_system_version:s' => \$coord_system_version,
            'logic_name:s' => \$logic_name,
            'input_id_type:s' => \$input_id_type,
            'help'            => \$help,
            'insert_analysis!' => \$insert_analysis,
            'date:s' => \$date,
            'patch_fix!' => \$patch_fix,
            'patch_novel!' => \$patch_novel,
            'all_patch_types!' => \$all_patch_types,
            'on_or_after:s' => \$on_or_after,
            'all_dates!' => \$all_dates,
            'slice_size:n' => \$slice_size,
           );


if ($help) {
    exec('perldoc', $0);
}

if(!$host || !$user || !$dbname){
  warn("Can't run without -dbhost $host -dbuser $user -dbname $dbname");
  $help = 1;
}

#which types of patch
throw "Must specify a patch type or types: -patch_novel -patch_fix or -all_patch_types\n" if !$patch_novel and !$patch_fix and !$all_patch_types;
throw "You've specified overlapping patch type options\n"if ($patch_novel or $patch_fix) and $all_patch_types;
my $code_patch_novel = 'patch_novel';
my $code_patch_fix = 'patch_fix';
my @pt;
if($all_patch_types){
  push @pt, $code_patch_novel, $code_patch_fix;
}
else{
  if($patch_novel){
    push @pt, $code_patch_novel;
  }
  if($patch_fix){
    push @pt, $code_patch_fix;
  }
}

#which date
throw "Must specify -on_or_after or -all_dates (format yyyy-mm-dd)\n" if !$on_or_after and !$all_dates;
throw "You've specified -on_or_after and -all_dates\n" if $on_or_after and $all_dates;
my $year;
my $month;
my $day;

if($on_or_after){
  ($year, $month, $day) = split /-/, $on_or_after;
  throw "Invalid year\n" if(length($year)!= 4);
  throw "Invalid month\n" if(length($month)!=2) or ($month > 12 or $month < 1);
  throw "Invalid day\n" if(length($day)!=2) or ($day >31 or $day < 1); 
}

my $db = new Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor(-host   => $host,
                                                      -user   => $user,
                                                      -pass   => $pass,
                                                      -port   => $port,
                                                      -dbname => $dbname);

if(!$logic_name){
  throw("Can't run with out a logic_name");
  exec('perldoc', $0);
}

my @input_id_slice_names;

my $sa = $db->get_SliceAdaptor();
#get non-ref but not duplicates
my @slices = @{$sa->fetch_all('toplevel',undef, 1)};
print scalar(@slices)."\n";
foreach my $slice (@slices){
  foreach my $type (@pt){
    my @slice_attributes = @{$slice->get_all_Attributes($type)};
    if (scalar(@slice_attributes) > 0){
      if($on_or_after){
        foreach my $attrib (@slice_attributes){
          my $value = $attrib->value;
          my ($patch_date, $patch_time) = split " ", $value;
          my ($patch_year, $patch_month, $patch_day) = split /-/, $patch_date;
          throw "Patch attrib value is not a date/timestamp\n" if (length($patch_year)!=4 or length($patch_month)!=2 or length($patch_day)!=2);
          throw "Patch month isn't a month\n" if $patch_month < 1 or $patch_month > 12;
          throw "Patch day isn't a day\n" if $patch_day < 1 or $patch_day > 31; 

          if($patch_year >= $year and $patch_month >= $month and $patch_day >= $day){
            store_slice($slice);
          }
        }
      }
      else{
        store_slice($slice);
      }
    }
  }
}

my $inputIDFactory = new Bio::EnsEMBL::Pipeline::Utils::InputIDFactory
  (
   -db => $db,
#   -top_level => $top_level,
#   -include_non_reference => $include_non_reference,
   -logic_name => $logic_name,
   -input_id_type => $input_id_type,
   -insert_analysis => $insert_analysis,
   -slice => 1,
   -coord_system => 'toplevel',
  );


$inputIDFactory->input_ids(\@input_id_slice_names);
$inputIDFactory->store_input_ids;

sub store_slice{
  my $slice = shift;

  #if slice size has been specified...
  if($slice_size){
    #...and slice is longer
    if($slice->length>$slice_size){
      my $start = 1;
      my $end = $slice_size;
      my $not_at_end = 1;
      while($not_at_end){
        my $sub_slice = $slice->sub_Slice($start, $end);
        push @input_id_slice_names, $sub_slice->name;
        print "Added ".$sub_slice->name."\n";

        #prepare for next section
        $start = $end + 1;
        $end = $end + $slice_size;

        #if next section will be the last one
        if($end>$slice->length){
          $end = $slice->length;
          $sub_slice = $slice->sub_Slice($start, $end);
          push @input_id_slice_names, $sub_slice->name;
          print "Added ".$sub_slice->name."\n";
          $not_at_end = 0;
        }
      }
    }
    else{
      push @input_id_slice_names, $slice->name;
      print "Added ".$slice->name."\n";
    }
  }
  #no slice size specified
  else{
    push @input_id_slice_names, $slice->name;
    print "Added ".$slice->name."\n";
  }

}
