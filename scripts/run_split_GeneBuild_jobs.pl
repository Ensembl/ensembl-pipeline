#!/usr/local/ensembl/bin/perl




use strict;
use Getopt::Long;

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::General qw(
							  GB_INPUTID_REGEX 
							 );
use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Databases qw (
					 GB_DBHOST
					 GB_DBNAME
					 GB_DBUSER
					 GB_DBPASS
                                         );
use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Scripts qw(GB_LENGTH_RUNNABLES);


use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBLoader;
my $dbtype = 'rdb';
my $host   = $GB_DBHOST;
my $port   = undef;
my $dbname = $GB_DBNAME;
my $dbuser = $GB_DBUSER;
my $dbpass = $GB_DBPASS;

my $runnable;
my $input_id;
my $analysis_logic_name;
my $write  = 0;
my $check  = 0;
my $params;
my $size;
my $pepfile;


&GetOptions( 
	     'input_id:s'  => \$input_id,
	     'runnable:s'  => \$runnable,
	     'analysis:n'  => \$analysis_logic_name,
             'write'       => \$write,
             'check'       => \$check,
             'split_size' => \$size,
             'parameters:s'=> \$params,
	     );

$| = 1;

die "No runnable entered" unless defined ($runnable);
(my $file = $runnable) =~ s/::/\//g;
require "$file.pm";

if ($check) {
   exit(0);
}

print STDERR "args: $host : $dbuser : $dbpass : $dbname : $input_id : $size";

if(!$size){
  print STDERR "Can't work as split_size isn't defined\n";
  exit(0);
}
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host             => $host,
    -user             => $dbuser,
    -dbname           => $dbname,
    -pass             => $dbpass,
    -perlonlyfeatures => 0,
);

die "No input id entered" unless defined ($input_id);

my $analysis = $db->get_AnalysisAdaptor->fetch_by_logic_name($analysis_logic_name);

my %hparams;
# eg -parameters param1=value1,param2=value2
if(defined $params){
  foreach my $p(split /,/, $params){
    my @sp = split /=/, $p;
    $sp[0] = '-' . $sp[0];
    $hparams{$sp[0]} =  $sp[1];
  }
}


my ($chr, $start, $end) = $input_id =~ /$GB_INPUTID_REGEX/;

my $new_start = $start;
my $new_end   = $start + ($size -1);

while ( $new_end <= $end && $new_start < $end ){
  
   my $split_input_id  = "$chr.$new_start-$new_end";

   my $runobj = "$runnable"->new(-db    => $db,
				 -input_id => $split_input_id,
				 -analysis => $analysis,
				 %hparams,
				);
   
   $runobj->fetch_input;
   
   $runobj->run;
   
 
   my @out = $runobj->output;

   if ($write) {
     $runobj->write_output;
   }
   
   $new_start = $new_start + $size;
   $new_end   = $new_end   + $size;
   if ( $new_end > $end ){
     $new_end = $end;
   }
   
 }


my $sql = "insert into input_id_analysis(input_id, analysis_id, created) values(?, ?, now())";

my $sth = $db->prepare($sql);

$sth->execute($input_id, $analysis->logic_name); 
