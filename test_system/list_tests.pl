#!/usr/local/ensembl/bin/perl -w

use strict;
use Bio::EnsEMBL::Utils::Exception qw(throw warning verbose);
use Getopt::Long;
use File::Path;
my $species = 'homo_sapiens';
my $conf_file;
my $help;

&GetOptions(
            'species:s' => \$species,
            'conf_file:s' => \$conf_file,
            'help!' => \$help,
           ) or ($help = 1);


$conf_file = 'TestDB.conf' if(!$conf_file);



if($help){
  print "list_tests.pl is a script for displaying what tests can be run ".
    "for a given species this is what the commandline should look like\n ".
      "./list_tests.pl -species homo_sapiens -conf_file TestDB.conf\n";
  exit;
}


my $conf = do $conf_file;

my $data_dir = $conf->{'data_dir'};


if(!-d $data_dir){
  throw("Can't unzip data if data directory $data_dir doesn't exist");
}
my $curr_dir = $ENV{'PWD'};
my $zip_file = $species.".zip";
my $path = $data_dir."/".$zip_file;
my $dest_dir = $curr_dir."/".$species;

if(!-e $path){
  throw("Can't unzip ".$zip_file." if ".$path." doesn't exist");
}

if(! -d $dest_dir){
  mkdir($dest_dir);
  my $cmd = "unzip -q $path -d $dest_dir ";
  system($cmd) == 0 or throw("Error running ".$cmd);
}elsif(-d $dest_dir){
  print "Will use existing unzipped directory ".$dest_dir."\n";
}



my $analysis_file = $dest_dir."/analysis.sql";

open(FH, $analysis_file) or throw("Can't open ".$analysis_file);
printf("%-15s %-15s\n", 'logic_name', 'module');
printf("%-15s %-15s\n", '----------', '------');
while(<FH>){
  chomp;
  my @values = split /\t/, $_;
  my $logic_name = $values[2];
  if($logic_name =~ /Submit/){
    next;
  }
  my $module = $values[10];
  printf("%-15s %-15s\n", $logic_name, $module);
}


close(FH) or throw("Couldn't close ".$analysis_file);
