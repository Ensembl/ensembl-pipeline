#!/usr/local/bin/perl
#
#This version store all data from dbSNP as a hash, keyed by refsnpid.
#go through all data to find the same version first.
#go through second time to find diff version and to write them out.
#then using convert_diff.pl to convert contig coordinates to clone coordinates.
#then use find_diff_entries.pl to get a set unique entries to be read into dbSNP


use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB::CrossSNPMap;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::ExternalData::SNPSQL::DBAdaptor; 
use Bio::EnsEMBL::DBLoader;
use Getopt::Long;

$| = 1;


my $dbtype     = 'rdb';
my $host       = 'ecs1d';
my $host_snp   = 'ecs1d';
my $port       = '';
my $user       = 'ensro';
my $dbname     = 'homo_sapiens_core_130';
my $dbname_snp = 'homo_sapiens_snp_130';
my $module_snp = 'Bio::EnsEMBL::ExternalData::SNPSQL::DBAdaptor';

my $locator_snp = "$module_snp/host=$host_snp;port=$port;dbname=$dbname_snp;user=$user";

my $debug=1;
my $score =80;
my $masklevel = 80;
my $minmatch = 20;

my $start_refnum = $ARGV[0];
my $end_refnum = $ARGV[1];

my $file="/work4/yuan/SNP_MAP/SNP_EXIST_$start_refnum";
open SAME, ">$file" || die "$file: $1";
my $file1="/work4/yuan/SNP_MAP/SNP_DIFF_$start_refnum";
open DIFF, ">$file1" || die "$file1: $1";
my $file2="/work4/yuan/SNP_MAP/NO_CLONE_$start_refnum";
open NO_CLONE, ">$file2" || die "$file2: $1";


print STDERR "Using $locator_snp for snp db\n";

my $snpdb =  Bio::EnsEMBL::DBLoader->new($locator_snp); 
my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
					     -host => $host,
					     -user => $user,
					     -dbname => $dbname);

my (%no_clone, %diff);
my $first=1;
my $crossmap = Bio::EnsEMBL::Pipeline::RunnableDB::CrossSNPMap->new(-snpdb=>$snpdb, 
								    -db   =>$db,
								    -score =>$score, 
								    -minmatch =>$minmatch,
								    -masklevel =>$masklevel,
								    -debug => $debug,
								    -first =>$first
								   );
print STDERR "Fetching input from snpdb\n";


$crossmap->fetch_input($start_refnum,$end_refnum);

my (%all_version, %same_version, %diff_version, %diff_real);
my $same_version = $crossmap->same_version;
my $all_version = $crossmap->all_version;

%same_version = %$same_version;
%all_version = %$all_version;

foreach my $snp (values %same_version) {
  print SAME $snp->id,"\t",$snp->acc,"\t",$snp->version,"\t",$snp->start,"\t",$snp->end,"\t\tsnp_map\t",$snp->strand,"\n";
}

foreach my $refsnpid (keys %all_version) {
  if (!$same_version{$refsnpid}) {
    $diff_version{$refsnpid}=1;
  }
}

foreach my $refsnpid (keys %diff_version) {
  my $crossmap1 = Bio::EnsEMBL::Pipeline::RunnableDB::CrossSNPMap->new(-snpdb=>$snpdb,
								    -db   =>$db,
								    -score =>$score,
								    -minmatch =>$minmatch,
								    -masklevel =>$masklevel,
								    -debug => $debug,
								   );

  
  my ($noclone) = $crossmap1->fetch_input($refsnpid);
  if ($noclone) {
    $no_clone{$refsnpid}=$noclone;
  }
  else {
    print STDERR "Running mapping for diff clone $refsnpid\n";
    $crossmap1->run;
    print STDERR "Writing output for diff clone $refsnpid\n";
    my ($last_line) = $crossmap1->write_output;
    if ($last_line) {
      my @all = split /\t/, $last_line;
      print DIFF "$last_line\n";
      $diff_real{$all[0]}=1;
    }
  }
}

my ($no_clone_count, $same_version_count, $diff_version_count, $no_flank_count, $tot_count, $no_cross_count);
foreach my $refsnpid (keys %all_version) {
  $tot_count++;
  if ($no_clone{$refsnpid}) {
    print NO_CLONE "$refsnpid ",$no_clone{$refsnpid}->acc,"\n";
    $no_clone_count++;
  }
  elsif ($same_version{$refsnpid}) {
    $same_version_count++;
  }
  elsif ($diff_version{$refsnpid}) {
    if ($diff_real{$refsnpid}) {
      $diff_version_count++;
    }
    else {
      $no_cross_count++;
    }
  }
}

close(SAME)     || die $!;
close(DIFF)     || die $!;
close(NO_CLONE) || die $!;

print STDERR "total count is $tot_count, $same_version_count have same versions, $diff_version_count have diff clones and $no_clone_count have no clones, $no_cross_count have no_cross_match results\n";










