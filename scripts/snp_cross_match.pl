#!/usr/local/bin/perl
#
#This version store all data from dbSNP as a hash, keyed by refsnpid.
#go through all data to find the same version first.
#go through second time to find diff version and to write them out.
#then using convert_diff.pl to convert contig coordinates to clone coordinates.
#then use find_diff_entries.pl to get a set unique entries to be read into dbSNP


use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB::CrossSNPMap;
use Bio::EnsEMBL::DBSQL::CrossMatchDBAdaptor;
use Bio::EnsEMBL::ExternalData::SNPSQL::DBAdaptor; 
use Bio::EnsEMBL::DBLoader;
use Getopt::Long;

$| = 1;


my $dbtype     = 'rdb';
my $host_snp   = 'ecs1a';
my $port       = '';
my $dbname_snp = 'snp120';
my $module     = 'Bio::EnsEMBL::DBSQL::CrossMatchDBAdaptor';
my $module_snp = 'Bio::EnsEMBL::ExternalData::SNPSQL::DBAdaptor';

my $locator_cross = "$module/host=$host_snp;port=$port;dbname=$dbname_snp;user=ensro";
my $locator_snp = "$module_snp/host=$host_snp;port=$port;dbname=$dbname_snp;user=ensro";

my $debug=2;
my $score =45;
my $masklevel = 80;
my $minmatch = 20;

my $start_refnum = $ARGV[0];
my $end_refnum = $ARGV[1];

my $file="/work4/yuan/SNP_MAP2/SNP_EXIST_$start_refnum";
open SAME, ">$file" || die "$file: $1";
my $file1="/work4/yuan/SNP_MAP2/SNP_DIFF_$start_refnum";
open DIFF, ">$file1" || die "$file1: $1";
my $file2="/work4/yuan/SNP_MAP2/NO_CLONE_$start_refnum";
open NO_CLONE, ">$file2" || die "$file2: $1";


print STDERR "Using $locator_snp for crossmatch db\n";
my $crossdb =  Bio::EnsEMBL::DBLoader->new($locator_cross);
my $snpdb =  Bio::EnsEMBL::DBLoader->new($locator_snp); 

print STDERR "Fetching input from snpdb\n";

my @infos = $snpdb->get_snp_info_between_two_refsnpid($start_refnum,$end_refnum);

my ($new_count, $old_count, $same_count, %id_acc, %acc_info, %done, @bases, %no_clone, %rec);

foreach my $info (@infos) {
  
  #my ($refsnpid,$snpclass,$snptype,$observed,$seq5,$seq3,$acc,$version,$start,$end,$strand) = @$info;
  if ($info) {
    my $refsnpid = $$info[0];
    my $mapweight = $$info[2];
    my $acc      = $$info[6];
    my $version  = $$info[7];
    if ($mapweight <=2) {
      push (@{$id_acc{$refsnpid}}, "$acc-$version"); 
      my $key = "$acc-$refsnpid-$version";
      $acc_info{$key} = $info;
    }
  }
}

REF : foreach my $refsnpid (sort numeric keys %id_acc) {
  if ($refsnpid) {
    print "refsnpid is $refsnpid and @{$id_acc{$refsnpid}}\n";
    #next CLONE;
    CLONE : foreach my $clone (@{$id_acc{$refsnpid}}) {
      if (!$done{$refsnpid}) {
	my ($acc, $version) = split /\-/, $clone;
	my $key = "$acc-$refsnpid-$version";
	if ($acc_info{$key}) {
	  my $first=1;
	  my $crossmap = Bio::EnsEMBL::Pipeline::RunnableDB::CrossSNPMap->new(-crossdb=>$crossdb, 
									      -score =>$score, 
									      -minmatch =>$minmatch,
									      -masklevel =>$masklevel,
									      -debug => $debug,
									      -first =>$first
									     );
	  my ($done,$start,$end,$strand) = $crossmap->check_clone_in_GP($acc_info{$key});
	  if ($done =~ /^done/) {
	    $done{$refsnpid}= "same version";
	    if ($debug>1) {
	      print STDERR "$refsnpid $clone has same version in first round\n";
	    }
	    print SAME "$refsnpid\t$acc\t$version\t$start\t$end\tsnp_map\t$strand\n";
	    next REF;
	  }
	  elsif ($done =~ /next/) {
	    next CLONE;
	  }
	}
      }
    }
  }
}

REF1 : foreach my $refsnpid (sort numeric keys %id_acc) {
  if ($refsnpid) {
    #print "refsnpid is $refsnpid and @{$id_acc{$refsnpid}}\n";
    #next CLONE;
    CLONE1 : foreach my $clone (@{$id_acc{$refsnpid}}) {
      if (!$done{$refsnpid}) {
	my ($acc, $version) = split /\-/, $clone;
	my $key = "$acc-$refsnpid-$version";
	if ($acc_info{$key}) {
	  my $crossmap = Bio::EnsEMBL::Pipeline::RunnableDB::CrossSNPMap->new(-crossdb=>$crossdb,
									      -score =>$score,
									      -minmatch =>$minmatch,
									      -masklevel =>$masklevel,
									      -debug => $debug,
									     );

	  my ($done,$start,$end,$strand) = $crossmap->check_clone_in_GP($acc_info{$key});
	  if ($done =~ /^done/) {
	    $done{$refsnpid}= "same version";
	    if ($debug>1) {
	      print STDERR "$refsnpid $clone has same version\n";
	    }
	    print SAME "$refsnpid\t$acc\t$version\t$start\t$end\tsnp_map\t$strand\n";
	    next REF1;
	  }
	  elsif ($done =~ /diff/) {
	    print STDERR "Running mapping for clone $refsnpid $acc\n";
	    $crossmap->run;
	    print STDERR "Writing output for clone $refsnpid $acc\n";
	    my ($last_line) = $crossmap->write_output;
	    if ($last_line) {
	      print DIFF "$last_line\n";
	      $done{$refsnpid}="diff";
	      next REF1;
	    }
	    else {
	      next CLONE1;
	    }
	  }
	  elsif ($done =~ /clone/) {
	    if ($debug>1) {
	      print STDERR "$refsnpid $acc is not in our database\n";
	    }
	    push (@{$no_clone{$refsnpid}}, "no clone $acc\n");
	    next CLONE1;
	  }
	  elsif ($done =~ /flank/) {
	    if ($debug>1) {
	      print STDERR "no flank seqs for $refsnpid $acc in our database\n";
	    }
	    push (@{$no_clone{$refsnpid}}, "no flank $acc-$version\n");
	    next CLONE1;
	  }
	  elsif ($done == "") {
	    if ($debug>1) {
	      print STDERR "Do not know the reseason for $refsnpid $acc\n";
	    }
	    push (@{$no_clone{$refsnpid}}, "no reason $acc-$version\n");
	    next CLONE1;
	  }
	}
	else {
	  print "There is no info for $acc-$refsnpid\n";
	  next CLONE1;
	}
      }
    }
  }
}

sub numeric {$a<=>$b}

my ($no_clone_count, $same_version_count, $diff_version_count, $no_flank_count, $tot_count, $no_cross_count);
foreach my $refsnpid (keys %id_acc) {
  $tot_count++;
  if (!$done{$refsnpid}) {
    if ($no_clone{$refsnpid}) {
      print NO_CLONE "$refsnpid @{$no_clone{$refsnpid}}\n";
      $no_clone_count++;
    }
    else {
      $no_cross_count++;
      print "no CROSS_MATCH results for $refsnpid\n";}
    
  }
  elsif ($done{$refsnpid} =~ /same/) {
    $same_version_count++;
  }
  elsif ($done{$refsnpid} =~ /diff/) {
    $diff_version_count++;
  }
  
}

close(SAME)     || die $!;
close(DIFF)     || die $!;
close(NO_CLONE) || die $!;

print STDERR "total count is $tot_count, $same_version_count have same versions, $diff_version_count have diff clones and $no_clone_count have no clones and $no_cross_count have no_cross_match results\n";










