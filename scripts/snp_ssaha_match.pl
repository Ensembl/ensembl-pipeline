#! /usr/local/bin/perl
#
#This script run after snp_cross_match.pl. snp_cross_match.pl generate a list of refsnpid for which none of the clones #are on GP. snp_ssaha_match.pl by using ssaha method to locate the contigs on GP and then use cross_match to find
#snp position
#

use strict;
use IPC::Open2; 
use Bio::EnsEMBL::Pipeline::RunnableDB::CrossSNPMap;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::ExternalData::SNPSQL::DBAdaptor; 
use Bio::EnsEMBL::DBLoader;
#use Getopt::Long;

my $host   = 'ecs1d';
my $host_snp = 'ecs1d';
my $port       = '';
my $user   = 'ensro';
#my $dbname_snp = 'mouse_snp';
my $dbname_snp = 'homo_sapiens_snp_120';
my $dbname     = 'homo_sapiens_core_120';

my $module_snp = 'Bio::EnsEMBL::ExternalData::SNPSQL::DBAdaptor';
my $locator_snp = "$module_snp/host=$host_snp;port=$port;dbname=$dbname_snp;user=$user";

print STDERR "Using $locator_snp for snp_db\n";

my $snpdb =  Bio::EnsEMBL::DBLoader->new($locator_snp);
my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(-host =>$host,
					     -user  =>$user,
					     -dbname=>$dbname);

my $score =50;
my $masklevel = 80;
my $minmatch = 20;
my $debug=2;

my ($file_num,$input,$diff,$seq_file,$ssaha_out,$mouse);

if ($dbname_snp =~ /mouse/i) {
  $mouse=1;
}

if ($ARGV[0] and $ARGV[1]) {
  $input = $ARGV[0];
  $diff = $ARGV[1];
  open FILE, "$input" || die "Can not open input file: $! \n";
  open DIFF, ">$diff"  || die "Can not open diff file: $! \n";
  $seq_file = "$input\_seq";
  $ssaha_out = "$input\_ssaha";
}

elsif ($ARGV[0] =~/^\d+|\d+$/ and !$ARGV[1]) {
  $file_num = $ARGV[0];
  open FILE, "/work4/yuan/MOUSE/$file_num" || die "Can not open file: $! \n";
  #open FILE, "/work4/yuan/OUT1/NO_CROSS_MATCH$file_num" || die "Can not open file: $! \n";
  open DIFF, ">/work4/yuan/MOUSE/DIFF\_$file_num" || die "Can not open FINAL: $! \n";
  #open DIFF, ">/work4/yuan/OUT1/SNP_SSAHA_FROM_NO_CROSS_MATCH$file_num" || die "Can not open FINAL: $! \n";
  
  #my $seq_file = "/nfs/acari/yuan/ensembl/test/SNP_SSAHA_FROM_NO_CLONE_seq_$file_num";
  $seq_file =  "/work4/yuan/MOUSE/SEQ\_$file_num";
  $ssaha_out = "/work4/yuan/MOUSE/ssaha\_$file_num";
}

my (%done, %no_clone, %find, %rec ,%snppos);

while (<FILE>) {
  if (/^\d+/ or /CROSS_MATCH/) { 
    s/^\s+//;
    chomp;
    my @all; my $refsnpid;
    @all = split;
    if (/^\d+/) {
      $refsnpid = $all[0];
    }
    elsif (/CROSS_MATCH/) { 
      $refsnpid = $all[4];
    }
    $no_clone{$refsnpid}=1;
  }
}
  
  
my $clone_count = keys %no_clone;
print "Total key is $clone_count\n";

my $snp_seq_all; 
foreach my $refsnpid (sort numeric keys %no_clone) {
  if ($debug>1) {
    print "refsnpid is $refsnpid\n";
  }
  my ($info) = $snpdb->get_snp_info_by_refsnpid($refsnpid,$mouse);
  if ($info) {
    if ($debug>1) {
      print "refsnpid is $refsnpid and info is $info\n";
    }
    print "THIS IS $info\n";
    my $observed = $info->alleles;
    $observed =~ s/\/.*//;
    my $seq5 = $info->upStreamSeq;
    my $seq3 = $info->dnStreamSeq;
    #print STDERR "seq5 is $seq5 and seq3 is $seq3\n";
    my $snp_pos = (length $seq5) +1;
    my $snp_seq = $seq5.$observed.$seq3;
    $snp_seq =~ s/\s+//g;
    #print "the length of query sequence is ", length($snp_seq), "\n";
    if (length $snp_seq <40 ) {
      #my $crossmap = Bio::EnsEMBL::Pipeline::RunnableDB::CrossCloneMap->new(-crossdb=>$crossdb, -score =>80, -debug=>$debug);
      #my $snp_clone_seq = $crossmap->get_seq_from_embl($acc,$version);
      #undef $crossmap;
      #my $snp_seq;
      #if ($snp_clone_seq) {
      #  $snp_seq = substr ($$snp_clone_seq, $start-100, 200);
      #}
      print "SHORT FLANK SEQ FOR $refsnpid\n";
      #$snp_pos = 100;
    }
    $snppos{$refsnpid} = $snp_pos;
    $rec{$refsnpid} = $info;
    $snp_seq_all .=">$refsnpid\n$snp_seq\n";
    #print "this is $snp_seq_all\n";
  }
}

open FASTA, ">$seq_file" || die "Can not open $seq_file: $! \n";
print FASTA $snp_seq_all;
close FASTA;

sub numeric {$a<=>$b}

if ($debug>1) {
  #print STDERR "The query seq is $snp_seq_all\n";
}

if (!(-e $ssaha_out) and -e $seq_file) {
  #### note port 50000 for human and 50006 is for mouse
  #open IN, "cat $seq_file | /usr/local/badger/bin/ssahaClient tcs1a 50006 20 0 0 0 DNA 10000 2 size |" ||  die "OPEN SSAHA failed :$!\n";
  open IN, "cat $seq_file | /usr/local/ensembl/bin/ssahaClient tcs1a 50000 20 0 0 0 DNA 10000 2 size |" ||  die "OPEN SSAHA failed :$!\n";
}
elsif (-e $ssaha_out) {
  open IN, "$ssaha_out" || die "$ssaha_out can not open :$!\n";
}

#open TEST, ">$ARGV[0]\_test";

while (<IN>) {
  chomp;
  my ($qname, $qstart, $qend, $sname, $sstart, $ssend, $dir, $num, $perc) = split (/\s+/, $_); 
  print "this is line $_ and snp_pos = $snppos{$qname}\n";
  #if ($snppos{$qname}>$qstart-20 and $snppos{$qname}<$qend+20) {
    push (@{$find{$qname}}, $sname);
    #if ($snppos{$qname}>$qstart and $snppos{$qname}<$qend) {#####added for testing
    #my $diff_po = $snppos{$qname}-$qstart;
    #my $snp_po = $sstart+$diff_po;
    #my ($clone, $version) = split /\./, $sname;
    #my $strand;
    #if ($dir eq 'F') {
    #  $strand=1;
    #}
    #else {
    #  $strand=-1;
    #}
    #print TEST "TEST $qname\t$sname\t$version\t$snp_po\t\tsnp_map\t$strand\n";
    #}
  #}
}
my $count;
foreach my $refsnpid (keys %find) {
  my (@contigs, @raw_contigs);
  my @contigs = @{$find{$refsnpid}};
  print "this is contigs @contigs\n";
  my $crossmap = Bio::EnsEMBL::Pipeline::RunnableDB::CrossSNPMap->new(-snpdb=>$snpdb, 
								      -db    =>$db,
								      -score =>$score, 
								      -minmatch =>$minmatch,
								      -masklevel =>$masklevel,
								      -debug => $debug,
								     );
  ##get rid off duplicated contigs
  my %uniq_contig;
  foreach my $contig (@contigs) {
    if (!$uniq_contig{$contig} and $contig) {
      $uniq_contig{$contig}=1;
      my $rc;
      eval {
	$rc = $db->get_Contig($contig);
      };
      if ($rc) {
	push (@raw_contigs, $rc);
      }
      else {
	print "$contig for $refsnpid is no in db $@\n";#incase ssaha searched wrong db
	if (!@raw_contigs and $contig =~/(\w+)\.\d+\.\d+\..*/) {
	  my $clone = $1;
	  my $clone_obj;
	  print "Try $clone if no raw_contigs...\n";
	  eval {
	    $clone_obj = $db->get_Clone($clone);
	  };
	  if ($clone_obj) {
	    my @contigs = $clone_obj->get_all_Contigs();
	    push (@raw_contigs, @contigs);
	  }
	  else {
	    print "$contig or $clone for $refsnpid is no in db $@\n";
	  }
	}
      }
    }
  }

  #print "this is raw_contigs @raw_contigs\n";
  
  print STDERR "fetching input for clone $refsnpid\n";
  $crossmap->fetch_input($refsnpid,"0",@raw_contigs);
  print STDERR "Running for clone $refsnpid\n";
  $crossmap->run;
  print STDERR "Writing output for clone $refsnpid\n";
  my ($last_line) = $crossmap->write_output ();
  if ($last_line) {
    print DIFF "$last_line\n";
    $done{$refsnpid}=1;
    $count++;
    print STDERR "$count $refsnpid has done\n";
  }
  else {
    print "no match at current criterion for $refsnpid\n";
  }
}


#open OUT, ">/work4/yuan/SNP_MAP/final_no_clone_from_no_clone_$file_num";
open OUT, ">/work4/yuan/SNP_MAP2/final_no_clone_from_cross_match$file_num";

my ($no_clone_count, $clone_mapped_count);
foreach my $refsnpid (keys %no_clone) {
  if ($done{$refsnpid}) {
    $clone_mapped_count++;
  }
  elsif (!$done{$refsnpid}) {
    print OUT "$refsnpid\n";
    $no_clone_count++;
  }
}

print STDERR "final $clone_mapped_count mapped and $no_clone_count do not mapped\n";









