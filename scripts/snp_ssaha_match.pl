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

my $host   = 'ecs1e';
my $host_snp = 'ecs1a';
my $port       = '';
my $user   = 'ensro';
#my $dbname_snp = 'mouse_snp';
my $dbname_snp = 'mouse_snp_105_NT29';
my $dbname     = 'mouse_whitehead_0401_denormalised';

my $module_snp = 'Bio::EnsEMBL::ExternalData::SNPSQL::DBAdaptor';
my $locator_snp = "$module_snp/host=$host_snp;port=$port;dbname=$dbname_snp;user=$user";

print STDERR "Using $locator_snp for snp_db\n";

my $snpdb =  Bio::EnsEMBL::DBLoader->new($locator_snp);

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(-host =>$host,
					     -user  =>$user,
					     -dbname=>$dbname);

my $score =20;
my $masklevel = 80;
my $minmatch = 20;
my $debug=2;

my ($file_num,$first,$input,$start_intnum,$diff,$seq_file,$ssaha_out,$mouse,%inseq);

if ($dbname_snp =~ /mouse|mus_musculus/i) {
  $mouse=1;
}

####if use $ARGV[0] eq "file", $input_seq has to be in same dir as $input###
if ($ARGV[0] eq "file") {
  $input = $ARGV[1];
  open FILE, "$input" || die "Can not open input file: $! \n";
  open DIFF, ">$input\_diff\_$score" || die "Can not open diff file: $! \n";
  $seq_file = "$input\_seq";
  $ssaha_out = "$input\_ssaha";
}
#####the first round check to see how many snps are already mapped in clones which is in GP
elsif ($ARGV[0] eq "first") {
  $input = $ARGV[1];
  open FILE, "$input" || die "Can not open input file: $! \n";
  $first=1;
  $input =~ /\_(\d+)/;
  $start_intnum = $1;
}
######read input file to mapped them
elsif ($ARGV[0] =~/^\d+|\d+$|all$/ and !$ARGV[1]) {
  $input = $ARGV[0];
  open FILE, "$input" || die "Can not open file: $! \n";
  open DIFF, ">DIFF\_$input" || die "Can not open FINAL: $! \n";
  $seq_file =  "$input\_seq";
  $ssaha_out = "$input\_ssaha";
}

my (%done, %input_id, %find, %rec ,%snppos);

while (<FILE>) {
  if (/^\d+/) { 
    s/^\s+//;
    chomp;
    my @all; my $refsnpid;
    @all = split;
    if (/^\d+/) {
      $refsnpid = $all[0];
    }
    $input_id{$refsnpid}=1;
  }
}
  
  
my $clone_count = keys %input_id;
print "Total key is $clone_count\n";

if ($first) {
  my $crossmap = Bio::EnsEMBL::Pipeline::RunnableDB::CrossSNPMap->new(-snpdb=>$snpdb, 
								      -db    =>$db,
								      -score =>$score, 
								      -minmatch =>$minmatch,
								      -masklevel =>$masklevel,
								      -debug => $debug,
								      -first => $first,
								     );
  $crossmap->fetch_input(
			 -inseq        => "",
			 -start_refnum =>$start_intnum,
			 -end_refnum   =>$start_intnum+500-1,
			 -contigs      =>"",
			 -mouse        =>$mouse);
  
  my $all_version = $crossmap->all_version;
  my $same_version = $crossmap->same_version;
  my %all_version = %$all_version;
  my %same_version = %$same_version;
  
  my %no_clones;
  open NO_CLONE, ">NO_CLONE_FIRST\_$input";
  open SAME_CLONE, ">SAME_CLONE_FIRST\_$input";
  
  my $all_version_count = keys %all_version;
  print "all_version_count is $all_version_count\n";

  foreach my $refsnpid (keys %input_id) {
    if (!$same_version{$refsnpid}) {
      print NO_CLONE "$refsnpid\n";
    }
    else {
      print SAME_CLONE "$refsnpid\n";
    }
  }
  exit;
}

my $snp_seq_all;

####generate snp query sequence file############
if ($ARGV[0] ne "file" and !(-e $seq_file)) {
  foreach my $refsnpid (sort numeric keys %input_id) {
    if ($debug>1) {
      print "refsnpid is $refsnpid\n";
    }
    
    ##get rid off duplicated contigs
    my ($info) = $snpdb->get_snp_info_by_refsnpid($refsnpid,$mouse);
    
    if ($info) {
      if ($debug>1) {
	print "refsnpid is $refsnpid and info is $info\n";
      }
      #print "THIS IS $info\n";
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
}
elsif ($ARGV[0] eq "file") { ###only useful for hgbase snp mapping###
  open SEQ, "$seq_file" || die "Can not open $seq_file: $! \n";
  my $name;
  while (<SEQ>) {
    if (/^\>/) {
      s/^\>|\n//g;
      $name = $_;
    }
    else {
      s/\n//;
      $inseq{$name}.=$_;
    }
  }
}

sub numeric {$a<=>$b}

if ($debug>1) {
  #print STDERR "The query seq is $snp_seq_all\n";
}

if (!(-e $ssaha_out) and -e $seq_file) {
  #### note port 50000 for human and 50006 is for mouse 50005 for denormalised coordinates
  #open IN, "cat $seq_file | /usr/local/badger/bin/ssahaClient tcs1h 50005 20 0 0 0 DNA 10000 2 size |" ||  die "OPEN SSAHA failed :$!\n";

  my $date1 = `date`;
  print "This is before ssaha $date1\n";

  system ("/usr/local/ensembl/bin/ssaha $seq_file /acari/work4/yuan/MOUSE/mouse_0402/ssahaTables/ssahaTable -qf fasta -sf hash -sm 3 -ms 100 -mg 5 -mi 5 -pf >$ssaha_out");# ||  die "OPEN SSAHA failed :$!\n";
  #`/usr/local/ensembl/bin/ssaha $seq_file /acari/scratch6/ensembl/mouse-105/ssahaTables/ssahaTable -qf fasta -sf hash -sm 3 -ms 100 -mg 5 -mi 5 -pf >$ssaha_out` ||  die "OPEN SSAHA failed :$!\n";
}

my $date2 = `date`;
print "This is after ssaha $date2\n";

if (-e $ssaha_out) {
  open IN, "$ssaha_out" || die "$ssaha_out can not open :$!\n";
}

open TEST, ">$input\_diff_ssaha";

while (<IN>) {
  chomp;
  my ($dir, $qname, $qstart, $qend, $sname, $sstart, $ssend, $num, $perc) = split (/\s+/, $_); 
  #print "this is line $_ and snp_pos = $snppos{$qname}\n";
  push (@{$find{$qname}}, $sname);
  #print "start is $qstart AND end is $qend and snp_pos is $snppos{$qname}\n";
  #if ($snppos{$qname}>$qstart and $snppos{$qname}<$qend) {#####added for testing
  #  #print "I AM HERE $qstart AND $qend\n";
  #  my $diff_po = $snppos{$qname}-$qstart;
  #  my $snp_po = $sstart+$diff_po;
  #  my ($clone, $version) = split /\./, $sname;
  #  my $strand;
  #  if ($dir eq 'FF') {
  #    $strand=1;
  #  }
  #  elsif($dir eq 'RF')  {
  #    $strand=-1;
  #  }
  #  print TEST "TEST $qname\t$sname\t$version\t$snp_po\t\tsnp_map\t$strand\n";
  #}
}
my $count;
FIND : foreach my $refsnpid (keys %find) {
  my (@contigs, @raw_contigs);
  @contigs = @{$find{$refsnpid}};
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
	if (@raw_contigs>2) {###get rid off hits more than 2 places
	  print STDEER "More than 2 hits found for this query---discarded\n";
	  next FIND;
	}
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
  
  $crossmap->fetch_input(
			 -inseq        => "$refsnpid\_$inseq{$refsnpid}",
			 -start_refnum =>$refsnpid,
			 -end_refnum   =>"",
			 -contigs      =>\@raw_contigs,
			 -mouse        =>$mouse);

  print STDERR "Running for snp $refsnpid\n";
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
open OUT, ">$input\_no_match";

my ($no_clone_count, $clone_mapped_count);
foreach my $refsnpid (keys %input_id) {
  if ($done{$refsnpid}) {
    $clone_mapped_count++;
  }
  elsif (!$done{$refsnpid}) {
    print OUT "$refsnpid\n";
    $no_clone_count++;
  }
}

print STDERR "final $clone_mapped_count mapped and $no_clone_count do not mapped\n";









