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
use Data::Dumper;

#use Getopt::Long;

my $host   = 'ecs1d';
my $host_snp = 'ecs1e';
my $port       = '3310';
my $user   = 'ensro';
my $dbname_snp = 'hum_snp_112';
#my $dbname     = 'otter_merged_end_jul';
#my $dbname     = 'human_chr10_comp';
my $dbname     = 'homo_sapiens_core_16_33_new';
#my $dbname     = 'anopheles_research_core_10_2';

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(-host =>$host,
					     -user  =>$user,
		                             #-port => $port,
					     -dbname=>$dbname);

$db->assembly_type('VEGA');

my $snp_db;

eval {
  $snp_db = Bio::EnsEMBL::ExternalData::SNPSQL::DBAdaptor->new(
								-host =>$host_snp,
								-user  =>$user,
								-dbname=>$dbname_snp)
};

my $snpdb;

eval {							        
  $snpdb = $snp_db->get_SNPAdaptor
};

print Dumper "db is $db\n";
print Dumper "snpdb is $snpdb\n";
 
my $score =20;
my $masklevel = 80;
my $minmatch = 20;
my $debug=1;
my $chr = 10;
my $flank_base = 500;

my ($file,$est,$first,$file_num,$first,$input,$search_file,$start_internal_id,$diff,$seq_file,$ssaha_out,$feature_out,$mouse,%inseq);

if ($dbname_snp =~ /mouse|mus_musculus/i) {
  $mouse=1;
}

####if use $ARGV[0] eq "file", $input_seq has to be in same dir as $input###
if ($ARGV[0] eq "file") {
  $file=1;
  $input = $ARGV[1];
  $search_file = $ARGV[2];
  $est   = $ARGV[3];
}
#####the first round check to see how many snps are already mapped in clones which is in GP
elsif ($ARGV[0] eq "first") {
  $input = $ARGV[1];
  $first=1;
  $input =~ /\_(\d+)/;
  $start_internal_id = $1;
}
######read input file to mapped them
elsif ($ARGV[0] =~/^\d+|\d+$|all$/ and !$ARGV[1]) {
  $input = $ARGV[0];
}

open FILE, "$input" || die "Can not open input file: $! \n";
if ($file) {
  $seq_file = $input;
}
else {
  $seq_file = "$input\_seq";
}

my @input_names = split /\//, $input;
$input = $input_names[-1];
open DIFF, ">DIFF\_$input" || die "Can not open diff file: $! \n";###added /tmp dir
$ssaha_out = "$input\_ssaha";
$feature_out = "$input\_features";


my $ssaha_command;

if ($input =~ /no\_match/ or $ARGV[4]==1) {
  $ssaha_command = "/usr/local/ensembl/bin/ssaha $seq_file $search_file -qf fasta -sf fasta -da 0 -mp 1 -sm 5 -ms 100 -mg 5 -mi 5 -pf > $ssaha_out";
}
else { 
  $ssaha_command = "/usr/local/ensembl/bin/ssaha_zn_old $seq_file $search_file -hit 10 -len 12 -cut 100 > $ssaha_out";
}



my (%done, %input_id, %find, %rec ,%snppos);

while (<FILE>) {
  my $refsnpid ;
  if (/^\d+/ or /^>/) {
    s/^>//;
    chomp;
    $refsnpid = $_;
    #print "refsnpid is $refsnpid\n";
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
								     );
  my $end_internal_id = $start_internal_id+500-1;
  my ($same_version) = $crossmap->check_data_for_acc_version($start_internal_id,$end_internal_id);
  
  my %same_version = %$same_version;
  
  my %no_clones;
  open NO_CLONE, ">NO_CLONE\_$input";
  open SAME_CLONE, ">SAME_CLONE\_$input";
  
  my $same_version_count = keys %same_version;
  print "same_version_count is $same_version_count\n";

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
    $refsnpid =~ s/\s+//g;
    ##get rid off duplicated contigs
    my ($info) = $snpdb->fetch_by_refsnpid($refsnpid,$mouse);
    
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
      my $snp_seq = lc($seq5).uc($observed).lc($seq3);
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
      $snp_seq_all .=">$refsnpid-$snp_pos\n$snp_seq\n";
      #print "this is $snp_seq_all\n";
    }
    else {
      print STDERR "snp info for $refsnpid is not found\n";
    }
  }
  
  open FASTA, ">$seq_file" || die "Can not open $seq_file: $! \n";
  print FASTA $snp_seq_all;
  close FASTA;
  #exit; ##if only want flanking sequences generated from refsnpids
}
elsif ($file) { ###it is useful for snp mapping###
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

#print STDERR "The query seq is $snp_seq_all\n";

####if total_ssaha_out is exist, always use it to save run ssaha many times####
if (!(-e "total_ssaha_out")) {
  if ((!(-e $ssaha_out) or -z $ssaha_out) and -e $seq_file) {
    my $date1 = `date`;
    print "This is before ssaha $date1\n";
    print "seq file is $seq_file\n";
    
    system ($ssaha_command);  ### ||  die "OPEN SSAHA failed :$!\n";
    
    my $date2 = `date`;
    print "This is after ssaha $date2\n";
  }
  if (-e $ssaha_out) {
    open IN, "$ssaha_out" || die "$ssaha_out can not open :$!\n";
  }
}
elsif (-e "total_ssaha_out") {
  open IN, "total_ssaha_out" || die "$ssaha_out can not open :$!\n";
}

#open TEST, ">$input\_diff_ssaha";

my (%rec_find, $qname, $zn_ssaha);

while (<IN>) {
  chomp;
  if (/^FF|RF/) {#####matching tony cox's ssaha
    my ($dir, $qname, $qstart, $qend, $sname, $sstart, $ssend, $num, $perc) = split (/\s+/, $_);
    #print "this is line $_ and $qname, $qstart, $qend, $sname, $sstart, $ssend, $num, $perc\n";
    if ($sname =~ /^\d+\-\d+$/) {
       if ($sname eq "1-4641698") {
          $chr = 25;
       }
       elsif ($sname eq "1-4754829") {
          $chr = 26;
       }
       elsif ($sname eq "1-170730988") {
          $chr = 6;
       }
       $sname = "$chr-$sstart-$ssend";
    }
    push (@{$find{$qname}}, [$sname,$sstart,$ssend]);
  }
  #if (/^\d+\s+\S+\s+\S+\s+\d+\s+\d+\s+\d+\s+\d+\s+[FC]{1}\s+\d+\s+\S+$/) {####matching ssaha2
    #my ($score,$qname,$sname,$qstart,$qend,$sstart, $ssend,$dir,$num,$perc) = split (/\s+/, $_); 
  elsif (/Input Query Sequence Name:\s+(\S+)/) {####match zn's ssaha
    $qname = $1;
    $zn_ssaha =1;
  }
  elsif (/^\S+\s+\d+\s+\d+\s+\d+\s+[+C]{1}\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+$/) {####match zn's ssaha
    my ($sname,$index,$length,$match,$rc,$qstart,$qend,$sstart, $ssend,$slength) = split (/\s+/, $_); 
    #print "this is line $_ and snp_pos = $snppos{$qname}\n";
    if ($sname =~ /^\d+\-\d+$|^\d+/) {
       if ($sname eq "1-4641698") {
          $chr = 25;
       }
       elsif ($sname eq "1-4754829") {
          $chr = 26;
       }
       elsif ($sname eq "1-170730988") {
          $chr = 6;
       }

      if ($rc eq "+") {
	$sname = "$chr-$sstart-$ssend";
      }
      elsif ($rc eq "C") {
	my $start = $sstart;
	$sstart = $slength-$ssend;
	$ssend = $slength-$start;
	$sname = "$chr-$sstart-$ssend";
      }
    }
    my $h={};
    $h->{'qname'} = $qname;
    $h->{'sname'} = $sname;
    $h->{'match'} = $match;
    $h->{'sstart'} = $sstart;
    $h->{'send'} = $send;
    #print "$qname, $sname, $match\n" if ($qname == 193385);
    push (@{$rec_find{$qname}}, $h);
  }
}

if ($zn_ssaha) {
  foreach my $qname (keys %rec_find) {
    my @h =  @{$rec_find{$qname}};
    @h = sort {$b->{'match'}<=>$a->{'match'}} @h;
    foreach my $h (@h[0..2]) {
      #print "qname is $qname and match is ",$h->{'match'},"\n" if ($qname == 193385);
      push @{$find{$qname}}, [$h->{'sname'},$h->{'sstart'},$h->{'send'}];
    }
  }
}

my ($count,%feature);

my $date = `date`;
print "This is after parsing ssaha file $date\n";

#exit;

FIND : foreach my $refsnpid (keys %input_id) { ##only table id from input file not ssaha output
  my (@contigs, @raw_contigs); print "this is refsnpid $refsnpid\n";
  @contigs = @{$find{$refsnpid}} if ($find{$refsnpid}); ##only ssaha output for the id is exist
  print "this is contigs @contigs\n";

  my $crossmap = Bio::EnsEMBL::Pipeline::RunnableDB::CrossSNPMap->new(-snpdb=>$snpdb, 
								      -db    =>$db,
								      -score =>$score, 
								      -minmatch =>$minmatch,
								      -masklevel =>$masklevel,
								      -est       =>$est,
								      -debug => $debug,
								     );
  ##get rid off duplicated contigs
  my (%uniq_contig,$slice,$coord_system_name,$seq_region_name,$start,$end,$strand,$version);

  CONTIGS : foreach my $contig_ref (@contigs) {
    $contig_name = $contig_ref->[0];
    $start = $contig_ref->[1];
    $end = $contig_ref->[2];

    $contig_name =~ s/^Contig\://;
    if (!$uniq_contig{$contig_name} and $contig_name) {
      $uniq_contig{$contig_name}=1;
      
      #print "contig_name is $contig_name\n";
      
      if ($contig_name =~/^\w+\.\d+\.\d+\.\d+$/) {
	$coord_system_name = "contig";
      }
      elsif ($contig_name =~/^\w+\.\d+$/) {
	$coord_system_name = "clone";
      }
      elsif ($contig_name =~/\S+\-\d+\-\d+/) {
	($contig_name, $start, $end) = split /\-/, $contig_name;
	$coord_system_name = "chromosome";
	$start = $start-$flank_base;
	$start = 0 if ($start <0);
	$end = $end + $flank_base;
      }

      eval {
	$slice = $db->get_SliceAdaptor->fetch_by_region($coord_system_name, $contig_name, $start, $end,$strand,$version);
      };
      if ($slice) {
	push (@raw_contigs, $slice);
	if (@raw_contigs>10) {###only take first 10 contigs
	  #print STDEER "More than 3 hits found for this query---discard the later contigs\n";
	  last CONTIGS;
	}
      }
    }
  }
  
  print "this is raw_contigs @raw_contigs\n";
  
  print STDERR "fetching input for clone $refsnpid\n";
  
  if ($file) {
    $crossmap->seq_map($refsnpid,\$inseq{$refsnpid},\@raw_contigs);
  }
  else {
    $crossmap->fetch_input($refsnpid, \@raw_contigs,$mouse);
  }
  
  print STDERR "Running for snp $refsnpid\n";
  $crossmap->run;
  print STDERR "Writing output for clone $refsnpid\n";
  
  #my (@fp)=$crossmap->output_features;
  #@{$feature{$refsnpid}} = @fp; 

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


#open FEA, ">$feature_out";
open OUT, ">$input\_no_match";

my ($no_clone_count, $clone_mapped_count);
foreach my $refsnpid (keys %input_id) {
  if ($done{$refsnpid}) {
    $clone_mapped_count++;
  }
  elsif ($refsnpid and !$done{$refsnpid}) {
    print OUT "$refsnpid\n";
    $no_clone_count++;
  }
  #foreach my $fp (@{$feature{$refsnpid}}) {
  #  print FEA $fp->seqname,"\t",$fp->start,"\t",$fp->end,"\t",$fp->strand,"\t", 
  #  $fp->hseqname,"\t",$fp->hstart,"\t",$fp->hend,"\t",$fp->hstrand,"\t",$fp->score,"\n";
  #}
}

print STDERR "final $clone_mapped_count mapped and $no_clone_count do not mapped\n";








