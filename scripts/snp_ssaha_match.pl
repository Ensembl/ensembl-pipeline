#! /usr/local/ensembl/bin/perl
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
use Bio::SeqIO;

#use Getopt::Long;

my $host   = 'ecs2f';
my $host_snp = 'ecs1e';
#my $port       = '3310';
my $user   = 'ensro';
my $dbname_snp = 'hum_snp_116';
#my $dbname     = 'otter_merged_end_jul';
#my $dbname     = 'fish_zv3_core';
#my $dbname     = 'human_ncbi34_raw';
#my $dbname     = 'homo_sapiens_core_18_34';
#my $dbname     = 'anopheles_research_core_10_2';
my $dbname     = 'rattus_norvegicus_core_20_3a';

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(-host =>$host,
					     -user  =>$user,
		                             #-port => $port,
					     -dbname=>$dbname);

#$db->assembly_type('VEGA');

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
my $chr = 7;
my $flank_base = 500;

my ($file,$id_file,$first,$file_num,$first,$input,$search_file,$non_ensembl,$start_internal_id,$diff,$seq_file,$ssaha_out,$feature_out,$mouse,%inseq);

###use following if a species only have RefSNP table and without mapping info (i.e. Hit, ContigHit tables)
###the initial mouse dbSNP was like that, so the following code
#if ($dbname_snp =~ /mouse|mus_musculus/i) {
#  $mouse=1;
#}

####if $ARGV[0] eq "seq_file", the file contains a list of flanking sequences to be mapped###
if ($ARGV[0] eq "seq_file") {
  $file=1;
  $input = $ARGV[1];
  ($chr) = $input =~ /\d+\-ch(\d+)\_query_seq$/;
  $search_file = $ARGV[2];
  $non_ensembl=1 if ($ARGV[3] eq "non_ensembl");
  $seq_file = $input;
}
#####the first round check to see how many snps are already mapped in clones which is in GP
elsif ($ARGV[0] eq "first") {
  $input = $ARGV[1];
  $first=1;
  $input =~ /\_(\d+)/;
  $start_internal_id = $1;
}
######read input file, the input file contains a list of refsnpids to be mapped
elsif ($ARGV[0] eq "id_file") {
  $id_file = 1;
  $input = $ARGV[1];
  $search_file = $ARGV[2];
  $seq_file = "$input\_seq";
}

my (%done, %input_id, %find, %rec ,%snppos);

open FILE, "$input" or die "Can not open input file: $! \n";

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

my @input_names = split /\//, $input;
$input = $input_names[-1];
open DIFF, ">DIFF\_$input" || die "Can not open diff file: $! \n";###added /tmp dir
$ssaha_out = "$input\_ssaha";
$feature_out = "$input\_features";


my $ssaha_command;

###this ssaha is slower and catch more hits with lower standard, I will limit score >0.5 for this
if ($input =~ /no\_match/ or $ARGV[3]==1) {
  $ssaha_command = "/usr/local/ensembl/bin/ssaha $seq_file $search_file -qf fasta -sf fasta -da 0 -mp 1 -sm 5 -ms 100 -mg 5 -mi 5 -pf > $ssaha_out";
}
###ssaha_zn is quicker and more strict, but will miss some hits. I use this first, then use ssaha with no_matched ones
else { 
  $ssaha_command = "/usr/local/ensembl/bin/ssaha_zn $seq_file $search_file -hit 10 -len 12 -cut 100 > $ssaha_out";
  #$ssaha_command = "/nfs/acari/yuan/ensembl/perl/ssaha_zn_old $seq_file $search_file -hit 10 -len 12 -cut 100 > $ssaha_out";
}


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

my ($snp_seq_all,@searchseqs);

####generate snp query sequence file############
if ($ARGV[0] eq "id_file" and !(-e $seq_file)) {
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

      my $snp_seq = lc($seq5).uc($observed).lc($seq3);
      $snp_seq =~ s/\s+//g;

      $inseq{$refsnpid} = $snp_seq;

      $snp_seq_all .=">$refsnpid\n$snp_seq\n";
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
elsif ($file or ($id_file and -e $seq_file)) { ###it is useful for snp mapping###
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

if ($non_ensembl) {

  my $in = Bio::SeqIO->new(
			   -FILE   => $search_file,
			   -FORMAT => 'FASTA',
			  );
  
  while (my $seq = $in->next_seq) {       
    print STDERR $seq->display_id,"\n"; 
    push @searchseqs, $seq;
  }
}


sub numeric {$a<=>$b}

#print STDERR "The query seq is $snp_seq_all\n";

my $ssaha_gz = $ssaha_out."\.gz";

print "my ssaha_gz is $ssaha_gz\n";

####if total_ssaha_out is exist, always use it to save run ssaha many times####
if (!(-e "total_ssaha_out")) {
  if (((!(-e $ssaha_out) and !(-e $ssaha_gz)) or (-z $ssaha_out and -z $ssaha_gz)) and -e $seq_file) {
    my $date1 = `date`;
    print "This is before ssaha $date1\n";
    print "seq file is $seq_file\n";
    
    system ($ssaha_command);  ### ||  die "OPEN SSAHA failed :$!\n";

    system ("gzip $ssaha_out") if (!(-z $ssaha_out));
    
    my $date2 = `date`;
    print "This is after ssaha $date2\n";
  }
  if (-e $ssaha_out) {
    open IN, "$ssaha_out" || die "$ssaha_out can not open :$!\n";
  }
  elsif (-e $ssaha_gz) {
    open IN, "gzip -dc $ssaha_gz |" || die "$ssaha_out can not open :$!\n";  }
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
    if ($sname =~ /^(\S+)\-\d+\-\d+$/ | $sname =~ /^\d+\-\d+$/) {
      $chr = $1 if ($1);
      $chr = give_chr_name($db,$sname) if (!$chr);

      ###for chr6 and haplotypes, vega use different chr No. I need map snps on them at the same time.
      #$chr = give_chr_name($db,$sname);
       #if ($sname eq "1-4641698") {
       #   $chr = 25;
       #}
       #elsif ($sname eq "1-4754829") {
       #   $chr = 26;
       #}
       #elsif ($sname eq "1-170730988") {
       #   $chr = 6;
       #}
      #print "chr_name is $chr\n" if ($chr !~ /^\d+/);
      $sname = "$chr-$sstart-$ssend"; ### if ($sstart >32600000 and $ssend < 32620000);###for HAP
      #print "$qname, $sname\n" if ($sname !~ /^\d+/);
    }
    push (@{$find{$qname}}, $sname) ;###if ($sname =~ /\S+\-\d+\-\d+/);
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
    if ($sname =~ /^(\S+)\-\d+\-\d+$/ | $sname =~ /^\d+\-\d+$/) {
      $chr = $1 if ($1);
      $chr = give_chr_name($db,$sname) if (!$chr);
      
       #if ($sname eq "1-4641698") {
       #   $chr = 25;
       #}
       #elsif ($sname eq "1-4754829") {
       #   $chr = 26;
       #}
       #elsif ($sname eq "1-170730988") {
       #   $chr = 6;
       #}
      #print "chr_name is $chr\n" if ($chr !~ /^\d+/);
      if ($rc eq "+") {
	$sname = "$chr-$sstart-$ssend"; ### if ($sstart > 32000000 and $ssend < 32200000); ###for HAP
      }
      elsif ($rc eq "C") {
	my $start = $sstart;
	$sstart = $slength-$ssend;
	$ssend = $slength-$start;
	$sname = "$chr-$sstart-$ssend"; ### if ($sstart > 32000000 and $ssend < 32200000);###for HAP
      }
    }
    my $h={};
    $h->{'qname'} = $qname;
    $h->{'sname'} = $sname;
    $h->{'match'} = $match;
    #print "$qname, $sname, $match\n" if ($sname !~ /^\d+/);
    push (@{$rec_find{$qname}}, $h); ### if ($sname =~ /\S+\-\d+\-\d+/);
  }
}

if ($zn_ssaha) {
  foreach my $qname (keys %rec_find) {
    my @h =  @{$rec_find{$qname}};
    @h = sort {$b->{'match'}<=>$a->{'match'}} @h;
    foreach my $h (@h[0..2]) {
      #print "qname is $qname and match is ",$h->{'match'},"\n" if ($qname == 193385);
      push @{$find{$qname}}, $h->{'sname'};
    }
  }
}

my ($count,%feature);

my $date = `date`;
print "This is after parsing ssaha file $date\n";

#exit;

FIND : foreach my $refsnpid (keys %input_id) { ##only table id from input file not ssaha output
  my (@contigs, @raw_contigs); 
  print "this is refsnpid $refsnpid\n";
  @contigs = @{$find{$refsnpid}} if ($find{$refsnpid}); ##only ssaha output for the id is exist
  print "this is contigs @contigs\n";

  my $crossmap = Bio::EnsEMBL::Pipeline::RunnableDB::CrossSNPMap->new(-snpdb=>$snpdb, 
								      -db    =>$db,
								      -score =>$score, 
								      -minmatch =>$minmatch,
								      -masklevel =>$masklevel,
								      #-est       =>$est,
								      -debug => $debug,
								     );
  ##get rid off duplicated contigs
  my %uniq_contig;
  CONTIGS : foreach my $contig (@contigs) {
    if ($non_ensembl) {
      foreach my $seq (@searchseqs) {
	if (!$uniq_contig{$contig} and $contig =~ /^(\S+)\-(\d+)\-(\d+)$/) {
	  $uniq_contig{$contig}=1;
	  my ($name,$start,$end) = ($1,$2,$3);
	  my $start_slice = $start-$flank_base;
	  $start_slice = 0 if ($start_slice <=0);
	  my $end_slice = $end+$flank_base;
	  my $display_id = $name."."."$start_slice"."-"."$end_slice";
	  my $sub_seq = $seq->subseq ($start_slice,$end_slice);
	  my $sub_seq_obj = Bio::PrimarySeq->new( -display_id => "$display_id", -seq => $sub_seq);
	  push (@raw_contigs, $sub_seq_obj);
	}
	elsif (!$uniq_contig{$contig} and $seq->display_id eq $contig) {
	  $uniq_contig{$contig}=1;
	  push (@raw_contigs, $seq);
	}
      }
    }
    else {
      $contig =~ s/^Contig\://;
      if (!$uniq_contig{$contig} and $contig) {
	$uniq_contig{$contig}=1;
	my $rc;
	#print "contig is $contig\n";
	eval {
	  $rc = $db->get_RawContigAdaptor->fetch_by_name($contig);
	  #$rc = $db->get_Contig($contig); ####using old_schema
	};
	if ($rc) {
	  push (@raw_contigs, $rc);
	  if (@raw_contigs>10) {###only take first 10 contigs
	    #print STDEER "More than 3 hits found for this query---discard the later contigs\n";
	    last CONTIGS;
	  }
	}
	else {
	  print "$contig for $refsnpid is not in db $@\n";#incase ssaha searched wrong db
	  if ($contig =~/(\w+)\.\d+\.\d+\..*/) {
	    my $clone = $1;
	    my $clone_obj;
	    print "Try $clone if no raw_contigs...\n";
	    eval {
	      $clone_obj = $db->get_CloneAdaptor->fetch_by_name($clone);
	    };
	    if ($clone_obj) {
	      my @contigs = $clone_obj->get_all_Contigs();
	      push (@raw_contigs, @contigs);
	    }
	    else {
	      print "$contig or $clone for $refsnpid is no in db $@\n";
	    }
	  }####this is the case sequence mapped to chr
	  elsif ($contig =~/^\S+\-\d+\-\d+/) {
	    my $slice;
	    my ($chr, $start, $end) = split /\-/, $contig;
	    my $start_slice = $start-$flank_base;
	    $start_slice = 1 if ($start_slice <0);
	    $slice = $db->get_SliceAdaptor->fetch_by_chr_start_end($chr, $start_slice, $end+$flank_base);
	    #$slice->name($contig); ####slice will have a name defined for me like : 13.31482005-31483556
	    push (@raw_contigs, $slice);
	  }
	}
      }
    }
  }
  print "this is raw_contigs @raw_contigs\n";
  
  print STDERR "fetching input for clone $refsnpid\n";
  
  $crossmap->seq_map($refsnpid,\$inseq{$refsnpid},\@raw_contigs);
    
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


sub give_chr_name {
  
  my ($db,$name) = @_;
  #print "name is $name\n";
  
  my ($start,$length) = split /\-/,$name;
  my $query = "select name from chromosome where length = $length";

  my $sth = $db->prepare($query);
  $sth->execute();
 
  while (my ($chr_name) = $sth->fetchrow) {
    return ($chr_name);
  }
}





