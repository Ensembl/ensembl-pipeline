#!/usr/local/bin/perl

use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::SeqIO;
use Getopt::Long;
use Bio::EnsEMBL::Pipeline::GeneConf;
use Bio::EnsEMBL::Pipeline::ESTConf;


my $dbhost = $EST_DBHOST; 
my $dbuser = $EST_DBUSER; 
my $dbname = $EST_DBNAME;
my $dbpass = $EST_DBPASS;

my $dnadbhost = $EST_REFDBHOST;
my $dnadbuser = $EST_REFDBUSER;
my $dnadbname = $EST_REFDBNAME;
my $dnadbpass = $EST_REFDBPASS; 

my $genetype = 'exonerate_e2g';


&GetOptions(
	    'dbhost:s'        => \$dbhost,
	    'dbname:s'        => \$dbname,
	    'dnadbhost:s'     => \$dnadbhost,
	    'dnadbname:s'     => \$dnadbname,
	    'genetype:s'      => \$genetype,
	    );

unless ( 1 ){
  print  "script to check the density of genes\n";
  
  print  "Usage: $0 [-dbname -dbhost -dnadbname -dnadbhost -genetype]\n";
  exit(0);
}


my $dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					       '-host'   => $dnadbhost,
					       '-user'   => $dnadbuser,
					       '-dbname' => $dnadbname,
					       '-pass'   => $dnadbpass,
					      );


my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $dbhost,
					    '-user'   => $dbuser,
					    '-dbname' => $dbname,
					    '-pass'   => $dbpass,
					    '-dnadb'  => $dnadb,
					   );


print STDERR "connected to $dbname : $dbhost\n";


my $sgp = $db->get_StaticGoldenPathAdaptor();

my @chrnames = &get_chrnames($db);


foreach my $chr ( @chrnames ){
  print "fetching density of genes $genetype in chromosome $chr...\n";
  my ($chr_start,$chr_end) = &get_start_end($db,$chr);
  
  my $start = 1;
  my $size = 1100000;
  my $end = $start + $size - 1;

  my ($contigs,$pos) = &get_contigs($db,$chr);
  my @contigs = @$contigs;
  my %position = %$pos;
  
  my $length = 0;
  my @exons;
  my $start;
  for(my $i=0;$i<scalar(@contigs);$i++){
    unless ($start){
      $start = $position{$contigs[$i]}{start};
    }
    my $this_length = ($position{$contigs[$i]}{end} - $position{$contigs[$i]}{start} + 1 );
    if ( ($length+$this_length) > $size ){
      my $number;
      if ( @exons ){
	$number = scalar(@exons);
      }
      else{
	$number = "NONE"; # we label it with something easy to 'grep' afterwards
	}
      print "$chr\t$start\t".($start + $length - 1)."\t".$number." exons\n";
      $length = 0;
      $start = $position{$contigs[$i]}{start};
      @exons = ();
    }
    push( @exons, &get_exons($db,$contigs[$i]) );
    $length += $this_length;
  }
}
  

sub get_start_end{
  my $db   = shift;
  my $chr  = shift;
  
  if (!$db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor')) {
    die "get_chrlengths should be passed a Bio::EnsEMBL::DBSQL::DBAdaptor\n";
  }
  
  
  my $q = qq( SELECT min(chr_start),max(chr_end)
	      FROM static_golden_path 
	      WHERE chr_name='$chr'
	    );
  
  my $sth = $db->prepare($q) || $db->throw("can't prepare: $q");
  my $res = $sth->execute || $db->throw("can't execute: $q");
  
  my ($chr_start,$chr_end) = $sth->fetchrow_array;
  return ($chr_start,$chr_end);
}


sub get_chrnames{
  my $db   = shift;
  
  if (!$db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor')) {
    die "get_chrlengths should be passed a Bio::EnsEMBL::DBSQL::DBAdaptor\n";
  }
  
  my @chrnames;
  
  my $q = qq( SELECT distinct(chr_name) 
	      FROM static_golden_path 
	    );
  
  my $sth = $db->prepare($q) || $db->throw("can't prepare: $q");
  my $res = $sth->execute || $db->throw("can't execute: $q");
  
  while( my ($chr) = $sth->fetchrow_array) {
    push (@chrnames, $chr);
  }
  return @chrnames;
}

sub get_contigs{
  my $db   = shift;
  my $chr = shift;
  
  if (!$db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor')) {
    die "get_chrlengths should be passed a Bio::EnsEMBL::DBSQL::DBAdaptor\n";
  }
  
  my %position;
  my %id;
  my @contigs;
  
  my $q = qq( SELECT s.chr_start,s.chr_end,s.raw_id 
	      FROM static_golden_path s
	      WHERE s.chr_name = '$chr'
	      ORDER BY s.chr_start
	    );
  
  my $sth = $db->prepare($q) || $db->throw("can't prepare: $q");
  my $res = $sth->execute || $db->throw("can't execute: $q");
  
  while( my ($chr_start, $chr_end,$raw_id) = $sth->fetchrow_array) {
    push (@contigs, $raw_id);
    $position{$raw_id}{start} = $chr_start;
    $position{$raw_id}{end}   = $chr_end;
  }
  return (\@contigs,\%position);
}

sub get_exons{
  my $db   = shift;
  my $contig = shift;
  
  if (!$db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor')) {
    die "get_chrlengths should be passed a Bio::EnsEMBL::DBSQL::DBAdaptor\n";
  }
  
  my @exon_ids;
  my @contigs;
  
  my $q = qq( SELECT exon_id 
	      FROM exon
	      WHERE contig_id=$contig
	    );
  
  my $sth = $db->prepare($q) || $db->throw("can't prepare: $q");
  my $res = $sth->execute || $db->throw("can't execute: $q");
  
  my $count = 0;
  while( my ($id) = $sth->fetchrow_array) {
    push ( @exon_ids,$id);
  }
  return @exon_ids;
}

