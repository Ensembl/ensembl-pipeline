#!/usr/local/ensembl/bin/perl -w

use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;

$| = 1;

my ($dbname, $dbhost, $dbuser, $dbpass, $contig, $slice, $size);

GetOptions(
	   'dbname=s' => \$dbname,
	   'dbpass=s' => \$dbpass,
	   'dbhost=s' => \$dbhost,
	   'dbuser=s' => \$dbuser,
	   'contig=s' => \$contig,
	   'slice=s' => \$slice,
	   'size=s' => \$size,
) 
or die ("couldn't get options :$! ");

foreach my $arg($dbname, $dbhost, $dbuser, $dbpass){
  if ($arg eq '' ){
   print STDERR "you must defined -dbname $dbname -dbhost $dbhost -dbuser $dbuser and -dbpass $dbpass -contig -chr -slice -size on the commandline \n"; 
    exit(1);
  }
}

if(!$contig && !$slice){
  print STDERR "One of -contig -chr -slice must be defined in order for this script to do anything\n";
  # this must refer to the dunmmy analysis logic name for analyses which depend on input_ids which can be contig names, chr names, or slice names
  exit(1);
}

if($slice && !$size){
  print STDERR "if you want to generate slice input ids you miust specifiy what size pieces you want with -size\n";
  exit(1);
}

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $dbhost,
					      -user   => $dbuser,
					      -pass   => $dbpass,
					      -dbname => $dbname,
					     );

my %input_ids;


if($contig){
  my $sql = 'select name from contig';
  my $sth = $db->prepare($sql);
  $sth->execute;
  my @input_ids;
  while(my ($name) = $sth->fetchrow){
    push(@input_ids, $name);
  }
  if(!$input_ids{$contig}){
    $input_ids{$contig} = [];
    push(@{$input_ids{$contig}}, @input_ids);
  }else{
    push(@{$input_ids{$contig}}, @input_ids);
  }
 
}


if($slice){
  my %chr_hash = %{&get_chrlengths($db)}; 
  my @input_ids;
  foreach my $chr(keys %chr_hash) {
    #print STDERR "making input ids for ".$chr."\n";
    my $length = $chr_hash{$chr};
    my $count = 1;
    while ($count < $length) {
      my $start = $count;
      my $end   = $count + $size -1;
      
      if ($end > $length) {
	$end = $length;
      }
      
      my $input_id = $chr . "." . $start . "-" .  $end;
      push(@input_ids, $input_id);
      $count = $count + $size;
    }
  }
  if(!$input_ids{$slice}){
    $input_ids{$slice} = [];
    push(@{$input_ids{$slice}}, @input_ids);
  }else{
    push(@{$input_ids{$slice}}, @input_ids);
  }
}


foreach my $logicname(keys(%input_ids)){

  my @input_ids = @{$input_ids{$logicname}};
  my $analysis = $db->get_AnalysisAdaptor->fetch_by_logic_name($logicname);
  &insert_input_ids($analysis, $db, @input_ids);
}
sub insert_input_ids{
  my ($analysis, $db, @input_ids) = @_;
  
  foreach my $input_id(@input_ids){
    
    my $time = time;
    my $sql = "insert into input_id_analysis(input_id, analysis_id, created) values('$input_id', '".$analysis->dbID."', now())";
    #print STDERR $sql."\n";
    my $sth = $db->prepare($sql);
    $sth->execute();
  }
}


sub get_chrlengths{
  my($db) = @_;

  my %chr_hash; 
  my $q = "SELECT c.name, max(a.chr_end) 
           FROM   chromosome c, assembly a
           WHERE  c.chromosome_id = a.chromosome_id
           GROUP BY c.name";
  #print STDERR "trying query ".$q."\n";
  my $sth = $db->prepare($q) || $db->throw("can't prepare: $q");
  my $res = $sth->execute || $db->throw("can't execute: $q");
  
  while( my ($chr, $length) = $sth->fetchrow_array) {
    $chr_hash{$chr} = $length;
  }
  return \%chr_hash;
  #print "have ".keys(%chr_hash)." chromosomes\n";
}
