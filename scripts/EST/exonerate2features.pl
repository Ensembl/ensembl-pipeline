#!/usr/local/bin/perl -w
use strict;
use Bio::EnsEMBL::Pipeline::ESTConf;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

#$| = 1;


# contig ids are of the form >AB015355.1.1.43999
# look for instance /data/blastdb/Ensembl/NCBI_28_dusted_masked_golden_contigs.fa

# in mouse they are of the form internal_id ( a number )  stable_id ( e.g. 1.77500001-78000000 )

### parameters ###
my $refdb_host  = $EST_REFDBHOST;
my $refdb_name  = $EST_REFDBNAME;
my $refdb_user  = $EST_REFDBUSER;
my $refdb_path  = $EST_GOLDEN_PATH;

my $count = 0;

my $skip_count = 0;
my $feat_count = 0;
my $bad_count  = 0;

my $feature_file;

&GetOptions( 
	    'feature_file:s'  => \$feature_file,
	   );

unless ($feature_file){
  print STDERR "Usage: $0 -feature_file estfeatures > & log_file\n";
  print STDERR "Run the script in the directory where the exonerate out is\n";
  print STDERR "It assumes the output is in gzipped files\n";
  exit(0);
}

my $refdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $refdb_host,
					       -user   => $refdb_user,
					       -dbname => $refdb_name,
					      );

$refdb->static_golden_path_type( $refdb_path );

# make sure you've copied the analysisprocess table of the previous EST release
# into the current one, then set the analysis_id to be the same as in that table for exonerate:

my $analysis = 1; 
my $name ='exonerate';

# get genomic regions
my $internal_id = get_chrlengths($refdb,$refdb_path);
my %internal_id = %$internal_id;

open LS, 'ls |' or die "Can't open pipe from ls : $!";

open( OUT,">$feature_file" ) or die("cannot open file $feature_file");

FILES:
while(<LS>){
  chomp;
  my $file = $_;
  my $newfile;
  next FILES if ( $file =~/(\.)+\// );
  next FILES if ( $file =~/estfeatures/ );
  
  #my $start = 6391;
  #my $where;
  #if ( $file =~/_chunk_(\d+)_(\d+)/ ){
  #  $where = $1;
  #}
  #if ( $where < $start ){
  #  next FILES;
  #}

  if ( $file =~ /(\S+).gz/){
    $newfile = $1;
    print STDERR "UNCompressing $file\t";
    system("gunzip $file");
  }
  else{
    $newfile = $file;
  }
  
  print STDERR "PROCessing $newfile\n";
  #system("/nfs/acari/eae/ensembl-scripts/ESTs/exonerate2mysqltab_human.pl < $newfile >> /scratch2/ensembl/eae/est_NCBI_28/exonerate_est/results/estfeatures");
    
  open FILE, "<$newfile";
  
  while(<FILE>){
    
    # AC004073.1.1.79612   75445  76005 2293.00  89  1  gi|12779912|emb|AL516419.1|AL516419     296     861     1
    
    my ($c_id, 
	$c_start, 
	$c_end, 
	$score, 
	$perc_id, 
	$c_strand, 
	$est_id, 
	$est_start, 
	$est_end, 
	$est_strand) = split();

    # Note: Exonerate puts the beginning of the feature at 0, not at 1!!
    if ( $est_start == 0 ){
      $est_start = 1;
    }

    # print out to test:
    #print STDERR "est_start: $est_start - est_end: $est_end\n";
    
    # only take features above 90% percentage identity
    #if ( $perc_id < 90 ){
    #  print STDERR "rejected $est_id for LOW perc_id $perc_id < 90\n";
    #  $skip_count++;
    #  next;
    #}
    
    # calculate the est_id
    #if ( $est_id =~ /\S+\|\S+\|\S+\|(\S+)\|\S+/ ){
    #  $est_id =~ s/\S+\|\S+\|\S+\|(\S+)\|\S+/$1/;
    #}
    #if ( $est_id =~ /PROTEIN/ 
#	 || $est_id =~ /gene product/ 
	# || $est_id eq "58" 
#	 || $est_id eq "1"
	# || $est_id eq "4-ISOMERASE"
     #  ){
     # print STDERR "skipping est $est_id for dodgy id\n";
     # next;
    #}

    
    # in which internal_id contig this EST falls
    my $int_id = $internal_id{ $c_id };
    
    if ( !( $int_id) ){
      print STDERR "NOT defined int_id, skipping feature in $c_id\n";
      $bad_count++;
      next;
    }
    
    ######################
    # strand conventions #
    ######################
    #
    # if both contig and est strands are the same, convention is to set both to be 1
    # if they differ, convention is to set contig strand to -1, est strand to 1
    if($c_strand == $est_strand){
      $c_strand = 1;
      $est_strand    = 1;
    }
    else{
      $c_strand = -1;
      $est_strand = 1;
    }
    

    # need to print out tab delimited line suitable for mysql table 'feature' in ensembl database
    # Autogenerate feature internal_id
    $feat_count++;
    print OUT "\\N\t$int_id\t$c_start\t$c_end\t$score\t$c_strand\t$analysis\t$name\t$est_start\t$est_end\t$est_id\tNULL\t$perc_id\t0\t0\n";
    
    ### recall:
    #mysql> describe feature;
    #+-----------+------------------+------+-----+---------+----------------+---------------------------------+
    #| Field     | Type             | Null | Key | Default | Extra          | Privileges                      |
    #+-----------+------------------+------+-----+---------+----------------+---------------------------------+
    #| id        | int(10) unsigned |      | PRI | NULL    | auto_increment | select,insert,update,references |
    #| contig    | int(10) unsigned |      | MUL | 0       |                | select,insert,update,references |
    #| seq_start | int(10)          |      |     | 0       |                | select,insert,update,references |
    #| seq_end   | int(10)          |      |     | 0       |                | select,insert,update,references |
    #| score     | double(16,4)     |      |     | 0.0000  |                | select,insert,update,references |
    #| strand    | int(1)           |      |     | 1       |                | select,insert,update,references |
    #| analysis  | int(10) unsigned |      | MUL | 0       |                | select,insert,update,references |
    #| name      | varchar(40)      | YES  |     | NULL    |                | select,insert,update,references |
    #| hstart    | int(11)          |      |     | 0       |                | select,insert,update,references |
    #| hend      | int(11)          |      |     | 0       |                | select,insert,update,references |
    #| hid       | varchar(40)      |      | MUL |         |                | select,insert,update,references |
    #| evalue    | varchar(20)      | YES  |     | NULL    |                | select,insert,update,references |
    #| perc_id   | tinyint(10)      | YES  |     | NULL    |                | select,insert,update,references |
    #| phase     | tinyint(1)       | YES  |     | NULL    |                | select,insert,update,references |
    #| end_phase | tinyint(1)       | YES  |     | NULL    |                | select,insert,update,references |
    #+-----------+------------------+------+-----+---------+----------------+---------------------------------+
        
    
  }  # close loop over FILE

  close(FILE);
  system("gzip $newfile");

}    # close loop over LS


close(OUT);
close(LS);



print STDERR 
  "features: $feat_count\n".
  "skipped for low similarity     : $skip_count\n".
  "skipped for undefined contig id: $bad_count\n";



sub get_chrlengths{
  my $db   = shift;
  my $type = shift;

  if (!$db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor')) {
    die "get_chrlengths should be passed a Bio::EnsEMBL::DBSQL::DBAdaptor\n";
  }

  my %internal_id;

  my $q = qq(  SELECT raw_id,id 
	       FROM   static_golden_path sgp, contig c 
	       WHERE  sgp.raw_id = c.internal_id
	    );

  my $sth = $db->prepare($q) || $db->throw("can't prepare: $q");
  my $res = $sth->execute || $db->throw("can't execute: $q");

  while( my ($internal_id, $stable_id) = $sth->fetchrow_array) {
    $internal_id{$stable_id} = $internal_id;
  }
  return \%internal_id;
}
