#!/usr/local/bin/perl -w

=head1 NAME

supp_evidence2feature.pl

=head1 SYNOPSIS
 
supp_evidence2feature.pl

=head1 DESCRIPTION

converts supporting features to features

 +-----------------------+------------------+------+-----+---------+----------------+------------+
| Field                 | Type             | Null | Key | Default | Extra          | Privileges |
+-----------------------+------------------+------+-----+---------+----------------+------------+
| supporting_feature_id | int(10) unsigned |      | PRI | NULL    | auto_increment | select     |
| exon_id               | int(10) unsigned |      | MUL | 0       |                | select     |
| contig_id             | int(10) unsigned |      |     | 0       |                | select     |
| seq_start             | int(10)          |      |     | 0       |                | select     |
| seq_end               | int(10)          |      |     | 0       |                | select     |
| score                 | int(10)          |      |     | 0       |                | select     |
| strand                | int(1)           |      |     | 1       |                | select     |
| analysis              | int(10) unsigned |      | MUL | 0       |                | select     |
| name                  | varchar(40)      |      | MUL |         |                | select     |
| hstart                | int(11)          |      |     | 0       |                | select     |
| hend                  | int(11)          |      |     | 0       |                | select     |
| hid                   | varchar(40)      |      | MUL |         |                | select     |
| evalue                | double(16,4)     | YES  |     | NULL    |                | select     |
| perc_id               | int(10)          | YES  |     | NULL    |                | select     |
| phase                 | tinyint(1)       | YES  |     | NULL    |                | select     |
| end_phase             | tinyint(1)       | YES  |     | NULL    |                | select     |
| hstrand               | tinyint(1)       | YES  |     | NULL    |                | select     |
+-----------------------+------------------+------+-----+---------+----------------+------------+


+-----------+------------------+------+-----+---------+----------------+------------+
| Field     | Type             | Null | Key | Default | Extra          | Privileges |
+-----------+------------------+------+-----+---------+----------------+------------+
| id        | int(10) unsigned |      | PRI | NULL    | auto_increment | select     |
| contig    | int(10) unsigned |      | MUL | 0       |                | select     |
| seq_start | int(10)          |      |     | 0       |                | select     |
| seq_end   | int(10)          |      |     | 0       |                | select     |
| score     | double(16,4)     |      |     | 0.0000  |                | select     |
| strand    | int(1)           |      |     | 1       |                | select     |
| analysis  | int(10) unsigned |      | MUL | 0       |                | select     |
| name      | varchar(40)      | YES  |     | NULL    |                | select     |
| hstart    | int(11)          |      |     | 0       |                | select     |
| hend      | int(11)          |      |     | 0       |                | select     |
| hid       | varchar(40)      |      | MUL |         |                | select     |
| evalue    | varchar(20)      | YES  |     | NULL    |                | select     |
| perc_id   | tinyint(10)      | YES  |     | NULL    |                | select     |
| phase     | tinyint(1)       | YES  |     | NULL    |                | select     |
| end_phase | tinyint(1)       | YES  |     | NULL    |                | select     |
+-----------+------------------+------+-----+---------+----------------+------------+



=head1 OPTIONS

  -refdb
  -refuser
  -refhost
  -refpass
  -file

=cut

use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::FeatureAdaptor;
use Bio::EnsEMBL::Pipeline::ESTConf;
use Getopt::Long;

# ref db holds the exonerate_e2g or genomewise gene/exon/supporting feature data
my $refdbname    = $EST_GENE_DBNAME;
my $refuser      = $EST_GENE_DBUSER;
my $refhost      = $EST_GENE_DBHOST;
my $refpass      = $EST_GENE_DBPASS;

# supp_evidence holds the supp_evidence ids to be converted
my $supp_evidence;

&get_variables();
die ("$supp_evidence does not exist") unless -e $supp_evidence;

print STDERR "database: $refdbname : $refhost\n";
my ($refdb) = get_database($refdbname, $refhost, $refuser, $refpass);
die ("failed to create reference database DBAdaptor") unless defined $refdb;

my $analysis        = get_analysis($refdb);
die ("failed to create analysis object") unless defined $analysis;

&process_features($analysis, $supp_evidence, $refdb);

=head2 process_features

  Function: gets supporting feature info out of refdb, 
            makes feature pairs out of them and writes them to STDOUT
  Returns : nothing
  Args    : $analysis_obj: Bio::EnsEMBL::Analysis; 
;            $file: name of file containing exonerate gff
            $ref_db: Bio::EnsEMBL::DBSQL::DBAdaptor;   

=cut

sub process_features{
  my ($analysis_obj, $file, $refdb) = @_;
  die ("no analysis object")    unless defined $analysis_obj && $analysis_obj->isa("Bio::EnsEMBL::Analysis");
  die ("no exon id file")       unless defined $file && -e $file;
  die ("no reference database") unless defined $refdb && $refdb->isa("Bio::EnsEMBL::DBSQL::DBAdaptor");

  #my $source_tag  = $analysis_obj->program;
  #my $primary_tag = $analysis_obj->program;
  #my $analysis_id = 1;
  #my $name = 'est'; 

  # foreach exon, get out relevant bits of info from $refdb, make a feature pair and store it in $targetdb
  open(EVIDENCE, "<$supp_evidence") or die "Can't open $supp_evidence\n";

  while(<EVIDENCE>){
    my ($id, $contig, $seq_start, $seq_end, $score, $strand, $phase, $end_phase, $hstart, $hend, $hid, $evalue, $perc_id);
    #my $score = 100; # no real score

    chomp;
    my $evi_id = $_;
    my $query = qq(
		   SELECT evi.supporting_feature_id, evi.contig_id, evi.seq_start, evi.seq_end, evi.score, evi.strand, evi.hstart, evi.hend, evi.hid, evi.evalue, evi.perc_id, evi.phase, evi.end_phase
		   FROM  supporting_feature evi
		   WHERE evi.supporting_feature_id = '$evi_id'
		  );

    my $sth = $refdb->prepare($query) or die "can't prepeare query\n";
    my $res = $sth->execute;
    
    $sth->bind_columns(\$id, \$contig, \$seq_start, \$seq_end, \$score, \$strand, \$hstart, \$hend, \$hid, \$evalue, \$perc_id, \$phase, \$end_phase);
    


    # we create a new feature name
    my $name = 'Exonerate';

    # we give the same id as the analysis for the genomewise genes
    my $analysis_id = $analysis_obj->dbID;
    
    my $count = 0;

  FETCH:
    while($sth->fetch){
      if($count) { 
	print STDERR "more than one entry for $evi_id!\n";
	next FETCH;
      }
      $count++;
      
      $evalue    = "NULL" unless $evalue;
      $perc_id   = $score unless $perc_id;
      $phase     = 0 unless $phase;
      $end_phase = 0 unless $end_phase;
      
      $analysis_id = 2;

      #print STDERR "analysis_id = $analysis_id\n";
      #print STDERR "contig_id: $contig\n";
      #print STDERR "id       : $id\n";

print "\\N\t$contig\t$seq_start\t$seq_end\t$score\t$strand\t$analysis_id\t$name\t$hstart\t$hend\t$hid\t$evalue\t$perc_id\t$phase\t$end_phase\n";

    }
  }

  close EVIDENCE or die "error closing $supp_evidence\n";

}


=head2 get_variables

  Title   : get_variables
  Usage   : get_variables
  Function: initialiases global variables according to input parameters. If required parameters 
            are not provided, prints usgae statement and exits script.
  Returns : none - uses globals
  Args    : none - uses globals

=cut

sub get_variables{
&GetOptions( 
	    'refdb:s'      => \$refdbname,
	    'refuser:s'    => \$refuser,
	    'refhost:s'    => \$refhost,
	    'refpass:s'    => \$refpass,
	    'file:s'   => \$supp_evidence,
	   );
  if(!(defined $refhost   && 
       defined $refdbname && 
       defined $refuser && 
       defined $refpass && 
       defined $supp_evidence)){
    print STDERR "Usage: supp_evidence2feature.pl -refdb -refhost -refuser -refpass -file\n";
    exit (1);
  }
}

=head2 get_database

  Title   : get_database
  Usage   : get_database
  Function: initialiases global refdb according to input parameters 
  Returns : none - uses globals
  Args    : none - uses globals

=cut

sub get_database {
  my ($rname, $rhost, $ruser, $rpass) = @_;
  my $rdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $rhost,
					       -user   => $ruser,
					       -dbname => $rname,
					       -pass   => $rpass,
					      );
  return ($rdb);
}

=head2 get_analysis

  Title   : get_analysis
  Usage   : get_analysis
  Function: checks estdb for a pre-existing analysis to attach to features, and 
            makes a new one if necessary
  Returns : Bio::EnsEMBL::Analysis
  Args    : $refdb: Bio::EnsEMBL::DBSQL::DBAdaptor

=cut

sub get_analysis{
  my ($db) = @_;
  die("no ref db provided\n") unless defined $db;
  die("$db is not a Bio::EnsEMBL::DBSQL::DBAdaptor\n") unless $db->isa("Bio::EnsEMBL::DBSQL::DBAdaptor");

  my $logicname  = 'exonerate_e2g';
  my $anaAdaptor = $db->get_AnalysisAdaptor;
  my @analyses   = $anaAdaptor->fetch_by_logic_name($logicname);
  my $analysis;
  
  if(scalar(@analyses) > 1){
    die("panic! > 1 analysis for $logicname\n");
  }
  elsif(scalar(@analyses) == 1){
    $analysis = $analyses[0];
  }
  else{
    # only need to insert ONCE.
    $analysis = new Bio::EnsEMBL::Analysis(
					   -db              => 'dbEST',
					   -db_version      => 1,
					   -program         => 'exonerate_e2g',
					   -program_version => 1,
					   -gff_source      => 'exonerate_e2g',
					   -gff_feature     => 'similarity',
					   -logic_name      => $logicname,
					   -module          => 'FilterESTs_and_E2G',
					  );
  }

  return $analysis;
  
}
