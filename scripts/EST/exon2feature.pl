#!/usr/local/bin/perl -w

=head1 NAME

  exon2feature.pl

=head1 SYNOPSIS
 
  exon2feature.pl

=head1 DESCRIPTION

  exon2feature.pl converts exons and their supporting feature data into simple feature table entries.

=head1 OPTIONS

  -refdb
  -refuser
  -refhost
  -refpass
  -exonfile

=cut

use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;

# ref db holds the exonerate_e2g gene/exon/supporting feature data
my $refdbname    = 'mouse_sanger_Oct01_est';
my $refuser      = 'ensro';
my $refhost      = 'ecs1f';
my $refpass      = '';

# exonfile holds a list of exons to be converted
# select exon_id from exon into outfile
my $exonfile;


&get_variables();
die ("$exonfile does not exist") unless -e $exonfile;

my ($refdb) = get_database($refdbname, $refhost, $refuser, $refpass);
die ("failed to create reference database DBAdaptor") unless defined $refdb;

my $analysis        = get_analysis($refdb);
die ("failed to create analysis object") unless defined $analysis;

&process_exons($analysis, $exonfile, $refdb);

=head2 process_exons

  Title   : process_exons
  Usage   : process_exons
  Function: gets exon and supporting feature info out of refdb, 
            makes feature pairs out of them and writes them to STDOUT
  Returns : nothing
  Args    : $analysis_obj: Bio::EnsEMBL::Analysis; 
            $file: name of file containing exonerate gff
            $ref_db: Bio::EnsEMBL::DBSQL::DBAdaptor;   

=cut

sub process_exons{
  my ($analysis_obj, $file, $refdb) = @_;
  die ("no analysis object")    unless defined $analysis_obj && $analysis_obj->isa("Bio::EnsEMBL::Analysis");
  die ("no exon id file")       unless defined $file && -e $file;
  die ("no reference database") unless defined $refdb && $refdb->isa("Bio::EnsEMBL::DBSQL::DBAdaptor");

  my $source_tag  = $analysis_obj->program;
  my $primary_tag = $analysis_obj->program;
  my $analysis_id = 1;
  my $name = 'est'; 
  my %contig_features;

  # foreach exon, get out relevant bits of info from $refdb, make a feature pair and store it in $targetdb
  open(EXON, "<$exonfile") or die "Can't open $exonfile\n";

  while(<EXON>){
    my ($id, $contig, $seq_start, $seq_end, $strand, $phase, $end_phase, $hstart, $hend, $hid);
    my $score = 100; # no real score

    chomp;
    my $exon_id = $_;

    #    print STDERR "$exon_id\n";

    my $query = qq(
		   SELECT exon.exon_id, exon.contig_id, exon.seq_start, exon.seq_end, exon.strand, exon.phase, exon.end_phase, min(hstart), max(hend), hid
		   FROM exon, supporting_feature
		   WHERE exon.exon_id = supporting_feature.exon_id
		   AND exon.exon_id = '$exon_id'
		   GROUP BY hid
		  );

    my $sth = $refdb->prepare($query) or die "can't prepeare query\n";
    my $res = $sth->execute;
    
    # bind the columns
    $sth->bind_columns(undef, \$id, \$contig, \$seq_start, \$seq_end, \$strand, \$phase, \$end_phase, \$hstart, \$hend, \$hid);

    my $count = 0;
  FETCH:
    while($sth->fetch){
      if($count) { 
	print STDERR "more than one entry for $exon_id!\n";
	next FETCH;
      }
      $count++;

print "\\N\t$contig\t$seq_start\t$seq_end\t$score\t$strand\t$analysis_id\t$name\t$hstart\t$hend\t$hid\t'NULL'\t'NULL'\t$phase\t$end_phase\n";

    }

    print STDERR "$exon_id\n" unless $count;
  }

  close EXON or die "error closing $exonfile\n";

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
	    'exonfile:s'   => \$exonfile,
	   );
  if(!(defined $refhost   && 
       defined $refdbname && 
       defined $refuser && 
       defined $refpass && 
       defined $exonfile)){
    print "Usage: exon2feature.pl -refdb -refhost -refuser -refpass -exonfile\n";
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

  my $logicname  = 'ex_e2g_feat';
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
					   -program_version => 3,
					   -gff_source      => 'est',
					   -gff_feature     => 'similarity',
					   -logic_name      => $logicname,
					   -module          => 'Filter_ESTs_and_E2G',
					  );
  }

  return $analysis;
  
}
