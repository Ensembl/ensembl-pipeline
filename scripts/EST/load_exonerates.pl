#!/usr/local/bin/perl -w
=head1 NAME

  load_exonerates.pl

=head1 SYNOPSIS
 
  load_exonerates.pl
  loads exonerate results into the fature table os a given database using the contig 
  internal ids of a reference database.

=head1 DESCRIPTION


=head1 OPTIONS

    -esthost      host name for est database (gets put as host= in locator)
    -refhost      host name for ref database (gets put as host= in locator)

    -estport      For RDBs, what port to connect to (port= in locator)
    -refport      For RDBs, what port to connect to (port= in locator)

    -estdbname    For RDBs, what name to connect to (dbname= in locator)
    -refdbname    For RDBs, what name to connect to (dbname= in locator)

    -estdbuser    For RDBs, what username to connect as (dbuser= in locator)
    -refdbuser    For RDBs, what username to connect as (dbuser= in locator)

    -estdbpass    For RDBs, what password to use (dbpass= in locator)
    -refdbpass    For RDBs, what password to use (dbpass= in locator)

    -estfile      File containing est gff to be written in.

=cut

use strict;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::FeatureAdaptor;
use Bio::EnsEMBL::Analysis;

my $esthost   = undef;
my $refhost   = undef;
my $estport   = '410000';
my $refport   = '410000';
my $estdbname = undef;
my $refdbname = undef;
my $estdbuser = undef;
my $refdbuser = undef;
my $estpass   = undef;
my $refpass   = undef;
my $estfile   = undef;

&get_variables();
die ("$estfile does not exist") unless -e $estfile;

my ($refdb, $estdb) = get_databases($refdbname, $refhost, $refdbuser, $estdbname, $esthost, $estdbuser);
die ("failed to create reference database DBAdaptor") unless defined $refdb;
die ("failed to create est database DBAdaptor") unless defined $estdb;

my $analysis        = get_analysis($estdb);
die ("failed to create analysis object") unless defined $analysis;

&process_ests($analysis, $estfile, $refdb, $estdb);


=head2 process_ests

  Title   : process_ests
  Usage   : process_ests
  Function: parses lines in estfile, makes feature pairs out of them and write them to the est database
            they are attached to a contig with internal_id derived from the reference database.
  Returns : nothing
  Args    : $analysis_obj: Bio::EnsEMBL::Analysis; $file: name of file containing exonerate gff
            $ref_db: Bio::EnsEMBL::DBSQL::DBAdaptor;   $est_db: Bio::EnsEMBL::DBSQL::DBAdaptor

=cut

sub process_ests {
  my ($analysis_obj, $file, $ref_db, $est_db) = @_;
  die ("no analysis object") unless defined $analysis_obj && $analysis_obj->isa("Bio::EnsEMBL::Analysis");
  die ("no exonerate resutls file") unless defined $file && -e $file;
  die ("no reference database") unless defined $refdb && $refdb->isa("Bio::EnsEMBL::DBSQL::DBAdaptor");
  die ("no est database") unless defined $estdb && $estdb->isa("Bio::EnsEMBL::DBSQL::DBAdaptor");

  my $source_tag  = $analysis_obj->program;
  my $primary_tag = $analysis_obj->program;
  my $feat_adaptor = $est_db->get_FeatureAdaptor;
  my %contig_features;

  open(ESTFILE, "<$file") or die "Can't open $file:$!\n";
 ESTHIT:
  while(<ESTFILE>){
    next ESTHIT if (/runnable/i);
     # AC004073.1.1.79612      75445   76005   2293.00 89      1       gi|12779912|emb|AL516419.1|AL516419     296     861     1
    my ($contig_id, $contig_start, $contig_end, $score, $percent_id, $contig_strand, 
	$est_id, $est_start, $est_end, $est_strand) = split();
   
    $est_id =~ s/\S+\|\S+\|\S+\|(\S+)\|\S+/$1/;

    # if both contig and est strands are the same, convention is to set both to be 1
    # if they differ, convention is to set contig strand to -1, est strand to 1
    if($contig_strand == $est_strand){
      $contig_strand = 1;
      $est_strand    = 1;
    }
    else{
      $contig_strand = -1;
      $est_strand = 1;
    }

    my $feat1 = new Bio::EnsEMBL::SeqFeature  (-start       =>   $contig_start,
					       -end         =>   $contig_end,
					       -seqname     =>   $contig_id,
					       -strand      =>   $contig_strand,
					       -score       =>   $score,
					       -percent_id  =>   $percent_id, 
					       -source_tag  =>   $source_tag,
					       -primary_tag =>   $primary_tag,
					       -analysis    =>   $analysis_obj );
    
    my $feat2 = new Bio::EnsEMBL::SeqFeature  (-start       =>   $est_start,
					       -end         =>   $est_end,
					       -seqname     =>   $est_id,
					       -strand      =>   $est_strand,
					       -score       =>   $score,
					       -percent_id  =>   $percent_id, 
					       -source_tag  =>   $source_tag,
					       -primary_tag =>   $primary_tag,
					       -analysis    =>   $analysis_obj );
    
    #create featurepair
    my $fp = new Bio::EnsEMBL::FeaturePair  (-feature1 => $feat1,
                                             -feature2 => $feat2) ;

    push(@{$contig_features{$contig_id}}, $fp);
  }
  
  foreach my $contig_id( keys %contig_features){
    my $contig;
    eval{
      $contig = $ref_db->get_Contig($contig_id);
    };
    if($@){
      print STDERR "No contig for $contig_id part 1\n$@\n";
      next ESTHIT;
    }

    # lock db only once per contig - may still be too slow.
    my @features = @{$contig_features{$contig_id}};
    $feat_adaptor->store($contig, @features);
  }
  
  close ESTFILE;
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
	      'esthost:s'      => \$esthost,
	      'refhost:s'      => \$refhost,
	      'estport:n'      => \$estport,
	      'refport:n'      => \$refport,
	      'estdbname:s'    => \$estdbname,
	      'refdbname:s'    => \$refdbname,

	      'estdbuser:s'    => \$estdbuser,
	      'refdbuser:s'    => \$refdbuser,
	      'estpass:s'      => \$estpass,
	      'refpass:s'      => \$refpass,
	      'estfile:s'      => \$estfile,
	     );
  
  if(!(defined $esthost   && defined $refhost &&
       defined $estdbname && defined $refdbname &&
       defined $estdbuser && defined $refdbuser &&
       defined $estfile)){
    print "Usage: load_exonerates.pl -estdbname -refdbname -esthost -refhost -estdbuser -refdbuser -estfile\n" .
      "optional arguments: -estdbport -refdbport -estdbpass -refdbpass\n";
    exit (1);
  }
}

=head2 get_databases

  Title   : get_databases
  Usage   : get_databases
  Function: initialiases global refdb and estdb according to input parameters 
  Returns : none - uses globals
  Args    : none - uses globals

=cut

sub get_databases{
  my ($rname, $rhost, $ruser, $ename, $ehost, $euser) = @_;
  my $rdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $rhost,
					       -user   => $ruser,
					       -dbname => $rname,
					      );
  my $edb = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $ehost,
					       -user   => $euser,
					       -dbname => $ename,
					      );
  return ($rdb, $edb);
}

=head2 get_analysis

  Title   : get_analysis
  Usage   : get_analysis
  Function: checks estdb for a pre-existing analysis to attach to features, and 
            makes a new one if necessary
  Returns : Bio::EnsEMBL::Analysis
  Args    : $estdb: Bio::EnsEMBL::DBSQL::DBAdaptor

=cut

sub get_analysis{
  my ($db) = @_;
  die("no est db provided\n") unless defined $db;
  die("$db is not a Bio::EnsEMBL::DBSQL::DBAdaptor\n") unless $db->isa("Bio::EnsEMBL::DBSQL::DBAdaptor");

  my $logicname  = 'est_exonerate';
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
#					   -db              => 'dbEST',
					   -db              => 'MGC',
					   -db_version      => 1,
					   -program         => 'exonerate',
					   -program_version => 3,
					   -gff_source      => 'exonerate',
					   -gff_feature     => 'similarity',
					   -logic_name      => $logicname,
					   -module          => 'ExonerateESTs',
					  );
  }

  return $analysis;
  
}
