#!/usr/local/ensembl/bin/perl

=head1 NAME

predict ncRNA

=head1 SYNOPSIS

predict_ncRNAs.pl -dbhost ecs1a -dbuser ensadmin -dbpass **** -dbname pipeline_db 

=head1 DESCRIPTION

gets the best scoring 2000 blast hits fo each rfam familly and flags them for Infernal to 
run on

=head1 OPTIONS

   -dbhost    host name for database (gets put as host= in locator)

    -dbport    what port to connect to (port= in locator)

    -dbname    what name to connect to (dbname= in locator)

    -dbuser    what username to connect as (user= in locator)

    -dbpass    what password to use (pass= in locator)

    -chunk     Chunk size for input ids (default 50)

    -aln       Analysis Logic Name - name of analysis you want to run on predictions

    -dln       Dummy Logic Name - Dummy analysis to load into input_id_analysis

=cut

use strict;
use Getopt::Long;

use Bio::EnsEMBL::Utils::Exception qw(stack_trace);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::FlagAdaptor;
use Bio::EnsEMBL::Pipeline::Flag;
use Bio::EnsEMBL::Analysis::Config::Databases;
use Bio::EnsEMBL::Pipeline::Utils::InputIDFactory;
use Bio::SeqIO;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbname;
my $dbport;
my $chunk = 50;
my $dln;
my $aln;
my $help;
my %RFAM;
my @flags;
my @dbID;
my $count =0;
my $num=0;
my $start;
my $analysis_logic_name;
my $dummy_logic_name;
my @input_ids;
my @multiexon_files;

GetOptions(
	   '-dbname=s'  => \$dbname,
	   '-dbhost=s'  => \$dbhost,
	   '-dbuser=s'  => \$dbuser,
	   '-dbpass=s'  => \$dbpass,
	   '-dbport=s'  => \$dbport,
	   '-chunk=s'   => \$chunk,
           '-h|help'    => \$help,
	   '-aln=s'  => \$aln,
	   '-dln=s'     => \$dln,
           ) or &usage();

if(!$dbhost || !$dbuser || !$dbname || !$dbpass || $help){
  warn("Can't run without -dbhost $dbhost -dbuser $dbuser -dbname $dbname -dbpass $dbpass");
  $help = 1;
}

if ($help) {
    exec('perldoc', $0);
}


my $sdb = new Bio::EnsEMBL::DBSQL::DBAdaptor
  (
   -host   => $dbhost,
   -user   => $dbuser,
   -pass   => $dbpass,
   -port   => $dbport,
   -dbname => $dbname,
  );
my $dna_db = new Bio::EnsEMBL::DBSQL::DBAdaptor
  (
   '-host'   => $GB_DBHOST,
   '-user'   => $GB_DBUSER,
   '-dbname' => $GB_DBNAME,
   '-pass'   => $GB_DBPASS,
   '-port'   => $GB_DBPORT,
  );
#add dna_db
$sdb->dnadb($dna_db);

my $flag_dba = new Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
  (
   '-host'   => $dbhost,
   '-user'   => $dbuser,
   '-dbname' => $dbname,
   '-pass'   => $dbpass,
   '-port'   => $dbport,
  );


my $fa = Bio::EnsEMBL::Pipeline::DBSQL::FlagAdaptor->new($sdb);
my $aa = $sdb->get_AnalysisAdaptor;
my $sa = $sdb->get_SliceAdaptor;
my $analysis = $aa->fetch_by_logic_name($aln);
my $limit = 0;
my %rfam_blasts;
my %rfam_threshold;
my $dafa = $sdb->get_DnaAlignFeatureAdaptor;
my @dafs;
# fetch them by familly!!!!!!!!!
my %thr;
open( T, "/data/blastdb/Rfam/Rfam.thr" ) or die("can't file the Rfam.thr file");
while(<T>) {
  if( /^(RF\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\S+)\s*$/ ) {
    $thr{ $1 } = { 'id' => $2, 'thr' => $3, 'win' => $4, 'mode' => $5 };
  }
}
close T;
my $total = scalar(keys (%thr));
foreach my $domain (keys %thr){
  $count ++;
  print "Looking at domain $domain\t";
  my @dafs = @{$dafa->generic_fetch("left(hit_name,7) = \"$domain\"")};
  print int($count/$total*100)."% completed\n";
  @dafs = sort {$a->p_value <=> $b->p_value} @dafs if (scalar(@dafs) > 2000 );
 DAF:  foreach my $daf(@dafs){
    next if ($daf->score < 20);
    #    print $daf->p_value." ";
    last DAF if ($rfam_blasts{$domain} &&  scalar(@{$rfam_blasts{$domain}}) >= 2000 );
    push @{$rfam_blasts{$domain}},$daf;
  }
}

print "Generating Flags\n";

foreach my $domain(keys %rfam_blasts){
  my @hits = @{$rfam_blasts{$domain}};
  foreach my $hit (@hits){
    my $flag = Bio::EnsEMBL::Pipeline::Flag->new
      (
       '-type'         => 'dna_align_feature',
       '-ensembl_id'   => $hit->dbID,
       '-goalAnalysis' => $analysis,
      );
    push @flags,$flag;
  }
}

print "Storing and randomizing flags so the input ids dont run in famillies\n";
while (scalar(@flags >= 1)){
  my $random =  rand($#flags);
  # make sure you have an integer
  $random = int $random;
  $fa->store($flags[$random]);
  splice(@flags,$random,1);
}

$count =0;
print "Making input ids for $aln\n";

die("analysis object not found $analysis_logic_name\n") unless ($analysis);
my @ids = @{$fa->fetch_by_analysis($analysis)};
@ids = sort {$a->dbID <=> $b->dbID} @ids;
foreach my $id (@ids){
  if ($count==0){
    $start = $id->dbID; 
  }
  $count++;
  if ($count == $chunk){
    $num++;
    push @input_ids,"$start:".$id->dbID;
    $count=0;
  }
}
if ($count >0){
  push @input_ids,"$start:".$ids[$#ids]->dbID;
}

my $inputIDFactory = new Bio::EnsEMBL::Pipeline::Utils::InputIDFactory
  (
   -db => $flag_dba,
   -top_level => 'top_level',
   -slice     => 1,
   -logic_name => $dln,
    );
$inputIDFactory->input_ids(\@input_ids);
$inputIDFactory->store_input_ids;

print "Finished\n";



exit;


