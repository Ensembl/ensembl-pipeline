#!/usr/local/ensembl/bin/perl

use strict;
use Getopt::Long;

use Bio::EnsEMBL::Utils::Exception qw(stack_trace);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis::Config::Databases;
use Bio::EnsEMBL::Analysis::Config::Pseudogene;
use Bio::EnsEMBL::Pipeline::DBSQL::FlagAdaptor;
use Bio::EnsEMBL::Pipeline::Utils::InputIDFactory;
use Bio::SeqIO;


my $usage = "prepare_SplicedElsewhere.pl -logic_name <dummy analysis which will run spliced elsewhere>
Loads input ids into pseudo_db and gets output files from pseudogene and makes them into a blast db";
my @dbID;
my $count =0;
my $num=0;
my $start;
my $logic_name;
my @input_ids;
my @multiexon_files;

&GetOptions('-logic_name:s' => \$logic_name);

die $usage unless($logic_name);

my $ref_db = new Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
  (
   '-host'   => $GB_DBHOST,
   '-user'   => $GB_DBUSER,
   '-dbname' => $GB_DBNAME,
   '-pass'   => $GB_DBPASS,
   '-port'   => $GB_DBPORT,
  );

if ($SINGLE_EXON) {
  print "Making input ids for single exon genes\n";

  my $fa = Bio::EnsEMBL::Pipeline::DBSQL::FlagAdaptor->new($ref_db);
  my $aa = $ref_db->get_AnalysisAdaptor;
  my $analysis = $aa->fetch_by_logic_name($SINGLE_EXON);
  my $multifile = $PS_MULTI_EXON_DIR."/all_multi_exon_genes.fasta";
  die("analysis object not found $SINGLE_EXON\n") unless ($analysis);
  my @ids = @{$fa->fetch_by_analysis($analysis)};
  @ids = sort {$a->dbID <=> $b->dbID} @ids;
  foreach my $id (@ids){
    if ($count==0){
       $start = $id->dbID; 
     }
    $count++;
    if ($count == $PS_CHUNK){
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
     -db => $ref_db,
     -top_level => 'top_level',
     -slice     => 'ignore this warning its just stupid',,
     -logic_name => $logic_name,
    );
  $inputIDFactory->input_ids(\@input_ids);
  $inputIDFactory->store_input_ids;

  print STDERR "Pooling multiexon genes into single blastDB .\n";

  my $db_output = Bio::SeqIO->new(
				  -file   => ">$multifile",
				  -format => 'fasta'
				 );

  unless (opendir(DIR, $PS_MULTI_EXON_DIR)) {
    closedir(DIR);
    die "cannot read files from $PS_MULTI_EXON_DIR.";
  }
  foreach (readdir(DIR)) {
    my $file = "$_";
    if($file =~ m/^multi_exon_seq.*\.fasta$/){
      my $bioseq = Bio::SeqIO->new(
				   -file   => $PS_MULTI_EXON_DIR."/".$file,
				   -format => 'fasta'
				  );
      while (my $seq = $bioseq->next_seq) {
	$db_output->write_seq($seq);
      }
    }
  }
 #   system ("rm $file");
  system ("xdformat -p $multifile");
}

print "Finished\n";

exit 0;
