#! /usr/local/bin/perl -w

use strict;
use lib 't';
use Test;
use Bio::EnsEMBL::Pipeline::Runnable::MinimalBlast;
use Bio::EnsEMBL::Pipeline::Runnable::BlastDB;
use Bio::SeqIO;

BEGIN { $| = 1; plan test => 18;}

ok(1);

my $verbose;
$verbose = 1 if @ARGV;

my $seqio = Bio::SeqIO->new('-file'   => 't/data/relaxins.fa',
			    '-format' => 'fasta');

ok($seqio);

my $input_seq = $seqio->next_seq;

ok($input_seq->isa("Bio::Seq"));

my $report;

{
  my $blastdb 
    = Bio::EnsEMBL::Pipeline::Runnable::BlastDB->new(
        -dbfile     => 't/data/relaxins.fa',
        -index_type => 'wu_old',
        -type       => 'DNA',
	-workdir    => '/tmp',
        -copy       => 1);

  ok($blastdb->isa("Bio::EnsEMBL::Pipeline::Runnable::BlastDB"));

  ok($blastdb->run);

  my $blast 
    = Bio::EnsEMBL::Pipeline::Runnable::MinimalBlast->new(
        -blastdb  => $blastdb,
        -queryseq => $input_seq,
        -program  => 'wublastn');

  ok($blast->isa("Bio::EnsEMBL::Pipeline::Runnable::MinimalBlast"));

  ok($report = $blast->run);

  ok($report->isa("Bio::Tools::BPlite"));
}

{
  my $blastdb 
    = Bio::EnsEMBL::Pipeline::Runnable::BlastDB->new(
        -dbfile     => 't/data/relaxins.fa',
        -index_type => 'wu_new',
        -type       => 'DNA',
	-workdir    => '/tmp',
        -copy       => 1);

  ok($blastdb->isa("Bio::EnsEMBL::Pipeline::Runnable::BlastDB"));

  ok($blastdb->run);

  my $blast 
    = Bio::EnsEMBL::Pipeline::Runnable::MinimalBlast->new(
        -blastdb  => $blastdb,
        -queryseq => $input_seq,
        -program  => 'wublastn');

  ok($blast->isa("Bio::EnsEMBL::Pipeline::Runnable::MinimalBlast"));

  ok($report = $blast->run);

  ok($report->isa("Bio::Tools::BPlite"));
}

{
  my $blastdb 
    = Bio::EnsEMBL::Pipeline::Runnable::BlastDB->new(
        -dbfile     => 't/data/relaxins.fa',
        -index_type => 'ncbi',
        -type       => 'DNA',
	-workdir    => '/tmp',
        -copy       => 1);

  ok($blastdb->isa("Bio::EnsEMBL::Pipeline::Runnable::BlastDB"));

  ok($blastdb->run);

  my $blast 
    = Bio::EnsEMBL::Pipeline::Runnable::MinimalBlast->new(
        -blastdb  => $blastdb,
        -queryseq => $input_seq,
        -program  => 'blastn');

  ok($blast->isa("Bio::EnsEMBL::Pipeline::Runnable::MinimalBlast"));

  ok($report = $blast->run);

  ok($report->isa("Bio::Tools::BPlite"));
}


if ($verbose){
  print "Your run produced a : [" . $report . "]\n";
  print "Query      : [" . $report->query    . "]\n";
  print "Database   : [" . $report->database . "]\n";
  my $hits;
  while(my $sbjct = $report->nextSbjct) {
    
    print "Hit       : [" . $sbjct->name . "]\n";
    $hits++;
  }
  print "Total Hits : [$hits]\n";
}
