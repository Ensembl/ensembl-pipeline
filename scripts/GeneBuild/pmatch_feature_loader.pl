#!/usr/local/bin/perl -w

=head1 NAME

  pmatch_feature_loader.pl -  eventually to be incorporated into genebuild pmatch scripts

=head1 SYNOPSIS
 
  pmatch_feature_loader.pl

=head1 DESCRIPTION

 loads pmatch features into database

=head1 OPTIONS

  -dbname
  -dbuser
  -host
  -pass
  -pmfile

=cut

use strict;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::PmatchFeatureAdaptor;
require "Bio/EnsEMBL/Pipeline/GB_conf.pl";

my $dbname = $::db_conf{'dbname'};
my $dbuser = $::db_conf{'dbuser'};
my $host   = $::db_conf{'dbhost'};
my $pass   = $::db_conf{'dbpass'};
my $pmfile = $::scripts_conf{'pm_output'};
if(defined $pmfile && $pmfile ne ''){
  $pmfile .= "/pm_best.out";
}
my $cdnafile = $::scripts_conf{'cdna_pairs'};
my $create = 0;
my $path = $::db_conf{'golden_path'};
$path = 'UCSC' unless (defined $path && $path ne '');

&GetOptions( 
	    'create_tables' => \$create,
	   );

# usage
if(!defined $dbname  ||
   !defined $dbuser  ||
   !defined $host    ||
   !defined $pass    ||
   !defined $pmfile
  ){
  print  "USAGE: pmatch_feature_loader.pl [-create_tables]\nVarious parameters must be set in GB_conf.pl - here are your relevant current settings:\n " .
    "db_conf:\n" .
    "\tdbname      = $::db_conf{'dbname'}\n" .
    "\tdbuser      = $::db_conf{'dbuser'}\n" .
    "\tdbhost      = $::db_conf{'dbhost'}\n" .
    "\tdbpass      = $::db_conf{'dbpass'}\n" . 
    "\tgolden_path = $::db_conf{'golden_path'}\n" . 
    "scripts_conf:\n" .
    "\tpm_output   = $::scripts_conf{'pm_output'}\n" .
    "\tcdna_pairs  = $::scripts_conf{'cdna_pairs'}\n";

  exit(1);
}

# global stuff
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host             => $host,
    -user             => $dbuser,
    -dbname           => $dbname,
    -pass             => $pass,
);

$db->static_golden_path_type($path);
my $sgpa = $db->get_StaticGoldenPathAdaptor;
my $pmfa = new Bio::EnsEMBL::Pipeline::DBSQL::PmatchFeatureAdaptor($db);

if($create) { 
  &create_tables ;
}

#&process_cdnas;
&process_proteins;

# SUBROUTINES #

=head2 create_tables

 Title   : create_tables
 Usage   :
 Function: 
 Example :
 Returns : 
 Args    :


=cut

sub create_tables{
  my @sql = split /;/, ($pmfa->create_sql);
  my $create_sql= join "\n", @sql;
  
  my $sth = $db->prepare($sql[1]);
  $sth->execute;

  $sth = $db->prepare($sql[2]);
  $sth->execute;
}

=head2 process_proteins

 Title   : process_proteins
 Usage   :
 Function: 
 Example :
 Returns : 
 Args    :


=cut

sub process_proteins {
  my @features;
  open(PMATCHES, "<$pmfile") or die "Cannot open [$pmfile] to read pmatches\nClean out protein and pmatch_feature tables before rerunning script!\n";

  while(<PMATCHES>){

    # chr22:10602496,10603128:Q9UGV6:99.4
    if(!/(chr\S+):(\d+),(\d+):(\S+):(\S+)/){
      die "Cannot parse [$_]\nClean out protein and pmatch_feature tables before rerunning script!\n";
    }

    my $chr = $1;
    my $start = $2;
    my $end = $3;
    my $protein = $4;
    my $coverage = $5;

    if($start > $end){
      $start = $3;
      $end = $2;
    }

    my $cdna_id = $pmfa->get_cdna_id($protein);

    my $pmf = new Bio::EnsEMBL::Pipeline::PmatchFeature(-protein_id  => $protein,
							-start       => $start,
							-end         => $end,
							-chr_name    => $chr,
							-cdna_id     => $cdna_id,
							-coverage    => $coverage,
						       );
    push(@features, $pmf);

  }
  
  close PMATCHES;

  $pmfa->write_PmatchFeatures(@features);
}

=head2 process_cdnas

 Title   : process_cdnas
 Usage   :
 Function: 
 Example :
 Returns : 
 Args    :


=cut

sub process_cdnas {
  open(CDNA, "<$cdnafile") or die "Cannot open [$cdnafile] to read cdnas\nClean out protein and pmatch_feature tables before rerunning script!\n";

  while(<CDNA>){
    chomp;
#    NP_116272.1 : 
#    NP_115616.1 : AK027172
    if(!/(\S+) : (\S*)/){
      die "Cannot parse protein cdna pairs from [$_]\nClean out protein and pmatch_feature tables before rerunning script!\n";
    }

    my $protein = $1;
    my $cdna = $2;

    $pmfa->write_protein($protein, $cdna);
  }

  close CDNA;
}
