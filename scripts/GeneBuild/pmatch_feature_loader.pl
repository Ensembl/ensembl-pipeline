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
use  Bio::EnsEMBL::Pipeline::GeneConf qw (
					  GB_DBNAME
					  GB_DBUSER
					  GB_DBHOST
					  GB_DBPASS
					  GB_PM_OUTPUT
					 );

my $dbname = $GB_DBNAME;
my $dbuser = $GB_DBUSER;
my $host   = $GB_DBHOST;
my $pass   = $GB_DBPASS;
my $pmfile = $GB_PM_OUTPUT;
if(defined $pmfile && $pmfile ne ''){
  $pmfile .= "/pm_best.out";
}

# usage
if(!defined $dbname  ||
   !defined $dbuser  ||
   !defined $host    ||
   !defined $pass    ||
   !defined $pmfile
  ){
  print  "USAGE: pmatch_feature_loader.pl\nVarious parameters must be set in GeneConf.pm - here are your relevant current settings:\n " .
    "\tGB_DBNAME      = $GB_DBNAME\n" .
    "\tGB_DBUSER      = $GB_DBUSER\n" .
    "\tGB_DBHOST      = $GB_DBHOST\n" .
    "\tGB_DBPASS      = $GB_DBPASS\n" . 
    "\tGB_PM_OUTPUT   = $GB_PM_OUTPUT\n";

  exit(1);
}

# global stuff
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host             => $host,
    -user             => $dbuser,
    -dbname           => $dbname,
    -pass             => $pass,
);

#$db->assembly_type($path);
my $sgpa = $db->get_SliceAdaptor;
my $pmfa = new Bio::EnsEMBL::Pipeline::DBSQL::PmatchFeatureAdaptor($db);

# warn that pm_best.out has chr names in the form: chr_name.chr_start-chr_end
# and that the script will only use chr_name to store it in table pmatch_features
print STDERR "Note: pm_best.out contains chr_name.chr_start-chr_end\n";
print STDERR "pmatch_feature_loafer.pl will only store chr_name in pmatch_feature table\n";

&process_proteins;

# SUBROUTINES #

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

    # eg chr22:10602496,10603128:Q9UGV6:99.4
    if(!/(\S+):(\d+),(\d+):(\S+):(\S+)/){
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

    ## get chr_name from entry chr_name.chr_start-chr_end if necessary
    my $chr_name;
    my $chr_start;
    my $chr_end;
    if ( $chr =~/(\S+)\.(\d+)-(\d+)/ ){
      $chr_name  = $1;
      $chr_start = $2;
      $chr_end   = $3;
      $chr = $chr_name;
 }   

    # not used anymore, we deal with cDNAs separately
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

