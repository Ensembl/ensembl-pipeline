#!/usr/local/bin/perl -w
use strict;
use Getopt::Long;

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Databases;
use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Scripts qw(GB_PROTEOME_FILES);
use Bio::EnsEMBL::DBSQL::DBAdaptor;

# script to identify which NM goes with which NP and insert that info into
# the protein table of the pipeline DB
# The script assumes the early stages (pmatch, targetted) of the build have
# already been run ie the protein table already has protein_id information in it.
# options can either be passed in on the command line or read from config files.

my $refseqfile;
my $dbhost;
my $dbport;
my $dbname;
my $dbuser;
my $dbpass;

&GetOptions(
	    'refseq:s' => \$refseqfile,
	    'dbhost:s' => \$dbhost,
	    'dbport:n' => \$dbport,
	    'dbname:s' => \$dbname,
	    'dbuser:s' => \$dbuser,
	    'dbpass:s' => \$dbpass,
	   );

&check_parameters;

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $dbhost,
					    -user   => $dbuser,
					    -dbname => $dbname,
					    -pass   => $dbpass,
					    -port   => $dbport,
					   );

# update the protein table in the pipeline database, adding in cdna_id for the refseq entries. 
my $update_string = "update protein set cdna_id=? where protein_id=?";
my $sth = $db->prepare($update_string) || $db->throw("can't prepare: $update_string");

open(REFSEQ, "<$refseqfile") or die "Can't open $refseqfile : $!";


while(<REFSEQ>){
#  next unless  /^>\w+\|\w+\|\w+\|(NP\w+)\|\s.\((NM\S+)\)\|/;
  next unless /^>/;
# we do NOT want NGs or NCs
  next unless /^>\w+\|\w+\|\w+\|(NP\S+)\|.*(NM\S+)\)\s+/;

  my $cdna_id    = $2;
  my $protein_id = $1;

#  print STDERR "$protein_id:$cdna_id\n";

  my $res = $sth->execute($cdna_id, $protein_id) || $db->throw("can't execute: $update_string with $cdna_id, $protein_id");
}

close REFSEQ;

sub check_parameters{
  # refseqfile first
  if(!defined $refseqfile){
    print STDERR "Checking for refseq file in GB_PROTEOME_FILES\n";
    # check through entries in GB_PROTEOME FILES to try to identify the refseq file. Give up if we don't find it.
    foreach my $file( @$GB_PROTEOME_FILES ){
      my $file_name = $file->{file_path};
      open (IN, "<$file_name") or die "Can't open $file_name : $!";
      while(<IN>){
	next unless /^>/;
	last unless /^>\w+\|\w+\|\w+\|NP\w+\|\s.\(NM\S+\)\|/;
	print STDERR "Found refseqfile: $file_name\n";
	$refseqfile = $file_name;
      }
      close IN;
    }
  }

  if (!defined $refseqfile){
    print STDERR "You need to define a refseq filename either on the command line or in Config::GeneBuild::Scripts::GB_PROTEOME_FILES\n";
    &usage;
    exit;
  }

  if(!-e $refseqfile){
    print STDERR "[$refseqfile] does not exist\n";
    &usage;
    exit;
  }

  # now check the db; dbpass and dbport can be undefined
  if(!defined $dbhost){
    $dbhost = $GB_DBHOST;
  }

  if(!defined $dbname){
    $dbname = $GB_DBNAME;
  }

  if(!defined $dbuser){
    $dbuser = $GB_DBUSER;
    if(!defined $dbpass){
      $dbpass = $GB_DBPASS;
    }
  }

  if(!defined $dbport){
    $dbport = $GB_DBPORT;
  }


  if(!defined $dbhost || !defined $dbname || !defined $dbuser ){
    print STDERR "Database details insufficient ($dbhost, $dbname, $dbuser)\n";
    &usage;
    exit;
  }
}

