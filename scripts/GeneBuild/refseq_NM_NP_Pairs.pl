#!/usr/local/bin/perl -w
use strict;
use Getopt::Long;

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Databases;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

# script to identify which NM goes with which NP and insert that info into
# the protein table of the pipeline DB
# The script assumes the early stages (pmatch, targetted) of the build have
# already been run ie the protein table already has protein_id information in it.
# options can either be passed in on the command line or read from config files.

# we now have to use the refseq gpff (full genbank) entry file as the fasta headers 
# have changed and no longer pair up NM and NP. Grrrrr.

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

my $cdna_id;
my $protein_id;

while(<REFSEQ>){
  next unless /^VERSION|DBSOURCE/;
#  print $_;
  if(/VERSION/){
    if(/(XP\S+)/){
      print STDERR "skipping [$1]\n";
    }

    next unless /(NP\S+)/;

    if(defined $protein_id){
      die("previous protein_id [$protein_id] has not been cleared out\n");
    }
    if(defined $cdna_id){
      die("previous cdna_id [$cdna_id] has not been cleared out ...\n");
    }

    $protein_id = $1;
  }

  if(/DBSOURCE/){

    # don't want NCs or NGs
    if(/(NC\_\S+)|(NG\_\S+)|(XM\S+)/){
      print STDERR "Skipping [$_] - matches NC or NG or XM\n";
      $cdna_id = undef;
      $protein_id = undef;
      next;
    }

    if(!defined $protein_id){
      die("something very wrong - no protein_id for $_\n");
    }
    if (defined $cdna_id){
      die("previous cdna_id [$cdna_id] has not been cleared out ...\n");
    }

    next unless /(NM\S+)/;

    $cdna_id = $1;

    print  "$protein_id\t$cdna_id\n";
    my $res = $sth->execute($cdna_id, $protein_id) || $db->throw("can't execute: $update_string with $cdna_id, $protein_id");

    $cdna_id    = undef;
    $protein_id = undef;
  }
}

close REFSEQ;

sub check_parameters{
  # refseqfile first
  if (!defined $refseqfile){
    print STDERR "You need to define a refseq gpff filename on the command line\n";
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

sub usage {
  print STDERR "USAGE: refseq_NM_NP_pairs -refseq refseq.protein.gpff -dbhost -dbport -dbname -dbuser -dbpass\n";
}
