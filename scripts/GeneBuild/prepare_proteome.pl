#!/usr/local/bin/perl -w
use strict;

=head1 NAME

  prepare_proteome.pl

=head1 SYNOPSIS
 
  prepare_proteome.pl

=head1 DESCRIPTION

  prepare_proteome.pl prepares a fasta file of protein sequences from swissprot and refseq 
  input files (also in fasta format). This file is needed for pmatch comparisons and its 
  creation is the first part of the GeneBuild.

  The file has a description line consisting solely of the accession number after the leading >
  All U are replaced by X to prevent pmatch complaining.

  The final part of the script does a tiny pmatch test run to reveal any problems persisting in 
  the file that would prevent later full scale pmatches from running smoothly.

=head1 OPTIONS
  
  Options are to be set in GeneBuild config files
  The important ones for this script are:
     GB_REFSEQ      location of refseq file in fasta format
     GB_SPTR        location of swissprot file in fasta format
     GB_PFASTA      where to write the clean fasta file
     GB_PMATCH      location of the pmatch executable
     GB_KILL_LIST   location of text file listing Swissprot IDs to ignore

=cut

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Scripts qw (
                                                           GB_REFSEQ
                                                           GB_SPTR
                                                           GB_PFASTA
					                   GB_KILL_LIST
                                                           GB_PMATCH
                                                          );

my $refseq    = $GB_REFSEQ;
my $sptr      = $GB_SPTR;
my $protfile  = $GB_PFASTA;
my $kill_list = $GB_KILL_LIST;
my $pmatch    = $GB_PMATCH;

if( defined $refseq && -e $refseq ) { &parse_refseq; }
if( defined $sptr   && -e $sptr )   { &parse_sptr; }

&test_protfile;

### END MAIN

sub get_kill_list {
  my @kill_list = ();
  open KILL_LIST_FH, $kill_list or die "can't open $kill_list";
  while (<KILL_LIST_FH>) {
    my @element = split;
    if (scalar(@element) == 0) {	# blank or empty line
      next;
    }
    push @kill_list, $element[0];
  }
  close KILL_LIST_FH or die "file error for $kill_list";
  return \@kill_list;
}

sub parse_sptr {
print STDERR "here\n";

  my $kill_list = get_kill_list();
  open (IN, "<$sptr") or die "Can't open $sptr\n";
  open (OUT, ">>$protfile") or die "Can't open $protfile\n";
  
  my $killing = 0;	# nonzero if we're killing the current entry
  while(<IN>){
    # eg >143G_HUMAN (Q9UN99) 14-3-3 protein gamma
    if(/^>\S+\s+\((\S+)\)/){
      $killing = 0;
      foreach my $id_to_kill (@$kill_list) {
        if ($1 eq $id_to_kill){
          print STDERR "INFO: ignoring kill list entry $1\n";
	  $killing = 1;
	}
      }
      if(! $killing){
        if($1 eq 'P17013'){
  	  die("DYING: $sptr still contains P17013. \nThis will probably cause problems with pmatch.\nYou should REMOVE IT AND RERUN prepare_proteome!\n");
        }
        if($1 eq 'Q99784'){
	  die("DYING: $sptr still contains Q99784. \nThis will probably cause problems with pmatch.\nYou should REMOVE IT AND RERUN prepare_proteome!\n");
        }
        print OUT ">$1\n";
      }
    }
    else{	# not a header line
      if (! $killing) {
        print OUT $_;
      }
    }
  }
  
  close IN;
  close OUT;

}

sub parse_refseq {

  open (IN, "<$refseq") or die "Can't open $refseq\n";
  open (OUT, ">$protfile") or die "Can't open $protfile\n";

  while(<IN>){
    # eg >gi|4501893|ref|NP_001094.1| actinin, alpha 2 [Homo sapiens]
    if(/^>/){
      if(/^>\w+\|\w+\|\w+\|(\S+)\|/){
	print OUT ">$1\n";
      }
      else {
	print OUT $_;
      }
    }
    else {
      # sequence - sub U by X
      s/U/X/g;
      print OUT $_;
    }
  }
  close IN;
  close OUT;

}

sub test_protfile {

  # set up a temporary file
  my $time = time;
  chomp ($time);
  my $tmpfile = "cup.$$.$time.fa";
  open (SEQ, ">$tmpfile") or die "can't open $tmpfile\n";
  print SEQ ">test_seq\n";
  print SEQ 'cctgggctgcctggggaagcacccagggccagggagtgtgaccctgcaggctccacacaggactgccagaggcacac';
  close SEQ;

  # do a pmatch test run
  print "starting pmatch test ... \n";
  open(PM, "$pmatch -D $protfile $tmpfile | ") or die "Can't run $pmatch\n";
  while(<PM>) {
    print $_;
  }
  close PM;

  # tidy up
  unlink $tmpfile;

}
