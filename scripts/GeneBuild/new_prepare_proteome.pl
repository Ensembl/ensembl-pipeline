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

     GB_KILL_LIST   location of text file listing Swissprot IDs to ignore

=cut


use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Scripts;
use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Pmatch;




my @file_info = @$GB_PROTEOME_FILES;
my $protfile = $GB_PFASTA;
my $pmatch    = $GB_PMATCH;

my %kill_list;
%kill_list = %{&get_kill_list($GB_KILL_LIST)} if($GB_KILL_LIST);
my %ids;
foreach my $file(@file_info){
  my $file_name = $file->{file_path};
  open (IN, "<$file_name") or die "Can't open $file_name : $!";
  open (OUT, ">>$protfile") or die "Can't open $protfile : $!";
  my $killing = 0;
  while(<IN>){
    #print STDERR "matching ".$_." with ".$file->{header_regex}."\n";
    if(/$file->{header_regex}/){
      #print "have matched ".$_." with ".$file->{header_regex}."\n";
      $killing = 0;
      if($kill_list{$1}){
	$killing = 1;
      }
      if(!$ids{$1}){
	$ids{$1} = 1;
      }else{
	print STDERR "skipping ".$1." it has already appeared\n";
	$killing = 1;
      }
      if(!$killing){
	print OUT ">$1\n";
      }
    }else{
      if(!$killing){
	s/U/X/g;
	print OUT $_;
      }
    }
  }
  close IN;
  close OUT;
}


&test_protfile;


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
  open(PM, "$GB_PMATCH -D $protfile $tmpfile | ") or die "Can't run $GB_PMATCH\n";
  while(<PM>) {
    print $_;
  }
  close PM;

  # tidy up
  unlink $tmpfile;

}




sub get_kill_list {
  my ($kill_list) = @_;
  my %kill_list = ();
  open KILL_LIST_FH, $kill_list or die "can't open $kill_list";
  while (<KILL_LIST_FH>) {
    my @element = split;
    if (scalar(@element) == 0) {	# blank or empty line
      next;
    }
    $kill_list{$element[0]} = 1;
  }
  close KILL_LIST_FH or die "file error for $kill_list";
  return \%kill_list;
}



