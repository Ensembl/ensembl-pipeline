#!/usr/local/bin/perl -w
use strict;

open LS, 'ls |' or die "Can't open pipe from ls : $!";


FILES:
while(<LS>){
  chomp;
  my $file = $_;
  next FILES if ( $file =~/(\.)+\// );
  open(FILE, "<$file") or die "Can't open [$file]\n";
  
 THISFILE:
  while (<FILE>){
    chomp;
    my $line = $_;
    if ( $line =~ /glib/i     || 
	 $line =~ /ERROR/i    || 
	 $line =~ /exiting/i  ||
	 $line =~ /aborting/i || 
	 $line =~ /exit/i     || 
	 $line =~ /exception/i#|| 
	 #$line =~ /cannot/i 
       ){
      print $file."\n";
      print $line."\n\n";
      last THISFILE;
    }
  }
  close( FILE );
  #system("grep -i 'glib' $file > /work6a/eae.tmp/Mouse/jobs/error_file");
}
