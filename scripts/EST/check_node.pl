#!/usr/local/bin/perl -w

=head1 NAME

  check_node.pl

=head1 SYNOPSIS
 
  check_node.pl
  Checks to make sure the estfile and input_id file (if it does not look like chr_name.start-end))
  can be found.
  Intended to be run as the -E command in a bsub

=head1 DESCRIPTION


=head1 OPTIONS

=cut

use strict;
use Bio::EnsEMBL::Pipeline::Config::cDNAs_ESTs::Exonerate qw (
							      EST_CHUNKDIR
							      EST_GENOMIC
							      
							     );

my $chunkname = $ARGV[0];

if(!defined $chunkname){
  print STDERR "Usage: check_node.pl chunkname\n";
  exit(1);
}

my $estfiledir = $EST_CHUNKDIR;
my $estfile = $estfiledir . "/" . $chunkname;

my $input_id   = $EST_GENOMIC;

# check to see if can find $estfile
if (! -e $estfile){
    print STDERR "can't find $estfile\n";
    exit(1);
}

# check to see if $input_id is likely to be a file, and if it is, can we find it?
if(!($input_id =~/^\S+\.\d+-\d+$/)){
  if (! -e $input_id){
    print STDERR "can't find $input_id\n";
    exit(1);
  }
}
exit (0);


