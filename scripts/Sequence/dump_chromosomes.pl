#!/usr/local/bin/perl -w

=head1 NAME

  dump_chromosomes.pl

=head1 SYNOPSIS
 
=head1 DESCRIPTION

=head1 OPTIONS

  -repeatmask (off by default)
  -dbname
  -dbport
  -dbhost
  -dbuser
  -dbpass
  -outdir (the directory where the chromosome sequences will 
	   be written, one file per chromosome).
  -chr (specify a chromosome name to dump just one chromosome)

=cut

use strict;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::SeqIO;

my $dbname;
my $dbport;
my $dbhost;
my $dbuser;
my $dbpass;
my $repeatmask;
my $dust;
my $outdir;
my $chr;

&GetOptions('dbname:s'   => \$dbname,
	    'dbport:s'   => \$dbport,
	    'dbhost:s'   => \$dbhost,
	    'dbuser:s'   => \$dbuser,
	    'dbpass:s'   => \$dbpass,
	    'repeatmask' => \$repeatmask,
	    'dust'       => \$dust,
	    'outdir:s'   => \$outdir,
	    'chr:s'      => \$chr);

# Whack out some help if lost...
if(!defined $dbname || !defined $dbhost || !defined $dbuser || !defined $outdir){

  print "Usage: dump_chromosomes.pl -dbname [name] -dbhost [host] -dbuser [user]\n" .
    "-dbport [port] -dbpass [pass] -outdir [directory for chr files] -repeatmask\n" .
      "-dust -chr [chr name]\n";

  exit(0);
}


# Connect to db and stuff

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(-dbname => $dbname,
					     -user   => $dbuser,
					     -host   => $dbhost,
					     -port   => $dbport,
					     -pass   => $dbpass);

my $slice_adaptor = $db->get_SliceAdaptor;
my $chromosome_adaptor = $db->get_ChromosomeAdaptor;

# Figure out what kind of masking is wanted

my @mask_logicnames = ();
push (@mask_logicnames, 'RepeatMask') if $repeatmask;
push (@mask_logicnames, 'Dust')       if $dust;

# Loop through each chromosome in our database and dump each sequence 
# into its own eponymous file.

foreach my $chromosome (@{$chromosome_adaptor->fetch_all}){

  print "Dumping chromosome " . $chromosome->chr_name;

  my $outfile = $outdir . '/' . $chromosome->chr_name . '.fa';

  if (-e $outfile){
    print " ... output file exists already - skipping\n";
    next
  }

  if ($chr && $chr ne $chromosome->chr_name){
    print " ... skipping\n";
    next
  }

  print "\n";

  # Open a Bio::SeqIO output stream

  my $seqio = Bio::SeqIO->new(-file   => ">$outfile",
			      -format => 'Fasta');


  my $slice = $slice_adaptor->fetch_by_chr_name($chromosome->chr_name);

  unless ($repeatmask){
    $seqio->write_seq($slice);
    next
  }

  my $rm_seq = $slice->get_repeatmasked_seq(\@mask_logicnames,1);
  $rm_seq->display_id($slice->display_id . " @mask_logicnames");


  $seqio->write_seq($rm_seq);

}
