#!/usr/local/bin/perl -w

BEGIN {
  # oooh this is not nice
  my $script_dir = $0;
  $script_dir =~ s/(\S+\/)\S+/$1/;
  unshift (@INC, $script_dir);
}

=head1 NAME

  pm_filter.pl

=head1 SYNOPSIS
 
  pm_filter.pl

=head1 DESCRIPTION

  pm_filter.pl runs and filters the output from pmatch when running human 
  proteome vs genomic sequence - pmatch produces reams of info that can be usefully 
  combined and input to eg TargettedGenewise

  pmatch compares a DNA fasta file with a protein fasta file; this script 
  needs to be given the location of the fasta file contain the human proteome, 
  the directory at the head of the tree containing all the fpc contig fasta 
  files and the chromosome on which to run. The first two options are common 
  with other scripts and should be set in GeneConf.pm

  pm_filter.pl should be run once for each chromosome. The results are written out to 
  chrname.pm.dat in the GB_PM_OUTPUT directory specified in GeneConf.pm

  The script analyses the whole chromosome by running pmatch for the protein data 
  against each fpc contig file in turn, and then combining the results to find 
  the "best" hit for each protein ie the one with the greatest percentage 
  coverage of the protein sequence, as long as the coverage exceeds the lower threshold 
  of 25%. It also returns any hits that have coverage within 2% of this best hit, again 
  as long as coverage >= 25%. We cannot just take the longest hit for each 
  protein as we expect (at least) one hit per exon, not one hit per protein. 
  Percentage coverage is based on protein sequence length.

  Typical usage:

  ./pm_filter.pl -chr chr1
  

=head1 OPTIONS
  
  Options are to be set in GeneConf.pm
  The important ones for this script are:
     GB_PFASTA      clean file of proteins in fasta format
     GB_PMATCH      location of the pmatch executable
     GB_PM_OUTPUT   directory to write filtered output files
     GB_FPCDIR      top of the hierarchy where the fpcctg fasta files are
     GB_TMPDIR      scratch area

     eg.
	    GB_PFASTA      => /work2/vac/GeneBuild/human_proteome.fa,
	    GB_PMATCH      => /work2/vac/rd-utils/pmatch,
	    GB_PM_OUTPUT   => /work2/vac/GeneBuild/,
	    GB_FPCDIR      => /work2/vac/data/humangenome,
	    GB_TMPDIR      => /work5/scratch/vac,
  
  If pm_output directory is not provided, output file is written to current directory
    
=cut


use strict;
use pmatch_modules; 
use Bio::Seq;
use File::Find;
use Getopt::Long;

use Bio::EnsEMBL::Pipeline::GeneConf qw (
                                          GB_PFASTA
	                                  GB_PMATCH
	                                  GB_PM_OUTPUT
	                                  GB_FPCDIR
	                                  GB_TMPDIR  
                                        );

# global vars
my %plengths; # stores lengths of all proteins in $protfile 
my @hits;     # stores the hits from pmatch runs
my $protfile = $GB_PFASTA;
my $outdir   = $GB_PM_OUTPUT;
my $fpcdir   = $GB_FPCDIR;
my $pmatch   = $GB_PMATCH;
my $tmpdir   = $GB_TMPDIR;
my $check    = 0;
my $outfile  = '';
my $chrname;

&GetOptions( 
	    'check'      => \$check,
	    'chr:s'      => \$chrname,
	    'protfile:s' => \$protfile
	   );

if ($check) {
  exit(0);
}

# final check to make sure all parameters are set before we go ahead
if(!defined($protfile) || !defined($fpcdir) || !defined($pmatch) || 
   $protfile eq '' || $fpcdir eq '' || $pmatch eq ''){
  print "You must set GB_PFASTA, GB_FPCDIR and GB_PMATCH in the config file, GeneConf.pm: $protfile : $fpcdir : $pmatch\n";
  exit (0);
} 

if(!defined($chrname)) {
  print "You must specify a chromosome to be pmatched\n";
  exit (0);
} 

#$outfile = "$chrname.$protfile.pm.out";
$outfile = "$chrname.pm.out";
if (defined ($outdir) && $outdir ne '') {
  $outfile = $outdir . "/" . $outfile;
}

# HACK
#$protfile = '/work2/vac//Mouse/Dec01/mouse_proteome/' . $protfile;
#open (OUT, ">$outfile") or die "Can't open $outfile for output: $!\n";
open (OUT, ">>$outfile") or die "Can't open $outfile for output: $!\n";

### MAIN

# get a listing of all the fpc contig files
$fpcdir = $fpcdir . "/" . $chrname;
my @dir = ($fpcdir);
my @files;
find sub { push(@files, $File::Find::name) unless -d  }, @dir;

foreach my $file(@files){
  print STDERR "$file\n";
}

# populate %plengths
   &make_protlist;
   
 # run pmatch on each fpc contig
 foreach my $fpcfile(@files){
   next if $fpcfile =~ /contigFa.zip/;

   print STDERR "processing $fpcfile\n";

  # now make the first round filterers  
  my $pmf1 = new First_PMF(-plengths => \%plengths,
			   -protfile => $protfile,
			   -pmatch => $pmatch,
			   -tmpdir => $tmpdir,
			   -fpcfile  => $fpcfile);
  
  $pmf1->run;
  foreach my $hit($pmf1->output) {
   print OUT $hit->query  . ":" . 
     $hit->qstart . "," . 
     $hit->qend   . ":" . 
     $hit->target . ":" .
     $hit->tstart . "," .
     $hit->tend   . " " .
     $hit->coverage . "\n";
 }
#  push (@hits, $pmf1->output);
}




close (OUT) or die "Can't close $outfile: $!\n";

### END MAIN

### SUBROUTINES

=head2 make_protlist

 Title   : make_protlist
 Usage   :
 Function: gets lengths for proteins in $protfile & stores them in %plengths. For some 
           reason Bio::SeqIO just returns me the first sequence ...
 Example :
 Returns : 
 Args    :


=cut

sub make_protlist{
  print STDERR "making prot_length list from $protfile\n";
  open (INFILE, "<$protfile") or die "Can't open $protfile\n";
  $/ = '>'; # looks for fields separated by > instead of "\n"
 SEQFETCH:
  while(<INFILE>) {
    my @lines = split /\n/;
    next SEQFETCH unless scalar(@lines) > 1;
    my $description = shift (@lines);
    $description =~ /(\S+)/; # find the first block of non-whitespace

    next SEQFETCH unless defined $description and $description ne '';
    my $bs = new Bio::Seq;
    $bs->display_id($1);
    
    my $seq;
    foreach my $line(@lines) {
      chomp $line;
      $seq .= $line;
    }
    
    $bs->seq($seq);
    
    $plengths{$bs->display_id} = $bs->length;
  }

  # note double quotes ...    
  $/ = "\n";
  close INFILE;
}
