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
  with other scripts and should be set in GB_conf.pl

  pm_filter.pl should be run once for each chromosome. The results are written out to 
  chrname.pm.dat in the pm_output directory specified in GB_conf.pl

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
  
  Options are to be set in GB_conf.pl
  The important ones for this script are:
     pfasta      clean file of proteins in fasta format
     pmatch      location of the pmatch executable
     pm_output   directory to write filtered output files
     fpcdir      top of the hierarchy where the fpcctg fasta files are
     tmpdir      scratch area

     eg.
	    'pfasta'      => '/work2/vac/GeneBuild/human_proteome.fa',
	    'pmatch'      => '/work2/vac/rd-utils/pmatch',
	    'pm_output'   => '/work2/vac/GeneBuild/',
	    'fpcdir'      => '/work2/vac/data/humangenome',
	    'tmpdir'      => '/work5/scratch/vac',
  
  If pm_output directory is not provided, output file is written to current directory
    
=cut


use strict;
use pmatch_modules; 
use Bio::Seq;
use File::Find;
use Getopt::Long;
require "Bio/EnsEMBL/Pipeline/GB_conf.pl";

my %conf;

{
  # yes I know I only use it once. It isn't a typo!
  local $^W = 0;
  %conf =  %::scripts_conf; # configuration options
}

# global vars
my %plengths; # stores lengths of all proteins in $protfile 
my @hits;     # stores the hits from pmatch runs
my $protfile = $conf{'pfasta'};
my $outdir   = $conf{'pm_output'};
my $fpcdir      = $conf{'fpcdir'};
my $pmatch   = $conf{'pmatch'};
my $tmpdir   = $conf{'tmpdir'};
my $check    = 0;
my $outfile  = '';
my $chrname;

&GetOptions( 
	    'check' => \$check,
	    'chr:s' => \$chrname
	   );

if ($check) {
  exit(0);
}

# final check to make sure all parameters are set before we go ahead
if(!defined($protfile) || !defined($fpcdir) || !defined($pmatch) || 
   $protfile eq '' || $fpcdir eq '' || $pmatch eq ''){
  print "You must set pfasta, fpcdir and pmatch in the config file, GB_conf.pl\n";
  exit (0);
} 

if(!defined($chrname)) {
  print "You must specify a chromosome to be pmatched\n";
  exit (0);
} 

$outfile = "$chrname.pm.out";
if (defined ($outdir) && $outdir ne '') {
  $outfile = $outdir . "/" . $outfile;
}

open (OUT, ">$outfile") or die "Can't open $outfile for output: $!\n";

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
  print STDERR "making prot_length list\n";
  open (INFILE, "<$protfile") or die "Can't open $protfile\n";
  $/ = '>'; # looks for fields separated by > instead of "\n"
 SEQFETCH:
  while(<INFILE>) {
    my @lines = split /\n/;
    next SEQFETCH unless scalar(@lines) > 1;
    my $description = shift (@lines);
    $description =~ /(\S+)/; # find the first block of non-whitespace
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
