#!/usr/local/bin/perl -w

BEGIN {
  # oooh this is not nice
  my $script_dir = $0;
  $script_dir =~ s/(\S+\/)\S+/$1/;
  unshift (@INC, $script_dir);
  require "GB_conf.pl";
}

use strict;
use pmatch_modules; 
use Bio::Seq;
use File::Find;
use Getopt::Long;

my %conf =  %::GB_conf; # configuration options

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

  push (@hits, $pmf1->output);
}


 foreach my $hit(@hits) {
   print OUT $hit->query  . ":" . 
     $hit->qstart . "," . 
     $hit->qend   . ":" . 
     $hit->target . ":" .
     $hit->tstart . "," .
     $hit->tend   . " " .
     $hit->coverage . "\n";
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




