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
use FileHandle;

use Bio::EnsEMBL::Pipeline::GeneConf qw (
					 GB_PFASTA
					 GB_PMATCH
					 GB_PM_OUTPUT
					 GB_FPCDIR
					 GB_TMPDIR  
                                        );
$| = 1;

# global vars
my %plengths; # stores lengths of all proteins in $protfile 

my @hits;     # stores the hits from pmatch runs
my $protfile = $GB_PFASTA;
my $outdir   = $GB_PM_OUTPUT;
my $fpcdir   = $GB_FPCDIR;
my $pmatch   = $GB_PMATCH;
my $tmpdir   = $GB_TMPDIR;
my $check    = 0;
my $outfile;
my $chrname;
my $run_pmatch = 1; # 1 to run pmatch, 0 if pmatch already run
my $pmatchfile; # can provide a pmatch results file - will be sorted by both protein id and chr_start
                # gets ignored if run_pmatch is set

&GetOptions( 
	    'check'         => \$check,
	    'chr:s'         => \$chrname,
	    'protfile:s'    => \$protfile,
	    'run_pmatch'    => \$run_pmatch,
	    'pmatchfile:s'  => \$pmatchfile,
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

$outfile = "$chrname.pm.out";
if (defined ($outdir) && $outdir ne '') {
  $outfile = $outdir . "/" . $outfile;
}

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

foreach my $fpcfile(@files){
  next if $fpcfile =~ /contigFa.zip/;
  
  print STDERR "processing $fpcfile\n";

  if($run_pmatch){
    # generate pmatchfile
    $pmatchfile = "/tmp/$chrname.pmatch";
    my $pmatch = $GB_PMATCH;

    # run pmatch -D
    system("$pmatch -D $protfile $fpcfile > $pmatchfile");
  }
  elsif(defined $pmatchfile){
    # use a different outfile
    $outfile = $outdir . "/" . $pmatchfile . ".pm.out";
  }
  elsif(!defined $pmatchfile){
    die "you need to either tell me to run_pmatch, or you need to give me a pmatch results file\n";
  }

  &process_pmatches($fpcfile);

  # tidy up files we have created
  unlink $pmatchfile unless !$run_pmatch;

}

### END MAIN
### SUBROUTINES

sub process_pmatches{
  my($fpcfile) = @_;
  # sort pmatch results file
  system("sort -k6,6 -k3,3n $pmatchfile > $pmatchfile.tmp");
  rename "$pmatchfile.tmp", $pmatchfile;

  open (OUT, ">>$outfile") or die "Can't open $outfile for output: $!\n";
  OUT->autoflush(1);

  open(PMATCHES, "<$pmatchfile") or die "can't open [$pmatchfile]\n";  
  my $prot_id;
  my $current_pmf = new First_PMF(
				  -plengths => \%plengths,
				  -protfile => $protfile,
				  -pmatch   => $pmatch,
				  -tmpdir   => $tmpdir,
				  -fpcfile  => $fpcfile,
				  -maxintronlen => 2500000,
				 );
  
 PMATCH:  
  while(<PMATCHES>){
    my @cols = split;
    # dump out line to file just in case this turns into a problem
    
    if(!defined $prot_id || $cols[5] ne $prot_id){
	# process current_pmf
	foreach my $hit($current_pmf->merge_hits) {
	  print  OUT $hit->query  . ":" . 
	         $hit->qstart . "," . 
	         $hit->qend   . ":" . 
	         $hit->target . ":" .
	         $hit->tstart . "," .
	         $hit->tend   . " " .
	         $hit->coverage . "\n";
	}

      # start a new PMF
      $current_pmf = new First_PMF(
				   -plengths => \%plengths,
				   -protfile => $protfile,
				   -pmatch   => $pmatch,
				   -tmpdir   => $tmpdir,
				   -fpcfile  => $fpcfile,
				   -maxintronlen => 2500000,
				  );
      $prot_id = $cols[5];
      $current_pmf->make_coord_pair($_);
    }
    else{
      # add this hit into current PMF
      $current_pmf->make_coord_pair($_);
    }
  }

  # make sure we at least try to proces the last one!
  foreach my $hit($current_pmf->merge_hits) {
    print  OUT $hit->query  . ":" . 
               $hit->qstart . "," . 
	       $hit->qend   . ":" . 
	       $hit->target . ":" .
	       $hit->tstart . "," .
	       $hit->tend   . " " .
	       $hit->coverage . "\n";
  }



  close (OUT) or die "Can't close $outfile: $!\n";
  close (PMATCHES) or die "can't close pmatchfile: $!\n";

}



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

