#!/usr/local/bin/perl -w

# make sure we can find pmatch_modules.pm
# change this or run with perl -I ...
BEGIN {
#  unshift (@INC,"/work1/vac/TGW/pmatch_filtering/new_pmatch_filter");
  unshift (@INC,"/work2/vac/ensembl-pipeline/scripts");
}

=head1 NAME

  pmatch_filter.pl

=head1 SYNOPSIS
 
  pmatch_filter.pl

=head1 DESCRIPTION

  pmatch_filter.pl filters the output from pmatch when running human proteome vs
  genomic sequence - pmatch produces reams of info that can be usefully 
  combined and input to eg TargettedGenewise

  pmatch compares a DNA fasta file with a protein fasta file; this script 
  needs to be run twice, once with the output from pmatch piped into it, and 
  once taking its input from the file containing the matches output from the 
  first run.

  The first run gives the pmatch hits for each protein vs each DNA sequence
  split into their component fragments; the second run finds the "best" hit 
  for each protein ie the one with the greatest percentage coverage 
  of the protein sequence. We cannot just take the longest hit for each protein
  as we expect (at least) one hit per exon, not one hit per protein. Percentage coverage 
  is based on protein sequence length - so we need either the protein fasta file 
  which will be indexed by this script, or (better) an index for the protein 
  fasta file created by Bio::Index::Fasta

  For whole genome analysis we run pmatch on each fpc contig file individually and 
  then combine the results because of memory concerns. This might be unfounded but 
  it seems to work. 
  pmatch_filter.pl could easily be modified to just need running once - just a 
  case of modifying the call to print_MergedHit. TargettedGenewise currently can take
  either a single hit or a list of hit fragments as input.

  Typical usage:

  touch jobs.dat
  foreach i(/work2/gs2/data/humangenome/oct07/oo23/*/*/*)
  pmatch -D ./prot.fa $i | ./pmatch_filter.pl -indexfile=./prot.inx -first >> jobs.dat
  end

  then 

  ./pmatch_filter.pl -indexfile=./prot.inx -second < jobs.dat
  

=head1 OPTIONS

  User needs to provide either protfile or indexfile. If both are provided, indexfile 
  will be used
    -protfile  Name of the fasta file containing the protein sequences. Must be the 
               full path - relative path is not good enough!
    -indexfile Path to indexfile for protein sequences already indexed with Bio::Index::Fasta

  User needs to tell script whether this is the first or second round of processing
    -first     Input piped in from pmatch directly
    -second    Parsing output from -first run

=cut

use strict;
use pmatch_modules;
use Getopt::Long;
use Bio::Seq;
use Bio::Index::Fasta;
use Bio::SeqIO;

my $protfile;
my $indexfile;
my $first;
my $second;

&GetOptions( 
	    'protfile:s'     => \$protfile,
	    'indexfile:s'    => \$indexfile,
	    'first'          => \$first,
	    'second'         => \$second,
	   );

if(!defined($protfile) && !defined($indexfile)){
  print "You must provide a protein fasta file, or an indexfile\n";
  print "Usage: pmatch_filter -protfile=fasta_protein_file | -indexfile=fasta_protein_index -first|-second\n";
  exit;
}

if((defined($first) && defined($second)) || (!defined($first) && !defined($second))){
  print "Please specify whether this is the first or the second stage processing\n";
  print "Usage: pmatch_filter -protfile=fasta_protein_file | -indexfile=fasta_protein_index -first|-second\n";
  exit;
}

### MAIN ###
my %proteins;      # hash relating protein ID to ProteinHit objectX
my @merged;        # list of merged hits

LOOP: while (<>) { 
  next LOOP unless /\S+/;
  if(defined $first){ 
    make_coord_pair($_); 
  }
  
  elsif(defined $second){ 
    make_contig_hit($_); 
  }
  
  else { 
    die ("what do you want me to do?\n"); 
  }
}

# only proceed if we have some pmatch results to process
if(scalar(values(%proteins))){
  if(!defined $indexfile){
    # we need to make a temporary index from the protein fasta file
    $indexfile = "/tmp/pmatch_inx.$$";
    index_proteins();
    merge_hits();
    unlink $indexfile;
  }
  
  else {
    # we are given an indexfile and rely on the user to unlink it
    merge_hits();
  }
  
  # output the hits
  foreach my $mh(@merged) { print_MergedHit($mh,$first); }
}
### END MAIN ###

=head2 make_coord_pair

 Title   : make_coord_pair
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub make_coord_pair {
  my ($line) = @_;
  my @cols = split();
  # clean up the refseq IDs
  my $protid = $cols[5];
  # alter this as necessary 
  if($protid =~ /\w+\|\w+\|\w+\|([\w\.]+)\|/) {
#  if($protid =~ /\w+\|\w+\|\w+\|(\S+/) {
    $protid = $1;
  }
  
  # sort out strand
  my $strand = 1;
  if ( $cols[3] < $cols[2]  ) { $strand=-1; }
  
  my $cp = new CoordPair(
			 -query  => $cols[1],
			 -target => $protid,
			 -qstart => $cols[2],
			 -qend   => $cols[3],
			 -tstart => $cols[6],
			 -tend   => $cols[7],
			 -percent=> $cols[8],
			 -strand => $strand,
			);
  # where to put the CoordPair?
  # find the relevant ProteinHit, or make a new one
  my $ph = $proteins{$protid};
  
  if (!defined $ph) {
     #make a new one and add it into %proteins
    $ph = new ProteinHit(-id=>$protid);
    $proteins{$protid} = $ph;
  }
  
  # now find the relevant ContigHit, or make a new one
  my $ch = $ph->get_ContigHit($cols[1]);
  if (!defined $ch) {
     # make a new one and add it into $ph
    $ch = new ContigHit(-id => $cols[1]);
    $ph->add_ContigHit($ch);
  }
  
  # now add the CoordPair
  $ch->add_CoordPair($cp);
  
}

=head2 make_contig_hit

 Title   : make_contig_hit
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub make_contig_hit {
  my ($line) = @_;
  my @coord_pairs = split /\s+/, $line;

  # use the first element to make a ContigHit
  my $fpc;
  my $protid;
  my $first = $coord_pairs[0];
  if($first =~ /(\w+):\d+,\d+;(\S+):\d+,\d+/) {
    $fpc = $1;
    $protid = $2;
  }
  else {
    die "can't sort out this contighit: $first\n";
  }
  
  my $ph = $proteins{$protid};
  if (!defined $ph) {
    #make a new one and add it into %proteins
    $ph = new ProteinHit(-id=>$protid);
    $proteins{$protid} = $ph;
  }

  my $ch = new ContigHit(-id => $fpc);
  $ph->add_ContigHit($ch);

  foreach my $cp(@coord_pairs) {
   # eg ctg13631:3262,3023;O00361:1,80 
    if($cp =~ /(\w+):(\d+),(\d+);(\S+):(\d+),(\d+)/ ) {
      my $strand = 1;
      if ($2 > $3) { $strand = -1; }
      my $pair = new CoordPair(
			 -query  => $1,
			 -target => $4,
			 -qstart => $2,
			 -qend   => $3,
			 -tstart => $5,
			 -tend   => $6,
			 -percent=> 100,
			 -strand => $strand,
			);
      $ch->add_CoordPair($pair);
    }
    else {
      warn "can't process $cp\n";
    }
  }
  
}


=head2 print_CoordPair

 Title   : print_CoordPair
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub print_CoordPair {
  my ($cp) = @_;
  print "cp: $cp\n";
  print $cp->query()  . "\t" . 
        $cp->target() . "\t" . 
        $cp->qstart() . "\t" . 
	$cp->qend()   . "\t" . 
	$cp->tstart() . "\t" .
	$cp->tend()   . "\t" .  
	$cp->percent()   . "\n"; 
}

=head2 print_MergedHit

 Title   : print_MergedHit
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub print_MergedHit {
  my ($mh,$first) = @_;
  my @cps = $mh->each_CoordPair();
  if(defined $first) {
  # print the details of all the constituent coord pairs separated by white space 
    foreach my $cp(@cps) {
      print $cp->query()  . ":" .
	    $cp->qstart() . "," . 
	    $cp->qend()   . ";" .
	    $cp->target() . ":" .
	    $cp->tstart() . "," .
	    $cp->tend() . " ";
    }
  }
  else {
    # "second" analysis - just one all encompassing hit please!
    my $qstart = $cps[0]->qstart;
    my $qend = $cps[$#cps]->qend;
    my $tstart = $cps[0]->tstart;
    my $tend = $cps[$#cps]->tend;
    
    print $mh->query . ":$qstart,$qend:" . $mh->target . ":$tstart,$tend\n";
  }
  
  print "\n";
}

=head2 new_merged_hit

 Title   : new_merged_hit
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub new_merged_hit {
  my ($cp) = @_;
  my $coverage = $cp->tend - $cp->tstart + 1;
  my $mh = new MergedHit( -query    =>  $cp->query(),
			  -target   =>  $cp->target(),
			  -strand   =>  $cp->strand(),
			  -coverage =>  $coverage,
			);
  $mh->add_CoordPair($cp);
  return $mh;
}

=head2 merge_hits

 Title   : merge_hits
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub merge_hits {
  # merge the hits together.
  foreach my $p_hit(values %proteins) {
    my @allhits = make_mergelist($p_hit);
    my @chosen = prune(@allhits);
    print STDERR "\nNo hits good enough for " . $p_hit->id() . "\n"
      unless scalar(@chosen);
    push(@merged,@chosen);
  }
}


=head2 make_mergelist

 Title   : make_mergelist
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub make_mergelist {
  my ($ph) = @_;
  my @hits = ();

  foreach my $ch($ph->each_ContigHit()) {
    # forward strand
    my @cps = $ch->each_ForwardPair();
    # sort!
    @cps = sort {$a->qstart <=> $b->qstart} @cps;
    #reverse strand
    my @mps = $ch->each_ReversePair();
    @mps = sort {$b->qstart <=> $a->qstart} @mps;
    push (@cps,@mps);
    
    my $first = shift(@cps);
    my $mh = new_merged_hit($first);
    
    push(@hits,$mh);
      
    CP: 
    foreach my $cp(@cps) {
      # need to compare with the last entry in @merged
      my $prev = $hits[$#hits];
      my @prev_cps = $prev->each_CoordPair();

      # first the strand
      my $strand = 1;
      if ($cp->qend() < $cp->qstart()) { $strand = -1; }

      # ignore hits that don't extend the current hit
      next CP if( $cp->tend <= $prev_cps[$#prev_cps]->tend );

      # need a fudge factor - pmatch could overlap them by 1 or 2 ... or 3
      if( $strand == $prev->strand &&
	 ( ($cp->tstart >=  $prev_cps[$#prev_cps]->tend) ||
	 (($prev_cps[$#prev_cps]->tend -  $cp->tstart) <= 3)))
	{
	  #extend existing MergedHit
	  my $coverage = $cp->tend - $cp->tstart + 1;
	  $coverage += $prev->coverage();

	  # compensate for any overlap 
	  my $diff = $prev_cps[$#prev_cps]->tend -  $cp->tstart;
	  if($diff >=0) {
	    $diff++;
	    $coverage -= $diff;
	  }

	  $prev->coverage($coverage);

	  $prev->add_CoordPair($cp);
	}
      else {
	# make a new MergedHit
	my $mh = new_merged_hit($cp);
	push(@hits,$mh);	
      }
    }
  }

  # try to merge overlapping hits?

  my $protlen = read_protlength($ph->id);
  warn "No length for " . $ph->id . "\n" unless $protlen;
  return unless $protlen; 
  foreach my $hit(@hits) {
    my $percent = $hit->coverage;
    $percent *= 100;
    $percent /= $protlen;
    $percent=sprintf("%.1f",$percent);
    $hit->coverage($percent);
  }
  return @hits;
}

=head2 prune

 Title   : prune
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub prune {
  my (@all) = @_;
  my @chosen = ();

  # sort by descending order of coverage
  @all = sort {$b->coverage <=> $a->coverage} @all;
  
  my $first = shift(@all);
  push (@chosen,$first);
  # don't select any hits that have coverage less than that of the first hit, be it 100 or 99.9 or ...
  my $curr_pc = $first->coverage(); 

 PRUNE:
  foreach my $hit(@all) {
    last PRUNE if $hit->coverage() < $curr_pc;
    push (@chosen,$hit);
  }

  return @chosen;

}

=head2 index_proteins

 Title   : index_proteins
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub index_proteins{
  #set up the index 

  my $inx = Bio::Index::Fasta->new(
				   -filename=>$indexfile, 
				   -write_flag=>1);
  $inx->id_parser(\&get_accession);
  $inx->make_index($protfile);

}

=head2 read_protlength

 Title   : read_protlength
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub read_protlength{
  my ($id) = @_;

  my $inx  = Bio::Index::Fasta->new(
				    -filename   => $indexfile, 
				   );
  # really ought to catch seq fetching exceptions here ...
  my $seq = $inx->fetch($id);
  return $seq->length;
}

=head2 get_accession

 Title   : get_accession
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_accession{
  my $header = shift;
  my $acc = $header;

  chomp $acc;
  if($acc =~ /^>(\w+)\s+\S+/){
    $acc = $1;
    $acc =~ s/>//;
  }  
  elsif($acc =~ /^>(\w+)/){
    $acc = $1;
    $acc =~ s/>//;
  }

  elsif($acc =~ /^>\S+\|\S+\|\S+\|(\S+)\|/){
    $acc = $1;
  }
  else { $acc=""; }

  return $acc;
}
