#!/usr/local/bin/perl -w
use strict;
use Bio::SeqIO;

=head1 NAME

  protein2cdna.pl

=head1 SYNOPSIS
 
  protein2cdna.pl

=head1 DESCRIPTION

  protein2cdna.pl finds the corresponding cDNA (if any) for each of the proteins in our
  human proteome set; this information is necessary for our UTR predictions.

  Currently the EMBL accessions are parsed from the gnp and swiss format protein files, though
  once these are fully indexed in efetch all we will need are the protein names and we grab
  the records from efetch.

  There are 2 stages:
     1. for each protein, identify the EMBL accessions associated with it.
     2. for each EMBL entry, check to see if it is a cDNA or a clone (various simple checks) 
        and if we have more than one EMBL link for a protein, keep the longest one.
 
  Stage 2 is the same for both refseq and swissprot proteins. Stage 1 differs in the parsing.

=head1 OPTIONS
  
  Options are to be set in GB_conf.pl
  The important ones for this script are:
     rf_gnp      location of refseq file in gnp format 
     sp_swiss    location of swissprot file in swiss format
     efetch      location of the efetch executable
     cdna_pairs  where to write the protein:cdna pairs file

  The need for rf_gnp and sp_swiss will be removed once efetch indexing includes 
  the full descriptions of all the proteins therein.

     eg.
	    'rf.gnp'      => '/work2/vac/GeneBuild/rf.gnp',
	    'sp.swiss'        => '/work2/vac/GeneBuild/sptr.swiss',
	    'efetch'      => '/usr/local/ensembl/bin/efetch.new', 
	    'cdna_pairs'  => '/work2/vac/GeneBuild/cpp.out',
  
=cut


BEGIN {
  # oooh this is not nice
  my $script_dir = $0;
  $script_dir =~ s/(\S+\/)\S+/$1/;
  unshift (@INC, $script_dir);
  require "GB_conf.pl";
}

# global variables
my %pairs;

&parse_refseq();
&parse_sptr();
&write_output();

### END MAIN

### SUBROUTINES

sub parse_refseq {
  my %rfe_pairs = identify_rf_embl_pairs();
  
  while( my ($prot, $value) = each %rfe_pairs){
    foreach my $emid(@$value){
      check_embl($emid, $prot);
    }
  }  
}

sub parse_sptr {
  my %spe_pairs = identify_sp_embl_pairs();
  
  while( my ($prot, $value) = each %spe_pairs){
    foreach my $emid(@$value){
      check_embl($emid, $prot);
    }
  }
}

sub identify_rf_embl_pairs {
  my $refseq   = $::GB_conf{'rf_gnp'};  

  my %proteins;
  $/ = "//\n"; # record separator

  open (IN, "<$refseq") or die "Can't open [$refseq]\n";

  my $pid;

  while(<IN>){
    if (/VERSION\s+(\S+)/){
      $pid = $1;
    }
    
    if(/was derived from(.*)\n/ ){
      my $l = $1;
      
      # strip commas
      $l =~ s/,//g;
      # strip leading whitespace
      $l =~ s/\s+//;
      # strip trailing .
      $l =~s /.$//;
      
      my @ids = split /\s+/, $l;
      foreach my $id(@ids) {
	# lose the version number, for now; efetch will barf
	$id =~ s/(\w+)\.\d+/$1/;
      }
      $proteins{$pid} = \@ids;
    }  
  }
  
  $/ = "\n"; # reset before all hell breaks loose
  close IN;
  return %proteins;
  
}

sub identify_sp_embl_pairs {
  my $swiss   = $::GB_conf{'sp_swiss'};  

  my %proteins;
  $/ = "//\n"; # record separator

  open (IN, "<$swiss") or die "Can't open [$swiss]\n";

  my $pid;
 RECORD: while(<IN>){
    my @l = split /\n/, $_;
    $pid = undef;
    foreach my $line(@l){
      # record may have >1 AC line
      if ($line =~/^AC\s+(\S+);/){
	if(defined ($pid)){
	  print STDERR "$pid already set; ignoring $1\n";
	}
	else{
	  $pid = $1;
	}

	if(!defined($proteins{$pid})){
	  $proteins{$pid} = ();
	}
	
      }

      elsif($line =~ /DR\s+EMBL;\s+(\w+)/ ){
	if (!defined ($pid)) { die "gargh: $line\n"; }
	push (@{$proteins{$pid}}, $1);
      }
    }
  }
    
  $/ = "\n"; # reset before all hell breaks loose
  close IN;
  return %proteins;
}

sub check_embl {
  my ($id, $prot) = @_;
  my $efetch  = $::GB_conf{'efetch'};

  my $seq;

  open(EFETCH, "$efetch $id |") || die("Error running $efetch for id [$id]");
  eval {
    my $fh = Bio::SeqIO->new(-fh   => \*EFETCH, "-format"=>'embl');
    $seq = $fh->next_seq();
  };
  close EFETCH;

  if($@) {
    print STDERR "problem with $id: [$@]\n";
    return;
  }

  my $validseq = 1;
  my $seen_cds = 0;
  open(EFETCH, "$efetch $id |") || die("Error running $efetch for id [$id]");
  while(<EFETCH>){
    last unless $validseq;
    next unless /FT\s+CDS/;
    if ($seen_cds){
      # cDNAs won't have > 1 CDS in them ...
      $validseq = 0;
      last;
    }
    $seen_cds = 1;
    next if (/join/);
  }

  # is $seq the longest cDNA so far for $prot?
  if($validseq) {
    $pairs{$prot} = $seq unless defined($pairs{$prot});
    if(($pairs{$prot}->length) < ($seq->length)) { 
      $pairs{$prot} = $seq;
    }
  }

}

sub write_output {
  my $outfile = $::GB_conf{'cdna_pairs'};
  open (OUT, ">$outfile") or die "Can't open [$outfile]\n";
  foreach my $prot(keys %pairs) {
     print OUT "$prot : " . $pairs{$prot}->accession . "\n";
  }
  close OUT;
}
