#!/usr/local/bin/perl -w

=head1 NAME

  get_human_ests.pl

=head1 SYNOPSIS
 
  get_human_ests.pl

=head1 DESCRIPTION

  gets human ESTs from dbEST or cDNAs from embl/vertRNA and polyT/polyA clips them

=head1 OPTIONS

  -estfile
  -outfile

=cut


use strict; 
use Getopt::Long;
use Bio::Seq;
use Bio::SeqIO;
use Bio::EnsEMBL::Utils::PolyA;

#$| = 1; # disable buffering
#local $/ = '>';

my $estfile;
my $seqoutfile;
my $clip;
my $softmask;

############################################################
# we usually clip 20bp on either end of the EST to eliminate low quality sequence
my $clip_ends = 0;

############################################################
# we don't want any EST which is shorter than 60bp
my $min_length = 60;


&GetOptions( 
	    'estfile:s'     => \$estfile,
	    'outfile:s'     => \$seqoutfile,
	    'clip'          => \$clip,
	    'clip_ends:n'   => \$clip_ends,
	    'softmask'      => \$softmask,
	    'min_length'    => \$min_length,
	   );

# usage
if(!defined $estfile    ||
   !defined $seqoutfile 
  ){
  print STDERR "script to collect human ESTs.\n";
  print STDERR "It rejects ESTs labelled as pseudogenes, non-coding RNAs or cancer genes\n";
  print STDERR "\n";
  print STDERR "USAGE: get_human_ests.pl -estfile estfile -outfile outfile\n";
  print STDERR "                         -clip (clips polyA/T) -clip_ends n (clips n bases from both ends)\n";
  print STDERR "                                                             default = 0, 20 seems to be ok\n";
  print STDERR "                         -softmask ( softmask the polyA/T )\n";
  print STDERR "                         -min_length ( min_est_length )\n";
  exit(0);
}


my $seqin  = new Bio::SeqIO(-file   => "<$estfile",
			    -format => "Fasta",
			  );

my $seqout = new Bio::SeqIO(-file   => ">$seqoutfile", 
			    -format => "Fasta"
			   );

if ( $clip_ends ){
  print STDERR "clipping $clip_ends from both ends of ESTs\n";
}

SEQFETCH:
while( my $cdna = $seqin->next_seq ){

  next unless $cdna->length > $min_length;
  
  my $display_id  = $cdna->display_id;
  my $description = $cdna->desc;

  # First select the species:
  next SEQFETCH unless (   $description =~ /Homo sapiens/
			   || $description =~ /DNA.*coding.*human/
		       );
  
  if(  $description =~ /similar to/ || $description =~ /homolog/i ){
    
    next SEQFETCH unless ( $description =~ /Homo sapiens.*similar to/ 
			   || $description =~ /Homo sapiens.*homolog/i 
			 );
  }

  ############################################################
  # reject pseudogenes
  if ( $description =~ /pseudogene/i ){
    print STDERR "rejecting potential  pseudogene: $description\n";
    next SEQFETCH;
  }

  ############################################################
  # reject non-coding RNAs
  if ( $description =~/tRNA/i 
       && 
       !( $description =~/synthetase/i
	  ||
	  $description =~/protein/i
	)
     ){
    print STDERR "rejecting potential rRNA: $description\n";
    next SEQFETCH;
  }
  
  ############################################################
  # reject cancer ESTs
  if ( $description =~/similar to/ ){
    
    next SEQFETCH if ( $description =~/carcinoma.*similar to/i 
		       ||
		       $description =~/cancer.*similar to/i
		       ||
		       $description =~/tumor.*similar to/i
		     );
  }
  else{
    next SEQFETCH if ( $description =~/carcinoma/i 
		       ||
		       $description =~/cancer/i
		       ||
		       $description =~/tumor/i
		     );
  }

  #print STDERR "keeping $description\n";

  ############################################################
  # parse description to get the id

  # GenBank
  if ( $display_id =~/gi\|\S+\|\S+\|(\S+\.\d+)\|/ || $description =~/gi\|\S+\|\S+\|(\S+\.\d+)\|/ ){
    $display_id = $1;
  }
  # EMBL vert-RNA
  else{
    my @labels = split /\s+/, $description;
    $display_id = $labels[0];
  }
  
  $cdna->display_id($display_id);
  $cdna->desc("");
  
  if($@){
    warn("can't parse sequence for [$description]:\n$@\n");
    next SEQFETCH;
  }
  
  ############################################################
  # clipping? 
  my $polyA_clipper = Bio::EnsEMBL::Utils::PolyA->new();
  my $new_cdna;
  
  if ( $clip_ends ){
    my $seq = $cdna->seq;
    my $seq_length = length( $seq );
    
    # skip it if you are going to clip more than the actual length of the EST
    if ( $clip_ends > 2*$seq_length ){
      next SEQFETCH;
    }
    my $new_seq = substr( $seq, $clip_ends, $seq_length - 2*$clip_ends );
    
    # skip it if you are left with an EST of less than 100bp
    if ( length( $new_seq ) < 100 ){
      next SEQFETCH;
    }
    $new_cdna = new Bio::Seq;
    $new_cdna->display_id( $cdna->display_id );
    $new_cdna->seq($new_seq);
  }
  else{ 
    $new_cdna = $cdna;
  }
  
  my $new_new_cdna;
  if ($clip){
    #print STDERR "going to pass ".$new_cdna->display_id."\n";
    $new_new_cdna = $polyA_clipper->clip($new_cdna);
  }
  elsif( $softmask ){
    $new_new_cdna = $polyA_clipper->mask($new_cdna, 'soft');
  }
  else{
    $new_new_cdna = $new_cdna;
  }
  
  unless( $new_new_cdna ){
    next SEQFETCH;
  }

  # skip it if you are left with an EST of less than 100bp
  if ( length( $new_new_cdna->seq ) < 100 ){
    next SEQFETCH;
  }


  # write sequence
  $seqout->write_seq($new_new_cdna);
}

