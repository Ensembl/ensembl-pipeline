#!/usr/local/ensembl/bin/perl

=head1 NAME

   post_GeneBuild_checks.pl

=head1 DESCRIPTION

prototype for a suite of checks after the gene build. The checks are based on virtual
contigs of 5Mb (this size should be the same one as the one used during the genebuild). It checks so far:

1.- all exons in a gene are in the same strand
2.- it checks that phases are consistent between exons
3.- also checks for folded transcripts
4.- it flags also single exon genes and from these, the ones which are longer than 50000 bases


For mouse denormalised contigs: in order to know the size of the chromosomes, 
you must have a file with the internal_ids and raw_contig ids.
For standard contigs this should be changed. 

It ought to read the database parameters from the GeneBuild config
files but it all depends on whether your final genes are in the same
database as the one you put in the config files.

=head1 OPTIONS

These are to be set in the GeneBuild config files:
GeneBuild::Databases::GB_FINALDBHOST
GeneBuild::Databases::GB_FINALDBNAME
GeneBuild::Databases::GB_DBHOST
GeneBuild::Databases::GB_DBNAME

GeneBuild::GeneBuilder::GB_FINAL_GENETYPE

=cut

use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::SeqIO;
use Getopt::Long;
use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Databases qw (
							     GB_FINALDBHOST
							     GB_FINALDBNAME
							     GB_DBHOST
							     GB_DBNAME
							    );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::GeneBuilder qw (
							       GB_FINAL_GENETYPE
							      );




use Bio::EnsEMBL::Utils::Eprof('eprof_start','eprof_end','eprof_dump');

my $dbhost;
my $dbuser    = 'ensro';
my $dbname;
my $dbpass    = undef;

my $dnadbhost;
my $dnadbuser = 'ensro';
my $dnadbname;
my $dnadbpass = undef;

my $genetype = "ensembl"; # default genetype


$dbuser = "ensro";
&GetOptions(
	    'dbname:s'    => \$dbname,
	    'dbhost:s'    => \$dbhost,
	    'dnadbname:s' => \$dnadbname,
	    'dnadbhost:s' => \$dnadbhost,
	    'genetype:s'  => \$genetype,
);

unless ( $dbname && $dbhost ){
  print STDERR "script to check the sanity of genes and transcripts after a build\n";
  print STDERR "Usage: $0 -dbname -dbhost\n";
  exit(0);
}

#my $dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
#					       '-host'   => $dnadbhost,
#					       '-user'   => $dnadbuser,
#					       '-dbname' => $dnadbname,
#					       '-pass'   => $dnadbpass,
#					      );


my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $dbhost,
					    '-user'   => $dbuser,
					    '-dbname' => $dbname,
					    '-pass'   => $dbpass,
#					    '-dnadb'  => $dnadb,
					   );


print STDERR "connected to $dbname : $dbhost\n";

print STDERR "checking genes of type $genetype\n";


my @gene_ids = &get_gene_ids($db,$genetype);

GENE:
foreach my $gene_id ( @gene_ids){
  print STDERR "checking gene dbID: ".$gene_id."\n";
  
  my @transcript_ids = &get_transcript_ids($db,$gene_id);

 TRANSCRIPT:
  foreach my $tran_id ( @transcript_ids ){
   
    print STDERR "checking transcript dbID: ".$tran_id."\n";
    my ($exons,$info) = &check_transcript($db,$tran_id);
    #&print_transcript($exons,$info);
    
    my @exons = @$exons;
    my %info = %$info;
    
    my $strand;
    my $end_phase;
    my $previous_exon;
    my $exon_count = 0;
    
    # mark single exon transcripts
    if ( scalar(@exons) == 1){
      my $exon_id = shift @exons;
      my $length = $info{end}{$exon_id} - $info{start}{$exon_id} + 1;
      print STDERR "single exon transcript $tran_id of length: $length\n";      
    }

  EXON:
    foreach my $exon_id (@exons){
      $exon_count++;

      # is strand defined?
      if ( !defined $info{strand}{$exon_id} ){
	print STDERR "PROBLEM: exon $exon_id has no strand\n";
      }
      
      # is strand consistent among exons?
      unless ( $strand ){
	$strand = $info{strand}{$exon_id};
      }      
      if ( $info{strand}{$exon_id} != $strand ){
	print STDERR "PROBLEM: Problem with the strands in transcript $tran_id:\n";
	&print_transcript( $exons,$info);
      }
     
      # are phases consistent?
      if ($exon_count>1){
	unless ( $info{sticky}{$exon_id} == 1 && $info{sticky}{$previous_exon} == $info{sticky}{$exon_id} ){
	  if ( $end_phase != $info{start_phase}{$exon_id} ){
	    print STDERR "PROBLEM: Inconsistent phases in transcript $tran_id:\n";
	    &print_transcript( $exons,$info);
	  }
	}
      }
      $end_phase = $info{end_phase}{$exon_id};


      # folded transcripts?
      if ($exon_count>1){
	
	#unless both exons are sticky
	unless ( $info{sticky}{$exon_id} == 1 && $info{sticky}{$previous_exon} == $info{sticky}{$exon_id} ){
	  if ( $strand == 1){
	    if ( $info{start}{$exon_id} < $info{end}{$previous_exon} ){
	      print STDERR "PROBLEM: Transcript $tran_id folds back on itself\n";
	      &print_transcript( $exons,$info);
	    }
	  }
	  if ($strand == -1){
	    if  ( $info{end}{$exon_id} > $info{start}{$previous_exon} ){
	      print STDERR "PROBLEM: Transcript $tran_id folds back on itself\n";
	      &print_transcript( $exons,$info);
	    }
	  }
	}
      }
      $previous_exon = $exon_id;

    } # end of EXON
  }   # end of TRANSCRIPT
}     # end of GENE



########################################################


sub get_chrlengths{
  my $db   = shift;
  my $type = shift;

  if (!$db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor')) {
    die "get_chrlengths should be passed a Bio::EnsEMBL::DBSQL::DBAdaptor\n";
  }

  my %chrhash;

  my $q = qq( SELECT chr_name,max(chr_end) FROM static_golden_path as sgp
               WHERE sgp.type = '$type' GROUP BY chr_name
            );

  my $sth = $db->prepare($q) || $db->throw("can't prepare: $q");
  my $res = $sth->execute || $db->throw("can't execute: $q");

  while( my ($chr, $length) = $sth->fetchrow_array) {
    $chrhash{$chr} = $length;
  }
  return \%chrhash;
}


sub check_transcript{
  my $db = shift;
  my $t_id = shift;
  my $path = $db->assembly_type;

  my $q = qq( SELECT e.exon_id,  
	      if(ass.contig_ori=1,(e.contig_start-ass.contig_start+ass.chr_start),
		 (ass.chr_start+ass.contig_end-e.contig_end)) as start,
	      if(ass.contig_ori=1,(e.contig_end-ass.contig_start+ass.chr_start),
		 (ass.chr_start+ass.contig_end-e.contig_start)) as end,
	      if (ass.contig_ori=1,e.contig_strand,(-e.contig_strand)) as strand,
	      c.name,        
	      abs(e.contig_end-e.contig_start)+1 as length,
	      e.phase,
	      e.end_phase,
	      et.rank, 
	      if(e.exon_id=tl.start_exon_id,
		 (concat(tl.seq_start," (start)",
			 if(e.exon_id=tl.end_exon_id,(concat(" ",tl.seq_end," (end)")),("")))),
		 if (e.exon_id=tl.end_exon_id,(concat(tl.seq_end," (end)")),("")))   
	      as transcoord,
	      if(e.sticky_rank>1,(concat("sticky (rank = ",e.sticky_rank,")")),
		 ("")) as sticky
	      FROM  translation tl, exon e, transcript tr, exon_transcript et, 
	      assembly ass, chromosome c
	      WHERE e.exon_id=et.exon_id AND 
	      et.transcript_id=tr.transcript_id AND 
	      ass.chromosome_id=c.chromosome_id AND                              
	      ass.contig_id=e.contig_id AND 
	      ass.type = '$path' AND
	      tr.transcript_id = '$t_id' AND  
	      tr.translation_id=tl.translation_id
	      ORDER BY et.rank
	    );
  #print $q."\n";
  my $sth = $db->prepare($q) || $db->throw("can't prepare: $q");
  my $res = $sth->execute || $db->throw("can't execute: $q");
  
  my %exon;
  my @exons;
  
  while( my ($exon_id,$start,$end,$strand,$chr,$length,$phase,$end_phase,$rank,$transl,$sticky) = $sth->fetchrow_array) {
    #print STDERR "$exon_id\t$start\t$end\t$phase\t$end_phase\t$strand\t$chr\t$length\t$rank\t$transl\t$sticky\n";
    push (@exons,$exon_id);
    $exon{start}{$exon_id}       = $start;
    $exon{end}{$exon_id}         = $end;
    $exon{start_phase}{$exon_id} = $phase;
    $exon{end_phase}{$exon_id}   = $end_phase;
    $exon{strand}{$exon_id}      = $strand;
    $exon{chr}{$exon_id}         = $chr;
    if ($sticky){
      $exon{sticky}{$exon_id} = 1;
    }
    else{
      $exon{sticky}{$exon_id} = 0;
    }
  }
  return (\@exons,\%exon);
}

sub get_transcript_ids{
  my $db   = shift;
  my $gene_id = shift;

  if (!$db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor')) {
    die "get_chrlengths should be passed a Bio::EnsEMBL::DBSQL::DBAdaptor\n";
  }

  
  my $q = qq( SELECT transcript_id
	      FROM transcript
	      WHERE gene_id=$gene_id
            );

  my $sth = $db->prepare($q) || $db->throw("can't prepare: $q");
  my $res = $sth->execute || $db->throw("can't execute: $q");

  my @transcript_ids;
  while( my $id = $sth->fetchrow_array) {
    push(@transcript_ids,$id);
  }
  return @transcript_ids;
}


sub get_gene_ids{
  my $db   = shift;
  my $type = shift;

  if (!$db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor')) {
    die "get_chrlengths should be passed a Bio::EnsEMBL::DBSQL::DBAdaptor\n";
  }

  my %chrhash;

  my $q = qq( SELECT gene_id
	      FROM gene
	      WHERE type='$type'
            );

  my $sth = $db->prepare($q) || $db->throw("can't prepare: $q");
  my $res = $sth->execute || $db->throw("can't execute: $q");

  my @gene_ids;
  while( my $id = $sth->fetchrow_array) {
    push( @gene_ids, $id);
  }
  return @gene_ids;
}


sub print_transcript{
  my $exons = shift;
  my $info  = shift;

  print STDERR "printing $exons $info\n";

  my @exons = @$exons;
  print STDERR "exons: @exons\n";
  my %info = %$info;

  my $strand = $info{strand}{$exons[0]};
  if ($strand == 1){
    @exons = sort{ $info{start}{$a} <=> $info{start}{$b} } @exons;
  }
  elsif($strand == -1){
    @exons = sort{ $info{start}{$b} <=> $info{start}{$a} } @exons;
  }
  
  foreach my $exon_id ( @exons ){
    my $sticky_label ='';
    if ( $info{sticky}{$exon_id} == 1){
      $sticky_label = 'sticky';
    }
    
    print STDERR $exon_id." ".$info{start}{$exon_id}."-".$info{end}{$exon_id}.
      " phase: ".$info{start_phase}{$exon_id}.
	" end_phase: ".$info{end_phase}{$exon_id}.
	  " strand: ".$info{strand}{$exon_id}.
	    " chr: ".$info{chr}{$exon_id}.
	      " ".$sticky_label."\n";
  }
}
