#!/usr/local/ensembl/bin/perl

use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::SeqIO;
use Getopt::Long;

my $file;

my $dbhost = 'ecs2d';
my $dbuser = 'ensro';
my $dbname = 'homo_sapiens_core_9_30';
my $dbpass = undef;

my $dnadbhost = 'ecs2d';
my $dnadbuser = 'ensro';
my $dnadbname = 'homo_sapiens_core_9_30';
my $dnadbpass = undef;


my $genetype = 'ensembl';

#my $dbhost = ecs2d;
#my $dbuser = 'ensro';
#my $dbname;
#my $dbpass    = undef;

#my $dnadbhost;
#my $dnadbuser = 'ensro';
#my $dnadbname;
#my $dnadbpass = undef;

my $id_file;
my $gstable_id;
my $tstable_id;
my $g_id;
my $t_id;

&GetOptions(
	    'tstable_id:s' => \$tstable_id,
	    't_id:s' => \$t_id,
	    'tstable_id:s' => \$tstable_id,
	    'dbhost:s'        => \$dbhost,
	    'dbname:s'        => \$dbname,
	    'dbuser:s'        => \$dbuser,
	    'genetype:s'      => \$genetype,
            );

unless ( $t_id || $tstable_id){
  print STDERR "script to print out all the evidence info from a transcript. useful for sanity checks.\n";
 
  print STDERR "Usage: $0 -dbname -dbhost [ -genetype] -t_id -tstable_id\n";
  exit(0);
}

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $dbhost,
					    '-user'   => $dbuser,
					    '-dbname' => $dbname,
					    '-pass'   => $dbpass,
					   );



print "connected to $dbname : $dbhost\n";
my $path = $db->assembly_type;


#my $sa = $db->get_SliceAdaptor();


#open (OUT,">$file") or die("unable to open file $file");

my $seqio = Bio::SeqIO->new('-format' => 'Fasta' , -fh => \*STDOUT ) ;


if ( $t_id){

  my $tran = $db->get_TranscriptAdaptor->fetch_by_dbID($t_id);
  my @exons = @{$tran->get_all_Exons};
  print "internal_id\tcontig_id\tcontig_name\tstart-end\tphase\tend_ph\tstrand\tlength\thit_name\tstart-end\tlength\tscore\tperc_id\n";
  foreach my $exon (@exons){
    
    # if exon is sticky print each component
    if ( $exon->isa('Bio::EnsEMBL::StickyExon') ){
      
      foreach my $exon_c ( @{$exon->get_all_component_Exons} ){
	my $length = $exon_c->end - $exon_c->start +1;
	print "Exon: ".$exon_c->dbID."\t".$exon_c->contig->dbID."\t".$exon_c->contig->name."\t".
	  $exon_c->start."-".$exon_c->end."\t".$exon_c->phase."\t".
	    $exon_c->end_phase."\t".$exon_c->strand."\t".$length."\n";
	&print_evidence($exon_c);
	print "\n";
      }
    }
    else{
      my $length = $exon->end - $exon->start +1;
      print "Exon: ".$exon->dbID."\t".$exon->contig->dbID."\t".$exon->contig->name."\t".
	$exon->start."-".$exon->end."\t".$exon->phase."\t".
	  $exon->end_phase."\t".$exon->strand."\t".$length."\n";
      &print_evidence($exon);
      print "\n";
    }
  }
  
 # print STDERR "Transcript from slice:\n";
#  my @genes = $slice->get_all_Genes_by_Type($genetype,'evidence');
#  my $this_transcript;
# GENE:
#  foreach my $gene (@genes){
#    foreach my $transcript ( $gene->each_Transcript ){
#      if ( $transcript->dbID == $t_id ){
#	$this_transcript = $transcript;
#	last GENE;
#      }
#    }
#  }
#  unless ($this_transcript){
#    print STDERR "cannot find transcript $t_id\n";
#    exit(0);
#  }
#  my @exons = $this_transcript->get_all_Exons;
#  @exons = sort { $a->start <=> $b->start } @exons;
#  foreach my $exon (@exons){
#    print "Exon ".$exon->start."-".$exon->end." phase: ".$exon->phase." end_phase: ".$exon->end_phase." strand: ".$exon->strand."\n";
#  }


}
elsif ( $tstable_id){
  #print "\nFrom database: $tstable_id\n";
  #&check_transcript_without_translation_from_stable_id($db,$tstable_id);
  #print "\n";

  print  "From TranscriptAdaptor: $t_id\n";
  my $tran = $db->get_TranscriptAdaptor->fetch_by_stable_id($tstable_id);
  my @exons = @{$tran->get_all_Exons};
  print "internal_id\tcontig_id\tcontig_name\tstart-end\tphase\tend_ph\strand\tlength\thit_name\tstart-end\tlength\tscore\tperc_idn";
  foreach my $exon (@exons){
    
    # if exon is sticky print each component
    if ( $exon->isa('Bio::EnsEMBL::StickyExon') ){
      
      foreach my $exon_c ( @{$exon->get_all_component_Exons} ){
	my $length = $exon_c->end - $exon_c->start +1;
	print "Exon: ".$exon_c->dbID."\t".$exon_c->contig->dbID."\t".$exon_c->contig->name."\t".
	  $exon_c->start."-".$exon_c->end."\t".$exon_c->phase."\t".
	    $exon_c->end_phase."\t".$exon_c->strand."\t".$length."\n";
	&print_evidence($exon_c);
	print "\n";
      }
    }
    else{
      my $length = $exon->end - $exon->start +1;
      print "Exon: ".$exon->dbID."\t".$exon->contig->dbID."\t".$exon->contig->name."\t".
	$exon->start."-".$exon->end."\t".$exon->phase."\t".
	  $exon->end_phase."\t".$exon->strand."\t".$length."\n";
      &print_evidence($exon);
      print "\n";
    }
  }



}
else{
  print STDERR "No id entered\n";
}
  
############################################################


sub check_transcript{
  my $db = shift;
  my $t_id = shift;

  
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
	      ass.type = "$path" AND
	      tr.transcript_id = $t_id AND  
	      tr.translation_id=tl.translation_id
	      ORDER BY et.rank
	    );

  #print STDERR $q."\n";
  my $sth = $db->prepare($q) || $db->throw("can't prepare: $q");
  my $res = $sth->execute || $db->throw("can't execute: $q");
  
  print STDERR "exon_id\tstart\tend\tphase\tend_phase\tstrand\tchr\tlength\trank\ttransl\tsticky\n";
  while( my ($exon_id,$start,$end,$strand,$chr,$length,$phase,$end_phase,$rank,$transl,$sticky) = $sth->fetchrow_array) {
    print STDERR "$exon_id\t$start\t$end\t$phase\t$end_phase\t$strand\t$chr\t$length\t$rank\t$transl\t$sticky\n";
  }
}
  
############################################################

sub check_transcript_without_translation{
  my $db = shift;
  my $t_id = shift;

  
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
	      if(e.sticky_rank>1,(concat("sticky (rank = ",e.sticky_rank,")")),
		 ("")) as sticky
	      FROM  exon e, transcript tr, exon_transcript et, 
	      assembly ass, chromosome c
	      WHERE e.exon_id=et.exon_id AND 
	      et.transcript_id=tr.transcript_id AND 
	      ass.chromosome_id=c.chromosome_id AND                              
	      ass.contig_id=e.contig_id AND 
	      ass.type = "$path" AND
	      tr.transcript_id = $t_id  
	      ORDER BY et.rank
	    );

  #print STDERR $q."\n";
  my $sth = $db->prepare($q) || $db->throw("can't prepare: $q");
  my $res = $sth->execute || $db->throw("can't execute: $q");
  
  print STDERR "exon_id\tstart\tend\tphase\tend_phase\tstrand\tchr\tlength\trank\tsticky\n";
  while( my ($exon_id,$start,$end,$strand,$chr,$length,$phase,$end_phase,$rank,$transl,$sticky) = $sth->fetchrow_array) {
    print STDERR "$exon_id\t$start\t$end\t$phase\t$end_phase\t$strand\t$chr\t$length\t$rank\t$transl\t$sticky\n";
  }
}

############################################################

sub check_transcript_without_translation_from_stable_id{
  my $db = shift;
  my $stable_id = shift;

  
  my $q = qq( SELECT es.stable_id, 
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
	      if(e.sticky_rank>1,(concat("sticky (rank = ",e.sticky_rank,")")),
		 ("")) as sticky
	      FROM  exon e, transcript tr, exon_transcript et, 
	      assembly ass, chromosome c, transcript_stable_id ts, exon_stable_id es
	      WHERE e.exon_id        = et.exon_id       AND 
	            es.exon_id       = e.exon_id        AND
	            et.transcript_id = tr.transcript_id AND
	            ts.transcript_id = tr.transcript_id AND
	            ts.stable_id     = "$stable_id"     AND
	            ass.chromosome_id= c.chromosome_id  AND                              
	            ass.contig_id    = e.contig_id      AND 
	            ass.type         = "$path"          
	      ORDER BY et.rank
	    );

  #print STDERR $q."\n";
  my $sth = $db->prepare($q) || $db->throw("can't prepare: $q");
  my $res = $sth->execute || $db->throw("can't execute: $q");
  
  print STDERR "exon_id\tstart\tend\tphase\tend_phase\tstrand\tchr\tlength\trank\tsticky\n";
  while( my ($exon_id,$start,$end,$strand,$chr,$length,$phase,$end_phase,$rank,$transl,$sticky) = $sth->fetchrow_array) {
    print STDERR "$exon_id\t$start\t$end\t$phase\t$end_phase\t$strand\t$chr\t$length\t$rank\t$transl\t$sticky\n";
  }
}


############################################################

sub print_evidence{
  my $exon = shift;
  my @evidence = @{$exon->get_all_supporting_features};
  if ( @evidence ){
    foreach my $evi ( @evidence ){
      my $length = $evi->end - $evi->start + 1;
      my $hlength = $evi->hend - $evi->hstart + 1;
      print "Evidence ".$evi->dbID."\t".
	$evi->contig->dbID."\t".$evi->contig->name."\t".
	  $evi->start."-".$evi->end."\t".$evi->phase."\t".
	    $evi->end_phase."\t".$evi->strand."\t".$length."\t".
	      $evi->hseqname."\t".
		$evi->hstart."-".$evi->hend."\t".$hlength."\t".
		  $evi->score."\t".$evi->percent_id."\n";
    }
  }
  else{
    print  "No evidence\n";
  }
}
