#!/usr/local/ensembl/bin/perl -w

use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::SeqIO;
use Getopt::Long;

my $file;

my $dbhost = 'ecs2d';
my $dbuser = 'ensro';
my $dbname = 'homo_sapiens_core_9_30';
my $dbpass = undef;

my $genetype = 'ensembl';
my $path;

my $id_file;
my $gstable_id;
my $tstable_id;
my $g_id;
my $t_id;

&GetOptions(
	    'tstable_id:s' => \$tstable_id,
	    't_id:s' => \$t_id,
	    'dbhost:s'        => \$dbhost,
	    'dbname:s'        => \$dbname,
	    'genetype:s'      => \$genetype,
            'path:s'          => \$path,
);

unless ( $t_id || $tstable_id){
  print STDERR "script to print out all the info from a transcript. useful for sanity checks.\n";
 
  print STDERR "Usage: $0 -dbname -dbhost [-genetype -path] -t_id -tstable_id\n";
  exit(0);
}



my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $dbhost,
					    '-user'   => $dbuser,
					    '-dbname' => $dbname,
					    '-pass'   => $dbpass,
					   );



print STDERR "connected to $dbname : $dbhost\n";

unless( $path ){
  $path = $db->assembly_type;
}
print STDERR "path = $path\n";

my $seqio = Bio::SeqIO->new('-format' => 'Fasta' , -fh => \*STDOUT ) ;


if ( $t_id){

  print STDERR "\nFrom database: $t_id\n";
  &check_transcript($db,$t_id);
  print STDERR "\n";
  
  print STDERR "From TranscriptAdaptor: $t_id\n";
  my $tran = $db->get_TranscriptAdaptor->fetch_by_dbID($t_id);
  my @exons = @{$tran->get_all_Exons};
  print "contig_id\tcontig_name\texon_id\tstart\tend\tphase\tend_phase\tstrand\tlength\n";
  foreach my $exon (@exons){
    if ( $exon->isa('Bio::EnsEMBL::StickyExon') ){
      foreach my $exon_c ( @{$exon->get_all_component_Exons} ){
	my $length = $exon_c->end - $exon_c->start +1;
	print $exon_c->contig->dbID."\t".$exon_c->contig->name."\t".
	  $exon_c->dbID."\t".$exon_c->start."\t".$exon_c->end."\t".$exon_c->phase."\t".
	    $exon_c->end_phase."\t".$exon_c->strand."\t".$length."\n";
      }
    }
    else{
      my $length = $exon->end - $exon->start +1;
      print $exon->contig->dbID."\t".$exon->contig->name."\t".
	$exon->dbID."\t".$exon->start."\t".$exon->end."\t".$exon->phase."\t".
	  $exon->end_phase."\t".$exon->strand."\t".$length."\n";
    }
  }
}
elsif( $tstable_id ){
  
  print STDERR "\nFrom database: $tstable_id\n";
  &check_transcript_from_stable_id($db,$tstable_id);
  print STDERR "\n";
}
else{
  print STDERR "No id entered\n";
}



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

  my $sth = $db->prepare($q) || $db->throw("can't prepare: $q");
  my $res = $sth->execute || $db->throw("can't execute: $q");
  
  print STDERR "exon_id\tstart\tend\tphase\tend_phase\tstrand\tchr\tlength\trank\ttransl\tsticky\n";
  while( my ($exon_id,$start,$end,$strand,$chr,$length,$phase,$end_phase,$rank,$transl,$sticky) = $sth->fetchrow_array) {
    print STDERR "$exon_id\t$start\t$end\t$phase\t$end_phase\t$strand\t$chr\t$length\t$rank\t$transl\t$sticky\n";
  }
}
  

sub check_transcript_from_stable_id{
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
	      if(e.exon_id=tl.start_exon_id,
		 (concat(tl.seq_start," (start)",
			 if(e.exon_id=tl.end_exon_id,(concat(" ",tl.seq_end," (end)")),("")))),
		 if (e.exon_id=tl.end_exon_id,(concat(tl.seq_end," (end)")),("")))   
	      as transcoord,
	      if(e.sticky_rank>1,(concat("sticky (rank = ",e.sticky_rank,")")),
		 ("")) as sticky
	      FROM  translation tl, exon e, transcript tr, exon_transcript et, 
	      assembly ass, chromosome c, transcript_stable_id ts, exon_stable_id es
	      WHERE e.exon_id         = et.exon_id       AND 
	            es.exon_id        = e.exon_id        AND
	            et.transcript_id  = tr.transcript_id AND 
	            ts.transcript_id  = tr.transcript_id AND
	            ts.stable_id      = "$stable_id"     AND
	            ass.chromosome_id = c.chromosome_id  AND                              
	            ass.contig_id     = e.contig_id      AND 
	            ass.type          = "$path"          AND
	            tr.translation_id = tl.translation_id
	      ORDER BY et.rank
	    );

  my $sth = $db->prepare($q) || $db->throw("can't prepare: $q");
  my $res = $sth->execute || $db->throw("can't execute: $q");
  
  print STDERR "exon_id\tstart\tend\tphase\tend_phase\tstrand\tchr\tlength\trank\ttransl\tsticky\n";
  while( my ($exon_id,$start,$end,$strand,$chr,$length,$phase,$end_phase,$rank,$transl,$sticky) = $sth->fetchrow_array) {
    print STDERR "$exon_id\t$start\t$end\t$phase\t$end_phase\t$strand\t$chr\t$length\t$rank\t$transl\t$sticky\n";
  }
}
  
