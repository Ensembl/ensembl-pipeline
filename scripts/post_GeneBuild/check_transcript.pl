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


my $path = 'NCBI_30';
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
	    'dbhost:s'        => \$dbhost,
	    'dbname:s'        => \$dbname,
	    'dnadbhost:s'     => \$dnadbhost,
	    'dnadbname:s'     => \$dnadbname,
	    'genetype:s'      => \$genetype,
            'path:s'          => \$path,
);

unless ( $t_id || $tstable_id){
  print STDERR "script to print out all the info from a transcript. useful for sanity checks.\n";
 
  print STDERR "Usage: $0 [-dbname -dbhost -dnadbname -dnadbhost -genetype] -t_id -tstable_id\n";
  exit(0);
}


my $dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					       '-host'   => $dnadbhost,
					       '-user'   => $dnadbuser,
					       '-dbname' => $dnadbname,
					       '-pass'   => $dnadbpass,
					      );


my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $dbhost,
					    '-user'   => $dbuser,
					    '-dbname' => $dbname,
					    '-pass'   => $dbpass,
					    '-dnadb'  => $dnadb,
					   );



print STDERR "connected to $dbname : $dbhost\n";
my $path = $db->assembly_type;

#my $sa = $db->get_SliceAdaptor();


#open (OUT,">$file") or die("unable to open file $file");

my $seqio = Bio::SeqIO->new('-format' => 'Fasta' , -fh => \*STDOUT ) ;


if ( $t_id){

  print STDERR "\nFrom database: $t_id\n";
  &check_transcript($db,$t_id);
  print STDERR "\n";


  #my $slice = $sa->fetch_by_transcript_id($t_id);
  

  #my @new_exons;
  #foreach my $exon ( @exons ){
  #  my $new_exon = $exon->transform($slice);
  #  push (@new_exons, $exon);
  #}
  #@new_exons = sort { $a->start <=> $b->start } @new_exons;

  print STDERR "From TranscriptAdaptor: $t_id\n";
  my $tran = $db->get_TranscriptAdaptor->fetch_by_dbID($t_id);
  my @exons = @{$tran->get_all_Exons};
  print "contig_id\tcontig_name\texon_id\tstart\tend\tphase\tend_phase\strand\tlength\n";
  foreach my $exon (@exons){
    my $length = $exon->end - $exon->start +1;
    print $exon->contig->dbID."\t".$exon->contig->name."\t".
      $exon->dbID."\t".$exon->start."\t".$exon->end."\t".$exon->phase."\t".
	$exon->end_phase."\t".$exon->strand."\t".$length."\n";
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


sub translate_from_gene_stable_id{
  my $g_id = shift;
  eval {
    my $gene = $db->get_GeneAdaptor->fetch_by_dbID($g_id);
    print STDERR "gene id $g_id, type: ".$gene->type."\n";
    foreach my $trans ( $gene->each_Transcript ) {
      # get out first exon. Tag it to clone and gene on this basis
           
      print STDERR "trying to get ->translate from transcript ".$trans->dbID."\n";
      my $tseq = $trans->translate();

      print STDERR "trying to get ->translation from transcript ".$trans->dbID."\n";
      my $translation = $trans->translation();

      #print STDERR "translation->translate is a : $tseq\n";
      #print STDERR "translation: ".$tseq->seq()."\n";
      if ( $tseq->seq =~ /\*/ ) {
	print STDERR "translation of transcript: ".$trans->dbID." in chr 20 has stop codons. Skipping! (in clone)\n";
      }
      $tseq->desc("Gene:$g_id transcript: ".$trans->dbID);
      my $result = $seqio->write_seq($tseq);
    }
  };
  
   
  if( $@ ) {
    print STDERR "unable to process $g_id, due to \n$@\n";
  }
}


close (OUT);


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
  
