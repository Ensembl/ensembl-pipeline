#!/usr/local/ensembl/bin/perl

use strict;  
use Bio::EnsEMBL::DBSQL::DBAdaptor;
#use Bio::EnsEMBL::Pipeline::GeneComparison::GeneComparison;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Runnable::Blast;
use Getopt::Long;


# human_db
my $human_dbname = 'homo_sapiens_core_13_31';
my $human_dbhost = 'ecs2f';
my $human = 'Homo sapiens';

# human_dnadb
my $human_dnadbname = 'homo_sapiens_core_13_31';
my $human_dnadbhost = 'ecs2f';

# mouse_db
my $mouse_dbname = 'mus_musculus_core_13_30';
my $mouse_dbhost = 'ecs2f';
my $mouse = 'Mus musculus';

# mouse_dnadb
my $mouse_dnadbname = 'mus_musculus_core_13_30';
my $mouse_dnadbhost = 'ecs2f';

# compara_db
my $compara_dbname = 'ensembl_compara_13_1';
my $compara_dbhost = 'ecs2f';


my $from_file;

# options
&GetOptions( 
	    'from_file'  => \$from_file,
	   );


my $human_dnadb;
my $mouse_dnadb;
my $mouse_db;
my $human_db;
my $compara_db;

my %gene_pair;

############################################################
# connect to the databases 
$human_dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $human_dnadbhost,
						  -user  => 'ensro',
						  -dbname=> $human_dnadbname,
						 );

$human_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $human_dbhost,
					       -user  => 'ensro',
					       -dbname=> $human_dbname,
					       -dnadb => $human_dnadb,
					      );


$mouse_dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $mouse_dnadbhost,
						  -user  => 'ensro',
						  -dbname=> $mouse_dnadbname,
						 );

$mouse_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $mouse_dbhost,
					       -user  => 'ensro',
					       -dbname=> $mouse_dbname,
					       -dnadb => $mouse_dnadb,
					      );


my $human_adaptor = $human_db->get_GeneAdaptor;
my $mouse_adaptor = $mouse_db->get_GeneAdaptor;



unless( $from_file ){

  print STDERR "Using compara_db\n";
  #my $compara_config =  '/nfs/acari/eae/ensembl/ensembl-compara/modules/Bio/EnsEMBL/Compara/Compara.conf';
  $compara_db = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(
							     -user      => 'ensro',
							     -dbname    => $compara_dbname,
							     -host      => $compara_dbhost,
							     #-conf_file => $compara_config,
							    );
  
  $compara_db->add_db_adaptor( $human_db);
  $compara_db->add_db_adaptor( $mouse_db);
  
  #my $mouse_db = $compara_db->get_db_adaptor($target_species,'NCBIM30');
  #my $human_db  = $compara_db->get_db_adaptor($focus_species ,'NCBI31');
  
  
  my $homol_adaptor = $compara_db->get_HomologyAdaptor;
  my @human_ids     = $homol_adaptor->list_stable_ids_from_species('Homo_sapiens');
  
  foreach my $human_id ( @human_ids ){
    
    my @mouse_homologs = 
      $homol_adaptor->fetch_homologues_of_gene_in_species('Homo_sapiens',$human_id,'Mus musculus');
    
    foreach my $homology ( @mouse_homologs ){
      
      print STDERR "comparing $human_id and ".$homology->stable_id."\n";
      &compare( $human_id, $homology->stable_id );
      push (@{$gene_pair{$human_id}},$homology->stable_id);
    }  
  }
}

if ( $from_file ){
  while(<>){
    chomp;
    my @entries = split;
    print STDERR "comparing $entries[0] and $entries[1]\n";
    &compare( $entries[0], $entries[1]);
    push( @{$gene_pair{$entries[0]}}, $entries[1] );
  }   
}


sub compare{
  my ($human_id,$mouse_id) = @_;
  
  my $human_gene        = $human_adaptor->fetch_by_stable_id($human_id);
  my @human_transcripts = @{$human_gene->get_all_Transcripts};

  my $mouse_gene = $mouse_adaptor->fetch_by_stable_id( $mouse_id);
  my @mouse_transcripts = @{$mouse_gene->get_all_Transcripts};
  
  
  ############################################################
  # we make a pair only if the transcripts align with gaps no longer than some threshold
  # e.g.
  #ENST00000318538	335-915	        581	1	ENSMUST00000047641	3-588	        586	1	score:1089	perc_id:71.4
  #ENST00000318538	726-1594	869	1	ENSMUST00000047641	270-1130	861	1	score:3035	perc_id:85.5
  #ENST00000318538	1536-1876	341	1	ENSMUST00000047641	1557-1900	344	1	score:1153	perc_id:82.9
  #ENST00000318538	1859-2820	962	1	ENSMUST00000047641	1977-2943	967	1	score:3330	perc_id:84.4
  # 
  # there is a big skip between line 2 and line 3, that could be an extra exon ( exon-slipping )

  foreach my $human_t ( @human_transcripts ){
    foreach my $mouse_t ( @mouse_transcripts ){
      print STDERR "blasting isoforms\n";
      &blast_isoforms( $human_t, $mouse_t );
      
      print STDERR "comparing exons:\n";
      my %score;
      &compare_Exons( $human_t, $mouse_t );
    }
  }
}

############################################################


sub blast_isoforms{
  my ( $tran1,$tran2 ) = @_;

  my $id1;
  if ( $tran1->dbID ){
    $id1 = $tran1->stable_id || $tran2->dbID;
  }
  else{
    $id1 = "no id";
  }
  
  my $id2;
  if ( $tran2->dbID ){
    $id2 = $tran2->stable_id || $tran2->dbID;
  }
  else{
    $id2 = "no id";
  }

  print STDERR "comparing $id1 and $id2\n";

  my $seq1    = $tran1->seq;
  my $length1 = $seq1->length;
  unless ( $seq1->display_id ){
    $seq1->display_id($id1);
  }

  my $seq2    = $tran2->seq;
  my $length2 = $seq2->length;
  unless ( $seq2->display_id ){
    $seq2->display_id($id2);
  }
  
  
  ############################################################
  # create database
  my $file = 'seq_'.$$.'.fa';
  my $database = "/tmp/".$file;
  open( DB_SEQ,">$database") || die("Could not open $database $!");
  
  my $seqout = Bio::SeqIO->new('-format' => 'Fasta',
			       '-fh'     => \*DB_SEQ);
  
  $seqout->write_seq($seq2);
  close( DB_SEQ );
  
  system("pressdb $database");#/tmp/db_seqs.$$");
  
  ############################################################
  my $blast =  Bio::EnsEMBL::Pipeline::Runnable::Blast->new ('-query'     => $seq1,
							     '-program'   => 'wublastn',
							     '-database'  => $database,
							     -threshold_type => "PVALUE",
							     '-threshold' => 1e-10,
							     #'-filter'    => $filter,
							     '-options'   => 'V=1000000'
							    );
  
  $blast->add_regex($file,'(\S+)');
  $blast->run();
  
  my @featurepairs = $blast->output();
  
  foreach my $fp (sort {$a->hstart <=> $b->hstart} @featurepairs) {
    &print_Feature($fp);
    #print $fp->gffstring . "\n";
  }
  print STDERR "$id1 length = $length1\n";
  print STDERR "$id2 length = $length2\n";
}



############################################################

sub print_Feature{
  my $f = shift;
  print 
    $f->seqname."\t".
      $f->start."-".$f->end."\t".
	  ($f->end - $f->start + 1)."\t".
	    $f->strand."\t".
	      $f->hseqname."\t".
		$f->hstart."-".$f->hend."\t".
		    ($f->hend - $f->hstart + 1 )."\t".
		      $f->strand."\t".
			"score:".$f->score."\t".
			  "perc_id:".$f->percent_id."\n";
}
	    




############################################################

sub compare_Exons{
  my ($human_t, $mouse_t ) = @_;
  
  my %score_matrix;
  foreach my $exon1 ( @{$human_t->get_all_Exons} ){
    foreach my $exon2 ( @{$mouse_t->get_all_Exons} ){
      my $id1;
      if ( $exon1->dbID ){
	$id1 = $exon1->stable_id || $exon2->dbID;
      }
      else{
	$id1 = $exon1;
      }
      
      my $id2;
      if ( $exon2->dbID ){
	$id2 = $exon2->stable_id || $exon2->dbID;
      }
      else{
	$id2 = $exon2;
      }
      
      $score_matrix{$id1}{$id2} = &blast_Exons($exon1,$exon2);
      $score_matrix{$id2}{$id1} = $score_matrix{$id1}{$id2};
      print "score_matrix[$id1][$id2] = ".$score_matrix{$id1}{$id2}."\n";
    }
  }
}


############################################################

sub blast_Exons{
  my ($exon1, $exon2) =@_;
  
  my $id1;
  if ( $exon1->dbID ){
    $id1 = $exon1->stable_id || $exon2->dbID;
  }
  else{
    $id1 = $exon1;
  }
  
  my $id2;
  if ( $exon2->dbID ){
    $id2 = $exon2->stable_id || $exon2->dbID;
  }
  else{
    $id2 = $exon2;
  }

  #print STDERR "comparing $id1 and $id2\n";

  my $seq1    = $exon1->seq;
  my $length1 = $seq1->length;
  unless ( $seq1->display_id ){
    $seq1->display_id($id1);
  }

  my $seq2    = $exon2->seq;
  my $length2 = $seq2->length;
  unless ( $seq2->display_id ){
    $seq2->display_id($id2);
  }
  
  ############################################################
  # create database
  my $file = 'seq_'.$$.'.fa';
  my $database = "/tmp/".$file;
  open( DB_SEQ,">$database") || die("Could not open $database $!");
  
  my $seqout = Bio::SeqIO->new('-format' => 'Fasta',
			       '-fh'     => \*DB_SEQ);
  
  $seqout->write_seq($seq2);
  close( DB_SEQ );
  
  system("pressdb $database");#/tmp/db_seqs.$$");


  ############################################################
  my $blast =  Bio::EnsEMBL::Pipeline::Runnable::Blast->new ('-query'     => $seq1,
							     '-program'   => 'wutblastx',
							     '-database'  => $database,
							     -threshold_type => "PVALUE",
							     '-threshold' => 1e-10,
							     #'-filter'    => $filter,
							     '-options'   => 'V=1000000'
							    );
  
  
  $blast->add_regex($file,'(\S+)');
  $blast->run();
  
  my @featurepairs = $blast->output();

  if ( @featurepairs ){
    #my @pos_strand = grep { $_->strand == 1} @featurepairs;  
    #my @neg_strand = grep { $_->strand == -1} @featurepairs;  
    #foreach my $fp (sort{ $a->hstart <=> $b->hstart} @pos_strand) {
    #  print $fp->gffstring . "\n";
    #}
    #foreach my $fp (sort{ $a->hstart <=> $b->hstart} @neg_strand) {
    #  print $fp->gffstring . "\n";
    #}
    
    my @feat_by_score = sort { $b->score <=> $a->score } @featurepairs;
    return $feat_by_score[0]->score;
  }
  else{
    return 0;
  }
}
  




































