#!/usr/local/ensembl/bin/perl -w

# dump_seq_into_fastA.pl
# it reads a bit of sequence and dump it into a fasA file, to eb able to view it in Apollo

use strict;
use diagnostics;

use Bio::SeqIO;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Runnable::GenewiseHmm;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Bio::EnsEMBL::PepDnaAlignFeature;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Exon;
use Getopt::Long;

# get a contig with a piece-of/entire  chromosome

my $input_id;
my $hmm;
my $check = 0;

### rat ###
#my $dbhost    = 'ecs2f';
#my $dbuser    = 'ensro';
#my $dbname    = 'genewisedb_rat';
#my $dnadbname = 'rat_Nov02';
#my $dnadbhost = 'ecs2a';

### mouse ###
my $dbhost    = 'ecs2f';
my $dbuser    = 'ensro';
my $dbname    = 'genewisedb_mouse';
my $dnadbname = 'mus_musculus_core_12_3';
my $dnadbhost = 'ecs2d';


&GetOptions(
	    'input_id:s'     => \$input_id,
	    'hmm:s'          => \$hmm,
	    'check'          => \$check,
	   );

if ( $check ){
  exit(0);
}

unless( $input_id && $hmm ){
  print STDERR "Usage: $0 -input_id -hmm\n";
  exit(0);
}

# connect to the database
my $dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					       -host  => $dnadbhost,
					       -user  => $dbuser,
					       -dbname=> $dnadbname
					      );

my $db =  new Bio::EnsEMBL::DBSQL::DBAdaptor(
					     -host  => $dbhost,
					     -user  => 'ensadmin',
					     -pass  => 'ensembl',
					     -dbname=> $dbname,
					     -dnadb => $dnadb,
					    );

my ($chr,$start,$end);
if ($input_id =~ /(\S+)\.(\d+)-(\d+)/){
  $chr   = $1;
  $start = $2;
  $end   = $3;
}

# a Slice object:
my $slice   = $db->get_SliceAdaptor->fetch_by_chr_start_end($chr,$start,$end);


#a Bio::PrimarySeq:
my $genomic = $slice->get_repeatmasked_seq(['RepeatMask']);



my $analysis = $db->get_AnalysisAdaptor->fetch_by_logic_name('rat_HMM');



my $genewise = Bio::EnsEMBL::Pipeline::Runnable::GenewiseHmm->new(-genomic  => $genomic,
								  -hmmfile  => $hmm,
								  #-memory   => $memory,
								 );


$genewise->run;

my @genes = $genewise->output;
print STDERR scalar(@genes)." genes found\n";

my $count = 1;

my @final_genes;

foreach my $feature (@genes){
  print STDERR "gene $count\n";
  $count++;
  my $gene = Bio::EnsEMBL::Gene->new();
  my $transcript  = Bio::EnsEMBL::Transcript->new();
  my $translation = Bio::EnsEMBL::Translation->new();
  $gene->add_Transcript($transcript);
  $transcript->translation( $translation );
  
  my @features = $feature->sub_SeqFeature;
  print STDERR scalar(@features)." exons found\n";
  next unless ( @features );
  if ($features[0]->strand == +1 ){
    @features = sort { $a->start <=> $b->start } @features;
  }
  else{
    @features = sort { $b->start <=> $a->start } @features;
  }
  
  my $exoncount = 1;
  my $previous_exon;
  foreach my $sub_feature ( @features ){
    my $exon = Bio::EnsEMBL::Exon->new();
    $transcript->add_Exon($exon);
    #print STDERR $sub_feature->gffstring."\n";
    $exon->start  ($sub_feature->start);
    $exon->end    ($sub_feature->end);
    $exon->strand ($sub_feature->strand);
    $exon->contig($slice);
    $exon->phase( $sub_feature->phase );
    $exon->end_phase( ($exon->phase + $exon->end - $exon->start + 1)%3 );
    if ( $exoncount == 1 ){
      #$exon->phase(0);
      $translation->start_Exon($exon);
      $translation->start(1);      
    }
    if ( $exoncount == scalar( @features ) ){
      $translation->end_Exon($exon);
      $translation->end( $exon->end - $exon->start + 1);      
    }
    #if ($previous_exon ){
    #  $exon->phase( $previous_exon->end_phase );
    #}
    $previous_exon = $exon;
    $exoncount++;
    
    if ( $sub_feature->sub_SeqFeature ){
      my @supp_features = $sub_feature->sub_SeqFeature;
      my $supp_feature;
      #eval{
      #	$supp_feature = Bio::EnsEMBL::DnaPepAlignFeature->new( -features => \@supp_features);
      #};
      #if ( $@ || !defined $supp_feature ){
      #	warn("could not create supporting feature:\n$@");
      #      }
      foreach my $sup_feature ( @supp_features ){
	my @list = ( $sup_feature );
	$supp_feature = Bio::EnsEMBL::DnaPepAlignFeature->new( -features => \@list );
	$supp_feature->contig     ($exon->contig);
	#$supp_feature->seqname    ($sub_sub_feature->feature1->seqname);
	#$supp_feature->hseqname   ($sub_sub_feature->feature2->seqname);
	#$supp_feature->score      ($supp_features[0]->score);
	#$supp_feature->percent_id ($sub_sub_feature->feature2->percent_id);
	$exon->add_supporting_features($supp_feature);
	#print STDERR "supp_feature: ".$supp_feature->gffstring."\t".$supp_feature->score."\n";
      }
    }
    else{
      print STDERR "no supporting features\n";
    }
  }
  $gene->type('genewiseHMM');
  $gene->analysis($analysis);
  #foreach my $transcript ( @{$gene->get_all_Transcripts} ){
  #  Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Evidence( $transcript );
  #}
  my $transformed_gene;
  eval{
    $transformed_gene = $gene->transform;
  };
  if ($@){
    warn("Unable to map gene to rawcontig coordinates!!");
    if ( $@ =~/lies on a gap/ ){
      warn("Exon list on a gap");
    }
    else{
      print STDERR "$@\n";
    }
  }
  else{
    push(@final_genes,$transformed_gene);
  }
}


my $gene_adaptor = $db->get_GeneAdaptor;
foreach my $gene (@final_genes){
  eval{
    $gene_adaptor->store($gene);
  };
  if ($@){
    warn("Unable to store gene!! $@");
    foreach my $tran (@{$gene->get_all_Transcripts}){
      Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript( $tran );
    }
    print STDERR "Error message:\n$@";
  }
  else{
    print STDERR "stored gene ".$gene->type." ".$gene->dbID."\n";
    foreach my $transcript ( @{$gene->get_all_Transcripts} ){
      Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript( $transcript );
      #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_TranscriptEvidence(   $transcript );
    }
  }

}
