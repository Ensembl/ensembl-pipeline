=head1 NAME

ExonPair

=head1 SYNOPSIS

my $exon_pair = Bio::EnsEMBL::Pipeline::GeneComparison::ExonPair->new();

my $score = $exon_pair->blast_Exons( $exon1, $exon2 );

$gene_pair->compare_isoforms();

=head1 DESCRIPTION

Class to operate on potential orthologous exons.

=head1 CONTACT

eae@sanger.ac.uk

=cut

# Let the code begin ...

package Bio::EnsEMBL::Pipeline::GeneComparison::ExonPair;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::Runnable::Blast;
use Bio::SeqIO;

@ISA = qw(Bio::EnsEMBL::Root);

=head1 METHODS

=cut

#########################################################################

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);
  
    my( $gap_penalty ) = $self->_rearrange([qw(
					       GAP_PENALTY
					       )], 
					   @args);
 
    if ( $gap_penalty ){
	$self->gap_penalty( $gap_penalty );
    }
    else{
	$self->gap_penalty( -100 );
    }
    
    return $self;
}

############################################################


sub blast_Exons{
  my ($self,$exon1, $exon2) =@_;
  
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
  
  my $min_length = $length1;
  if ( $length2 < $length1 ){
    $min_length = $length2;
  }
  my $word = 5;
  #if ( $word*3 >= $min_length ){
  #  $word = int($min_length/3) - 1;
  #}
  if ( $word >= $min_length ){
      $word = $min_length - 1;
  }
  if ( $word < 1 ){
    return 0;
  }

  #print STDERR "word: $word\n";
  #print STDERR "comparing $id1:\t".$seq1->seq."\n";
  #print STDERR "      and $id2:\t".$seq2->seq."\n";


  ############################################################
  # create database
  my $file = 'seq_'.$$.'.fa';
  my $database = "/tmp/".$file;
  open( DB_SEQ,">$database") || die("Could not open $database $!");
  
  my $seqout = Bio::SeqIO->new('-format' => 'Fasta',
			       '-fh'     => \*DB_SEQ);
  
  $seqout->write_seq($seq2);
  close( DB_SEQ );
  
  system("pressdb $database > /dev/null 2>&1");
  
  ############################################################
  # Ian's parameters:
  #my $options = "W=$word M=1 N=-1 Q=3 R=3 S2=8"; 
  #my $options = "W=5";
  my $options = "W=$word -warnings";
  
  # tblastx options:
  #my $options = 'altscore="* any na" altscore="any * na" S2=12';
  #$options .= " V=200 B=200 ";
  #$options .= " -nogap ";
  #$options .= " W=$word ";
  
  #print STDERR "options: $options\n";
  
  #my $options = 'V=200 B=200 altscore="* any na" altscore="any * na" W=4 E=0.01 E2=0.01 -nogap';
  #my $options = 'V=200 B=200 W=9 E=0.01 E2=0.01';
  my $blast =  
    Bio::EnsEMBL::Pipeline::Runnable::Blast->new ('-query'          => $seq1,
						  #'-program'        => 'wutblastx',
						  '-program'        => 'wublastn',
						  '-database'       => $database,
						  '-threshold_type' => "PVALUE",
						  #'-threshold'      => 1e-10,
						  '-options'        => $options,
						 );
  
  
  $blast->add_regex($file,'(\S+)');

  eval{
    $blast->run();
  };
  if ( $@ ){
    print STDERR "Blast failed for\n";
    print STDERR "exon1: ".$seq1->seq."\n";
    print STDERR "exon2: ".$seq2->seq."\n";
    return (0);
  }    
  
  unlink( $database );
  
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
      my $score = 0;
      foreach my $fp ( @featurepairs ){
	  $score += $fp->score;
      }
      $score = $score/scalar(@featurepairs);
      return ($score,\@featurepairs);
      
      #my @feat_by_score = sort { $b->score <=> $a->score } @featurepairs;
      #return $feat_by_score[0]->score;
  }
  else{
      #print STDERR "no hits for:\n";
      #print STDERR "$id1:\t".$seq1->seq."\n";
      #print STDERR "$id2:\t".$seq2->seq."\n";
      return (0);
  }
}
  
############################################################

sub gap_penalty{
    my ($self,$value) = @_;
    if ( $value ){
	$self->{_gap_penalty} = $value;
    }
    return $self->{_gap_penalty};
}

############################################################

1;
