package SimilarityInputIdConfig;

use vars qw(@ISA);
use strict;


sub new{
  my ($class,@args) = @_;
  my $self = bless {},$class;
  return $self;
}

sub PAF_LOGICNAMES{
  my ($self, $arg) = @_;
  if($arg){
    $self->{PAF_LOGICNAMES} = $arg;
  }
  return $self->{PAF_LOGICNAMES}
}

sub PAF_MIN_SCORE_THRESHOLD{
  my ($self, $arg) = @_;
  if($arg){
    $self->{PAF_MIN_SCORE_THRESHOLD} = $arg;
  }
  return $self->{PAF_MIN_SCORE_THRESHOLD}
}

sub PAF_UPPER_SCORE_THRESHOLD{
  my ($self, $arg) = @_;
  if($arg){
    $self->{PAF_UPPER_SCORE_THRESHOLD} = $arg;
  }
  return $self->{PAF_UPPER_SCORE_THRESHOLD}
}



sub PAF_SOURCE_DB{
  my ($self, $arg) = @_;
  if($arg){
    $self->{PAF_SOURCE_DB} = $arg;
  }
  return $self->{PAF_SOURCE_DB}
}

sub GENE_SOURCE_DB{
  my ($self, $arg) = @_;
  if($arg){
    $self->{GENE_SOURCE_DB} = $arg;
  }
  return $self->{GENE_SOURCE_DB}
}




sub OUTPUT_DB{
  my ($self, $arg) = @_;
  if($arg){
    $self->{OUTPUT_DB} = $arg;
  }
  return $self->{OUTPUT_DB}
}

sub OUTPUT_BIOTYPE{
  my ($self, $arg) = @_;
  if($arg){
    $self->{OUTPUT_BIOTYPE} = $arg;
  }
  return $self->{OUTPUT_BIOTYPE}
}

sub GENEWISE_PARAMETERS{
  my ($self, $arg) = @_;
  if($arg){
    $self->{GENEWISE_PARAMETERS} = $arg;
  }
  return $self->{GENEWISE_PARAMETERS}
}

sub MINIGENEWISE_PARAMETERS{
  my ($self, $arg) = @_;
  if($arg){
    $self->{MINIGENEWISE_PARAMETERS} = $arg;
  }
  return $self->{MINIGENEWISE_PARAMETERS}
}

sub MULTIMINIGENEWISE_PARAMETERS{
  my ($self, $arg) = @_;
  if($arg){
    $self->{MULTIMINIGENEWISE_PARAMETERS} = $arg;
  }
  return $self->{MULTIMINIGENEWISE_PARAMETERS}
}

sub BLASTMINIGENEWISE_PARAMETERS{
  my ($self, $arg) = @_;
  if($arg){
    $self->{BLASTMINIGENEWISE_PARAMETERS} = $arg;
  }
  return $self->{BLASTMINIGENEWISE_PARAMETERS}
}


sub EXONERATE_PARAMETERS{
  my ($self, $arg) = @_;
  if($arg){
    $self->{EXONERATE_PARAMETERS} = $arg;
  }
  return $self->{EXONERATE_PARAMETERS}
}


sub FILTER_PARAMS{
  my ($self, $arg) = @_;
  if($arg){
    $self->{FILTER_PARAMETERS} = $arg;
  }
  return $self->{FILTER_PARAMETERS}
}



sub FILTER_OBJECT{
  my ($self, $arg) = @_;
  if($arg){
    $self->{FILTER_OBJECT} = $arg;
  }
  return $self->{FILTER_OBJECT}
}


sub BIOTYPES_TO_MASK{
  my ($self, $arg) = @_;
  if($arg){
    $self->{BIOTYPES_TO_MASK} = $arg;
  }
  return $self->{BIOTYPES_TO_MASK}
}


sub EXON_BASED_MASKING{
  my ($self, $arg) = @_;
  if($arg){
    $self->{EXON_BASED_MASKING} = $arg;
  }
  return $self->{EXON_BASED_MASKING}
}


sub GENE_BASED_MASKING{
  my ($self, $arg) = @_;
  if($arg){
    $self->{GENE_BASED_MASKING} = $arg;
  }
  return $self->{GENE_BASED_MASKING}
}


sub POST_GENEWISE_MASK{
  my ($self, $arg) = @_;
  if($arg){
    $self->{POST_GENEWISE_MASK} = $arg;
  }
  return $self->{POST_GENEWISE_MASK}
}

sub PRE_GENEWISE_MASK{
  my ($self, $arg) = @_;
  if($arg){
    $self->{PRE_GENEWISE_MASK} = $arg;
  }
  return $self->{PRE_GENEWISE_MASK}
}

sub REPEATMASKING{
  my ($self, $arg) = @_;
  if($arg){
    $self->{REPEATMASKING} = $arg;
  }
  return $self->{REPEATMASKING}
}

sub SEQFETCHER_OBJECT{
  my ($self, $arg) = @_;
  if($arg){
    $self->{SEQFETCHER_OBJECT} = $arg;
  }
  return $self->{SEQFETCHER_OBJECT}
}

sub SEQFETCHER_PARAMS{
  my ($self, $arg) = @_;
  if($arg){
    $self->{SEQFETCHER_PARAMS} = $arg;
  }
  return $self->{SEQFETCHER_PARAMS}
}

sub USE_KILL_LIST{
  my ($self, $arg) = @_;
  if($arg){
    $self->{USE_KILL_LIST} = $arg;
  }
  return $self->{USE_KILL_LIST}
}


sub LIMIT_TO_FEATURE_RANGE{
  my ($self, $arg) = @_;
  if($arg){
    $self->{LIMIT_TO_FEATURE_RANGE} = $arg;
  }
  return $self->{LIMIT_TO_FEATURE_RANGE}
}


sub FEATURE_RANGE_PADDING{
  my ($self, $arg) = @_;
  if($arg){
    $self->{FEATURE_RANGE_PADDING} = $arg;
  }
  return $self->{FEATURE_RANGE_PADDING}
}

sub WRITE_REJECTED{
  my ($self, $arg) = @_;
  if(defined($arg)){
    $self->{WRITE_REJECTED} = $arg;
  }
  return $self->{WRITE_REJECTED};
}

sub REJECTED_BIOTYPE{
  my ($self, $arg) = @_;
  if($arg){
    $self->{REJECTED_BIOTYPE} = $arg;
  }
  return $self->{REJECTED_BIOTYPE};
}


sub SOFTMASKING{
  my ($self, $arg) = @_;
  if($arg){
    $self->{SOFTMASKING} = $arg;
  }
  return $self->{SOFTMASKING}
}

1;


