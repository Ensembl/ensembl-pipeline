package Bio::EnsEMBL::Pipeline::Runnable::Fathom;

use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Pipeline::RunnableI;


@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);



sub new{
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  $self->{'_query'} = undef;
  $self->{'_genes'} = undef;
  $self->{'_fathom'} = undef;
  $self->{'_workdir'} = undef;
  $self->{'_filename'}  = undef; # file to store Bio::Seq object
  $self->{'_results'}   = undef; # file to store results of genscan
  $self->{'_protected'} = [];    # a list of file suffixes protected from 
  #deletion
  $self->{'_parameters'} = undef; #location of parameters for genscan
  $self->{'_hmmfile'} = undef;
  $self->{'_output'} = undef;
  my ($query, $genes, $fathom, $parameters, $hmmfile) = 
    $self->_rearrange([qw(QUERY GENES FATHOM PARAM HMMFILE)], @args);
  print STDERR "@args\n";
  $fathom = '/usr/local/ensembl/Zoe/bin/fathom' unless $fathom;

  $self->throw("A query sequence must be provided") unless($query);
  $self->query($query);
  $self->throw("You must provide a ref to an array of genes") 
    unless($genes);
  if(ref($genes) eq 'ARRAY'){
    $self->genes($genes);
  }else{
    $self->throw("[$genes] is not an array ref.");
  }
  $self->fathom($fathom);
  $self->parameters($parameters);
  $self->throw("you must provide a hmmfile") unless($hmmfile);
  $self->hmmfile($hmmfile);

  return $self;
}


sub query {
    my ($self, $seq) = @_;
    if ($seq){
      unless ($seq->isa("Bio::PrimarySeqI") || $seq->isa("Bio::SeqI")){
        $self->throw("Input isn't a Bio::Seq or Bio::PrimarySeq");
      }
      $self->{'_query'} = $seq ;
      $self->filename($seq->id.".$$.seq");
      $self->results($self->filename.".fathom");
      $self->genefile($self->filename.".genes");
      $self->file($self->results);
    }
    return $self->{'_query'};
}


sub fathom {
    my ($self, $location) = @_;
    if ($location){
      $self->throw("fathom not found at $location: $!")
        unless (-e $location);
      $self->{'_fathom'} = $location ;
    }
    return $self->{'_fathom'};
}


sub genes{
  my ($self, $genes) = @_;

  if($genes){
    $self->{'_genes'} = $genes;
  }

  return $self->{'_genes'};
}


sub genefile{
  my ($self, $genefile) = @_;
  if($genefile){
    $self->{'_genefile'} = $genefile;
  }

  return $self->{'_genefile'};
}

sub parameters {
    my ($self, $param) = @_;
   
    if ($param){
      $self->{'_parameters'} = $param;
    }
    return $self->{'_parameters'};
}

sub hmmfile {
  my ($self, $location) = @_;
  
  if ($location) {
    $self->{'_hmmfile'} = $location ;
  }
  
  return $self->{'_hmmfile'};
}


sub run{
  my ($self) = @_;
  $self->workdir('/tmp') unless $self->workdir();
  $self->checkdir();
  $self->writefile();
  $self->write_genefile();
  $self->run_fathom();
  $self->parse_results();
  $self->deletefiles();
  unlink ($self->genefile) or
    $self->throw ("Couldn't delete ".$self->genefile.":$!");

  return 1;
}



sub write_genefile{
  my ($self) =  @_;
  open (FH, ">".$self->genefile) or die("can't open ".$self->genefile);
  my $fh = \*FH;
  print $fh ">".$self->query->id."\n";
  foreach my $gene(@{$self->genes}){
  TRANS:foreach my $trans(@{$gene->get_all_Transcripts}){
      my $trans_name = $trans->stable_id;
      $trans_name = $trans->dbID unless($trans_name);
      $trans->sort;
      my @exons = @{$trans->get_all_Exons};
      if(@exons == 1){
        my $e = $exons[0];
        if($e->strand == -1){
          $self->single_exon($fh, $e->end, $e->start, $trans_name) 
        }else{
          $self->single_exon($fh, $e->start, $e->end, $trans_name);
        }
      }else{
        my $count = 1;
        foreach my $e(@exons){
          if($e->start < 1 || $e->end < 1){
            next TRANS;
          }
          if($count == 1){
            if($e->strand == -1){
              $self->start_exon($fh, $e->end, $e->start, $trans_name) 
            }else{
              $self->start_exon($fh, $e->start, $e->end, $trans_name);
            }
          }elsif($count == scalar(@exons)){
            if($e->strand == -1){
              $self->end_exon($fh, $e->end, $e->start, $trans_name) 
            }else{
              $self->end_exon($fh, $e->start, $e->end, $trans_name);
            }
          }else{
            if($e->strand == -1){
              $self->internal_exon($fh, $e->end, $e->start, $trans_name) 
            }else{
              $self->internal_exon($fh, $e->start, $e->end, $trans_name);
            }
          }
          $count++;
        }
      }
    }
  }
  close(FH) or die("Couldn't close ".$self->genefile);
}

sub single_exon{
  my ($self, $fh, $start, $end, $trans_name) = @_;
  
  print $fh "Esngl\t".$start."\t".$end."\t".$trans_name."\n";

}

sub start_exon{
  my ($self, $fh, $start, $end, $trans_name) = @_;

  print $fh "Einit\t".$start."\t".$end."\t".$trans_name."\n";
}
sub internal_exon{
  my ($self, $fh, $start, $end, $trans_name) = @_;

  print $fh "Exon\t".$start."\t".$end."\t".$trans_name."\n";
}
sub end_exon{
  my ($self, $fh, $start, $end, $trans_name) = @_;

  print $fh "Eterm\t".$start."\t".$end."\t".$trans_name."\n";
}


sub run_fathom{
  my ($self) = @_;

  $self->check_environment;
  my $command = $self->fathom." -score-genes ".$self->hmmfile." ".
    $self->genefile." ".$self->filename." > ".$self->results;
  print STDERR $command."\n";
  system($command);
  $self->throw($command." didn't run sucessfully no results exist") 
    unless(-e $self->results);

}


sub parse_results{
  my ($self) = @_;

  open(OUT, $self->results) or
    $self->throw("Couldn't open fathom_results.txt");
  my %score_hash;
  while(<OUT>){
    chomp;
    #ENST00000248853	score=-354
    /(\S+)\s+score\=(\S+)/;
    my $name = $1;
    my $score = $2;
    print "name ".$name." score ".$score."\n";
    $score_hash{$name} = $score;
  }
  print "have ".%score_hash."\n";
  print "have ".\%score_hash."\n";
  $self->output(\%score_hash);
}


sub output{
  print "args @_\n";
  my ($self, $score_hash) = @_;
  if($score_hash){
    if(ref($score_hash) eq 'HASH'){
      $self->{'_output'} = $score_hash;
    }else{
      $self->throw("[$score_hash] is not a hashref");
    }
  }
  return $self->{'_output'};
}




sub check_environment {
  my($self) = @_;
  if (! -d $ENV{ZOE}) {
    $self->throw("No ZOE ["  . $ENV{ZOE} . "]");
  }
}


1;
