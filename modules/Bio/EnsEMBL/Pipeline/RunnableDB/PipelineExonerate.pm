package Bio::EnsEMBL::Pipeline::RunnableDB::PipelineExonerate;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::RunnableDB::ExonerateToGenes;
use Bio::EnsEMBL::Pipeline::Config::cDNAs_ESTs::Exonerate;

use vars qw(@ISA);

@ISA = qw (Bio::EnsEMBL::Pipeline::RunnableDB);



sub fetch_input{
  my ($self) = @_;
  
  my @sequences;
  my $file = $EST_CHUNKDIR."/".$self->input_id;

  my $seqio= Bio::SeqIO->new(
			     -format => "Fasta",
			     -file   => $file,
			    );
  
  ############################################################
  # read the query and create Bio::Seqs
  while( my $seq = $seqio->next_seq() ){
    if ($seq && !($seq->seq eq '') 
	&& !($seq->display_id eq '') ){
      push( @sequences, $seq );
    }
    else{
      print STDERR "problems getting sequence ".$seq->display_id."\n".
	$seq->seq."\n";
    }
  }
  #print STDERR "got ".scalar(@sequences)." sequence objects\n";

  my $runobj = Bio::EnsEMBL::Pipeline::RunnableDB::ExonerateToGenes
    ->new(-db         => $self->db,
	  -input_id   => $self->input_id,
	  -rna_seqs   => \@sequences,
	  -analysis   => $self->analysis,
	  -database   => $EST_GENOMIC,
	  -query_type => 'dna',
	  -target_type=> 'dna',
	  -options    => $EST_EXONERATE->{OPTIONS},
	 );

  $self->exonerate_runnabledb($runobj);

}


sub exonerate_runnabledb{
  my ($self, $arg) = @_;

  if($arg){
    $self->{'exonerate2genes'} = $arg;
  }

  return $self->{'exonerate2genes'}
}


sub run{
  my ($self) =@_;

  my $runnabledb = $self->exonerate_runnabledb;

  
  $runnabledb->fetch_input;
  $runnabledb->run; 

  # the following is only necessary for testing this runnableDB, because
  # the write_output method delegates to the child runnable db
  $self->output($runnabledb->output);
}

sub write_output{
  my ($self) = @_;

  $self->exonerate_runnabledb->write_output;
}


############################################################

# must override RunnableDB::output()

sub output {
  my ($self, @output) = @_;
  if (@output){
    push( @{$self->{_output} }, @output);
  }
  
  my @ret_output;
  if($self->{_output}){
    @ret_output = @{$self->{_output}};
  }

  return @ret_output;
}

1;
