

package Bio::EnsEMBL::Pipeline::RunnableDB::STS_GSS;

use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::Root::RootI;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Getseqs;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Pipeline::Runnable::STS_GSS;

require "Bio/EnsEMBL/Pipeline/GB_conf.pl";

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);


sub new {
    my ($new,@args) = @_;
    my $self = $new->SUPER::new(@args);    
           
    # dbobj, input_id, seqfetcher, and analysis objects are all set in
    # in superclass constructor (RunnableDB.pm)

    $self->{'_fplist'} = []; #create key to an array of feature pairs
    my $seqfetcher = $self->make_seqfetcher();
    #print $seqfetcher."\n"; 
    $self->seqfetcher($seqfetcher);
	
    return $self;
}



sub fetch_input {
  my( $self) = @_;
  #print "running fetch input\n";  
  my @fps;
  my %ests;
  my @estseqs;
  $self->throw("No input id") unless defined($self->input_id);
  
  my $contigid  = $self->input_id;
  my $contig    = $self->dbobj->get_Contig($contigid);
  #print "got contig\n";
  my $genseq   = $contig->primary_seq;
  my $repeatmasked_seq = $contig->get_repeatmasked_seq;
  #print "got dnaseq\n";
  

  #print "got all ests\n";
  #print "logic name = ".$self->analysis->logic_name."\n";
  my $runnable  = Bio::EnsEMBL::Pipeline::Runnable::STS_GSS->new(-query => $repeatmasked_seq,
								 -unmasked =>$genseq,
								 -database =>$self->analysis->db,
								 -program =>$self->analysis->program,
								 -options => $self->analysis->parameters,
								 -seqfetcher => $self->seqfetcher);
  #print "created STS_GSS Runnable\n";  
  $self->runnable($runnable);
  #print "finshed fetching input\n";
}    
      
  
    
    
    
    



sub runnable {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->throw("[$arg] is not a Bio::EnsEMBL::Pipeline::RunnableI") unless $arg->isa("Bio::EnsEMBL::Pipeline::RunnableI");
	
	$self->{_runnable} = $arg;
    }

    return $self->{_runnable};
}

sub run {
    my ($self) = @_;

    my $runnable = $self->runnable;
    $runnable || $self->throw("Can't run - no runnable object");

    $runnable->run;
    push (@{$self->{'_output'}}, $runnable->output);
    #foreach my $f(@{$self->{'_output'}}){
    #  print $f->source_tag."\n";
    #}
}

sub output {
    my ($self) = @_;
    return @{$self->{_output}};
}



sub make_seqfetcher {
  my ( $self ) = @_;
  
  my $seqfetcher = new Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;
  

  return $seqfetcher;

}


1;
