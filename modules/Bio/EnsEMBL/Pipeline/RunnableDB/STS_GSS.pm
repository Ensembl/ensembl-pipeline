

package Bio::EnsEMBL::Pipeline::RunnableDB::STS_GSS;

use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Finished_Pfetch;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Getseqs;
use Bio::EnsEMBL::Pipeline::Config::General;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Pipeline::Runnable::STS_GSS;


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

=head2 fetch_input

  Arg [1]   : none
  Function  : fetches contig and creates sts_gss runnable
  Returntype: none
  Exceptions: throws if not given a contig inputId
  Caller    : 
  Example   : 

=cut

sub fetch_input {
  my( $self) = @_;
 
  my @fps;
  my %ests;
  my @estseqs;
  $self->throw("No input id") unless defined($self->input_id);
    
  my $contigid  = $self->input_id;
  my $rawContigAdaptor = $self->db->get_RawContigAdaptor();
  my $contig   = $rawContigAdaptor->fetch_by_name($contigid);

  my $genseq = $contig;
  my $repeatmasked_seq = $contig->get_repeatmasked_seq($PIPELINE_REPEAT_MASKING);
 
  my $percent_id = 95;

  
  my $runnable  = Bio::EnsEMBL::Pipeline::Runnable::STS_GSS->new(-query => $repeatmasked_seq,
								 -unmasked =>$genseq,
								 -database =>$self->analysis->db,
								 -program =>$self->analysis->program,
								 -options => $self->analysis->parameters,
								 -seqfetcher => $self->seqfetcher,
								-percent_id => $percent_id
								-percent_filter => 1);
  
  $self->runnable($runnable);
}    
      
  
    
    
    
=head2 runnable

  Arg [1]   : runnable object
  Function  : sets runnable varible to runnable passed
  Returntype: runnable object
  Exceptions: if arg passed isn't a runnableI object'
  Caller    : 
  Example   : 

=cut    



sub runnable {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->throw("[$arg] is not a Bio::EnsEMBL::Pipeline::RunnableI") unless $arg->isa("Bio::EnsEMBL::Pipeline::RunnableI");
	
	$self->{_runnable} = $arg;
    }

    return $self->{_runnable};
}

=head2 run

  Arg [1]   : none 
  Function  : runs runnable and pushes output onto output array
  Returntype: none
  Exceptions: throws if hasn't got a runnable'
  Caller    : 
  Example   : 

=cut

sub run {
    my ($self) = @_;
    
    my $runnable = $self->runnable;
    $runnable || $self->throw("Can't run - no runnable object");
    eval{
        $runnable->run;
    };
    if(my $err = $@){
        chomp $err;
        $self->failing_job_status($1) 
            if $err =~ /^\"([A-Z_]{1,40})\"$/i; # only match '"ABC_DEFGH"' and not all possible throws
        $self->throw("$@");
    }
    push (@{$self->{'_output'}}, $runnable->output);
    foreach my $f(@{$self->{'_output'}}){
      #$f->source_tag($self->analysis->db);
      #print $f->source_tag."\n";
    }
}


=head2 output

  Arg [1]   : none
  Function  : returns output array
  Returntype: 
  Exceptions: none
  Caller    : 
  Example   : 

=cut

sub output {
    my ($self) = @_;
    return @{$self->{_output}};
}

=head2 make_seqfetcher

  Arg [1]   :  none
  Function  : makes a pfetch obj 
  Returntype: seqfetcher object
  Exceptions: none
  Caller    : 
  Example   : 

=cut


sub make_seqfetcher {
  my ( $self ) = @_;
  
  my $seqfetcher = new Bio::EnsEMBL::Pipeline::SeqFetcher::Finished_Pfetch;
  

  return $seqfetcher;

}


1;
