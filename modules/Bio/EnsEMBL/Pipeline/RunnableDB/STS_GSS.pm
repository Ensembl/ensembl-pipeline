

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
  my $contig    = $self->dbobj->get_Contig($contigid);
  
  my $genseq   = $contig->primary_seq;
  my $repeatmasked_seq = $contig->get_repeatmasked_seq;


    print STDERR "Set genseq to " . $repeatmasked_seq. "\n";
    
    # input sequence needs to contain at least 3 consecutive nucleotides
    my $seq = $repeatmasked_seq->seq;
    if ($seq =~ /[CATG]{3}/) {
        $self->input_is_void(0);
    }
    else {
        $self->input_is_void(1);
        $self->warn("Need at least 3 nucleotides");
    }

    #extract parameters into a hash
    my ( $parameter_string ) = $self->analysis->parameters();
    my ( $arguments, $thresh_type, $thresh );
    
    if ($parameter_string)
    {
        $parameter_string =~ s/\s+//g;
        my @pairs = split ( /,/, $parameter_string );
        foreach my $pair ( @pairs )
        {
            my ($key, $value) = split (/=>/, $pair);
            if ($key eq '-threshold_type' && $value) {
                $thresh_type = $value;
            }
            elsif ($key eq '-threshold' && $value) {
                $thresh = $value;
            }
            else
	    # remaining arguments not of '=>' form
	    # are simple flags (like -p1)
            {
                $arguments .= " $key ";
            }
        }
    }


  
    my $percent_id = 95;
    
    my $runnable  = Bio::EnsEMBL::Pipeline::Runnable::STS_GSS->new(-query => $repeatmasked_seq,
							         -unmasked =>$genseq,
							         -database =>$self->analysis->db,
							         -program =>$self->analysis->program,
							         -options => $arguments,
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
    
    $runnable->run;
    
    push (@{$self->{'_output'}}, $runnable->output);
    foreach my $f(@{$self->{'_output'}}){
      $f->source_tag($self->analysis->db);
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
  
  my $seqfetcher = new Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;
  

  return $seqfetcher;

}


1;
