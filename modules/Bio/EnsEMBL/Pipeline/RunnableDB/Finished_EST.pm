
### Bio::EnsEMBL::Pipeline::RunnableDB::Finished_EST

package Bio::EnsEMBL::Pipeline::RunnableDB::Finished_EST;

use strict;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Getseqs;

sub new {
    my ($new,@args) = @_;
    my $self = $new->SUPER::new(@args);    
           
    # dbobj, input_id, seqfetcher, and analysis objects are all set in
    # in superclass constructor (RunnableDB.pm)
  	
    return $self;
}

sub fetch_input {
    my( $self) = @_;

    my @fps;

    $self->throw("No input id") unless defined($self->input_id);

    my $contigid  = $self->input_id;
    my $contig    = $self->dbobj->get_Contig($contigid);
    my $genseq    = $contig->primary_seq;
    my $masked    = $contig->get_repeatmasked_seq->seq;

    # Make seqfetcher
    my $seqfetcher = $self->make_seqfetcher;
    $self->seqfetcher($seqfetcher);

    my $runnable = Bio::EnsEMBL::Pipeline::Runnable::Finished_EST->new(
        '-query'        => $masked,
        '-unmasked'     => $genseq,
        '-seqfetcher'   => $seqfetcher,
        '-analysis'     => $self->analysis,
        );

    $self->runnable($runnable);
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
}

sub output {
    my ($self) = @_;

    my @runnable = $self->runnable;
    my @results;
    
    foreach my $runnable (@runnable){
      print STDERR "runnable = ".$runnable[0]."\n";
      push(@results, $runnable->output);
    }
    return @results;  
}


sub make_seqfetcher {
    my ( $self, $index ) = @_;

    my( $seqfetcher );
    if (my $dbf = $self->analysis->db_file) {
        my $index = "$ENV{BLASTDB}/$dbf";
        $seqfetcher = Bio::EnsEMBL::Pipeline::SeqFetcher::Getseqs->new('-db' => [$index]);
    } else {
        $seqfetcher = Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch->new;
    }
    return $seqfetcher;
}



1;

__END__

=head1 NAME - Bio::EnsEMBL::Pipeline::RunnableDB::Finished_EST

=head1 AUTHOR

James Gilbert B<email> jgrg@sanger.ac.uk

