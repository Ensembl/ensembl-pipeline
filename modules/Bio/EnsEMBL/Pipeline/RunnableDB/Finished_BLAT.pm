# Copyright GRL/EBI 2002
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::Finished_BLAT;


=head1 SYNOPSIS

my $blat2genes = Bio::EnsEMBL::Pipeline::RunnableDB::BlatSSAHA->new(
    -dbobj      => $db,
    -input_id   => $input_id
    -analysis   => $analysis
);

$blat2genes->fetch_input();
$blat2genes->run();
$blat2genes->output();
$blat2genes->write_output(); #writes to DB

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Pipeline::Runnable::Blat
It is meant to provide the interface for mapping ESTs to the genome
sequence and writing the results as genes. By the way Blat is run
(similar to the way Exonerate is run) we do not cluster transcripts into
genes and only write one transcript per gene.
A Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor is required for database storage.

=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::Finished_BLAT;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::BLAT;
use Carp 'cluck';

use vars qw(@ISA);
@ISA = qw (Bio::EnsEMBL::Pipeline::RunnableDB);

sub new {
    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);

    $self->{'_fplist'}   = [];
    $self->{'_genseq'}   = undef;
    $self->{'_runnable'} = undef;
    return $self;
}

sub fetch_input {
    my ($self) = @_;

    $self->throw("No input id") unless defined( $self->input_id );

    my $contigid = $self->input_id;


    my $contig = $self->dbobj->get_Contig($contigid);
    my $genseq = $contig->get_repeatmasked_seq() or $self->throw("Unable to fetch contig");

    $self->genseq($genseq);
    print STDERR "Set genseq to " . $self->genseq . "\n";

    my $seq = $self->genseq->seq;

    if ( $seq =~ /[CATG]{3}/ ) {
        $self->input_is_void(0);
    }
    else {
        $self->input_is_void(1);
        $self->warn("Need at least 3 nucleotides");
    }

    # softmask repeats for blat
    $seq =~ s/N/n/g;
    $self->genseq->seq($seq);


    my $parameter_string = $self->analysis->parameters();

    $parameter_string =~ s/,/ /g;

    my $runnable = Bio::EnsEMBL::Pipeline::Runnable::BLAT->new(
        -database => $self->analysis->db,
        -query    => $self->genseq,
        -blat     => $self->analysis->program,
        -options  => $parameter_string
    );

    $self->runnable($runnable);

}


sub run {
    my ($self) = @_;
    my @results;

    # get the funnable
    $self->throw("Can't run - no runnable objects") unless defined( $self->runnable );
    my $runnable = $self->runnable;

    # run the runnable
    $runnable->run;
    $self->runnable($runnable);

}

sub output {
    my ($self) = @_;

    my @runnable = $self->runnable;
    my @results;

    foreach my $runnable (@runnable) {
        print STDERR "runnable = " . $runnable[0] . "\n";
        push ( @results, $runnable->output );
    }
    return @results;
}

#
# get/set methods
#

sub runnable {
    my ( $self, $runnable ) = @_;
    if ( defined($runnable) ) {
        unless ( $runnable->isa("Bio::EnsEMBL::Pipeline::RunnableI") ) {
            $self->throw("$runnable is not a Bio::EnsEMBL::Pipeline::RunnableI");
        }
        $self->{_runnable} = $runnable;
    }
    return $self->{_runnable};
}

1;
