#
# Written by Stephen Keenan
#
# Copyright GRL/EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::FindPrimer

=head1 SYNOPSIS

  # create and fill Bio::Seq object
  my $seqfile = '/path/to/seq.fa';
  my $seq = Bio::Seq->new();
  my $seqstream = Bio::SeqIO->new(
    -file => $seqfile,
    -fmt => 'Fasta'
  );
  $seq = $seqstream->next_seq();

  # create Bio::EnsEMBL::Pipeline::Runnable::FindPrimer object
  my $find_primer = Bio::EnsEMBL::Pipeline::Runnable::FindPrimer->new(
    -CLONE => $seq
  );
  $findprimer->run();

  # get results
  my @results = $trf->output();

=head1 DESCRIPTION

FindPrimer takes a Bio::Seq (or Bio::PrimarySeq) object and finds exacdt matches to primers conatined in 
fastfile. Asolute path must be set in analysisprocess table as db_file. 

=head2 Methods:

=over4

=item new($seq_obj)

=item run()

=item output()

=back

=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::Runnable::FindPrimer;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::RootI;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::Repeat;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::Analysis;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Root::RootI;
use FileHandle;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

=head2 new

=cut

sub new {
    my ( $class, @args ) = @_;

    my $self = bless {}, $class;
    my ( $clone, $analysis ) = $self->_rearrange(
        [
            qw{
            CLONE
            ANALYSIS
            }
        ],
        @args
    );

    die "No UNMASKED (genomic sequence) given" unless $clone;

    $self->clone($clone);
    $self->analysis($analysis);

    return $self;

}


sub clone {
    my ( $self, $seq ) = @_;

    if ($seq) {
        unless ( $seq->isa("Bio::PrimarySeqI") || $seq->isa("Bio::Seq") ) {
            $self->throw("Input isn't a Bio::Seq or Bio::PrimarySeq");
        }
        $self->{'_clone'} = $seq;
    }
    return $self->{'_clone'};
}

sub analysis {
    my ( $self, $analysis ) = @_;

    if ($analysis) {
        $self->{'_analysis'} = $analysis;
    }
    return $self->{'_analysis'};
}

###########
# Analysis methods
##########

=head2 run

=cut

sub run {
    my ($self) = @_;

    # check seq
    my $seq = $self->clone() || $self->throw("Seq required for FindPrimer\n");

    $self->find_primer();

}

=head2 parsefile

=cut

sub find_primer {
    my ($self) = @_;

    # options() set in clone() at same time as results()
    # both are a catenation of the parameters
    my $primer_file = $self->analysis->db_file;

    my $primer_in = Bio::SeqIO->new(
        -FILE   => $primer_file,
        -FORMAT => 'fasta',
    );

    my (@primer);
    while ( my $pri = $primer_in->next_seq ) {
        push ( @primer, $pri );
    }

    my $gen = $self->clone;

    my $gen_id = $gen->id;

    foreach my $strand ( 1, -1 ) {
        my ($seq);
        if ( $strand == 1 ) {
            $seq = $gen;
        }
        else {
            $seq = $gen->revcom;
        }
        my $gen_str = lc $seq->seq;
        foreach my $pri (@primer) {
            my $pri_str = lc $pri->seq;
            my $pri_id  = $pri->id;
            my $i       = 0;
            while ( ( $i = index( $gen_str, $pri_str, $i ) ) != -1 ) {
                my ( $start, $end );
                my $pri_length = length($pri_str);
                if ( $strand == 1 ) {
                    $start = $i + 1;
                    $end   = $i + $pri_length;
                }
                else {
                    my $gen_length = length($gen_str);
                    $end   = $gen_length - $i;
                    $start = $end - $pri_length + 1;
                }
               # printf "%-15s  %-15s  %+2d  %9d  %9d\n", $gen_id, $pri_id,
               #   $strand, $start, $end;

                my ( %feat1, %feat2 );

                $feat1{name}    = $gen_id;
                $feat1{start}   = $start;
                $feat1{end}     = $end;
                $feat1{strand}  = $strand;
                $feat1{score}   = '100';
                $feat1{source}  = 'FindPrimer';
                $feat1{primary} = 'primer';
                $feat1{percent} = '100';
                $feat1{p_value} = '';

                $feat2{name}    = $pri_id;
                $feat2{start}   = 1;
                $feat2{end}     = $pri_length;
                $feat2{strand}  = $strand;
                $feat2{score}   = '0';
                $feat2{source}  = 'FindPrimer';
                $feat2{primary} = 'primer';
                $feat2{percent} = '100';
                $feat2{p_value} = '';

                $feat2{db}         = undef;
                $feat2{db_version} = undef;
                $feat2{program}    = '';
                $feat2{p_version}  = '1';
                $feat2{source}     = 'primer';
                $feat2{primary}    = 'primer';

                $self->createfeaturepair( \%feat1, \%feat2 );
                $i++;
            }
        }
    }
}


sub output {
    my ($self) = @_;
    if (!defined($self->{'_fplist'})) {
    $self->{'_fplist'} = [];
    }
    return @{$self->{'_fplist'}};
}

__END__

=head1 NAME - Bio::EnsEMBL::Ace::Filter::Repeatmasker

=head1 AUTHOR

Stephen Keenan B<email> keenan@sanger.ac.uk

1;
