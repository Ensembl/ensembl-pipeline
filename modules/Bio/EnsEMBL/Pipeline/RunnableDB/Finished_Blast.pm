
package Bio::EnsEMBL::Pipeline::RunnableDB::Finished_Blast;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB::Blast;
use Bio::EnsEMBL::Pipeline::Runnable::Finished_Blast;

use Data::Dumper;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB::Blast);

=head2 run

    Title   :   run
    Usage   :   $self->run();
    Function:   Runs Bio::EnsEMBL::Pipeline::Runnable::Blast->run()
    Returns :   none
    Args    :   none

=cut

sub run {
    my ($self) = @_;

    #need to pass one peptide at a time
    $self->throw("Input must be fetched before run") unless ( $self->genseq );

    #extract parameters into a hash
    my ($parameter_string) = $self->analysis->parameters();
    my %parameters;
    my ( $thresh, $thresh_type, $arguments );
    my $split_gapped = 0;

    if ($parameter_string) {
        $parameter_string =~ s/\s+//g;
        my @pairs = split ( /,/, $parameter_string );
        foreach my $pair (@pairs) {
            my ( $key, $value ) = split ( /=>/, $pair );
            if ( $key eq '-threshold_type' && $value ) {
                $thresh_type = $value;
            }
            elsif ( $key eq '-threshold' && $value ) {
                $thresh = $value;
            }
            elsif ( $key eq '-split_gapped_alignments' ) {
                $split_gapped = 1;
            }
            else

              # remaining arguments not of '=>' form
              # are simple flags (like -p1)
            {
                $arguments .= " $key ";
            }
        }
    }
    print $arguments. "\n";
    $parameters{'-query'}      = $self->genseq;
    $parameters{'-database'}   = $self->analysis->db;
    $parameters{'-program'}    = $self->analysis->program;
    $parameters{'-options'}    = $arguments if $arguments;
    $parameters{'-seqfetcher'} = $self->make_seqfetcher;
    if ( $thresh && $thresh_type ) {
        $parameters{'-threshold'}      = $thresh;
        $parameters{'-threshold_type'} = $thresh_type;
    }
    else {
        $parameters{'-threshold'}      = 1e-4;
        $parameters{'-threshold_type'} = 'PVALUE';
    }

    my $runnable = Bio::EnsEMBL::Pipeline::Runnable::Finished_Blast->new(%parameters);
    $runnable->split_gapped_alignments($split_gapped);

    $runnable->run();
    $self->runnable($runnable);

    # Write the descrptions
    my @output = $runnable->output;
    my $dbobj = $self->dbobj;
    my $seqfetcher = $self->make_seqfetcher;
    my %ids = map { $_->hseqname, $_ } @output;
    $seqfetcher->write_descriptions($dbobj, keys(%ids) );
   

}

sub make_seqfetcher {
    my ( $self, $index ) = @_;

    my $seqfetcher;
    if ( my $dbf = $self->analysis->db_file ) {

        my @dbs = $dbf;
        $seqfetcher = Bio::EnsEMBL::Pipeline::SeqFetcher::OBDAIndexSeqFetcher->new(
            -db => \@dbs,
        );
    }
    else {
        $seqfetcher = Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch->new;
    }
    return $seqfetcher;
}


1;
