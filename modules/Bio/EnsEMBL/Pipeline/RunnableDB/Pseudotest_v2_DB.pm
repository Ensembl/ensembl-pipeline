# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::Pseudotest_v2_DB.pm

=head1 SYNOPSIS



=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Pipeline::Tools::Pseudotest_v2.pm to add
functionality to read from databases (so far).

=head1 CONTACT

Describe contact details here

=head1 APPENDIX



=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::Pseudotest_v2_DB;

use strict;

use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Tools::Pseudotest_v2;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data for Pseudotest_v2.pm from the database
    Returns :   none
    Args    :   none

=cut



sub fetch_input {
    my( $self) = @_;
    
    $self->throw("No input id") unless defined($self->input_id);

    $self->fetch_sequence;

    my %parameters = $self->parameter_hash;

    $parameters{'-query'} = $self->query;

    my $runname = "Bio::EnsEMBL::Pipeline::Tools::Pseudotest_v2";

    my $runnable = $runname->new
      ( '-query'  => $parameters{'-query'},
	'-max_intron_length' => $parameters{'-max_intron_length'},
	'-max_intron_coverage' => $parameters{'-max_intron_coverage'},
	'-max_exon_coverage' => $parameters{'-max_exon_coverage'},
      );

    $self->runnable($runnable);

    return 1;

}

1;
