# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

package Bio::EnsEMBL::Pipeline::RunnableDB::Accumulator;

use strict;

use Bio::EnsEMBL::Pipeline::RunnableDB;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Dummy method to comply to the interface
    Returns :   none
    Args    :   none

=cut

sub fetch_input {
    my( $self) = @_;
    
    $self->throw("No input id") unless defined($self->input_id);

    return 1;

}

sub run {
    my ($self) = @_;

    print "Dummy RunnableDB - no runnable to run\n";

}

sub write_output {
    my ($self) = @_;

    print "Dummy RunnableDB - no output to write\n";

    return 1;
}

1;
