# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

package Bio::EnsEMBL::Pipeline::RunnableDB::PlaceHolder;

use strict;

use Bio::EnsEMBL::Pipeline::RunnableDB;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data for CPG from the database
    Returns :   none
    Args    :   none

=cut

sub fetch_input {
    my( $self) = @_;
    
    $self->throw("No input id") unless defined($self->input_id);

    #my $runname = "Bio::EnsEMBL::Pipeline::Runnable::PlaceHolder";

    #my $runnable = $runname->new;

    #$self->runnable($runnable);

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
