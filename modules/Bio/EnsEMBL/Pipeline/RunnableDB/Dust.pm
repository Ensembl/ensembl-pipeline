#
# You may distribute this module under the same terms as perl itself
#

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::Dust

=head1 SYNOPSIS

my $db      = Bio::EnsEMBL::DBLoader->new($locator);
my $dust = Bio::EnsEMBL::Pipeline::RunnableDB::Dust->new(
    -dbobj    => $db,
    -input_id => $input_id,
    -analysis => $analysis
);
$dust->fetch_input;
$dust->run;
$dust->output;
$dust->write_output; #writes to DB

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Pipeline::Runnable::Dust to add
functionality to read and write to databases. The appropriate
Bio::EnsEMBL::Analysis object must be passed for extraction of
parameters. A Bio::EnsEMBL::Pipeline::DBSQL::Obj is required
for database access.

=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::Dust;

use strict;

use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::Dust;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

=head2 fetch_input

=cut

sub fetch_input {
    my( $self) = @_;
   
    $self->throw("No input id") unless defined($self->input_id);

    $self->fetch_sequence;
    my %parameters      = $self->parameter_hash;
    $parameters{-dust}  = $self->analysis->program_file || undef;
    $parameters{-query} = $self->query;

    my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::Dust(%parameters);

    $self->runnable($runnable);

    return 1;
}

1;
