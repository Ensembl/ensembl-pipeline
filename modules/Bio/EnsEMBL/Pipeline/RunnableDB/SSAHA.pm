#
# Written by Simon Potter
#
# Copyright GRL/EBI 2002
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::SSAHA

=head1 SYNOPSIS

my $db    = Bio::EnsEMBL::DBLoader->new($locator);
my $ssaha = Bio::EnsEMBL::Pipeline::RunnableDB::SSAHA->new(
    -dbobj      => $db,
    -input_id   => $input_id
    -analysis   => $analysis
);

$ssaha->fetch_input();
$ssaha->run();
$ssaha->output();
$ssaha->write_output(); #writes to DB

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Pipeline::Runnable::SSAHA to add
functionality for reading and writing to databases. The appropriate
Bio::EnsEMBL::Analysis object must be passed for extraction of
appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
is required for databse access.

=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::SSAHA;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::SSAHA;

use vars qw(@ISA);

@ISA = qw (Bio::EnsEMBL::Pipeline::RunnableDB);


sub fetch_input {
    my( $self) = @_;

    $self->throw("No input id") unless defined($self->input_id);

    $self->fetch_sequence;
    my %parameters = $self->parameter_hash;

    my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::SSAHA(
            -query  => $self->query,
            -ssaha  => $self->analysis->program_file,
            -db     => $self->analysis->db_file,
            -min_pc => $parameters{'-min_pc'},
            -length => $parameters{'-length'}
    );

    $self->runnable($runnable);
}


1;
