# Copyright GRL/EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::TRF

=head1 SYNOPSIS

my $db      = Bio::EnsEMBL::DBLoader->new($locator);
my $trf = Bio::EnsEMBL::Pipeline::RunnableDB::TRF->new(
    -dbobj      => $db,
    -input_id   => $input_id
    -analysis   => $analysis
);
$trf->fetch_input();
$trf->run();
$trf->output();
$trf->write_output(); #writes to DB

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Pipeline::Runnable::TRF to add
functionality to read and write to databases. The appropriate
Bio::EnsEMBL::Analysis object must be passed for extraction of
parameters. A Bio::EnsEMBL::Pipeline::DBSQL::Obj is required
for databse access.

=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::TRF;

use strict;

use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::TRF;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

=head2 fetch_input

=cut

sub fetch_input {
    my( $self) = @_;
   
    $self->throw("No input id") unless defined($self->input_id);

    $self->fetch_sequence;

    my %parameters      = $self->parameter_hash;
    $parameters{-trf}   = $self->analysis->program_file || undef;
    $parameters{-query} = $self->query;

    my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::TRF(%parameters);

    $self->runnable($runnable);

    return 1;
							     

}





1;
