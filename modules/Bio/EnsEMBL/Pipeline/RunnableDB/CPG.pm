# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::CPG

=head1 SYNOPSIS

my $db      = Bio::EnsEMBL::DBLoader->new($locator);

my $cpg = Bio::EnsEMBL::Pipeline::RunnableDB::CPG->new ( 
                                   -dbobj      => $db,
			           -input_id   => $input_id
                                   -analysis   => $analysis 
                                    );

$cpg->fetch_input();

$cpg->run();

$cpg->output();

$cpg->write_output(); #writes to DB

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Pipeline::Runnable::CPG to add
functionality to read and write to databases. The appropriate
Bio::EnsEMBL::Analysis object must be passed for extraction
of appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
is required for database access.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::CPG;

use strict;

use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::CPG;

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

    $self->fetch_sequence;

    my %parameters = $self->parameter_hash;

    $parameters{'-query'} = $self->query;
    $parameters{'-cpg'} = $self->analysis->program_file;

    my $runname = "Bio::EnsEMBL::Pipeline::Runnable::CPG";

    my $runnable = $runname->new
      ( '-query'  => $parameters{'-query'},
	'-length' => $parameters{'-length'},
	'-gc'     => $parameters{'-gc'},
	'-oe'     => $parameters{'-oe'},
	'-cpg'    => $parameters{'-cpg'},
      );

    $self->runnable($runnable);

    return 1;

}

1;
