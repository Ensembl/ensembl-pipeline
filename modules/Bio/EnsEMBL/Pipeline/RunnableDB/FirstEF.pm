#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB:FirstEF

=head1 SYNOPSIS

my $fef = Bio::EnsEMBL::Pipeline::RunnableDB::FirstEF->new(-dbobj     => $db,
			                                   -input_id  => $input_id,
                                                           -analysis  => $analysis );
$fef->fetch_input;
$fef->run;
my @output = $fef->output;
$fef->write_output;

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Pipeline::Runnable::FirstEF to add
functionality for reading and writing to databases. The appropriate
Bio::EnsEMBL::Analysis object must be passed for extraction
of appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
is required for database access.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::FirstEF;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::FirstEF;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data from the database
    Returns :   none
    Args    :   none

=cut

sub fetch_input {
    my( $self) = @_;    
    
    $self->throw("No input id") unless defined($self->input_id);
    
    my $contigid  = $self->input_id;
    my $contig    = $self->db->get_RawContigAdaptor->fetch_by_name($contigid);

    $self->query($contig);

    # For practical purposes the following information is hard-coded.
    # If this runnableDB ends up in more widespread use it would
    # make sense to construct a config file and place it in 
    # Bio::EnsEMBL::Pipeline::Config

    # Directory where firstef.* binaries and FirstEF_parser.pl 
    # are located.  At runtime, the runnable determines which
    # platform specific binary should be used.
    my $APPLICATION_DIR = '/usr/local/ensembl/firstef/';

    # Usually, this directory is a sub-directory of the firstef
    # application directory called 'parameters'.  This directory
    # contains all of the training data that firstef uses for
    # identifying first exons.
    my $PARAMETER_DIR   = '/usr/local/ensembl/firstef/parameters';


    my $runnable = Bio::EnsEMBL::Pipeline::Runnable::FirstEF->new(
		     -query       => $self->query,
		     -db          => $self->db,
		     -firstef_dir => $APPLICATION_DIR,
		     -param_dir   => $PARAMETER_DIR);

							   
    $self->runnable($runnable);
    
    return 1;
}

1;
