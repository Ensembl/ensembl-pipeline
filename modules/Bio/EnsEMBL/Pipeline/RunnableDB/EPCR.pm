# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::EPCR

=head1 SYNOPSIS

my $db      = Bio::EnsEMBL::DBLoader->new($locator);
my $epcr   = Bio::EnsEMBL::Pipeline::RunnableDB::EPCR->new ( 
                                                    -dbobj      => $db,
			                            -input_id   => $input_id
                                                    -analysis   => $analysis );
$epcr->fetch_input();
$epcr->run();
$epcr->output();
$epcr->write_output(); #writes to DB

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Pipeline::Runnable::EPCR to add
functionality for reading and writing to databases. The appropriate
Bio::EnsEMBL::Analysis object must be passed for extraction
of appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
is required for databse access.

=head1 CONTACT

For general Ensembl comments mail to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::EPCR;

use strict;

use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::EPCR;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data for epcr from the database
    Returns :   none
    Args    :   none

=cut

sub fetch_input {
    my($self) = @_;
    my $sts;
    
    $self->throw("No input id") unless defined($self->input_id);

    $self->fetch_sequence;
    my %parameters = $self->parameter_hash;

    my $db_file = $self->analysis->db_file;
    if ($db_file) {
	$sts = $db_file;
    }
    else {
        $sts = $self->db->get_MarkerAdaptor->fetch_all;
	$self->throw("No markers in database") unless @{$sts};
    }

    my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::EPCR(
	    -query => $self->query,
	    -sts   => $sts,	  # either a filename or an array reference
            -pcr   => $self->analysis->program_file,
            -nmin  => $parameters{'-NMIN'},
            -nmax  => $parameters{'-NMAX'},
            -w     => $parameters{'-W'},
            -m     => $parameters{'-M'}
    );

    $self->runnable($runnable);

    return 1;
}


1;
