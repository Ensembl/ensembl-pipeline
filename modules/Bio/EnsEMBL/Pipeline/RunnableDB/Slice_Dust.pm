# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::Slice_Dust

=head1 SYNOPSIS

my $db      = Bio::EnsEMBL::DBLoader->new($locator);

my $eponine = Bio::EnsEMBL::Pipeline::RunnableDB::Slice_Dust->new ( 
                                   -db          => $db,
			           -input_id   => $input_id
                                   -analysis   => $analysis 
                                    );

$eponine->fetch_input();

$eponine->run();

$eponine->output();

$eponine->write_output(); #writes to DB

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Pipeline::Runnable::Dust to add
functionality to read and write to databases. The appropriate
Bio::EnsEMBL::Analysis object must be passed for extraction
of appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
is required for database access.

=head1 CONTACT

For general Ensembl comments mail to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::Slice_Dust;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::Dust;
use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);


=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data for Dust from the database
    Returns :   none
    Args    :   none

=cut

sub fetch_input {
    my ($self) = @_;
    
    $self->throw("No input id") unless defined($self->input_id);

    my $slice_str = $self->input_id;
    my ($chr, $start, $end, $sgp) =
     $slice_str =~ m!(\S+)\.(\d+)\.(\d+):?([^:]*)!;

    $self->db->assembly_type($sgp) if $sgp;

    my $slice = $self->db->get_SliceAdaptor->
     fetch_by_chr_start_end($chr, $start, $end);

    $self->throw("Unable to fetch virtual contig") unless $slice;

    $self->query($slice);

    my %parameters = $self->parameter_hash;

    $self->runnable(Bio::EnsEMBL::Pipeline::Runnable::Dust->new(
        '-level' => $parameters{'-level'},
	'-dust'  => $self->analysis->program_file,
        '-query' => $self->query
    ));
    return 1;
}



sub write_output {
    my($self) = @_;
    my $contig;

    my $db  = $self->db;
    my $sfa = $db->get_SimpleFeatureAdaptor;

    my $slice = $self->query;
    my @mapped_features;

    foreach my $f ($self->output) {

	$f->analysis($self->analysis);
	$f->contig($slice);

	my (@mapped) = $f->transform;

	unless (@mapped == 1) {
	    print STDERR "Warning: can't map $f";
	    next;
	}

	push @mapped_features, $mapped[0];

	my $contig_id = $f->contig->dbID;
    }

    $sfa->store(@mapped_features) if @mapped_features;

    return 1;
}

1;
