# You may distribute this module under the same terms as perl itself

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::Slice_CPG

=head1 SYNOPSIS

my $db   = Bio::EnsEMBL::DBLoader->new($locator);
my $cpg = Bio::EnsEMBL::Pipeline::RunnableDB::Slice_CPG->new( 
    -dbobj      => $db,
    -input_id   => $input_id,   # chr.start.end
    -analysis   => $analysis
);
$cpg->fetch_input();
$cpg->run();
$cpg->output();
$cpg->write_output(); #writes to DB

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Pipeline::Runnable::CPG to add
functionality for reading and writing to databases. The appropriate
Bio::EnsEMBL::Analysis object must be passed for extraction
of appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
is required for database access.
Takes as input a Slice ("chr.start.end"). cpg is too quick to run on
RawContigs.

=head1 CONTACT

For general Ensembl comments mail to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::Slice_CPG;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::CPG;
use Bio::EnsEMBL::Pipeline::Config::General;
use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);


=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data for cpg from the database
    Returns :   none
    Args    :   none

=cut


sub fetch_input {
    my($self) = @_;
    
    $self->throw("No input id") unless defined($self->input_id);

    my $slice_str = $self->input_id;
    my ($chr, $start, $end, $sgp) = $slice_str =~ m!$SLICE_INPUT_ID_REGEX!;

    $self->db->assembly_type($sgp) if $sgp;

    my $slice = $self->db->get_SliceAdaptor->fetch_by_chr_start_end($chr, $start, $end);

    $self->throw("Unable to fetch slice") unless $slice;
    $self->query($slice);

    my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::CPG(
              -query   => $self->query,
              -cpg     => $self->analysis->program_file,
              -args    => $self->arguments
    );

    $self->runnable($runnable);

    return 1;
}


sub write_output {
    my($self) = @_;

    my $db  = $self->db;
    my $sfa = $self->db->get_SimpleFeatureAdaptor;
    
    my @mapped_features;
  
    my $slice = $self->query;

    foreach my $f ($self->output) {

	$f->analysis($self->analysis);
	$f->contig($slice);
	my @mapped = $f->transform;

        if (@mapped == 0) {
	    $self->warn("Couldn't map $f - skipping");
	    next;
        }
        if (@mapped == 1 && $mapped[0]->isa("Bio::EnsEMBL::Maper::Gap")) {
	    $self->warn("$f seems to be on a gap - something bad has happened ...");
	    next;
        }

	push @mapped_features, $mapped[0];

    }
    $sfa->store(@mapped_features) if @mapped_features;

    return 1;
}


=head2 fetch_output

    Title   :   fetch_output
    Usage   :   $self->fetch_output($file_name);
    Function:   Fetches output data from a frozen perl object
                stored in file $file_name
    Returns :   array of repeats (with start and end)
    Args    :   none

=cut


1;
