#
#
# Cared for by Michele Clamp  <michele@sanger.ac.uk>
#
# Copyright Michele Clamp
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code
#
# Modified 11.2001 by SCP to run on Virtual Contigs

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::Slice_EPCR

=head1 SYNOPSIS

my $db   = Bio::EnsEMBL::DBLoader->new($locator);
my $epcr = Bio::EnsEMBL::Pipeline::RunnableDB::Slice_EPCR->new( 
    -dbobj      => $db,
    -input_id   => $input_id,
    -analysis   => $analysis
);
$epcr->fetch_input();
$epcr->run();
$epcr->output();
$epcr->write_output(); #writes to DB

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Pipeline::Runnable::EPCR to add
functionality for reading and writing to databases. The appropriate
Bio::EnsEMBL::Analysis object must be passed for extraction
of appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
is required for database access.

=head1 CONTACT

For general Ensembl comments mail to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::Slice_EPCR;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::EPCR;
use Bio::EnsEMBL::Pipeline::Config::General;
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

    my $slice_str = $self->input_id;
    my ($chr, $start, $end, $sgp) = $slice_str =~ m!$SLICE_INPUT_ID_REGEX!;

    $self->db->assembly_type($sgp) if $sgp;

    my $slice = $self->db->get_SliceAdaptor->fetch_by_chr_start_end($chr, $start, $end);

    $self->throw("Unable to fetch contig") unless $slice;

    $self->query($slice);

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
        -epcr  => $self->analysis->program_file,
        -sts   => $sts,   # either a filename or an array reference
        -nmin  => $parameters{'-NMIN'},
        -nmax  => $parameters{'-NMAX'},
        -w     => $parameters{'-W'},
        -m     => $parameters{'-M'}
    );

    $self->runnable($runnable);

    return 1;
}


sub write_output {
    my($self) = @_;

    my $db  = $self->db;
    my $mfa = $self->db->get_MarkerFeatureAdaptor;
    
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

	# if a primer has N's at the 5' end, the reported start of
	# the STS may be in a gap region. ignoring these cases.
	# if this happens, the best solution is probably to edit
	# the primer

	push @mapped_features, $mapped[0];

    }

    # need to change adaptor to store > 1 feature at once
    foreach my $mf (@mapped_features) {
        $mfa->store($mf);
    }

    return 1;
}


1;
