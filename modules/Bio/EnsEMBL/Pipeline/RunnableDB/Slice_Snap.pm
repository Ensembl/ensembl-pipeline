#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code
#
# author: mongin@ebi.ac.uk

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::Slice_Snap

=head1 SYNOPSIS

  my $db      = Bio::EnsEMBL::DBLoader->new($locator);
  my $Snap = Bio::EnsEMBL::Pipeline::RunnableDB::Slice_Snap->new ( 
                                                    -dbobj      => $db,
			                            -input_id   => $input_id
                                                    -analysis   => $analysis );
  $snap->fetch_input();
  $snap->run();
  $snap->output();
  $snap->write_output(); #writes to DB

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Pipeline::Runnable::Snap to add
functionality to read and write to databases.
The appropriate Bio::EnsEMBL::Analysis object must be passed for
extraction of appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::Obj is
required for databse access.

snap is a gene predictor written by Ian Korf (ik1@sanger.ac.uk) part the Zoe software library.

=head1 CONTACT

For general queries please contact <ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::Slice_Snap;

use strict;

use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::Snap;
use Bio::EnsEMBL::Pipeline::Config::General;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);



=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data for Snap from the database
    Returns :   none
    Args    :   none

=cut

sub fetch_input {
    my ($self) = @_;
    
    $self->throw("No input id") unless defined($self->input_id);

    my $slice_str = $self->input_id;
    my ($chr, $start, $end, $sgp) = $slice_str =~ m!$SLICE_INPUT_ID_REGEX!;

    $self->db->assembly_type($sgp) if $sgp;

    my $slice = $self->db->get_SliceAdaptor->fetch_by_chr_start_end($chr, $start, $end) or $self->throw("Unable ot fetch Slice");;

    my $genseq;

    if (@$SNAP_MASKING) {
	$genseq = $slice->get_repeatmasked_seq($SNAP_MASKING) or $self->throw("Unable ot fetch seq");;
    }
    else {
	
	$genseq = $slice;
    }

    $self->throw("Unable to fetch virtual contig") unless $slice;
    #print STDERR $slice->seq."\n";
    $self->query($genseq);
    

    my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::Snap(
							      -query   => $self->query,
							      -snap => $self->analysis->program_file,
							      -hmmfile  => $self->analysis->db_file,
							      -args    => $self->arguments
							      );
    
    $self->runnable($runnable);

    return 1;
}



sub write_output {
    my($self) = @_;
    my $contig;

    my $db  = $self->db;

    my $pred_adp = $self->db->get_PredictionTranscriptAdaptor;

    my $slice = $self->query;

    
    my @mapped_features;

    GENE: foreach my $f ($self->output) {
	$f->analysis($self->analysis);
	
	eval {
	    $f->transform
	    };
	if ($@) {
	    print STDERR "Transcript does not map back, won't be loaded: $@\n";
	    next GENE;
	}
	    

	print STDERR "FE: ".$f."\n";
	push(@mapped_features,$f);

    }

    $pred_adp->store(@mapped_features) if @mapped_features;

    return 1;
}

1;






