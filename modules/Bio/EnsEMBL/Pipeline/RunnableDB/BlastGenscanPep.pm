#
#
# Cared for by Michele Clamp  <michele@sanger.ac.uk>
#
# Copyright Michele Clamp
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::BlastGenscanPep

=head1 SYNOPSIS

my $db          = Bio::EnsEMBL::DBLoader->new($locator);
my $genscan     = Bio::EnsEMBL::Pipeline::RunnableDB::BlastGenscanPep->new ( 
                                                    -dbobj      => $db,
			                            -input_id   => $input_id
                                                    -analysis   => $analysis );

$genscan->fetch_input();
$genscan->run();
$genscan->output();

=head1 DESCRIPTION

This object runs Bio::EnsEMBL::Pipeline::Runnable::Blast on peptides
constructed from assembling genscan predicted features to peptide
sequence. The resulting blast hits are written back as
DnaPepAlignFeature's.

The appropriate Bio::EnsEMBL::Analysis object must be passed for
extraction of appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::Obj is
required for database access.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _'

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::BlastGenscanPep;

use strict;

use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::BlastGenscanPep;
use Bio::PrimarySeq;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);


=head2 fetch_input

  Args       : none
  Example    : $runnable->fetch_input
  Description: Fetches input data for BlastGenscanPep and makes runnable
  Returntype : none
  Exceptions : $self->input_id is not defined
  Caller     : run_RunnableDB, Bio::EnsEMBL::Pipeline::Job

=cut

sub fetch_input {
    my($self) = @_;

    $self->throw("No input id") unless defined($self->input_id);

    my $contigid  = $self->input_id;
    my $contig    = $self->db->get_RawContigAdaptor->fetch_by_name($contigid)
        or $self->throw("Unable to find contig ($contigid)\n");

    $self->query($contig);

    my $pta = $self->db->get_PredictionTranscriptAdaptor;
    #need to get features predicted by genscan
    my @genscan_peps = @{$pta->fetch_all_by_RawContig($contig, 'Genscan')};
    $self->_transcripts(@genscan_peps);
    my @r;

    my ($thr, $thr_type);
    my %p = $self->parameter_hash;

    if (defined $p{-threshold} && defined $p{-threshold_type}) {
	$thr      = $p{-threshold};
	$thr_type = $p{-threshold_type};
    }
    else {
	$thr_type = 'PVALUE';
	$thr      = 0.001;
    }

    foreach my $t (@genscan_peps) {
        foreach my $db (split ',', ($self->analysis->db_file)) {
	    $self->runnable(Bio::EnsEMBL::Pipeline::Runnable::BlastGenscanPep->new(
                -genomic        => $self->query,
	        -peptide        => $t,
	        -database       => $db,
                -program        => $self->analysis->program_file,
                -options        => $self->arguments || undef,
		-threshold      => $thr,
		-threshold_type => $thr_type
            ));
        }
    }
    return 1;
}


=head2 _transcripts

  Args[1..]  : @Bio::EnsEMBL::PredictionTranscript
  Example    : $runnable->fetch_input
  Description: Internal method to store/retrieve transcripts
  Returntype : @Bio::EnsEMBL::PredictionTranscript
  Exceptions : arg is not a Bio::EnsEMBL::PredictionTranscript
  Caller     : Bio::EnsEMBL::Pipeline::RunnableDB::BlastGenscanPep

=cut

sub _transcripts {
    my ($self, @transcripts) = @_;

    $self->{'_transcripts'} ||= [];
    
    if (@transcripts)
    {
        foreach (@transcripts)
        {
            $self->throw("Input $_ is not a Bio::EnsEMBL::PredictionTranscript\n")
                unless $_->isa("Bio::EnsEMBL::PredictionTranscript");
        }
        push (@{$self->{'_transcripts'}}, @transcripts);
    }
    return @{$self->{'_transcripts'}};
}


1;
