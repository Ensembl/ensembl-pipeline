#
# Mickeymouse implementation of RunnableI
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

Bio::EnsEMBL::Pipeline::Runnable::CloneExonPair

=head1 SYNOPSIS

    my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::CloneExonPair(-clone => $clone);
    my $status   = $runnable->run;
    my @pairs    = $runnable->output;

=head1 DESCRIPTION

Implementation of RunnableI that creates exon pairs from
exons with supporting evidence.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::Runnable::CloneExonPair;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::RootI;

use Bio::EnsEMBL::Pipeline::RunnableI;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

sub new {
    my ($class,@args) = @_;

    my $self = $class->SUPER::new(@args);

    my ($clone) = $self->_rearrange([qw(QUERY)],@args);

    $clone                                  || $self->throw("No input clone object");
    $clone->isa("Bio::EnsEMBL::DB::CloneI") || $self->throw("Clone isn't a Bio::EnsEMBL::DB::CloneI");
    
    $self->clone($clone);
    $self->{'_pairs'} = [];
    return $self;
}

=head2 run
  Title   : run
  Usage   : $self->run
  Function: Runs a ps command and stores the
            output in an array
  Returns : nothing
  Args    : none

=cut

sub run {
    my ($self) = @_;

    my $clone = $self->clone;

    foreach my $contig ($clone->get_all_Contigs) {
	my @genes    = $contig->get_all_Genes;
	my @features = $contig->get_all_SimilarityFeatures;
	
	foreach my $gene (@genes) {
	    foreach my $exon ($gene->each_unique_Exon) {
		$exon->find_supporting_evidence(\@features);
	    }
	    $contig->make_ExonPairs($gene->each_unique_Exon);
	}
    }
	$self->{'_clone'} = undef;
}


=head2 output
  Title   : output
  Usage   : my @out = $self->output
  Function: Returns the output from the ps
            command in an array of hashes
            Each element of the array contains
            details of one process
  Returns : @HASH
  Args    : none

=cut

sub output {
    my ($self) = @_;
    
    if (!(defined($self->{'_pairs'}))) {
        $self->{'_pairs'} = [];
    }
    return @{$self->{'_pairs'}};

}


=head2 add_ExonPair
  Title   : add_ExonPair
  Usage   : $self->add_ExonPair
  Function: Adds an exon pair 
  Returns : nothing
  Args    : Bio::EnsEMBL::Pipeline::ExonPair

=cut

sub add_ExonPair {
    my ($self,$arg) = @_;

    push(@{$self->{'_pairs'}},$arg);
}


=head2 clone
  Title   : clone
  Usage   : $self->clone($clone)
  Function: Get/set method for the input clone
            that we want to generate exon pairs for
  Returns : Bio::EnsEMBL::DB::CloneI
  Args    : Bio::EnsEMBL::DB::CloneI

=cut


sub clone {
    my ($self,$arg) = @_;
    
    if (defined($arg)) {
	if ($arg->isa("Bio::EnsEMBL::DB::CloneI")) {
	    $self->{'_clone'} = $arg;
	} else {
	    
	    $self->throw("[$arg] is not a Bio::EnsEMBL::DB::CloneI");
	}
    }
    
    return $self->{'_clone'};
}

1;

