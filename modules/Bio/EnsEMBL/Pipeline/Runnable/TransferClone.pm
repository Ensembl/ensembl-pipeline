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

Bio::EnsEMBL::Pipeline::Runnable::TransferClone

=head1 SYNOPSIS

    my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::TransferClone
    ('-from_locator' => $locator1,
     '-to_locator'   => $locator2);
    my $status   = $runnable->run;
    my @pairs    = $runnable->output;

=head1 DESCRIPTION

Transfers a clones genes and exons from the core database to the
analysis database

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::Runnable::TransferClone;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::EnsEMBL::Root;

use Bio::EnsEMBL::Pipeline::RunnableI;

use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::DBSQL::Obj;

use Bio::EnsEMBL::Pipeline::DBSQL::Obj;
use Bio::EnsEMBL::Pipeline::DBSQL::Clone;
use Bio::EnsEMBL::Pipeline::DBSQL::Contig;

use Bio::EnsEMBL::Root;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI Bio::EnsEMBL::Root);

sub new {
    my ($class,@args) = @_;

    my $self = $class->SUPER::new(@args);

    my ($from_locator,$to_locator,$cloneid) = 
	$self->_rearrange([qw(FROM_LOCATOR
			      TO_LOCATOR
			      CLONEID
			      )],@args);

    $from_locator || $self->throw("No from locator string");
    $to_locator   || $self->throw("No to locator string");
    $cloneid      || $self->throw("No clone id input");

    $self->from_locator($from_locator);
    $self->to_locator  ($to_locator);
    $self->cloneid     ($cloneid);

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

    my $tdb    = new Bio::EnsEMBL::DBLoader($self->to_locator);
    my $fdb    = new Bio::EnsEMBL::DBLoader($self->from_locator);

    my $fclone = $fdb   ->get_Clone($self->cloneid);
    my $tclone = new Bio::EnsEMBL::Pipeline::DBSQL::Clone(-disk_id => "none",
							  -dbobj   => $tdb);
    my @contig;

    foreach my $contig ($fclone->get_all_Contigs) {

	my @features = $contig->get_all_SeqFeatures;
	my @genes    = $contig->get_all_Genes;


	my $tcontig  = new Bio::EnsEMBL::Pipeline::DBSQL::Contig(-id    => $contig->id,
								 -dbobj => $tdb);

	$tclone->add_Contig($tcontig);
	
	foreach my $f (@features) {
	    $tcontig->add_SimilarityFeature($f);
	}
	
	foreach my $gene (@genes) {
	    
	    foreach my $exon ($gene->each_unique_Exon) {
		$exon->find_supporting_evidence(\@features);
	    }
	    $tcontig->add_Gene($gene);
	}
    }
    $self->clone($tclone);
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

    return $self->clone;
}

sub to_locator {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{'_to_locator'} = $arg;
    }
    return $self->{'_to_locator'};
}

sub from_locator {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{'_from_locator'} = $arg;
    }

    return $self->{'_from_locator'};
}


sub cloneid {			
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{'_cloneid'} = $arg;
    }
    return $self->{'_cloneid'};
}

sub clone {
    my ($self,$arg) = @_;
    
    if (defined($arg)) {
	$self->{'_clone'} = $arg;
    }
    return $self->{'_clone'};
}

1;
