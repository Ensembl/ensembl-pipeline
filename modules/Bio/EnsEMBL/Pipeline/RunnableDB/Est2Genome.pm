#!/usr/local/bin/perl

#
#
# Cared for by EnsEMBL  <ensembl-dev@ebi.ac.uk>
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::Est2Genome

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::RunnableDB::Est2Genome->new(
					     -dbobj     => $db,
					     -input_id  => $id
                                             );
    $obj->fetch_input
    $obj->run

    my @newfeatures = $obj->output;


=head1 DESCRIPTION

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::RunnableDB::Est2Genome;

use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::RootI;
use Bio::EnsEMBL::Pipeline::Runnable::AlignFeature;


@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB Bio::Root::RootI);


sub new {
    my ($class,@args) = @_;
    my $self = bless {}, $class;

    $self->{'_fplist'} = []; #create key to an array of feature pairs
    
    my( $dbobj, $input_id, $seqfetcher ) = $self->_rearrange(['DBOBJ',
							      'INPUT_ID',
							      'SEQFETCHER'], @args);
       
    $self->throw("No database handle input")           unless defined($dbobj);
    $self->throw("[$dbobj] is not a Bio::EnsEMBL::Pipeline::DB::ObjI") unless $dbobj->isa("Bio::EnsEMBL::Pipeline::DB::ObjI");
    $self->dbobj($dbobj);

    $self->throw("No input id input") unless defined($input_id);
    $self->input_id($input_id);
    
    if(!defined $seqfetcher) {
      $seqfetcher = new Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;
    }

    $self->seqfetcher($seqfetcher);

    return $self; # success - we hope!
}

=head2 dbobj

Inherited from C<Bio::EnsEMBL::Pipeline::RunnableDBI>.

=cut 


=head2 fetch_output

    Title   :   fetch_output
    Usage   :   $self->fetch_output
    Function:   Fetchs output data from a frozen perl object
    Returns :   array of exons (with start and end)
    Args    :   none

=cut

sub fetch_output {
    my($self,$output) = @_;
    
    $output || $self->throw("Name of frozen object data file not given");
    
    my $object;
    open (IN,"<$output") || $self->throw("Could not open output data file '$output'");
    
    while (<IN>) {
	$_ =~ s/\[//;
	$_ =~ s/\]//;
	$object .= $_;
    }
    close(IN);
    my @out;
    my (@obj) = FreezeThaw::thaw($object);
    foreach my $array (@obj) {
	foreach my $object (@$array) {
	    push @out, $object;
	}
    }
    return @out;
}

=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output
    Function:   Writes output data to db
    Returns :   array of exons (with start and end)
    Args    :   none

=cut

sub write_output {
    my($self,@features) = @_;

    my $db=$self->dbobj();

    @features || $self->throw("No frozen object passed for the output");
    my $contig;
    eval {
	$contig = $db->get_Contig($self->input_id);
    };
    if ($@) {
	print STDERR "Contig not found, skipping writing output to db\n";
    }
    else {
	my $feat_Obj=Bio::EnsEMBL::DBSQL::Feature_Obj->new($db);
	$feat_Obj->write($contig,@features);
    }
    return 1;
}

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetchs input data for est2genome from the database
    Returns :   nothing
    Args    :   none

=cut

sub fetch_input {
    my( $self) = @_;
    

    $self->throw("No input id") unless defined($self->input_id);

    my $contigid  = $self->input_id;
    my $contig    = $self->dbobj->get_Contig($contigid);
    my $genseq   = $contig->primary_seq;
    my @features = $contig->get_all_SimilarityFeatures;
    my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::AlignFeature('-genomic'    => $genseq,
								      '-features'   => \@features,
								      '-seqfetcher' => $self->seqfetcher);

    $self->runnable($runnable);
}


sub runnable {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->throw("[$arg] is not a Bio::EnsEMBL::Pipeline::RunnableI") unless $arg->isa("Bio::EnsEMBL::Pipeline::RunnableI");
	
	$self->{_runnable} = $arg;
    }

    return $self->{_runnable};
}

sub run {
    my ($self) = @_;

    my $runnable = $self->runnable;
    $runnable || $self->throw("Can't run - no runnable object");

    return $runnable->minirun;
}

sub output {
    my ($self) = @_;

    my $runnable = $self->runnable;
    $runnable || $self->throw("Can't return output - no runnable object");

    return $runnable->output;
}

1;


