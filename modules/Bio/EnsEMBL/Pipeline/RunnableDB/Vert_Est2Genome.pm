#!/usr/local/bin/perl

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

Bio::EnsEMBL::Pipeline::RunnableDB::AlignFeature

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

package Bio::EnsEMBL::Pipeline::RunnableDB::Vert_Est2Genome;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object;
use Bio::EnsEMBL::Pipeline::Runnable::AlignFeature;


@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDBI Bio::Root::Object );

sub _initialize {
    my ($self,@args) = @_;
    my $make = $self->SUPER::_initialize(@_);    
           
    $self->{'_fplist'} = []; #create key to an array of feature pairs
    
    my( $dbobj, $db2,$input_id ) = $self->_rearrange(['DBOBJ',
						      'DB2',
						      'INPUT_ID'], @args);
       
    $self->throw("No database handle input")           unless defined($dbobj);
    $self->throw("[$dbobj] is not a Bio::EnsEMBL::Pipeline::DB::ObjI") unless $dbobj->isa("Bio::EnsEMBL::Pipeline::DB::ObjI");
    if (defined($db2)) {
    $self->warn("[$db2] is not a Bio::EnsEMBL::DB::ObjI") unless $db2->isa("Bio::EnsEMBL::DBSQL::Obj");
    }
    $self->dbobj($dbobj);
    $self->db2($db2);

    $self->throw("No input id input") unless defined($input_id);
    $self->input_id($input_id);
    
    return $self; # success - we hope!
}
sub input_id {
	my ($self,$arg) = @_;

   if (defined($arg)) {
      $self->{_input_id} = $arg;
   }

   return $self->{_input_id};
}

=head2 dbobj

    Title   :   dbobj
    Usage   :   $self->dbobj($db)
    Function:   Get/set method for database handle
    Returns :   Bio::EnsEMBL::Pipeline::DB::ObjI
    Args    :   

=cut 

sub dbobj {
    my( $self, $value ) = @_;    
    if ($value) {

        $value->isa("Bio::EnsEMBL::Pipeline::DB::ObjI") || $self->throw("Input [$value] isn't a Bio::EnsEMBL::Pipeline::DB::ObjI");
        $self->{'_dbobj'} = $value;
    }
    return $self->{'_dbobj'};
}
sub db2 {
    my ($self,$value) = @_;

    if (defined($value)) {
	$self->{_db2} = $value;
    }
    return $self->{_db2};
}
=head2 fetch_output

    Title   :   fetch_output
    Usage   :   $self->fetch_output
    Function:   Fetchs output data from a frozen perl object
    Returns :   array of exons (with start and end)
    Args    :   none

=cut

sub fetch_output {
    my($self,$output) = @_;
    
    $output || $self->throw("No frozen object passed for the output");
    
    my $object;
    open (IN,"<$output") || do {print STDERR ("Could not open output data file... skipping job\n"); next;};
    
    while (<IN>) {
	$_ =~ s/\[//;
	$_ =~ s/\]//;
	$object .= $_;
    }
    close(IN);
    my @out;
   
    if (defined($object)) {
    my (@obj) = FreezeThaw::thaw($object);
    foreach my $array (@obj) {
	foreach my $object (@$array) {
	    push @out, $object;
	}
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
	$contig = $db->Contig($self->input_id);
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
    my $contig    = $self->db2->get_Contig($contigid);
    my $genseq   = $contig->seq;
    my @features = $contig->get_all_SimilarityFeatures;

    my @mrnafeatures;

    foreach my $f (@features) {
	if (defined($f->analysis) && $f->analysis->db eq "vert"  && $f->score > 1000) {
	    my $organism = $self->get_organism($f->hid);
	    if ($organism eq "Homo sapiens (human)") {
		push(@mrnafeatures,$f);
	    } else {
		print STDERR "Invalid organism $organism\n";
	    }
	}
    }
    print STDERR "Number of features is " . scalar(@mrnafeatures) . "\n";

    my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::AlignFeature(-genomic  => $genseq,
							    -features => \@mrnafeatures);

    $self->runnable($runnable);
}
sub get_organism {
    my ($self,$hid) = @_;

    $self->throw("No hid input") unless defined($hid);

    eval {
	open (ORG,"efetch $hid |");

	my $org;

	while (<ORG>) {
	    chomp;
	    if (/^OS\s+(.*)/) {
		$org = $1;
		close(ORG);
		return $org;
	    }
	}
	close(ORG);
    }; 
    if ($@) {
	$self->warn("Efetch failed for $hid. Skipping. Error was [$@]\n");
    }
    return $org;
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

    $self->throw("Can't run - no runnable object") unless defined($self->runnable);

    $self->runnable->run;
}

sub output {
    my ($self) = @_;

    $self->throw("Can't return output  - no runnable object") unless defined($self->runnable);

    return $self->runnable->output;
}

1;


