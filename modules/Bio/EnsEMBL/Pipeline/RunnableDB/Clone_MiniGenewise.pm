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

    my $obj = Bio::EnsEMBL::Pipeline::RunnableDB::Clone_MiniGenewise->new(-dbobj     => $db,
									  -input_id  => $id);

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

package Bio::EnsEMBL::Pipeline::RunnableDB::Clone_MiniGenewise;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object;
use Bio::EnsEMBL::Pipeline::RunnableDBI;
use Bio::EnsEMBL::Pipeline::RunnableDB::MiniGenewise;
use Bio::EnsEMBL::Analysis::MSPcrunch;
use Bio::SeqIO;
use Bio::Tools::Blast;

use Data::Dumper;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDBI Bio::Root::Object );

sub _initialize {
    my ($self,@args) = @_;
    my $make = $self->SUPER::_initialize(@_);    
           
    $self->{'_fplist'} = []; #create key to an array of feature pairs
    
    my( $dbobj,$input_id ) = $self->_rearrange(['DBOBJ',
						'INPUT_ID'], @args);
       
    $self->throw("No database handle input")           unless defined($dbobj);
    $self->throw("[$dbobj] is not a Bio::EnsEMBL::DB::ObjI") unless $dbobj->isa("Bio::EnsEMBL::DB::ObjI");
    $self->dbobj($dbobj);

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
	
        $value->isa("Bio::EnsEMBL::DB::ObjI") || $self->throw("Input [$value] isn't a Bio::EnsEMBL::DB::ObjI");
        $self->{'_dbobj'} = $value;
    }
    return $self->{'_dbobj'};
}

=head2 fetch_output

    Title   :   fetch_output
    Usage   :   $self->fetch_output($file_name);
    Function:   Fetchs output data from a frozen perl object
                stored in file $file_name
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

    my $db       = $self->dbobj();

    @features || $self->throw("No frozen object passed for the output");
   
    my %contighash;

    foreach my $f (@features) {
        print STDERR "Writing feature " . $f->seqname . "\n";
	my $contigid = $f->seqname; 

        eval {
          if (! defined($contighash{$contigid})) {
  	  my $contig = $db->get_Contig($contigid);
          $contighash{$contigid} = $contig;
          $contighash{$contigid}{features} = [];
          }
        };

        if ($@) {
   	   $self->throw("Contig $contigid not found, skipping writing output to db [$@]");
        } else {
            print STDERR $contighash{$contigid}{features} . "\n";
           push(@{$contighash{$contigid}{features}},$f);
        }
    }

    my $feat_Obj=Bio::EnsEMBL::DBSQL::Feature_Obj->new($db);

    foreach my $con (keys %contighash) {
        if (scalar(@{$contighash{$con}{features}}) > 0) {
          print(STDERR "Number of features for $con is " . scalar({$contighash{$con}{features}}) . "\n");
          $feat_Obj->write($contighash{$con},@{$contighash{$con}{features}});
        }
    }
    return 1;
}

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data for est2genome from the database
    Returns :   nothing
    Args    :   none

=cut

sub fetch_input {
    my( $self) = @_;
    

    $self->throw("No input id") unless defined($self->input_id);

    my $cloneid  = $self->input_id;
    my $clone    = $self->dbobj->get_Clone($cloneid);
    my @contigs  = $clone->get_all_Contigs;

    foreach my $contig (@contigs) {
	my $runnable = new Bio::EnsEMBL::Pipeline::RunnableDB::MiniGenewise(-dbobj    => $self->dbobj,
									    -input_id => $contig->id);
	$runnable->fetch_input;
	$self->add_ContigRunnable($runnable);
    }
    
}


sub run {
    my ($self) = @_;

    foreach my $contigrunnable ($self->each_ContigRunnable) {
	$contigrunnable->run;
    }
    
}

sub add_ContigRunnable {
    my ($self,$runnable) = @_;

    if (defined($runnable)) {
	if (!($runnable->isa("Bio::EnsEMBL::Pipeline::RunnableDBI"))) {
	    $self->throw("Runnable must be Bio::EnsEMBL::Pipeline::RunnableDBI");
	} 
	if (!defined($self->{_contigrunnable})) {
	    $self->{_contigrunnable} = [];
	}
	print STDERR "Adding runnable to " . $self->{_contigrunnable} . "\n";
	push(@{$self->{_contigrunnable}},$runnable);
    } else {
	$self->throw("No runnable input to add_ContigRunnable");
    }
}

sub each_ContigRunnable {
    my ($self) = @_;

    if (!defined($self->{_contigrunnable})) {
	$self->{_contigrunnable} = [];
    }
    return @{$self->{_contigrunnable}};
}
    
sub output {
    my ($self) = @_;

    my @output;

    foreach my $runnable ($self->each_ContigRunnable) {
	eval {
	    push(@output,$runnable->output);
	};
	if ($@) {
	    print STDERR "No output for $runnable [$@]\n";
	}
    }
    return @output;
}


1;






