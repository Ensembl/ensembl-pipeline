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

Bio::EnsEMBL::Pipeline::RunnableDB::Clone_Vert_Est2Genome

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::RunnableDB::Clone_Vert_Est2Genome->new(
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

package Bio::EnsEMBL::Pipeline::RunnableDB::Clone_Vert_Est2Genome;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::RootI;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::AlignFeature;
use Bio::EnsEMBL::Pipeline::RunnableDB::Vert_Est2Genome;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;
use Bio::EnsEMBL::Analysis::MSPcrunch;
use Bio::SeqIO;
use Bio::Tools::Blast;
use Bio::Root::RootI;
use Data::Dumper;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

=head2 new

    Title   :   new
    Usage   :   $self->new(-DBOBJ       => $db
                           -INPUT_ID    => $id
                           -ANALYSIS    => $analysis);
                           
    Function:   creates a Bio::EnsEMBL::Pipeline::RunnableDB::Clone_Vert_Est2Genome object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::Clone_Vert_Est2Genome object
    Args    :   -dbobj:      A Bio::EnsEMBL::DB::Obj (required), 
                -input_id:   Contig input id (required), 
                -seqfetcher: A Sequence Fetcher Object (required)
                -analysis:   A Bio::EnsEMBL::Pipeline::Analysis (optional) 
=cut

sub new {
    my ($new,@args) = @_;
    my $self = $self->SUPER::new(@args);    
           
    # dbobj, input_id, seqfetcher, and analysis objects are all set in
    # in superclass constructor (RunnableDB.pm)

    $self->{'_fplist'} = []; #create key to an array of feature pairs
    return $self;
}

=head2 input_id

    Title   :   input_id
    Usage   :   $self->input_id($input_id);
    Function:   Gets or sets the value of input_id
    Returns :   valid input id for this analysis (if set) 
    Args    :   input id for this analysis 

=head2 dbobj

    Title   :   dbobj
    Usage   :   $self->dbobj($db)
    Function:   Get/set method for database handle
    Returns :   Bio::EnsEMBL::Pipeline::DB::ObjI
    Args    :   

=head2 seqfetcher

    Title   :   seqfetcher
    Usage   :   $self->seqfetcher($seqfetcher)
    Function:   Get/set method for SeqFetcher
    Returns :   Bio::DB::RandomAccessI object
    Args    :   Bio::DB::RandomAccessI object

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
    my $analysis = $db->get_OldAnalysis(8);

    @features || $self->throw("No frozen object passed for the output");
   
    my %contighash;

    foreach my $f (@features) {
        print STDERR "Writing feature " . $f->hseqname . "\n";
	my $contigid = $f->seqname; 
        $f->feature1->analysis($analysis);
        $f->feature1->source_tag("est2genome");
        $f->feature1->primary_tag("similarity");
        $f->feature2->analysis($analysis);
        $f->feature2->source_tag("est2genome");
        $f->feature2->primary_tag("similarity");

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
	my $runnable = new Bio::EnsEMBL::Pipeline::RunnableDB::Vert_Est2Genome
	    ('-dbobj'      => $self->dbobj,
	     '-input_id'   => $contig->id,
	     '-seqfetcher' => $self->seqfetcher);
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
	push(@output,$runnable->output);
    }
    return @output;
}

1;
