#
# BioPerl module for Blastp.pm
#
# Cared for by Emmanuel Mongin <mongin@ebi.ac.uk>
#
# Copyright Emmanuel Mongin
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Prints.pm - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Blastp;
use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::Blast;
use Bio::EnsEMBL::DBSQL::Protein_Adaptor;
use Bio::EnsEMBL::DBSQL::Protein_Feature_Adaptor;


@ISA = qw (Bio::EnsEMBL::Pipeline::RunnableDB);


=head2 new

    Title   :   new
    Usage   :   $self->new(-DBOBJ       => $db
                           -INPUT_ID    => $id
                           -ANALYSIS    => $analysis);
                           
    Function:   creates a Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Profile object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::Blast object
    Args    :   -dbobj:     A Bio::EnsEMBL::DB::Obj, 
                -input_id:   Contig input id , 
                -analysis:  A Bio::EnsEMBL::Analysis 

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    
    $self->{'_fplist'}      = [];
    $self->{'_pepseq'}      = undef;
    $self->{'_runnable'}    = undef;            
    return $self;
}


=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data for repeatmasker from the database
    Returns :   none
    Args    :   none

=cut

sub fetch_input {
    my($self) = @_;
    
    $self->throw("No input id") unless defined($self->input_id);

    

    my $translriptid  = $self->input_id;
    my $prot_adapt = $self->dbobj->get_Protein_Adaptor();
    
    my $prot = $prot_adapt->fetch_Protein_by_dbid($self->input_id);

    my $pepseq    = $prot->seq;


    my $peptide  =  Bio::PrimarySeq->new(  '-seq'         => $pepseq,
					   '-id'          => $self->input_id,
					   '-accession'   => $self->input_id,
					   '-moltype'     => 'protein');

    $self->genseq($peptide);

# input sequence needs to contain at least 3 consecutive nucleotides
    my $seq = $self->genseq;
}

#get/set for runnable and args
sub runnable {
    my ($self) = @_;
    
    if (!defined($self->{'_runnable'})) {
      my $run = Bio::EnsEMBL::Pipeline::Runnable::Blast->new(-query     => $self->genseq,
							     -database  => $self->analysis->db,
							     -program   => $self->analysis->program);
 
           
      $self->{'_runnable'} = $run;
    }
    
    return $self->{'_runnable'};
}


=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output
    Function:   Writes output data to db
    Returns :   array of repeats (with start and end)
    Args    :   none

=cut

sub write_output {
    my($self) = @_;

    my @features = $self->output();

     my $feat_Obj= $self->dbobj->get_Protfeat_Adaptor;

    foreach my $feat(@features) {
	
	$feat_Obj->write_Protein_feature($feat);
    }

    return 1;
}

sub output {
    my ($self) = @_;

    my $runnable = $self->runnable;
    $runnable || $self->throw("Can't return output - no runnable object");

    return $runnable->output;
}

1;










