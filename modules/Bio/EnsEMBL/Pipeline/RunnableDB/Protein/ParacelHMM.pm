#
# Ensembl module for ParacelHMM
#
# Cared for by Emmanuel Mongin <mongin@ebi.ac.uk>
#
# Copyright Emmanuel Mongin
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

ParacelHMM - DESCRIPTION of Object

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

    
package Bio::EnsEMBL::Pipeline::RunnableDB::Protein::ParacelHMM;

use strict;
use vars qw(@ISA);

use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::Protein::ParacelHMM;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);


=head2 new

 Title    : new
 Usage    : $self->new ( -dbobj       => $db
                         -input_id    => $id
                         -analysis    => $analysis,
                       );
 Function : creates a Bio::EnsEMBL::Pipeline::RunnableDB::Protein::ParacelHMM object
 Example  : 
 Returns  : a Bio::EnsEMBL::Pipeline::RunnableDB::Protein::ParacelHMM object
 Args     : -dbobj    :  a Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
            -input_id :  input id
            -analysis :  a Bio::EnsEMBL::Pipeline::Analysis
 Throws   :

=cut

sub new {
    my ($class, @args) = @_;

    # this new method also parses the @args arguments,
    # and verifies that -dbobj and -input_id have been assigned
    my $self = $class->SUPER::new(@args);
    $self->throw ("Analysis object required") unless ($self->analysis);
    $self->{'_genseq'}      = undef;
    $self->{'_runnable'}    = undef;
    
    # set up seg specific parameters,
    # my $params = $self->parameters;  # we don't have to read the parameters column from the database
                                       # in this case; no parameters are passed on to the Runnable

    return $self;
}

# IO methods

=head2 fetch_input

 Title    : fetch_input
 Usage    : $self->fetch_input
 Function : fetches the query sequence from the database
 Example  :
 Returns  :
 Args     :
 Throws   :

=cut

sub fetch_input {
    my ($self) = @_;
    my $proteinAdaptor = $self->dbobj->get_Protein_Adaptor;
    my $prot;
    my $peptide;

    eval {
	$prot = $proteinAdaptor->fetch_Protein_by_dbid ($self->input_id);
    };
    
    if (!$@) {
	#The id is a protein id, that's fine, create a PrimarySeq object
	#This module can work with single sequence object althought is completly useless to use the paracel box in that way. 
	my $pepseq    = $prot->seq;
	$peptide  =  Bio::PrimarySeq->new(  '-seq'         => $pepseq,
					    '-id'          => $self->input_id,
					    '-accession'   => $self->input_id,
					    '-moltype'     => 'protein');
    }

    else {
	#An error has been returned...2 solution, either the input is a peptide file and we can go on or its completly rubish and we throw an exeption.
	
	
	#Check if the file exists, if not throw an exeption 
	$self->throw ("The input_id given is neither a protein id nor an existing file") unless (-e $self->input_id);
	$peptide = $self->input_id;
    }

    
    $self->genseq($peptide);
}



=head2 write_output

 Title    : write_output
 Usage    : $self->write_output
 Function : writes the features to the database
 Example  :
 Returns  :
 Args     :
 Throws   :

=cut

sub write_output {
    my ($self) = @_;
    my $proteinFeatureAdaptor = $self->dbobj->get_Protfeat_Adaptor;
    my @features = $self->output;
    
    foreach my $feat(@features) {
	$proteinFeatureAdaptor->write_Protein_feature($feat);
    }

}

=head2 run

    Title   :   run
    Usage   :   $self->run();
    Function:   Runs Bio::EnsEMBL::Pipeline::Runnable::Protein::Seq->run()
    Returns :   none
    Args    :   none

=cut

sub run {
    my ($self,$dir) = @_;
    $self->throw("Runnable module not set") unless ($self->runnable());
    $self->throw("Input not fetched")      unless ($self->genseq());

    $self->runnable->run($dir);
}


# runnable method

=head2 runnable

 Title    :  runnable
 Usage    :  $self->runnable($arg)
 Function :  sets a runnable for this RunnableDB
 Example  :
 Returns  :  Bio::EnsEMBL::Pipeline::RunnableI
 Args     :  Bio::EnsEMBL::Pipeline::RunnableI
 Throws   :

=cut

sub runnable {
    my ($self) = @_;
    
    if (!defined($self->{'_runnable'})) {
	
	my $run = Bio::EnsEMBL::Pipeline::Runnable::Protein::ParacelHMM->new(-clone     => $self->genseq,
									     -analysis  => $self->analysis	);
	
	
	$self->{'_runnable'} = $run;
    }
    
    return $self->{'_runnable'};
}

1;


