# Author: Marc Sohrmann (ms2@sanger.ac.uk)
# Copyright (c) Marc Sohrmann, 2001
# You may distribute this code under the same terms as perl itself
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

  Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Seg

=head1 SYNOPSIS

  my $seg = Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Seg->new ( -dbobj      => $db,
	    	                                                    -input_id   => $input_id,
                                                                    -analysis   => $analysis,
                                                                  );
  $seg->fetch_input;  # gets sequence from DB
  $seg->run;
  $seg->output;
  $seg->write_output; # writes features to to DB

=head1 DESCRIPTION

  This object wraps Bio::EnsEMBL::Pipeline::Runnable::Protein::Seg
  to add functionality to read and write to databases.
  A Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor is required for database access (dbobj).
  The query sequence is provided through the input_id.
  The appropriate Bio::EnsEMBL::Pipeline::Analysis object
  must be passed for extraction of parameters.

=head1 CONTACT

  Marc Sohrmann: ms2@sanger.ac.uk

=head1 APPENDIX

  The rest of the documentation details each of the object methods. 
  Internal methods are usually preceded with a _.

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Seg;

use strict;
use vars qw(@ISA);

use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::Protein::Seg;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);


=head2 new

 Title    : new
 Usage    : $self->new ( -dbobj       => $db
                         -input_id    => $id
                         -analysis    => $analysis,
                       );
 Function : creates a Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Seg object
 Example  : 
 Returns  : a Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Seg object
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
	
	my $run = Bio::EnsEMBL::Pipeline::Runnable::Protein::Seg->new(-clone     => $self->genseq,
									  -analysis  => $self->analysis	);
	
	
	$self->{'_runnable'} = $run;
    }
    
    return $self->{'_runnable'};


}

1;
