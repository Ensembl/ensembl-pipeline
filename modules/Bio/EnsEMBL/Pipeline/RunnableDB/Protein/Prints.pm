
#
# BioPerl module for Prints.pm
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

 $self->new(-DBOBJ       => $db
                           -INPUT_ID    => $id
                           -ANALYSIS    => $analysis);

Where the analysis id can be either a translation internal id or the location of a file. That\'s feature is used only for the protein annotation. Some of the protein analysis are extremely fast to run (eg: seg) in that cases a full protein dataset will be given to the RunnableDB (running a protein at a time would be to expensive).

=head1 DESCRIPTION

 This object wraps Bio::EnsEMBL::Pipeline::Runnable::Protein::Seg
  to add functionality to read and write to databases.
  A Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor is required for database access (dbobj).
  The query sequence is provided through the input_id.
  The appropriate Bio::EnsEMBL::Pipeline::Analysis object
  must be passed for extraction of parameters.

=head1 CONTACT

mongin@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Prints;
use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::Protein::Prints;
use Bio::EnsEMBL::DBSQL::Protein_Adaptor;
use Bio::EnsEMBL::DBSQL::Protein_Feature_Adaptor;


@ISA = qw (Bio::EnsEMBL::Pipeline::RunnableDB);


=head2 new

    Title   :   new
    Usage   :   $self->new(-DBOBJ       => $db
                           -INPUT_ID    => $id
                           -ANALYSIS    => $analysis);
                           
    Function:   creates a Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Prints object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::Blast object
    Args    :   -dbobj:     A Bio::EnsEMBL::DB::Obj, 
                -input_id:   Contig input id , 
                -analysis:  A Bio::EnsEMBL::Pipeline::Analysis 

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

#get/set for runnable and args
sub runnable {
  my ($self) = @_;
    
    if (!defined($self->{'_runnable'})) {
	
	my $run = Bio::EnsEMBL::Pipeline::Runnable::Protein::Prints->new(-clone     => $self->genseq,
									  -analysis  => $self->analysis	);
	
	
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
    
     my $feat_Obj=$self->dbobj->get_Protfeat_Adaptor;  

    foreach my $feat(@features) {
	
	$feat_Obj->write_Protein_feature($feat);
    }

    return 1;
}



=head2 run

    Title   :   run
    Usage   :   $self->run();
    Function:   Runs Bio::EnsEMBL::Pipeline::Runnable::Protein::Prints->run()
    Returns :   none
    Args    :   none

=cut

sub run {
    my ($self,$dir) = @_;
    $self->throw("Runnable module not set") unless ($self->runnable());
    $self->throw("Input not fetched")      unless ($self->genseq());

    $self->runnable->run($dir);
}



sub output {
    my ($self) = @_;

    my $runnable = $self->runnable;
    $runnable || $self->throw("Can't return output - no runnable object");

    return $runnable->output;
}

1;










