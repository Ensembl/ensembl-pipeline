#
# BioPerl module for Protein_Annotation.pm
#
# Cared for by ensembl <ensembl-dev@ebi.ac.uk>
#
# Copyright ensembl
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Protein_Annotation.pm - DESCRIPTION of Object

=head1 SYNOPSIS

this is the base class for the Protein_annotation runnabledbs and shouldn't really be created on its own as it doesn't do anything other than reduce code duplication

=head1 DESCRIPTION

 This object wraps Bio::EnsEMBL::Pipeline::Runnable::Protein::Seg
  to add functionality to read and write to databases.
  A Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor is required for database access (db).
  The query sequence is provided through the input_id.
  The appropriate Bio::EnsEMBL::Analysis object
  must be passed for extraction of parameters.

=head1 CONTACT

ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Protein_Annotation;
use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::DBSQL::ProteinAdaptor;
use Bio::EnsEMBL::DBSQL::ProteinFeatureAdaptor;
use Bio::EnsEMBL::Pipeline::Config::Protein_Annotation::General;

@ISA = qw (Bio::EnsEMBL::Pipeline::RunnableDB);

=head2 new

    Title   :   new
    Usage   :   $self->new(-DB       => $db
                           -INPUT_ID    => $id
                           -ANALYSIS    => $analysis);
                           
    Function:   creates a Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Prints object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::Blast object
    Args    :   -db:     A Bio::EnsEMBL::DBSQL::DBAdaptor, 
                -input_id:   Contig input id , 
                -analysis:  A Bio::EnsEMBL::Analysis 

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    $self->throw ("Analysis object required") unless ($self->analysis);
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
  my $proteinAdaptor = $self->db->get_ProteinAdaptor;
  my $prot;
  my $peptide;
  
  my $input_id;
 
  if($self->input_id eq 'proteome'){
    $input_id = $PA_PEPTIDE_FILE;
    $self->throw($PA_PEPTIDE_FILE." doesn't exist\n") unless(-e $PA_PEPTIDE_FILE);
  }elsif($self->input_id =~ /chunk/){
    $input_id = $PA_CHUNKS_DIR."/".$self->input_id;
    $self->throw($input_id." doesn't exist\n") unless(-e $input_id);
  }else{

    eval {
	$prot = $proteinAdaptor->fetch_by_translation_id ($self->input_id);
    };
    if(($@) || (!$prot)){
      $self->throw($self->input_id." either isn't a transcript dbID or doesn't exist in the database : $@\n");
    }
   
	#The id is a protein id, that's fine, create a PrimarySeq object
    my $pepseq    = $prot->seq;
    $input_id  =  Bio::PrimarySeq->new(  '-seq'         => $pepseq,
					 '-id'          => $self->input_id,
					 '-accession'   => $self->input_id,
					 '-moltype'     => 'protein');
  }

    

    
  $self->query($input_id);
}


sub write_output {
    my($self) = @_;

    my @features = $self->output();
    
     my $feat_Obj=$self->db->get_ProteinFeatureAdaptor();  

    foreach my $feat(@features) {
	
	$feat_Obj->store($feat);
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
    $self->throw("Input not fetched")      unless ($self->query());

    $self->runnable->run($dir);
}



sub output {
    my ($self) = @_;

    my $runnable = $self->runnable;
    $runnable || $self->throw("Can't return output - no runnable object");

    return $runnable->output;
}

1;
