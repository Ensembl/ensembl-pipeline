#
#
# Cared for by Laura Clarke  <lec@sanger.ac.uk>
#
# Copyright Laura Clarke

#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB:HalfwiseHMM
=head1 SYNOPSIS

my $db          = Bio::EnsEMBL::DBLoader->new($locator);
my $halfwisehmm     = Bio::EnsEMBL::Pipeline::RunnableDB::HalfwiseHMM->new ( 
                                                    -dbobj      => $db,
			                                        -input_id   => $input_id
                                                    -analysis   => $analysis );
$halfwisehmm->fetch_input();
$halfwisehmm->run();
$halfwisehmm->output();
$halfwisehmm->write_output(); #writes to DB

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Pipeline::Runnable::HalfwiseHMM to add
functionality for reading and writing to databases. The appropriate
Bio::EnsEMBL::Analysis object must be passed for extraction
of appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
is required for database access.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::HalfwiseHMM;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::HalfwiseHMM;
use Data::Dumper;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

=head2 new

    Title   :   new
    Usage   :   $self->new(-DBOBJ       => $db
                           -INPUT_ID    => $id
                           -ANALYSIS    => $analysis);      
                          
    Function:   creates a Bio::EnsEMBL::Pipeline::RunnableDB::HalfwiseHMM object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::HalfwiseHMM object
    Args    :   -dbobj:     A Bio::EnsEMBL::DB::Obj, 
                input_id:   Contig input id , 
                -analysis:  A Bio::EnsEMBL::Analysis 

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    
    $self->{'_fplist'}      = [];
    $self->{'_genseq'}      = undef;
    $self->{'_runnable'}    = undef;
    $self->{'_parameters'}  = undef;
    
    $self->throw("Analysis object required") unless ($self->analysis);
    $self->init('Bio::EnsEMBL::Pipeline::Runnable::HalfwiseHMM');
    
    return $self;
}


=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data for HalfwiseHMM from the database
    Returns :   none
    Args    :   none

=cut

sub fetch_input {
    my($self) = @_;
    
    $self->throw("No input id") unless defined($self->input_id);

    my $contigid  = $self->input_id;
    my $contig    = $self->dbobj->get_Contig($contigid);
    my $genseq    = $contig->get_repeatmasked_seq() or $self->throw("Unable to fetch contig");
    $self->genseq($genseq);

}

=head2 init

    Title   :   init
    Usage   :   $self->init
    Function:   initializes the halfwisehmm runnable
    Returns :   nothing
    Args    :   name of runnable

=cut

#get/set for runnable and args
sub init {
    my ($self, $runnable) = @_;
    
    if ($runnable) {
      #extract parameters into a hash
      my %parameters;
      my ($parameter_string) = $self->parameters();
      if ($parameter_string) {
	my @pairs = split (/,/, $parameter_string);
	foreach my $pair (@pairs) {
	  
	  my ($key, $value) = split (/=>/, $pair);
	  $key =~ s/\s+//g;
	  $parameters{$key} = $value;
	}
	
      }
      $parameters {'-genewise'}  = $self->analysis->program_file;
      $parameters {'-hmmdb'}   = $self->analysis->db_file;
      #creates empty Bio::EnsEMBL::Runnable::HalfwiseHMM object
      my $runnable = $runnable->new(%parameters);
      $self->runnable($runnable);
    }
}


=head2 result_quality_tag

    Title   :   result_quality_tag
    Usage   :   $self->result_quality_tag
    Function:   Returns an indication of whether the data is suitable for 
                further analyses. Allows distinction between failed run and 
                no hits on a sequence.
    Returns :   'VALID' or 'INVALID'
    Args    :   none

=cut
#a method of writing back result quality
sub result_quality_tag {
    my ($self) = @_;
    
    if ($self->output)
    {
        return 'VALID';
    }
    else
    {
        return 'INVALID';
    }
}

1;
