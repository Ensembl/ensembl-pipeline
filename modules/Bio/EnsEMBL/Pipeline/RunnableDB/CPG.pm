#
#
# Cared for by Val Curwen  <vac@sanger.ac.uk>
#
# Copyright Val Curwen
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::CPG

=head1 SYNOPSIS

my $db      = Bio::EnsEMBL::DBLoader->new($locator);

my $cpg = Bio::EnsEMBL::Pipeline::RunnableDB::CPG->new ( 
                                   -dbobj      => $db,
			           -input_id   => $input_id
                                   -analysis   => $analysis 
                                    );

$cpg->fetch_input();

$cpg->run();

$cpg->output();

$cpg->write_output(); #writes to DB

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Pipeline::Runnable::CPG to add
functionality to read and write to databases. The appropriate
Bio::EnsEMBL::Analysis object must be passed for extraction
of appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
is required for database access.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::CPG;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::CPG;
use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

=head2 new

    Title   :   new
    Usage   :   $self->new(-DBOBJ       => $db
                           -INPUT_ID    => $id
                           -ANALYSIS    => $analysis);
                           
    Function:   creates a Bio::EnsEMBL::Pipeline::RunnableDB::CPG object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::CPG object
    Args    :   -dbobj:     A Bio::EnsEMBL::DBSQL::DBAdaptor, 
                -input_id:   Contig input id , 
                -analysis:  A Bio::EnsEMBL::Analysis

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);

    $self->{'_fplist'}      = []; # ???   
    $self->{'_genseq'}      = undef;
    $self->{'_runnable'}    = undef;
    
    $self->throw("Analysis object required") unless ($self->analysis);
    
    &Bio::EnsEMBL::Pipeline::RunnableDB::CPG::runnable($self,'Bio::EnsEMBL::Pipeline::Runnable::CPG');
    return $self;
}

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data for CPG from the database
    Returns :   none
    Args    :   none

=cut

sub fetch_input {
    my( $self) = @_;
    
    $self->throw("No input id") unless defined($self->input_id);
    my $contigid  = $self->input_id;
    my $contig    = $self->db->get_RawContigAdaptor->fetch_by_name($contigid);
    $self->genseq($contig);
}

#get/set for runnable and args
sub runnable {
    my ($self, $runnable) = @_;
    if ($runnable)
    {
        #extract parameters into a hash
        my ($parameter_string) = $self->parameters() ;
        my %parameters;
        if ($parameter_string)
        {
            $parameter_string =~ s/\s+//g;
            my @pairs = split (/,/, $parameter_string);
            
            foreach my $pair (@pairs)
            {
                my ($key, $value) = split (/=>/, $pair);
                $parameters{$key} = $value;
            }
        }
        $parameters{'-cpg'} = $self->analysis->program_file;
        #creates empty Bio::EnsEMBL::Runnable::CPG object
        $self->{'_runnable'} = $runnable->new
	    ( '-query' => $parameters{'-query'},
	      '-length' => $parameters{'-length'},
	      '-gc' => $parameters{'-gc'},
	      '-oe' => $parameters{'-oe'},
	      '-cpg' => $parameters{'-cpg'},
	      );
    }
    return $self->{'_runnable'};
}

sub write_output{
  my ($self) = @_;

  my @features = $self->output();
  my $simple_f_a = $self->db->get_SimpleFeatureAdaptor();
  my $contig;
  eval {
    $contig = $self->db->get_RawContigAdaptor->fetch_by_name($self->input_id);
  };

  if ($@) {
      print STDERR "Contig not found, skipping writing output to db: $@\n";
      return;
  }
  foreach my $f(@features){
    $f->analysis($self->analysis);
    $f->attach_seq($contig);
    $simple_f_a->store($f);
  }


}
1;
