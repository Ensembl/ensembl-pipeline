#
#
# Cared for by Thomas Down <td2@sanger.ac.uk>
#
# Based on CPG.pm by Val Curwen
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::EponineTSS

=head1 SYNOPSIS

my $db      = Bio::EnsEMBL::DBLoader->new($locator);

my $eponine = Bio::EnsEMBL::Pipeline::RunnableDB::EponineTSS->new ( 
                                   -dbobj      => $db,
			           -input_id   => $input_id
                                   -analysis   => $analysis 
                                    );

$eponine->fetch_input();

$eponine->run();

$eponine->output();

$eponine->write_output(); #writes to DB

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Pipeline::Runnable::EponineTSS to add
functionality to read and write to databases. The appropriate
Bio::EnsEMBL::Pipeline::Analysis object must be passed for extraction
of appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
is required for database access.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::EponineTSS;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::EponineTSS;
use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

=head2 new

    Title   :   new
    Usage   :   $self->new(-DBOBJ       => $db
                           -INPUT_ID    => $id
                           -ANALYSIS    => $analysis);
                           
    Function:   creates a Bio::EnsEMBL::Pipeline::RunnableDB::EponineTSS object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::EponineTSS object
    Args    :   -dbobj:     A Bio::EnsEMBL::DB::Obj, 
                -input_id:   Contig input id , 
                -analysis:  A Bio::EnsEMBL::Pipeline::Analysis

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);

    $self->{'_fplist'}      = []; # ???   
    $self->{'_genseq'}      = undef;
    $self->{'_runnable'}    = undef;
    
    $self->throw("Analysis object required") unless ($self->analysis);
    
    &Bio::EnsEMBL::Pipeline::RunnableDB::EponineTSS::runnable($self,'Bio::EnsEMBL::Pipeline::Runnable::EponineTSS');
    return $self;
}

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data for Eponine from the database
    Returns :   none
    Args    :   none

=cut

sub fetch_input {
    my( $self) = @_;
    
    $self->throw("No input id") unless defined($self->input_id);
    
    my $contigid  = $self->input_id;
    my $contig    = $self->dbobj->get_RawContigAdaptor->fetch_by_name($contigid);
    my $genseq    = $contig->primary_seq() or $self->throw("Unable to fetch contig");
   
    $self->genseq($genseq);
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
        $parameters{'-java'} = $self->analysis->program_file;
        #creates empty Bio::EnsEMBL::Runnable::EponineTSS object
        $self->{'_runnable'} = $runnable->new
	    ( '-threshold' => $parameters{'-threshold'},
	      '-epojar' => $parameters{'-epojar'},
	      '-java' => $parameters{'-java'},
	      );
    }
    return $self->{'_runnable'};
}

sub write_output{
  my ($self) = @_;

  my @features = $self->output();
  my $simple_f_a = $self->dbobj->get_SimpleFeatureAdaptor();
  my $contig;
  eval 
    {
      $contig = $self->dbobj->get_RawContigAdaptor->fetch_by_name($self->input_id);
    };

  if ($@) 
    {
      print STDERR "Contig not found, skipping writing output to db: $@\n";
    }
  foreach my $f(@features){
    $f->analysis($self->analysis);
    $simple_f_a->store($contig->dbID, $f);
  }


}


1;
