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

Bio::EnsEMBL::Pipeline::RunnableDB::tRNAscan_SE

=head1 SYNOPSIS

my $db      = Bio::EnsEMBL::DBLoader->new($locator);
my $trna = Bio::EnsEMBL::Pipeline::RunnableDB::tRNAscan_SE->new ( 
                                                    -dbobj      => $db,
			                            -input_id   => $input_id
                                                    -analysis   => $analysis );
$trna->fetch_input();
$trna->run();
$trna->output();
$trna->write_output(); #writes to DB

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Pipeline::Runnable::tRNAscan_SE to add
functionality to read and write to databases.
The appropriate Bio::EnsEMBL::Analysis object must be passed for
extraction of appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::Obj is
required for database access.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::tRNAscan_SE;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::tRNAscan_SE;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

=head2 new

    Title   :   new
    Usage   :   $self->new(-DBOBJ       => $db
                           -INPUT_ID    => $id
                           -ANALYSIS    => $analysis);
                           
    Function:   creates a Bio::EnsEMBL::Pipeline::RunnableDB::tRNAscan_SE object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::tRNAscan_SE object
    Args    :   -dbobj:     A Bio::EnsEMBL::DB::Obj, 
                -input_id:   Contig input id , 
                -analysis:  A Bio::EnsEMBL::Analysis

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);

    $self->{'_fplist'}      = [];  
    $self->{'_genseq'}      = undef;
    $self->{'_runnable'}    = undef;
    

    $self->throw("Analysis object required") unless ($self->analysis);
    
    # set up cpg specific parameters
    my $params = $self->parameters();
    if ($params ne "") { $params .= " , "; }
    $params .= "-tRNAscan_SE=>".$self->analysis->program_file();
    $self->parameters($params);

    $self->runnable('Bio::EnsEMBL::Pipeline::Runnable::tRNAscan_SE');
    return $self;
}

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data for tRNAscan_SE from the database
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

=head2 runnable

    Title   :   runnable
    Usage   :   $self->runnable($arg)
    Function:   Sets a runnable for this RunnableDB
    Returns :   Bio::EnsEMBL::Pipeline::RunnableI
    Args    :   Bio::EnsEMBL::Pipeline::RunnableI

=cut
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
#            $parameter_string =~ s/\s+//g;
            my @pairs = split (/,/, $parameter_string);
            foreach my $pair (@pairs)
            {
                my ($key, $value) = split (/=>/, $pair);
		$key =~ s/\s+//g;
                $parameters{$key} = $value;
            }
        }

        $self->{'_runnable'} = $runnable->new(%parameters);
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
