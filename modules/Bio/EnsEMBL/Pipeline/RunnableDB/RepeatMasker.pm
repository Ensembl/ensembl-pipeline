#
#
# Cared for by Michele Clamp  <michele@sanger.ac.uk>
#
# Copyright Michele Clamp
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::RepeatMasker

=head1 SYNOPSIS

my $db      = Bio::EnsEMBL::DBLoader->new($locator);
my $repmask = Bio::EnsEMBL::Pipeline::RunnableDB::RepeatMasker->new ( 
                                                    -dbobj      => $db,
			                                        -input_id   => $input_id
                                                    -analysis   => $analysis );
$repmask->fetch_input();
$repmask->run();
$repmask->output();
$repmask->write_output(); #writes to DB

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Pipeline::Runnable::RepeatMasker to add
functionality to read and write to databases.
The appropriate Bio::EnsEMBL::Analysis object must be passed for
extraction of appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::Obj is
required for databse access.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::RepeatMasker;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::RepeatMasker;
use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

=head2 new

    Title   :   new
    Usage   :   $self->new(-DBOBJ       => $db
                           -INPUT_ID    => $id
                           -ANALYSIS    => $analysis);
                           
    Function:   creates a Bio::EnsEMBL::Pipeline::RunnableDB::RepeatMasker object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::RepeatMasker object
    Args    :    -db:     A Bio::EnsEMBL::DB::Obj, 
                input_id:   Contig input id , 
                -analysis:  A Bio::EnsEMBL::Analysis

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    
    $self->{'_fplist'}      = [];
    $self->{'_genseq'}      = undef;
    $self->{'_runnable'}    = undef;
    
    # db input_id mandatory and read in by BlastableDB
    # anlaysis not mandatory for BlastableDB, so we check here 
    $self->throw("Analysis object required") unless ($self->analysis);
    
    &Bio::EnsEMBL::Pipeline::RunnableDB::RepeatMasker::runnable($self,'Bio::EnsEMBL::Pipeline::Runnable::RepeatMasker');
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
    my( $self) = @_;
    
    #my @times = times;
    #print STDERR "starting fetching input @times \n";
    $self->throw("No input id") unless defined($self->input_id);
    
    my $contigid  = $self->input_id;
    my $contig    = $self->db->get_RawContigAdaptor->fetch_by_name($contigid);
    #my $genseq    = $contig->seq() or $self->throw("Unable to fetch contig");
    #@times = times;
    #print STDERR "end fetching input @times\n";
    #print "have ".$genseq."\n";
    $self->genseq($contig);
}

#get/set for runnable and args
sub runnable {
    my ($self, $runnable) = @_;
    my $arguments = "";
    
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
		if ($key && $value) {
		    $parameters{$key} = $value;
		}
		else {
		    $arguments .= " $key ";
		}
            }
        }
        $parameters{'-repm'} = $self->analysis->program_file || undef;
        $parameters{'-args'} = $arguments;
        #creates empty Bio::EnsEMBL::Runnable::RepeatMasker object
        $self->{'_runnable'} = $runnable->new(%parameters);;
    }
    return $self->{'_runnable'};
}

sub write_output{
  my ($self) = @_;

  my @features = $self->output();
  my $repeat_f_a = $self->db->get_RepeatFeatureAdaptor();
  my $contig;
  eval 
    {
      $contig = $self->db->get_RawContigAdaptor->fetch_by_name($self->input_id);
    };

  if ($@) 
    {
      print STDERR "Contig not found, skipping writing output to db: $@\n";
    }
  foreach my $f(@features){
    $f->analysis($self->analysis);
    $repeat_f_a->store($contig->dbID, $f);
  }


}


1;
