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

Bio::EnsEMBL::Pipeline::RunnableDB::Blast

=head1 SYNOPSIS

my $db      = Bio::EnsEMBL::DBLoader->new($locator);
my $blast   = Bio::EnsEMBL::Pipeline::RunnableDB::Blast->new ( 
                                                    -dbobj      => $db,
			                                        -input_id   => $input_id
                                                    -analysis   => $analysis );
$blast->fetch_input();
$blast->run();
$blast->output();
$blast->write_output(); #writes to DB

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Pipeline::Runnable::Blast to add
functionality for reading and writing to databases.
The appropriate Bio::EnsEMBL::Analysis object must be passed for
extraction of appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::Obj is
required for databse access.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::Blast;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::Blast;

use vars qw(@ISA);

@ISA = qw (Bio::EnsEMBL::Pipeline::RunnableDB);

=head2 new

    Title   :   new
    Usage   :   $self->new(-DBOBJ       => $db
                           -INPUT_ID    => $id
                           -ANALYSIS    => $analysis);
                           
    Function:   creates a Bio::EnsEMBL::Pipeline::RunnableDB::Blast object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::Blast object
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

    my $contigid  = $self->input_id;
    my $contig    = $self->dbobj->get_Contig($contigid);
    my $genseq    = $contig->get_repeatmasked_seq() or $self->throw("Unable to fetch contig");

    print STDERR "Setting genseq to " . $genseq. "\n";

    $self->genseq($genseq);
    print STDERR "Set genseq to " . $self->genseq. "\n";
# input sequence needs to contain at least 3 consecutive nucleotides
    my $seq = $self->genseq->seq;
    if ($seq =~ /[CATG]{3}/) {
        $self->input_is_void(0);
    }
    else {
        $self->input_is_void(1);
        $self->warn("Need at least 3 nucleotides");
    }
}

#get/set for runnable and args

sub runnable {
    my ($self, @runnable) = @_;
    if (@runnable)
    {
        foreach my $runnable (@runnable)
        {
            $runnable->isa("Bio::EnsEMBL::Pipeline::RunnableI") or
                $self->throw("Input to runnable is not Bio::EnsEMBL::Pipeline::RunnableI");
        }
        push (@{$self->{'_runnable'}}, @runnable);
    }
    return @{$self->{'_runnable'}};
}

=head2 run

    Title   :   run
    Usage   :   $self->run();
    Function:   Runs Bio::EnsEMBL::Pipeline::Runnable::Blast->run()
    Returns :   none
    Args    :   none

=cut

sub run {
    my ($self) = @_;

    #need to pass one peptide at a time
    $self->throw("Input must be fetched before run") unless ($self->genseq);
    

    #extract parameters into a hash
    my ($parameter_string) = $self->analysis->parameters();
    my %parameters;
    my ($thresh, $thresh_type, $arguments);

    if ($parameter_string)
    {
        $parameter_string =~ s/\s+//g;
        my @pairs = split (/,/, $parameter_string);
        foreach my $pair (@pairs)
        {
            my ($key, $value) = split (/=>/, $pair);
            if ($key eq '-threshold_type' && $value) {
                $thresh_type = $value;
            }
            elsif ($key eq '-threshold' && $value) {
                $thresh = $value;
            }
            else
	    # remaining arguments not of '=>' form
	    # are simple flags (like -p1)
            {
                $arguments .= " $key ";
            }
        }
    }
    print $arguments."\n";
    $parameters{'-query'} = $self->genseq;
    $parameters{'-database'} = $self->analysis->db;
    $parameters{'-program'} = $self->analysis->program;
    $parameters{'-options'} = $arguments if $arguments;
    if ($thresh && $thresh_type) {
	$parameters{'-threshold'} = $thresh;
	$parameters{'-threshold_type'} = $thresh_type;
    }
    else {
	$parameters{'-threshold'} = 1e-4;
	$parameters{'-threshold_type'} = 'PVALUE';
    }

    
	my $runnable = Bio::EnsEMBL::Pipeline::Runnable::Blast->new(
	    %parameters
	);

	$runnable->run();
	$self->runnable($runnable);                                        

}

=head2 output

    Title   :   output
    Usage   :   $self->output();
    Function:   Runs Bio::EnsEMBL::Pipeline::Runnable->output()
    Returns :   An array of Bio::EnsEMBL::Repeat objects (FeaturePairs)
    Args    :   none

=cut

sub output {
    my ($self) = @_;

    my @runnable = $self->runnable;
    my @results;
    
    foreach my $runnable(@runnable){
      print STDERR "runnable = ".$runnable[0]."\n";
      push(@results, $runnable->output);
    }
    return @results;  
}

1;
