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

Bio::EnsEMBL::Pipeline::RunnableDB::Clone_CPG

=head1 SYNOPSIS

my $db      = Bio::EnsEMBL::DBLoader->new($locator);

my $cpg = Bio::EnsEMBL::Pipeline::RunnableDB::Clone_CPG->new ( 
                                            -dbobj      => $db,
         		                    -input_id   => $input_id
                                            -analysis   => $analysis 
                                            );

$cpg->fetch_input();

$cpg->run();

$cpg->output();

$cpg->write_output();

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Pipeline::Runnable::CPG to add
functionality to read and write to databases. 
This object takes clone ids, while Bio::EnsEMBL::Pipeline::RunnabdleDB::CPG 
acts on contigs. This allows us to submit one Job per clone rather than one 
per contig and should speed things up ...

The appropriate Bio::EnsEMBL::Analysis object must be passed for
extraction of appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::Obj is
required for databse access.

=head1 CONTACT

Val Curwen vac@sanger.ac.uk

=head1 APPENDIX


The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::Clone_CPG;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB::CPG;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB::CPG);

=head2 new

    Title   :   new
    Usage   :   $self->new(-DBOBJ       => $db
                           -INPUT_ID    => $id
                           -ANALYSIS    => $analysis);
                           
    Function:   creates a Bio::EnsEMBL::Pipeline::RunnableDB::Clone_CPG object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::Clone_CPG object
    Args    :   -dbobj:     A Bio::EnsEMBL::DB::Obj, 
                -input_id:   Clone input id , 
                -analysis:  A Bio::EnsEMBL::Analysis

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    
#    $self->{'_flist'}      = []; # list of features
    $self->{'_runnable'}    = []; # list of runnables, one per contig associated with input clone
    
    $self->throw("Analysis object required") unless ($self->analysis);
    return $self;
}

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data for cpg from the database
    Returns :   none
    Args    :   none

=cut

sub fetch_input {
    my( $self) = @_;
    
    $self->throw("No input id") unless defined($self->input_id);

    my $cloneid     = $self->input_id;
    my $clone       = $self->dbobj->get_Clone($cloneid);
    foreach my $contig  ($clone->get_all_Contigs())
    {       
        my $genseq
            = $contig->primary_seq() or $self->throw("Unable to fetch contig");
        $self->runnable('Bio::EnsEMBL::Pipeline::Runnable::CPG', $genseq);
    }
}

#get/set for runnable and args

sub runnable {
    my ($self, $runnable, $genseq) = @_;
    print "runnable is ", $self->{'_runnable'}, "\n";
    if ($runnable && $genseq)
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

        $parameters {'-clone'} = $genseq;
        #creates empty Bio::EnsEMBL::Runnable::CPG object
        push (@{$self->{'_runnable'}}, $runnable->new(
						      -clone => $parameters{'-clone'},
						      -length => $parameters{'-length'},
						      -gc => $parameters{'-gc'},
						      -oe => $parameters{'-oe'},
						      -cpg => $parameters{'-cpg'}
						     ));
    }
    return @{$self->{'_runnable'}};
}

=head2 run

    Title   :   run
    Usage   :   $self->run();
    Function:   Runs Bio::EnsEMBL::Pipeline::Runnable::CPG->run()
    Returns :   none
    Args    :   none

=cut

sub run {
    my ($self) = @_;
    $self->throw("Runnable modules not set") unless ($self->runnable());
    foreach my $runnable ($self->runnable())
    {
        $runnable->run();
    }
}

=head2 output

    Title   :   output
    Usage   :   $self->output();
    Function:   Runs Bio::EnsEMBL::Pipeline::Runnable::CPG->output() on each runnable
    Returns :   An array of SeqFeatures representing cpg islands
    Args    :   none

=cut

sub output {
    my ($self) = @_;
    
    my @output;
    foreach my $runnable ($self->runnable())
    {
        push (@output, $runnable->output());
    }
    return @output;
}

=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output
    Function:   Writes output data to db
    Returns :   array of SeqFeatures representing cpg islands
    Args    :   none

=cut

sub write_output {
    my($self) = @_;

    $self->throw("fetch_input must be called before write_output\n") 
        unless ($self->runnable());

    my $db=$self->dbobj();
    foreach my $runnable ($self->runnable())
    {
        my $contig;
        my @islands = $runnable->output();

        foreach my $is (@islands) {
	  $is->analysis($self->analysis);
        }
        eval 
        {
	    $contig = $db->get_Contig($runnable->clone->display_id());
        };
        if ($@) 
        {
	    print STDERR "Contig not found, skipping writing output to db\n";
        }
        elsif (@islands) 
	{
	    my $feat_adp=Bio::EnsEMBL::DBSQL::FeatureAdaptor->new($db);
	    $feat_adp->store($contig, @islands);
	}
    } 
}

1;
