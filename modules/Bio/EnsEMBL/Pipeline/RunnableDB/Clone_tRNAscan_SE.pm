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

Bio::EnsEMBL::Pipeline::RunnableDB::Clone_tRNAscan_SE

=head1 SYNOPSIS

my $db      = Bio::EnsEMBL::DBLoader->new($locator);
my $trna = Bio::EnsEMBL::Pipeline::RunnableDB::Clone_tRNAscan_SE->new ( 
                                                    -dbobj      => $db,
         		                            -input_id   => $input_id
                                                    -analysis   => $analysis );
$trna->fetch_input();
$trna->run();
$trna->output();
$trna->write_output();

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Pipeline::Runnable::tRNAscan_SE to add
functionality to read and write to databases. 
This object takes clone ids, while Bio::EnsEMBL::Pipeline::RunnabdleDB::tRNAscan_SE 
acts on contigs. This allows us to submit one Job per clone rather than one 
per contig and should speed things up ...

The appropriate Bio::EnsEMBL::Pipeline::Analysis object must be passed for
extraction of appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::Obj is
required for databse access.

=head1 CONTACT

Val Curwen vac@sanger.ac.uk

=head1 APPENDIX


The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::Clone_tRNAscan_SE;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB::tRNAscan_SE;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

=head2 new

    Title   :   new
    Usage   :   $self->new(-DBOBJ       => $db
                           -INPUT_ID    => $id
                           -ANALYSIS    => $analysis);
                           
    Function:   creates a Bio::EnsEMBL::Pipeline::RunnableDB::Clone_tRNAscan_SE object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::Clone_tRNAscan_SE object
    Args    :   -dbobj:     A Bio::EnsEMBL::DB::Obj, 
                -input_id:   Clone input id , 
                -analysis:  A Bio::EnsEMBL::Pipeline::Analysis

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    
    $self->{'_runnable'}    = []; # list of runnables, one per contig associated with input clone        
    # dbobj, analysis, input_id are set in RunnableDB parent
    $self->throw("Analysis object required") unless ($self->analysis);
    
    return $self;
}

=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output
    Function:   Writes output data to db
    Returns :   array of SeqFeatures representing 
    Args    :   none

=cut

#
# needs to be able to deal with genes, featurepairs and features.
#
sub write_output {
    my($self) = @_;
    
    $self->throw("fetch_input must be called before write_output\n") 
        unless ($self->runnable());

    my $db=$self->dbobj();
    foreach my $runnable ($self->runnable())
    {
        my $contig;
        my @tRNAs = $runnable->output();
        eval 
        {
	        $contig = $db->get_Contig($runnable->clone->display_id());
        };
        if ($@) 
        {
	        print STDERR "Contig not found, skipping writing output to db\n";
        }
        elsif (@tRNAs) 
	  {
	    my $feat_adp=Bio::EnsEMBL::DBSQL::FeatureAdaptor->new($db);
	    $feat_adp->store($contig, @tRNAs);
	  }
      } 
  }

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data for tRNAscan-SE from the database
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
        $self->runnable('Bio::EnsEMBL::Pipeline::Runnable::tRNAscan_SE', $genseq);
    }
}

#get/set for runnable and args

sub runnable {
    my ($self, $runnable, $genseq) = @_;
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
        #creates empty Bio::EnsEMBL::Runnable::tRNAscan_SE object
        push (@{$self->{'_runnable'}}, $runnable->new(%parameters));
    }
    return @{$self->{'_runnable'}};
}

=head2 run

    Title   :   run
    Usage   :   $self->run();
    Function:   Runs Bio::EnsEMBL::Pipeline::Runnable::tRNAscan_SE->run()
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
    Function:   Runs Bio::EnsEMBL::Pipeline::Runnable::tRNAscan_SE->output() 
                on each runnable
    Returns :   An array of SeqFeatures representing tRNA predictions
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
