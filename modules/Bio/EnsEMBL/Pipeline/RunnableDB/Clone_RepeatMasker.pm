#!/usr/local/bin/perl -w

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
The appropriate Bio::EnsEMBL::Pipeline::Analysis object must be passed for
extraction of appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::Obj is
required for databse access.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::Clone_RepeatMasker;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDBI;
use Bio::EnsEMBL::Pipeline::Runnable::RepeatMasker;
use Bio::EnsEMBL::Pipeline::RunnableDB::RepeatMasker;
use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB::RepeatMasker);

=head2 new

    Title   :   new
    Usage   :   $self->new(-DBOBJ       => $db
                           -INPUT_ID    => $id
                           -ANALYSIS    => $analysis);
                           
    Function:   creates a Bio::EnsEMBL::Pipeline::RunnableDB::RepeatMasker object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::RepeatMasker object
    Args    :    -dbobj:     A Bio::EnsEMBL::DB::Obj, 
                input_id:   Contig input id , 
                -analysis:  A Bio::EnsEMBL::Pipeline::Analysis

=cut

sub new {
    my ($class, @args) = @_;
    my $self = bless {}, $class;
    
    $self->{'_fplist'}      = [];
    $self->{'_runnable'}    = [];
    $self->{'_input_id'}    = undef;
    
    my ( $dbobj, $input_id, $analysis) = 
            $self->_rearrange (['DBOBJ', 'INPUT_ID', 'ANALYSIS'], @args);
    
    $self->throw('Need database handle') unless ($dbobj);
    $self->throw("[$dbobj] is not a Bio::EnsEMBL::DB::ObjI")  
                unless $dbobj->isa ('Bio::EnsEMBL::DB::ObjI');
    $self->dbobj($dbobj);
    
    $self->throw("No input id provided") unless ($input_id);
    $self->input_id($input_id);
    
    $self->throw("Analysis object required") unless ($analysis);
    $self->throw("Analysis object is not Bio::EnsEMBL::Pipeline::Analysis")
                unless ($analysis->isa("Bio::EnsEMBL::Pipeline::Analysis"));
    $self->set_parameters($analysis);
    
    return $self;
}

=head2 parameters

    Title   :   parameters
    Usage   :   $self->parameters($param);
    Function:   Gets or sets the value of parameters
    Returns :   A string containing parameters for Bio::EnsEMBL::Runnable run
    Args    :   A string containing parameters for Bio::EnsEMBL::Runnable run

=cut

=head2 dbobj

    Title   :   dbobj
    Usage   :   $self->dbobj($obj);
    Function:   Gets or sets the value of dbobj
    Returns :   A Bio::EnsEMBL::Pipeline::DB::ObjI compliant object
                (which extends Bio::EnsEMBL::DB::ObjI)
    Args    :   A Bio::EnsEMBL::Pipeline::DB::ObjI compliant object

=cut

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data for repeatmasker from the database
    Returns :   none
    Args    :   none

=cut

sub fetch_input {
    my( $self) = @_;
    
    $self->throw("No input id") unless defined($self->input_id);

    print STDERR "Input id " . $self->input_id . "\n";
    print STDERR "Dbobj " . $self->dbobj . "\n";
    
    my $cloneid     = $self->input_id;
    my $clone       = $self->dbobj->get_Clone($cloneid);
    foreach my $contig  ($clone->get_all_Contigs())
    {       
        my $genseq
            = $contig->primary_seq() or $self->throw("Unable to fetch contig");
        $self->runnable('Bio::EnsEMBL::Pipeline::Runnable::RepeatMasker', $genseq);
    }
}

#get/set for runnable and args

sub runnable {
    my ($self, $runnable, $genseq) = @_;
    if ($runnable && $genseq)
    {
        #extract parameters into a hash
        my ($parameter_string) = $self->parameters() ;
        $parameter_string =~ s/\s+//g;
        my @pairs = split (/,/, $parameter_string);
        my %parameters;
        foreach my $pair (@pairs)
        {
            my ($key, $value) = split (/=>/, $pair);
            $parameters{$key} = $value;
        }
        $parameters {'-clone'} = $genseq;
        #creates empty Bio::EnsEMBL::Runnable::RepeatMasker object
        push (@{$self->{'_runnable'}}, $runnable->new(%parameters));
    }
    return @{$self->{'_runnable'}};
}

=head2 run

    Title   :   run
    Usage   :   $self->run();
    Function:   Runs Bio::EnsEMBL::Pipeline::Runnable::RepeatMasker->run()
    Returns :   none
    Args    :   none

=cut

sub run {
    my ($self) = @_;
    $self->throw("Runnable modules not set") unless ($self->runnable());
    foreach my $runnable ($self->runnable)
    {
        $runnable->run();
    }
}

=head2 output

    Title   :   output
    Usage   :   $self->output();
    Function:   Runs Bio::EnsEMBL::Pipeline::Runnable::RepeatMasker->output()
    Returns :   An array of Bio::EnsEMBL::Repeat objects (FeaturePairs)
    Args    :   none

=cut

sub output {
    my ($self) = @_;
    
    my @output;
    foreach my $runnable ($self->runnable)
    {
        push (@output, $runnable->output());
    }
    return @output;
}

=head2 fetch_output

    Title   :   fetch_output
    Usage   :   $self->fetch_output($file_name);
    Function:   Fetchs output data from a frozen perl object
                stored in file $file_name
    Returns :   array of repeats (with start and end)
    Args    :   none

=cut

=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output
    Function:   Writes output data to db
    Returns :   array of repeats (with start and end)
    Args    :   none

=cut

sub write_output {
    my($self) = @_;
    
    $self->throw("fetch_input must be called before write_output\n") 
        unless ($self->runnable);

    my $db=$self->dbobj();
    foreach my $runnable ($self->runnable)
    {
        my $contig;
        my @repeats = $runnable->output();
        eval 
        {
	        $contig = $db->get_Contig($runnable->clone->display_id);
        };
        if ($@) 
        {
	        print STDERR "Contig not found, skipping writing output to db\n";
        }
        elsif (@repeats) 
        {
	        foreach my $repeat (@repeats)
            {
                print STDERR ($repeat->hseqname()."\t");
            }
            my $feat_Obj=Bio::EnsEMBL::DBSQL::Feature_Obj->new($db);
	        $feat_Obj->write($contig, @repeats);
        }
        return 1;
    } 
}
