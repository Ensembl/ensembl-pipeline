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
The appropriate Bio::EnsEMBL::Pipeline::Analysis object must be passed for
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
use Bio::EnsEMBL::Pipeline::RunnableDBI;
use Bio::EnsEMBL::Pipeline::Runnable::Blast;
use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Pipeline::Runnable::Blast);

=head2 new

    Title   :   new
    Usage   :   $self->new(-DBOBJ       => $db
                           -INPUT_ID    => $id
                           -ANALYSIS    => $analysis);
                           
    Function:   creates a Bio::EnsEMBL::Pipeline::RunnableDB::Blast object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::Blast object
    Args    :   -dbobj:     A Bio::EnsEMBL::DB::Obj, 
                input_id:   Contig input id , 
                -analysis:  A Bio::EnsEMBL::Pipeline::Analysis 

=cut

sub new {
    my ($class, @args) = @_;
    my $self = bless {}, $class;
    
    $self->{'_fplist'}      = [];
    $self->{'_genseq'}      = undef;
    $self->{'_runnable'}    = undef;
    $self->{'_input_id'}    = undef;
    $self->{'_parameters'}  = undef;
        
    my ( $dbobj, $input_id, $analysis) = 
            $self->_rearrange (['DBOBJ', 'INPUT_ID', 'ANALYSIS'], @args);
    
    $self->throw('Need database handle') unless ($dbobj);
    $self->throw("[$dbobj] is not a Bio::EnsEMBL::DB::ObjI")  
                unless ($dbobj->isa ('Bio::EnsEMBL::DB::ObjI'));
    $self->dbobj($dbobj);
    
    $self->throw("No input id provided") unless ($input_id);
    $self->input_id($input_id);
    
    $self->throw("Analysis object required") unless ($analysis);
    $self->throw("Analysis object is not Bio::EnsEMBL::Pipeline::Analysis")
                unless ($analysis->isa("Bio::EnsEMBL::Pipeline::Analysis"));
    $self->analysis($analysis);
    
    $self->runnable('Bio::EnsEMBL::Pipeline::Runnable::Blast');
    
    return $self;
}

sub analysis {
    my ($self, $analysis) = @_;
    
    if ($analysis)
    {
        $self->throw("Not a Bio::EnsEMBL::Pipeline::Analysis object")
            unless ($analysis->isa("Bio::EnsEMBL::Pipeline::Analysis"));
        $self->{'_analysis'} = $analysis;
        $self->parameters($analysis->parameters);
    }
    return $self->{'_analysis'}
}

=head2 parameters

    Title   :   parameters
    Usage   :   $self->parameters($param);
    Function:   Gets or sets the value of parameters
    Returns :   A string containing parameters for Bio::EnsEMBL::Runnable run
    Args    :   A string containing parameters for Bio::EnsEMBL::Runnable run

=cut

sub parameters {
    my ($self, $parameters) = @_;
    $self->analysis->parameters($parameters) if ($parameters);
    return $self->analysis->parameters();
}

=head2 dbobj

    Title   :   dbobj
    Usage   :   $self->dbobj($obj);
    Function:   Gets or sets the value of dbobj
    Returns :   A Bio::EnsEMBL::Pipeline::DB::ObjI compliant object
                (which extends Bio::EnsEMBL::DB::ObjI)
    Args    :   A Bio::EnsEMBL::Pipeline::DB::ObjI compliant object

=cut

sub dbobj {
    my( $self, $value ) = @_;
    
    if ($value) 
    {
        $value->isa("Bio::EnsEMBL::DB::ObjI")
            || $self->throw("Input [$value] isn't a Bio::EnsEMBL::DB::ObjI");
        $self->{'_dbobj'} = $value;
    }
    return $self->{'_dbobj'};
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
    $self->genseq($genseq);
}

sub input_id {
    my ($self, $input) = @_;
    if ($input)
    {
        $self->{'_input_id'} = $input;
    }
    return $self->{'_input_id'};
}

sub genseq {
    my ($self, $genseq) = @_;
    if ($genseq)
    {
        $self->{'_genseq'} = $genseq;
    }
    return $self->{'_genseq'}
}

#get/set for runnable and args
sub runnable {
    my ($self, $runnable) = @_;
    
    if ($runnable)
    {
        #extract parameters into a hash
        my ($parameter_string) = $self->analysis->parameters() ;
        my %parameters;
        if ($parameter_string)
        {
            my @pairs = split (/,/, $parameter_string);
            foreach my $pair (@pairs)
            {
                my ($key, $value) = split (/=>/, $pair);
                $key =~ s/\s+//g;
                $parameters{$key} = $value;
            }
        }
        $parameters {'-db'}      = $self->analysis->db();
        $parameters {'-blast'}   = $self->analysis->program();  
        #creates empty Bio::EnsEMBL::Runnable::Blast object
        $self->{'_runnable'} = $runnable->new(%parameters);
    }
    return $self->{'_runnable'};
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
    $self->throw("Runnable module not set") unless ($self->runnable());
    $self->throw("Input not fetched") unless ($self->genseq());
    $self->runnable->clone($self->genseq());
    $self->runnable->run();
}

=head2 output

    Title   :   output
    Usage   :   $self->output();
    Function:   Runs Bio::EnsEMBL::Pipeline::Runnable::Blast->output()
    Returns :   An array of Bio::EnsEMBL::Repeat objects (FeaturePairs)
    Args    :   none

=cut

sub output {
    my ($self) = @_;

    my $runnable = $self->runnable;
    $runnable || $self->throw("Can't return output - no runnable object");

    return $runnable->output;
}

=head2 fetch_output

    Title   :   fetch_output
    Usage   :   $self->fetch_output($file_name);
    Function:   Fetches output data from a frozen perl object
                stored in file $file_name
    Returns :   array of repeats (with start and end)
    Args    :   none

=cut

#cut & pasted from Vert_Est2genome
sub fetch_output {
    my($self,$output) = @_;
    
    $output || $self->throw("No frozen object passed for the output");
    
    my $object;
    open (IN,"<$output") || do {print STDERR ("Could not open output data file... skipping job\n"); next;};
    
    while (<IN>) 
    {
        $_ =~ s/\[//;
	    $_ =~ s/\]//;
	    $object .= $_;
    }
    close(IN);
    my @out;
   
    if (defined($object)) 
    {
        my (@obj) = FreezeThaw::thaw($object);
        foreach my $array (@obj) 
        {
	        foreach my $object (@$array) 
            {
	            push @out, $object;
	        }
        }
    }
    return @out;
}

=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output
    Function:   Writes output data to db
    Returns :   array of repeats (with start and end)
    Args    :   none

=cut

sub write_output {
    my($self) = @_;

    my $db=$self->dbobj();
    my @features = $self->output();
    
    my $contig;
    eval 
    {
        $contig = $db->get_Contig($self->input_id);
    };
    
    if ($@) 
    {
	    print STDERR "Contig not found, skipping writing output to db: $@\n";
    }
    elsif (@features) 
    {
        #should add conditional for evalue here
	    print STDERR "Writing features to database\n";
        my $feat_Obj=Bio::EnsEMBL::DBSQL::Feature_Obj->new($db);
	    $feat_Obj->write($contig, @features);
    }
    return 1;
}
