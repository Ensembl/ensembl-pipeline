#!/usr/local/bin/perl -w

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
functionality to read and write to databases.
The appropriate Bio::EnsEMBL::Pipeline::Analysis object must be passed for
extraction of appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::Obj is
required for database access.

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
    Args    :   -dbobj:     A Bio::EnsEMBL::DB::Obj, 
                -input_id:   Contig input id , 
                -analysis:  A Bio::EnsEMBL::Pipeline::Analysis

=cut

sub new {
    my ($class, @args) = @_;
    my $self = bless {}, $class;
 

    $self->{'_fplist'}      = []; # ???   
    $self->{'_genseq'}      = undef;
    $self->{'_runnable'}    = undef;
    $self->{'_input_id'}    = undef;
    $self->{'_analysis'}    = undef;
    
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
    $self->analysis($analysis);
    
    $self->runnable('Bio::EnsEMBL::Pipeline::Runnable::CPG');
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
    my $contig    = $self->dbobj->get_Contig($contigid);
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
        $parameters{'-cpg'} = $self->analysis->program_file;
        #creates empty Bio::EnsEMBL::Runnable::CPG object
        $self->{'_runnable'} = $runnable->new( 
					      -clone => $parameters{'-clone'},
					      -length => $parameters{'-length'},
					      -gc => $parameters{'-gc'},
					      -oe => $parameters{'-oe'},
					      -cpg => $parameters{'-cpg'},
					     );
    }
    return $self->{'_runnable'};
}

1;
