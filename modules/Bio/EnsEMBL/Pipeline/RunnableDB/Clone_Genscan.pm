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

Bio::EnsEMBL::Pipeline::RunnableDB::Clone_Genscan

=head1 SYNOPSIS

my $db          = Bio::EnsEMBL::DBLoader->new($locator);
my $genscan     = Bio::EnsEMBL::Pipeline::RunnableDB::Clone_Genscan->new ( 
                                                    -dbobj      => $db,
			                            -input_id   => $input_id
                                                    -analysis   => $analysis );
$genscan->fetch_input();
$genscan->run();
$genscan->output();
$genscan->write_output(); #writes to DB

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Pipeline::Runnable::Genscan to add
functionality for reading and writing to databases. This object takes
clone ids, while Bio::EnsEMBL::Pipeline::RunnableDB::Genscan acts on
contigs. This allows us to submit one Job per clone rather than one
per contig and should speed things up ...

The appropriate Bio::EnsEMBL::Pipeline::Analysis object must be passed for
extraction of appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::Obj is
required for databse access.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::Clone_Genscan;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::RunnableDB::Genscan;
use Data::Dumper;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB::Genscan);

=head2 new

    Title   :   new
    Usage   :   $self->new(-DBOBJ       => $db
                           -INPUT_ID    => $id
                           -ANALYSIS    => $analysis);
                           
    Function:   creates a Bio::EnsEMBL::Pipeline::RunnableDB::Clone_Genscan object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::Genscan object
    Args    :   -dbobj:     A Bio::EnsEMBL::DB::Obj, 
                input_id:   Clone input id , 
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
                unless ($dbobj->isa ('Bio::EnsEMBL::DB::ObjI'));
    $self->dbobj($dbobj);
    
    $self->throw("No input id provided") unless ($input_id);
    $self->input_id($input_id);
    
    $self->throw("Analysis object required") unless ($analysis);
    $self->throw("Analysis object is not Bio::EnsEMBL::Pipeline::Analysis")
                unless ($analysis->isa("Bio::EnsEMBL::Pipeline::Analysis"));
    $self->analysis($analysis);
    
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

    my $cloneid  = $self->input_id;
    my $clone    = $self->dbobj->get_Clone($cloneid);
    foreach my $contig ($clone->get_all_Contigs()) {
	my $genseq = $contig->get_repeatmasked_seq() or
	 $self->throw("Unable to fetch contig");
	$self->runnable($genseq);
    }
}

#get/set for runnable and args
sub runnable {
    my ($self, $genseq) = @_;
    
    if ($genseq)
    {
        #extract parameters into a hash
        my %parameters;
        my ($parameter_string) = $self->parameters();
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
        $parameters {'-genscan'}  = $self->analysis->program_file;
        $parameters {'-matrix'}   = $self->analysis->db_file;
        $parameters {'-clone'}    = $genseq;
        #creates empty Bio::EnsEMBL::Runnable::Genscan object
        push @{$self->{'_runnable'}}, Bio::EnsEMBL::Pipeline::Runnable::Genscan->new(%parameters);
    }
    return @{$self->{'_runnable'}};
}

=head2 run

    Title   :   run
    Usage   :   $self->run();
    Function:   Runs Bio::EnsEMBL::Pipeline::Runnable::Genscan->run()
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
    Function:   Runs Bio::EnsEMBL::Pipeline::Runnable::Genscan->output() on each runnable
    Returns :   An array of SeqFeatures representing exons
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
    Returns :   array of exons (with start and end)
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
	my @features = $self->output();
	eval                           
	{                              
	    $contig = $db->get_Contig($runnable->clone->display_id());
	};                                                          
	if ($@)                                                     
	{      
	    print STDERR "Contig not found, skipping writing output to db: $@\n";        
	}
	elsif (@features)
	{                
	    print STDERR "Writing features to database\n";
	    my $feat_Obj = Bio::EnsEMBL::DBSQL::Feature_Obj->new($db);
	    $feat_Obj->write($contig, @features);
	}
    }                                            
    return 1;
}            

=head2 result_quality_tag

    Title   :   result_quality_tag
    Usage   :   $self->result_quality_tag
    Function:   Returns an indication of whether the data is suitable for 
                further analyses. Allows distinction between failed run and 
                no hits on a sequence.
    Returns :   'VALID' or 'INVALID'
    Args    :   none

=cut
#a method of writing back result quality
sub result_quality_tag {
    my ($self) = @_;
    
    if ($self->output)
    {
        return 'VALID';
    }
    else
    {
        return 'INVALID';
    }
}
