#
#
# Cared for by Steve Searle  <michele@sanger.ac.uk>
#
# Copyright Steve Searle
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::HsperateGenscanPep

=head1 SYNOPSIS

my $db          = Bio::EnsEMBL::DBLoader->new($locator);
my $genscan     = Bio::EnsEMBL::Pipeline::RunnableDB::HsperateGenscanPep->new ( 
                                                    -dbobj      => $db,
			                            -input_id   => $input_id
                                                    -analysis   => $analysis );

$genscan->fetch_input();
$genscan->run();
$genscan->output();

=head1 DESCRIPTION

This object runs Bio::EnsEMBL::Pipeline::Runnable::Hsperate on peptides constructed from 
assembling genscan predicted features to peptide sequence. The resulting hsperate hits are
written back as FeaturePairs.

The appropriate Bio::EnsEMBL::Analysis object must be passed for
extraction of appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::Obj is
required for database access.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::HsperateGenscanPep;

use strict;

use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::HsperateGenscanPep;
use Bio::PrimarySeq;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

=head2 new

    Title   :   new
    Usage   :   $self->new(-DBOBJ       => $db
                           -INPUT_ID    => $id
                           -ANALYSIS    => $analysis);
                           
    Function:   creates a Bio::EnsEMBL::Pipeline::RunnableDB::HsperateGenscanPep object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::HsperateGenscanPep object
    Args    :   -dbobj:     A Bio::EnsEMBL::DB::Obj, 
                -input_id:  Contig input id , 
                -analysis:  A Bio::EnsEMBL::Analysis 

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);

    $self->{'_runnable'}    = [];    

    $self->{'_genseq'}      = undef;
    $self->{'_transcripts'} = [];
    $self->{'_parameters'}  = undef;

    $self->{'_featurepairs'}= [];

    $self->throw("Analysis object required") unless (defined($self->analysis));
    
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
    print STDERR "Fetching contig $contigid\n";
    my $contig    = $self->dbobj->get_Contig($contigid) 
        or $self->throw("Unable to find contig ($contigid)\n");
    my $genseq    = $contig->primary_seq() 
        or $self->throw("Unable to fetch contig sequence");

    $self->genseq($genseq);
    #need to get features predicted by genscan
    $self->transcripts($contig->get_genscan_peptides);
}

sub transcripts {
    my ($self, @transcripts) = @_;
    
    if (@transcripts)
    {
        foreach (@transcripts)
        {
            $self->throw("Input $_ is not a Bio::EnsEMBL::Transcript\n")
                unless $_->isa("Bio::EnsEMBL::Transcript");
        }
        push (@{$self->{'_transcripts'}}, @transcripts);
    }
    return @{$self->{'_transcripts'}};
}


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
    Function:   Runs Bio::EnsEMBL::Pipeline::Runnable::Hsperate->run()
    Returns :   none
    Args    :   none

=cut

sub run {
    my ($self) = @_;

    #need to pass one peptide at a time
    $self->throw("Input must be fetched before run") unless ($self->genseq);
    print STDERR "Running against ".scalar($self->transcripts)." predictions\n";

    #extract parameters into a hash
    my ($parameter_string) = $self->analysis->parameters();
    my %parameters;
    my ($thresh, $thresh_type, $arguments);

    if ($parameter_string)
    {
        #$parameter_string =~ s/\s+//g;
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
            elsif ($key eq '-hsperate_arguments' && $value) {
                $arguments .= $value;
            }
            else {
                $self->throw("Unable to parse parameters from $parameter_string");
            }
        }
    }

    $parameters{'-genomic'} = $self->genseq;
    $parameters{'-database'} = $self->analysis->db;
    $parameters{'-program'} = $self->analysis->program;
    $parameters{'-options'} = $arguments if $arguments;
    if ($thresh && $thresh_type) {
	$parameters{'-threshold'} = $thresh;
	$parameters{'-threshold_type'} = $thresh_type;
    }
    else {
	$parameters{'-threshold'} = 150;
	$parameters{'-threshold_type'} = 'SCORE';
    }


    #Temporary hack to ensure transcripts have unique IDs
    $self->_make_transcripts_unique();

    @{$parameters{'-transcripts'}} = $self->transcripts;

    my $runnable = Bio::EnsEMBL::Pipeline::Runnable::HsperateGenscanPep->new(
	    %parameters
    );

    $runnable->run();
    $self->runnable($runnable);                                        

}


sub _make_transcripts_unique {
    my ($self) = @_;

    my $transnum = 1;

    foreach my $transcript ($self->transcripts) {
      $transcript->temporary_id($transcript->temporary_id() . "." . $transnum++);
    }
}

=head2 output

    Title   :   output
    Usage   :   $self->output();
    Function:   Runs Bio::EnsEMBL::Pipeline::Runnable::Hsperate->output()
    Returns :   An array of Bio::EnsEMBL::FeaturePair objects (FeaturePairs)
    Args    :   none

=cut

sub output {
    my ($self) = @_;

    my @output;
    foreach my $run ($self->runnable) {
      my @tmp = $run->output;
      foreach my $f (@tmp) {
	$f->analysis($self->analysis);
      }
      push(@output,@tmp);
    }
    return @output;
}

1;
