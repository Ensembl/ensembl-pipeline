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

Bio::EnsEMBL::Pipeline::RunnableDB::BlastGenscanDNA

=head1 SYNOPSIS

my $db          = Bio::EnsEMBL::DBLoader->new($locator);
my $genscan     = Bio::EnsEMBL::Pipeline::RunnableDB::BlastGenscanDNA->new ( 
                                                    -dbobj      => $db,
			                            -input_id   => $input_id
                                                    -analysis   => $analysis );

$genscan->fetch_input();
$genscan->run();
$genscan->output();
$genscan->write_output(); 

=head1 DESCRIPTION

This object runs Bio::EnsEMBL::Pipeline::Runnable::Blast on peptides constructed from 
assembling genscan predicted features to peptide sequence. The resulting blast hits are
written back as FeaturePairs.

The appropriate Bio::EnsEMBL::Pipeline::Analysis object must be passed for
extraction of appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::Obj is
required for database access.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::BlastGenscanDNA;

use strict;

use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::BlastGenscanDNA;
use Bio::PrimarySeq;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

=head2 new

    Title   :   new
    Usage   :   $self->new(-DBOBJ       => $db
                           -INPUT_ID    => $id
                           -ANALYSIS    => $analysis);
                           
    Function:   creates a Bio::EnsEMBL::Pipeline::RunnableDB::BlastGenscanPep object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::BlastGenscanPep object
    Args    :   -dbobj:     A Bio::EnsEMBL::DB::Obj, 
                -input_id:  Contig input id , 
                -analysis:  A Bio::EnsEMBL::Pipeline::Analysis 

=cut

sub new {
    my ($class, @args) = @_;
    my $self = bless {}, $class;

    $self->{'_input_id'}    = undef;
    $self->{'_runnable'}    = [];
    

    $self->{'_genseq'}      = undef;
    $self->{'_transcripts'} = [];
    $self->{'_parameters'}  = undef;

    $self->{'_featurepairs'}= [];
    
    my ( $dbobj, $input_id, $analysis) = 
            $self->_rearrange (['DBOBJ', 'INPUT_ID', 'ANALYSIS'], @args);

    # Check the input parameters are valid

    $self->throw('Need database handle')     unless (defined($dbobj));
    $self->throw("No input id provided")     unless (defined($input_id));
    $self->throw("Analysis object required") unless (defined($analysis));

    $self->throw("[$dbobj] is not a Bio::EnsEMBL::DBSQL::Obj")  
                unless ($dbobj->isa ('Bio::EnsEMBL::DBSQL::Obj'));

    # $self->throw("[$dbobj] is not a Bio::EnsEMBL::Pipeline::DBSQL::Obj")  
    #             unless ($dbobj->isa ('Bio::EnsEMBL::Pipeline::DBSQL::Obj'));

    $self->throw("Analysis object is not Bio::EnsEMBL::Pipeline::Analysis")
                unless ($analysis->isa("Bio::EnsEMBL::Pipeline::Analysis"));

    # Finally set parameters in the object
    $self->dbobj   ($dbobj);
    $self->input_id($input_id);
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
    Function:   Runs Bio::EnsEMBL::Pipeline::Runnable::Blast->run()
    Returns :   none
    Args    :   none

=cut

sub run {
    my ($self) = @_;

    #need to pass one peptide at a time
    $self->throw("Input must be fetched before run") unless ($self->genseq);
    print STDERR "Running against ".scalar($self->transcripts)." predictions\n";

    foreach my $transcript ($self->transcripts) {

      my $runnable = Bio::EnsEMBL::Pipeline::Runnable::BlastGenscanDNA->new(-genomic   => $self->genseq,
									    -peptide   => $transcript,
									    -database  => $self->analysis->db,
									    -program   => $self->analysis->program,
									    -threshold => 1e-10);

      $runnable->run();
      $self->runnable($runnable);                                        

    }

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

    my @output;
    foreach my $run ($self->runnable) {
      push(@output,$run->output);
    }
    return @output;
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

    my $db       = $self->dbobj();
    my @features = $self->output();

    return if scalar(@features == 0);
      
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
        my $feat_Obj=Bio::EnsEMBL::DBSQL::Feature_Obj->new($db);
	    $feat_Obj->write($contig, @features);
    }
    return 1;
}


sub input_id  {
  my ($self,$arg) = @_;

  if (defined($arg)) {
    $self->{_input_id} = $arg;
  }
  return $self->{_input_id};
}

sub analysis {
  my ($self,$arg) = @_;

  if (defined($arg)) {
    $self->{_analysis} = $arg;
  }
  return $self->{_analysis};
}

sub genseq {
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    $self->{_genseq} = $arg;
  }
  return $self->{_genseq};
}
1;
